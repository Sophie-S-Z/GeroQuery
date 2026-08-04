"""Wrapping of real biolearn clocks.

biolearn is an optional dependency and does not install on every Python version
(its `ecos` dependency has no wheel for 3.14). So the wrapping logic is exercised
against a fake module shaped like biolearn's real API — model_definitions with
`output`/`year`/`species`/`tissue`, a gallery `.get()`, models with
`.methylation_sites()` and `.predict(GeoData)`. The `live` tests below run the
real thing wherever biolearn is genuinely importable.
"""

from __future__ import annotations

import sys
import types

import numpy as np
import pandas as pd
import pytest

from geroquery.clocks import library as lib
from geroquery.clocks.registry import ClockRegistry
from geroquery.exceptions import ClockInputError

# ---- outcome mapping (no biolearn needed) ---------------------------------


@pytest.mark.parametrize(
    "output,expected",
    [
        ("Age (Years)", "chronological_age"),
        ("Age (years)", "chronological_age"),
        ("Human Cortex Age (Years)", "chronological_age"),
        ("Mortality Risk", "mortality"),
        ("Mortality Adjusted Age (Years)", "mortality_adjusted_age"),
        ("Aging Rate (Years/Year)", "pace_of_aging"),
        ("Gestational Age", "gestational_age"),
    ],
)
def test_known_outputs_map_to_the_right_outcome(output, expected):
    assert lib.map_outcome(output) == expected


@pytest.mark.parametrize("output", ["Smoking Status", "BMI", "Total Cholesterol", "Sex"])
def test_non_age_predictors_are_never_labelled_chronological_age(output):
    """biolearn's gallery includes plenty of models that are not aging clocks.
    Defaulting an unmapped output to chronological_age would present a smoking
    predictor as an age clock — exactly the confusion this metadata exists to
    prevent."""
    assert lib.map_outcome(output) != "chronological_age"


def test_unmapped_output_is_preserved_not_discarded():
    assert lib.map_outcome("Depression Risk") == "depression_risk"
    assert lib.map_outcome(None) == "unknown"


# ---- fake biolearn ---------------------------------------------------------

SITES = ["cg00000001", "cg00000002", "cg00000003"]


class _FakeGeoData:
    def __init__(self, metadata, dnam=None, rna=None):
        self.metadata = metadata
        self.dnam = dnam
        self.rna = rna


class _FakeModel:
    """Mimics biolearn's LinearMethylationModel closely enough to test the seam."""

    def __init__(self, name, attr="dnam", n_out=None):
        self.name = name
        self._attr = attr
        self._n_out = n_out

    def methylation_sites(self):
        return list(SITES)

    def predict(self, geo_data):
        matrix = getattr(geo_data, self._attr)
        if matrix is None:
            raise ValueError(f"model expected {self._attr} data")
        missing = set(SITES) - set(matrix.index)
        if missing:
            raise ValueError(f"Missing required CpG sites: {sorted(missing)}")
        # biolearn returns a DataFrame indexed by sample.
        values = matrix.loc[SITES].sum(axis=0) * 10.0
        n = self._n_out if self._n_out is not None else len(values)
        return pd.DataFrame({"Predicted": values.to_numpy()[:n]}, index=values.index[:n])


DEFINITIONS = {
    "Horvathv1": {
        "year": 2013,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://example.invalid/horvath",
        "output": "Age (Years)",
        "model": {"type": "LinearMethylationModel", "file": "Horvath1.csv"},
    },
    "DNAmPhenoAgeMortality": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "output": "Mortality Risk",
        "model": {"type": "LinearMethylationModel", "file": "pheno.csv"},
    },
    "SmokingPredictor": {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "output": "Smoking Status",
        "model": {"type": "LinearMethylationModel", "file": "smoke.csv"},
    },
    "SomeFutureModelType": {
        "year": 2030,
        "species": "Human",
        "tissue": "Blood",
        "output": "Age (Years)",
        "model": {"type": "NotAKnownBuilder", "file": "x.csv"},
    },
}


class _FakeGallery:
    def __init__(self, definitions=None, broken: set[str] = frozenset()):
        self.model_definitions = definitions if definitions is not None else DEFINITIONS
        self._broken = broken

    def get(self, name, imputation_method=None):
        if name in self._broken:
            raise FileNotFoundError(f"coefficient file for {name} is missing")
        return _FakeModel(name)


@pytest.fixture
def fake_biolearn(monkeypatch):
    """Install a fake `biolearn` package for the duration of a test."""
    gallery_holder = {"gallery": _FakeGallery()}

    pkg = types.ModuleType("biolearn")
    gallery_mod = types.ModuleType("biolearn.model_gallery")
    gallery_mod.ModelGallery = lambda *a, **k: gallery_holder["gallery"]
    data_mod = types.ModuleType("biolearn.data_library")
    data_mod.GeoData = _FakeGeoData

    monkeypatch.setitem(sys.modules, "biolearn", pkg)
    monkeypatch.setitem(sys.modules, "biolearn.model_gallery", gallery_mod)
    monkeypatch.setitem(sys.modules, "biolearn.data_library", data_mod)
    monkeypatch.setattr(lib, "biolearn_available", lambda: True)

    # get_registry() memoizes a module-level registry. Without resetting it, a
    # ClockService built inside a faked test reuses whichever registry an
    # earlier test happened to construct first -- so these tests pass alone and
    # fail in a full run.
    from geroquery.clocks import registry as registry_module

    monkeypatch.setattr(registry_module, "_REGISTRY", None)
    return gallery_holder


def _matrix(n=4):
    rng = np.random.default_rng(3)
    return pd.DataFrame(
        rng.uniform(0.1, 0.9, size=(n, len(SITES))),
        columns=SITES,
        index=[f"S{i}" for i in range(n)],
    )


# ---- discovery -------------------------------------------------------------


def test_clocks_are_namespaced_and_carry_their_declared_outcome(fake_biolearn):
    clocks = lib.library_clocks()
    assert "biolearn:Horvathv1" in clocks
    assert clocks["biolearn:Horvathv1"].info.predicted_outcome == "chronological_age"
    assert clocks["biolearn:DNAmPhenoAgeMortality"].info.predicted_outcome == "mortality"
    assert clocks["biolearn:SmokingPredictor"].info.predicted_outcome == "smoking_status"


def test_unknown_model_types_are_skipped_not_registered_broken(fake_biolearn):
    """We cannot know which GeoData slot an unrecognized builder reads, so
    registering it would produce a clock that raises the moment it is called."""
    assert "biolearn:SomeFutureModelType" not in lib.library_clocks()


def test_a_model_that_fails_to_load_does_not_break_the_others(fake_biolearn):
    fake_biolearn["gallery"] = _FakeGallery(broken={"Horvathv1"})
    clocks = lib.library_clocks()
    assert "biolearn:Horvathv1" not in clocks
    assert "biolearn:DNAmPhenoAgeMortality" in clocks


def test_metadata_records_library_units_and_training_population(fake_biolearn):
    info = lib.library_clocks()["biolearn:Horvathv1"].info
    assert info.library == "biolearn"
    assert info.units == "years"
    assert info.input_type == "methylation"
    assert "Human" in info.training_population and "2013" in info.training_population
    assert "Age (Years)" in info.notes


def test_advertised_feature_list_is_capped(fake_biolearn, monkeypatch):
    """A methylation clock needs tens of thousands of CpGs; listing them all in
    every /v1/clocks response would dwarf the payload."""
    many = [f"cg{i:08d}" for i in range(500)]
    monkeypatch.setattr(_FakeModel, "methylation_sites", lambda self: list(many))

    clock = lib.library_clocks()["biolearn:Horvathv1"]
    assert len(clock.info.required_features) == lib.MAX_ADVERTISED_FEATURES
    # ...but validation still uses the complete set.
    assert len(clock.required_feature_set) == 500


# ---- prediction ------------------------------------------------------------


def test_predict_transposes_samples_to_columns(fake_biolearn):
    """GeroQuery passes samples as rows, biolearn wants them as columns. Getting
    this backwards silently produces one prediction per CpG."""
    clock = lib.library_clocks()["biolearn:Horvathv1"]
    matrix = _matrix(n=4)
    out = clock.predict(matrix)
    assert out.shape == (4,)
    expected = matrix[SITES].sum(axis=1).to_numpy() * 10.0
    np.testing.assert_allclose(out, expected)


def test_predict_reports_missing_features_as_a_clock_input_error(fake_biolearn):
    clock = lib.library_clocks()["biolearn:Horvathv1"]
    with pytest.raises(ClockInputError) as exc:
        clock.predict(_matrix()[SITES[:1]])
    assert exc.value.detail["n_missing"] == 2
    assert exc.value.detail["n_required"] == 3


def test_library_error_is_translated_not_leaked(fake_biolearn, monkeypatch):
    """biolearn raises bare ValueError. Callers of the API need the domain error
    so it maps to a 422 instead of a 500."""
    clock = lib.library_clocks()["biolearn:Horvathv1"]

    def _boom(self, geo_data):
        raise ValueError("something upstream went wrong")

    monkeypatch.setattr(_FakeModel, "predict", _boom)
    with pytest.raises(ClockInputError, match="could not run on this input"):
        clock.predict(_matrix())


def test_wrong_length_result_is_rejected(fake_biolearn, monkeypatch):
    """A silent length mismatch would misalign every prediction with its sample."""
    monkeypatch.setattr(_FakeModel, "predict", lambda self, g: pd.DataFrame({"P": [1.0, 2.0]}))
    clock = lib.library_clocks()["biolearn:Horvathv1"]
    with pytest.raises(ClockInputError, match="different number of predictions"):
        clock.predict(_matrix(n=4))


def test_series_return_is_accepted(fake_biolearn, monkeypatch):
    monkeypatch.setattr(_FakeModel, "predict", lambda self, g: pd.Series([1.0, 2.0, 3.0, 4.0]))
    clock = lib.library_clocks()["biolearn:Horvathv1"]
    np.testing.assert_allclose(clock.predict(_matrix(n=4)), [1.0, 2.0, 3.0, 4.0])


# ---- registry integration --------------------------------------------------


def test_registry_merges_library_clocks_with_reference_clocks(fake_biolearn):
    registry = ClockRegistry()
    ids = {c.clock_id for c in registry.list_clocks()}
    assert "clinical_phenoage_demo" in ids  # reference tier survives
    assert "biolearn:Horvathv1" in ids
    assert registry.get("biolearn:Horvathv1").info.library == "biolearn"


def test_every_clock_declares_a_predicted_outcome(fake_biolearn):
    for info in ClockRegistry().list_clocks():
        assert info.predicted_outcome
        assert info.training_population


def test_registry_loads_without_biolearn(monkeypatch):
    """The optional integrations must never be able to break the registry.

    Both library tiers are disabled here, not just biolearn: on a machine with
    pyaging installed the registry legitimately holds 170+ more clocks, and this
    test is about the floor the reference tier guarantees.
    """
    from geroquery.clocks import pyaging_clocks as pyac_mod

    monkeypatch.setattr(lib, "biolearn_available", lambda: False)
    monkeypatch.setattr(pyac_mod, "pyaging_available", lambda: False)
    assert lib.library_clocks() == {}
    ids = {c.clock_id for c in ClockRegistry().list_clocks()}
    # phenoage ships with the package: real published coefficients, no optional
    # dependency. It is present whether or not biolearn/pyaging are installed.
    assert ids == {
        "phenoage",
        "clinical_phenoage_demo",
        "clinical_mortality_demo",
        "transcriptomic_demo",
    }


def test_a_broken_gallery_degrades_to_reference_clocks_only(monkeypatch):
    monkeypatch.setattr(lib, "biolearn_available", lambda: True)
    broken = types.ModuleType("biolearn.model_gallery")

    def _explode(*a, **k):
        raise RuntimeError("gallery data files are corrupt")

    broken.ModelGallery = _explode
    monkeypatch.setitem(sys.modules, "biolearn.model_gallery", broken)
    assert lib.library_clocks() == {}


# ---- live ------------------------------------------------------------------


@pytest.mark.live
@pytest.mark.skipif(not lib.biolearn_available(), reason="biolearn not installed")
def test_live_real_biolearn_gallery_wraps_horvath():
    clocks = lib.library_clocks()
    assert clocks, "biolearn is installed but no clocks were wrapped"
    horvath = clocks.get("biolearn:Horvathv1")
    assert horvath is not None
    assert horvath.info.predicted_outcome == "chronological_age"
    assert len(horvath.required_feature_set) > 300  # Horvath uses 353 CpGs


@pytest.mark.live
@pytest.mark.skipif(not lib.biolearn_available(), reason="biolearn not installed")
def test_live_no_real_clock_is_mislabelled_as_chronological_age():
    for clock in lib.library_clocks().values():
        declared = clock.info.notes
        if "Smoking" in clock.info.clock_id or "BMI" in clock.info.clock_id:
            assert clock.info.predicted_outcome != "chronological_age", declared


# ---- status reporting ------------------------------------------------------


def test_status_names_the_reason_when_biolearn_is_absent(monkeypatch):
    monkeypatch.setattr(lib, "biolearn_available", lambda: False)
    clocks, status = lib.load_library_clocks()
    assert clocks == {}
    assert status.installed is False and status.usable is False
    assert "not installed" in status.reason


def test_status_surfaces_an_import_failure_instead_of_swallowing_it(monkeypatch):
    """Regression: this returned {} silently, so a missing transitive dependency
    (biolearn imports torch and seaborn at module scope without declaring them)
    was indistinguishable from 'biolearn is not installed'. Both produced an
    empty clock list and no explanation."""
    monkeypatch.setattr(lib, "biolearn_available", lambda: True)
    broken = types.ModuleType("biolearn.model_gallery")

    def _explode(*a, **k):
        raise ModuleNotFoundError("No module named 'seaborn'")

    broken.ModelGallery = _explode
    monkeypatch.setitem(sys.modules, "biolearn.model_gallery", broken)

    clocks, status = lib.load_library_clocks()
    assert clocks == {}
    assert status.installed is True and status.usable is False
    assert "seaborn" in status.reason


def test_status_reports_skipped_models(fake_biolearn):
    _clocks, status = lib.load_library_clocks()
    assert status.usable is True
    # SomeFutureModelType has an unknown builder and cannot be wrapped.
    assert status.n_skipped == 1
    assert "SomeFutureModelType" in status.reason


def test_registry_exposes_library_status(fake_biolearn):
    registry = ClockRegistry()
    assert registry.library_status.usable is True
    assert registry.library_status.n_clocks == len(DEFINITIONS) - 1


# ---- caller data is not mutated --------------------------------------------


def test_predict_does_not_mutate_the_callers_matrix(fake_biolearn, monkeypatch):
    """biolearn's LinearModel.predict does `matrix_data.loc["intercept"] = 1`.
    Without a defensive copy that lands in the caller's DataFrame."""

    def _mutating_predict(self, geo_data):
        geo_data.dnam.loc["intercept"] = 1.0
        return pd.DataFrame({"P": np.zeros(geo_data.dnam.shape[1])}, index=geo_data.dnam.columns)

    monkeypatch.setattr(_FakeModel, "predict", _mutating_predict)
    clock = lib.library_clocks()["biolearn:Horvathv1"]
    matrix = _matrix(n=4)
    before = matrix.copy()

    clock.predict(matrix)
    pd.testing.assert_frame_equal(matrix, before)


@pytest.mark.skipif(
    pd.__version__ >= "3",
    reason="pandas 3 makes copy-on-write mandatory and always returns read-only "
    "arrays; this is why the `clocks` extra pins pandas<3",
)
def test_predict_hands_the_library_a_writable_buffer(fake_biolearn, monkeypatch):
    """DunedinPACE quantile-normalizes in place. A read-only numpy buffer fails
    with 'assignment destination is read-only'."""
    seen = {}

    def _writes_in_place(self, geo_data):
        arr = geo_data.dnam.to_numpy()
        seen["writeable"] = arr.flags.writeable
        return pd.DataFrame({"P": np.zeros(geo_data.dnam.shape[1])}, index=geo_data.dnam.columns)

    monkeypatch.setattr(_FakeModel, "predict", _writes_in_place)
    lib.library_clocks()["biolearn:Horvathv1"].predict(_matrix(n=4))
    assert seen["writeable"] is True


# ---- covariate metadata ----------------------------------------------------


def test_metadata_is_forwarded_and_aligned_to_the_matrix(fake_biolearn, monkeypatch):
    """GrimAge is trained on age and sex alongside methylation and raises
    without them. The first version of this wrapper always passed an empty
    frame, which made every covariate-dependent clock unusable."""
    captured = {}

    def _needs_metadata(self, geo_data):
        captured["meta"] = geo_data.metadata
        if "sex" not in geo_data.metadata.columns:
            raise ValueError("GrimAge requires 'sex' column in metadata")
        return pd.DataFrame({"P": np.zeros(geo_data.dnam.shape[1])}, index=geo_data.dnam.columns)

    monkeypatch.setattr(_FakeModel, "predict", _needs_metadata)
    clock = lib.library_clocks()["biolearn:Horvathv1"]
    matrix = _matrix(n=4)

    with pytest.raises(ClockInputError, match="requires 'sex'"):
        clock.predict(matrix)

    # Deliberately out of order: it must be reindexed onto the matrix, not zipped.
    meta = pd.DataFrame(
        {"sex": [2, 1, 2, 1], "age": [40, 30, 60, 50]},
        index=["S1", "S0", "S3", "S2"],
    )
    clock.predict(matrix, metadata=meta)
    assert list(captured["meta"].index) == list(matrix.index)
    assert list(captured["meta"]["age"]) == [30, 40, 50, 60]


def test_service_folds_chronological_age_into_clock_metadata(fake_biolearn, monkeypatch):
    from geroquery.clocks.service import ClockService

    captured = {}

    def _capture(self, geo_data):
        captured["meta"] = geo_data.metadata
        return pd.DataFrame(
            {"P": np.arange(geo_data.dnam.shape[1], dtype=float)},
            index=geo_data.dnam.columns,
        )

    monkeypatch.setattr(_FakeModel, "predict", _capture)
    service = ClockService()
    matrix = _matrix(n=4)
    service.apply_clock(
        "biolearn:Horvathv1",
        matrix,
        chronological_age=[20.0, 30.0, 40.0, 50.0],
        sample_metadata=pd.DataFrame({"sex": [1, 2, 1, 2]}, index=matrix.index),
    )
    assert list(captured["meta"]["age"]) == [20.0, 30.0, 40.0, 50.0]
    assert list(captured["meta"]["sex"]) == [1, 2, 1, 2]


def test_reference_clocks_ignore_metadata(fake_biolearn):
    """Reference clocks take a matrix alone. Passing covariates through the
    service must not break them."""
    from geroquery.clocks.registry import CLINICAL_FEATURES
    from geroquery.clocks.service import ClockService

    matrix = pd.DataFrame({f: [1.0, 2.0] for f in CLINICAL_FEATURES})
    result = ClockService().apply_clock(
        "clinical_phenoage_demo",
        matrix,
        chronological_age=[50.0, 60.0],
        sample_metadata=pd.DataFrame({"sex": [1, 2]}, index=matrix.index),
    )
    assert result.n_samples == 2


# ---- pyaging tier ----------------------------------------------------------

from geroquery.clocks import pyaging_clocks as pyac  # noqa: E402


@pytest.mark.parametrize(
    "predicts,expected",
    [
        (["chronological age"], "chronological_age"),
        (["transcriptomic age"], "chronological_age"),
        (["mortality risk"], "mortality"),
        (["pace of aging"], "pace_of_aging"),
        (["gestational age"], "gestational_age"),
        (["mitotic age"], "mitotic_age"),
        ("chronological age", "chronological_age"),
    ],
)
def test_pyaging_outcomes_map(predicts, expected):
    assert pyac.map_outcome(predicts) == expected


@pytest.mark.parametrize(
    "predicts", [["smoking exposure"], ["grip strength"], ["VO2max"], ["frailty risk"]]
)
def test_pyaging_non_age_predictors_keep_their_own_label(predicts):
    """pyaging's catalog includes grip strength, VO2max and smoking exposure.
    None of those are aging clocks."""
    assert pyac.map_outcome(predicts) != "chronological_age"


def test_pyaging_outcome_of_empty_or_missing():
    assert pyac.map_outcome([]) == "unknown"
    assert pyac.map_outcome(None) == "unknown"


def test_pyaging_status_when_not_installed(monkeypatch):
    monkeypatch.setattr(pyac, "pyaging_available", lambda: False)
    clocks, status = pyac.load_pyaging_clocks()
    assert clocks == {} and status.installed is False
    assert "not installed" in status.reason


def test_pyaging_status_when_catalog_unavailable(monkeypatch):
    """Offline with a cold cache must degrade to 'no pyaging clocks', with the
    reason stated — not to an exception that takes the registry down."""
    monkeypatch.setattr(pyac, "pyaging_available", lambda: True)
    monkeypatch.setattr(
        pyac, "_catalog", lambda allow_network: (_ for _ in ()).throw(OSError("no cache"))
    )
    clocks, status = pyac.load_pyaging_clocks(allow_network=False)
    assert clocks == {}
    assert status.installed is True and status.usable is False
    assert "GEROQUERY_ALLOW_NETWORK" in status.reason


def test_pyaging_registration_is_metadata_only(monkeypatch):
    """Feature lists are not known until an artifact is downloaded, so clocks
    register with an empty required_features rather than triggering a
    multi-hundred-megabyte fetch just to answer GET /v1/clocks."""
    catalog = {
        "horvath2013": {
            "data_type": "DNA methylation",
            "species": "Homo sapiens",
            "year": 2013,
            "tissue": ["multi-tissue"],
            "predicts": ["chronological age"],
            "citation": "Horvath 2013",
        },
        "grimage2packyrs": {
            "data_type": "DNA methylation",
            "species": "Homo sapiens",
            "predicts": ["smoking exposure"],
        },
    }
    monkeypatch.setattr(pyac, "pyaging_available", lambda: True)
    monkeypatch.setattr(pyac, "_catalog", lambda allow_network: catalog)

    clocks, status = pyac.load_pyaging_clocks(allow_network=True)
    assert status.usable is True and status.n_clocks == 2

    horvath = clocks["pyaging:horvath2013"]
    assert horvath.info.predicted_outcome == "chronological_age"
    assert horvath.info.input_type == "methylation"
    assert horvath.info.units == "years"
    assert horvath.required_feature_set == ()
    assert "Homo sapiens" in horvath.info.training_population

    assert clocks["pyaging:grimage2packyrs"].info.predicted_outcome == "smoking_exposure"


def test_pyaging_and_biolearn_ids_cannot_collide(monkeypatch):
    """Both libraries ship a Horvath 2013. Namespaced ids keep them distinct so
    one silently overwrites the other in the registry."""
    monkeypatch.setattr(pyac, "pyaging_available", lambda: True)
    monkeypatch.setattr(
        pyac,
        "_catalog",
        lambda allow_network: {
            "horvath2013": {"data_type": "DNA methylation", "predicts": ["chronological age"]}
        },
    )
    py_clocks, _ = pyac.load_pyaging_clocks(allow_network=True)
    assert set(py_clocks) == {"pyaging:horvath2013"}
    assert not set(py_clocks) & set(lib.library_clocks())


@pytest.mark.live
@pytest.mark.skipif(
    not (lib.biolearn_available() and pyac.pyaging_available()),
    reason="needs both biolearn and pyaging installed",
)
def test_live_biolearn_and_pyaging_agree_on_horvath_2013():
    """The strongest validation available for these wrappers.

    Two independently-maintained libraries implement the same published clock.
    Fed byte-identical input through two different code paths — biolearn wants
    features-as-rows, pyaging wants samples-as-rows — they must land on the same
    number. A transpose bug, a feature-ordering bug, or a wrong output column in
    either wrapper breaks this and essentially nothing else catches it.
    """
    bl_clocks, _ = lib.load_library_clocks()
    py_clocks, _ = pyac.load_pyaging_clocks(allow_network=True)

    sites = sorted(bl_clocks["biolearn:Horvathv1"].required_feature_set)
    rng = np.random.default_rng(0)
    matrix = pd.DataFrame(
        rng.beta(2, 2, size=(5, len(sites))),
        columns=sites,
        index=[f"S{i}" for i in range(5)],
    )

    from_biolearn = bl_clocks["biolearn:Horvathv1"].predict(matrix)
    from_pyaging = py_clocks["pyaging:horvath2013"].predict(matrix)

    np.testing.assert_allclose(from_biolearn, from_pyaging, rtol=0, atol=0.05)


@pytest.mark.live
@pytest.mark.skipif(not pyac.pyaging_available(), reason="pyaging not installed")
def test_live_pyaging_refuses_when_no_features_match():
    """pyaging raises a bare, message-less NameError when every feature is
    missing. Unwrapped that reads as an internal bug rather than 'your columns
    are not CpG ids'."""
    py_clocks, _ = pyac.load_pyaging_clocks(allow_network=True)
    rng = np.random.default_rng(1)
    bogus = pd.DataFrame(
        rng.beta(2, 2, size=(3, 40)),
        columns=[f"not_a_cpg_{i}" for i in range(40)],
        index=["a", "b", "c"],
    )
    with pytest.raises(ClockInputError, match="none of its required features"):
        py_clocks["pyaging:horvath2013"].predict(bogus)


@pytest.mark.live
@pytest.mark.skipif(
    not (lib.biolearn_available() and pyac.pyaging_available()),
    reason="needs both biolearn and pyaging installed",
)
def test_live_pyaging_refuses_a_mostly_imputed_prediction():
    """pyaging imputes absent features from reference values instead of failing.
    Past a threshold the output describes pyaging's reference cohort, not these
    samples, so it must not be returned as a measurement."""
    bl_clocks, _ = lib.load_library_clocks()
    py_clocks, _ = pyac.load_pyaging_clocks(allow_network=True)

    sites = sorted(bl_clocks["biolearn:Horvathv1"].required_feature_set)[:20]
    rng = np.random.default_rng(2)
    sparse = pd.DataFrame(
        rng.beta(2, 2, size=(3, len(sites))), columns=sites, index=["a", "b", "c"]
    )
    with pytest.raises(ClockInputError, match="impute"):
        py_clocks["pyaging:horvath2013"].predict(sparse)
