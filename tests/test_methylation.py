"""GEO methylation series parsing and the PhenoAge clock.

Offline: the series matrix is a gzipped fixture with the real file's structure.
What is under test is the preparation, because that is where the failures were —
a clock that "requires features not present in the input" turned out to be a
data-preparation bug, not an incompatible clock.

``ruff: noqa: E501`` for the fixture below. A GEO series matrix puts one physical
line per characteristic, so those rows are long by construction; wrapping them
would make the fixture stop resembling the file it imitates.
"""

# ruff: noqa: E501

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

from geroquery.clocks import phenoage
from geroquery.clocks.registry import PHENOAGE_INFO, PhenoAgeClock, get_registry
from geroquery.exceptions import ClockInputError, SourceError
from geroquery.sources import methylation as meth
from geroquery.sources.manifest import METHYLATION_PANEL

SERIES_MATRIX = """!Series_title\t"Test methylomes"
!Series_platform_id\t"GPL13534"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\t"GSM5"\t"GSM6"\t"GSM7"
!Sample_characteristics_ch1\t"age: 25.5"\t"age: 31.2"\t"age: 44.0"\t"age: 71.8"\t"age: Newborn"\t"age: 80 years"\t"age: NA"
!Sample_characteristics_ch1\t"Sex: female"\t"Sex: male"\t"Sex: female"\t"Sex: male"\t"Sex: female"\t"Sex: male"\t"Sex: female"
!Sample_characteristics_ch1\t"tissue: whole blood"\t"tissue: whole blood"\t"tissue: whole blood"\t"tissue: whole blood"\t"tissue: whole blood"\t"tissue: whole blood"\t"tissue: whole blood"
!Sample_characteristics_ch1\t"disease status: Control"\t"disease status: Control"\t"disease status: Control"\t"disease status: Control"\t"disease status: Control"\t"disease status: Case"\t"disease status: Control"
!Sample_characteristics_ch1\t"dna methylation age: 26.1"\t"dna methylation age: 30.4"\t"dna methylation age: 45.9"\t"dna methylation age: 69.2"\t"dna methylation age: 0.4"\t"dna methylation age: 78.0"\t"dna methylation age: NA"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\t"GSM5"\t"GSM6"\t"GSM7"
"cg00000029"\t0.10\t0.12\t0.20\t0.31\t0.05\t0.33\t0.11
"cg00000108"\t0.80\t0.79\t0.70\t0.61\t0.90\t0.60\t0.81
"cg00000109"\t0.50\tnull\t0.52\t0.55\t0.48\t0.56\t0.51
"cg99999999"\tnull\tnull\tnull\tnull\tnull\tnull\tnull
!series_matrix_table_end
"""

TEST_SERIES = meth.MethylationSeries(
    accession="GSE00000",
    title="Test methylomes",
    platform=meth.PLATFORM_450K,
    n_samples=7,
    tissue="whole blood",
    characteristics={
        "age": "age",
        "Sex": "sex",
        "tissue": "tissue",
        "disease status": "disease_status",
        "dna methylation age": "published_horvath_age",
    },
    control_filter=("disease_status", "Control"),
    reference_age_column="published_horvath_age",
)


@pytest.fixture
def matrix_path(tmp_path):
    path = tmp_path / "GSE00000_series_matrix.txt.gz"
    path.write_bytes(gzip.compress(SERIES_MATRIX.encode("utf-8")))
    return path


# ---- age parsing ----------------------------------------------------------


@pytest.mark.parametrize(
    ("label", "years"),
    [
        ("25.5", 25.5),
        ("80 years", 80.0),
        ("97 yrs", 97.0),
        ("Newborn", 0.0),
        ("cord blood", 0.0),
    ],
)
def test_parse_age(label, years):
    assert meth.parse_age(label) == pytest.approx(years)


@pytest.mark.parametrize("label", ["NA", "", "unknown", None, "adult"])
def test_unparseable_age_is_none_not_zero(label):
    """A zero would put a newborn's worth of error into a clock's MAE for every
    sample whose age was simply not recorded."""
    assert meth.parse_age(label) is None


# ---- parsing --------------------------------------------------------------


def test_read_series_restricts_to_controls_and_parseable_ages(matrix_path):
    ds = meth.read_series(matrix_path, TEST_SERIES)
    # GSM6 is a case; GSM7 has no age.
    assert list(ds.metadata.index) == ["GSM1", "GSM2", "GSM3", "GSM4", "GSM5"]
    assert ds.dropped["samples_not_control"] == 1
    assert ds.dropped["samples_no_age"] == 1
    assert ds.metadata["age"].tolist() == [25.5, 31.2, 44.0, 71.8, 0.0]
    assert ds.metadata["published_horvath_age"].iloc[0] == pytest.approx(26.1)


def test_all_missing_cpgs_are_dropped_but_partial_ones_are_imputed(matrix_path):
    """Dropping every CpG with a gap deletes probes the clocks are built on:
    three of Horvath's coefficients are among GSE30870's incomplete CpGs."""
    ds = meth.read_series(matrix_path, TEST_SERIES)
    assert "cg99999999" not in ds.betas.index  # nothing to impute from
    assert "cg00000109" in ds.betas.index  # one gap, imputed
    assert ds.dropped["cpgs_all_missing"] == 1
    assert ds.dropped["cpgs_with_any_gap"] == 1
    assert not np.isnan(ds.betas.to_numpy()).any()
    # Imputed with the CpG's mean across the retained samples.
    retained = [0.50, 0.52, 0.55, 0.48]
    assert ds.betas.loc["cg00000109", "GSM2"] == pytest.approx(np.mean(retained), abs=1e-5)


def test_imputation_can_be_switched_off(matrix_path):
    ds = meth.read_series(matrix_path, TEST_SERIES, impute=False)
    assert np.isnan(ds.betas.loc["cg00000109", "GSM2"])


def test_matrix_is_samples_as_rows(matrix_path):
    """The clock wrappers take samples as rows. Handing them the stored
    CpG-by-sample orientation is what made the first full panel run fail 236
    times with 'requires features not present'."""
    ds = meth.read_series(matrix_path, TEST_SERIES)
    assert list(ds.matrix().index) == list(ds.metadata.index)
    assert "cg00000029" in ds.matrix().columns


def test_cpg_allow_list_narrows_the_read(matrix_path):
    ds = meth.read_series(matrix_path, TEST_SERIES, cpgs=["cg00000029"])
    assert list(ds.betas.index) == ["cg00000029"]


def test_a_renamed_characteristic_fails_loudly(matrix_path, tmp_path):
    """A silently-absent key would yield a frame of NaN ages that looks fine
    until every clock reports a nonsense error."""
    broken = tmp_path / "broken.txt.gz"
    broken.write_bytes(
        gzip.compress(SERIES_MATRIX.replace("dna methylation age:", "epigenetic age:").encode())
    )
    with pytest.raises(SourceError, match="no longer in the series matrix"):
        meth.read_series(broken, TEST_SERIES)


def test_an_unexpected_platform_is_refused(matrix_path, tmp_path):
    other = tmp_path / "other.txt.gz"
    other.write_bytes(gzip.compress(SERIES_MATRIX.replace("GPL13534", "GPL8490").encode()))
    with pytest.raises(SourceError, match="expected platform"):
        meth.read_series(other, TEST_SERIES)


def test_panel_series_are_declared_and_pinned():
    """Every series the adapter knows how to read must be a pinned artifact, and
    the URLs must agree, or a verified download is a different file."""
    for accession, series in meth.SERIES.items():
        assert accession in METHYLATION_PANEL
        assert series.url == METHYLATION_PANEL[accession].url


# ---- PhenoAge -------------------------------------------------------------


def _phenoage_frame(n: int = 6) -> pd.DataFrame:
    rng = np.random.default_rng(3)
    return pd.DataFrame(
        {
            "albumin": rng.normal(4.3, 0.2, n),
            "creatinine": rng.normal(0.9, 0.1, n),
            "glucose": rng.normal(95.0, 8.0, n),
            "crp": rng.normal(2.0, 0.5, n),
            "lymphocyte_pct": rng.normal(30.0, 4.0, n),
            "mcv": rng.normal(89.0, 3.0, n),
            "rdw": rng.normal(13.2, 0.5, n),
            "alp": rng.normal(70.0, 10.0, n),
            "wbc": rng.normal(7.0, 1.0, n),
            "age": np.linspace(25, 80, n),
        }
    )


def test_phenoage_accepts_geroquery_column_names():
    """The clinical frame names columns after the marker; the coefficients are
    named after the unit. Neither side should have to rename to talk."""
    frame = _phenoage_frame()
    renamed = phenoage.normalize_columns(frame)
    assert set(phenoage.REQUIRED_FEATURES) <= set(renamed.columns)
    assert np.isfinite(phenoage.phenotypic_age(frame)).all()


def test_phenoage_rises_with_chronological_age():
    predicted = phenoage.phenotypic_age(_phenoage_frame(12))
    ages = np.linspace(25, 80, 12)
    assert np.corrcoef(predicted, ages)[0, 1] > 0.9


def test_phenoage_yields_a_mortality_risk_from_the_same_model():
    """Phenotypic age and 10-year mortality risk are two readouts of one fit.

    Checked by *rank* correlation, not Pearson: the risk-to-age mapping is a
    double logarithm, so the relationship is exactly monotone and markedly
    non-linear. A Pearson threshold here would be measuring the curvature, not
    the property being asserted.
    """
    from scipy import stats

    frame = _phenoage_frame(8)
    risk = phenoage.mortality_risk_10yr(frame)
    assert ((risk >= 0) & (risk <= 1)).all()
    rho, _ = stats.spearmanr(risk, phenoage.phenotypic_age(frame))
    assert rho == pytest.approx(1.0)


def test_phenoage_is_registered_without_any_optional_dependency():
    clock = get_registry().get("phenoage")
    assert isinstance(clock, PhenoAgeClock)
    assert clock.info.predicted_outcome == "chronological_age"
    assert clock.info.units == "years"
    assert len(PHENOAGE_INFO.required_features) == 10


def test_phenoage_names_the_features_it_is_missing():
    clock = get_registry().get("phenoage")
    with pytest.raises(ClockInputError) as excinfo:
        clock.predict(_phenoage_frame().drop(columns=["alp", "wbc"]))
    detail = excinfo.value.detail
    assert set(detail["missing_features"]) == {"alp_ul", "wbc_1000ul"}


def test_phenoage_refuses_missing_values_rather_than_imputing():
    frame = _phenoage_frame()
    frame.loc[0, "albumin"] = np.nan
    with pytest.raises(ClockInputError, match="missing"):
        get_registry().get("phenoage").predict(frame)
