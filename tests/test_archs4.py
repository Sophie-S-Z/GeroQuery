"""ARCHS4 selection rule and count handling, offline against planted metadata."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from geroquery.sources import archs4


def _meta(rows):
    """rows: (gsm, series, characteristics, source[, single-cell probability])."""
    return pd.DataFrame(
        [(*r, 0.1) if len(r) == 4 else r for r in rows],
        columns=[
            "geo_accession",
            "series_id",
            "characteristics_ch1",
            "source_name_ch1",
            "singlecellprobability",
        ],
    )


def _arm(series, age, n, start, extra="", source="liver"):
    return [
        (f"GSM{start + i}", series, f"strain: C57BL/6,tissue: liver,age: {age}{extra}", source)
        for i in range(n)
    ]


def test_characteristics_and_age_units():
    fields = archs4.characteristics("Age (months): 24,Tissue: Liver,sex: M")
    assert fields == {"age (months)": "24", "tissue": "Liver", "sex": "M"}
    assert archs4.sample_age(fields, "mouse") == "24 months"  # unit carried from the key
    # A bare number is years for human, unreadable for mouse.
    assert archs4.sample_age({"age": "67"}, "human") == "67 years"
    assert archs4.sample_age({"age": "24"}, "mouse") == "24"
    assert archs4.sample_age({"age": "8wks"}, "mouse") == "8wks"


def test_a_clean_young_vs_old_series_becomes_one_contrast():
    meta = _meta(_arm("GSE1", "3 months", 3, 0) + _arm("GSE1", "24 months", 4, 10))
    (contrast,) = archs4.select_contrasts(meta, "mouse")
    assert contrast.study_id == "ARCHS4:GSE1:tissue_liver"
    assert contrast.tissue == "liver"
    assert contrast.young_samples == ("GSM0", "GSM1", "GSM2")
    assert len(contrast.old_samples) == 4
    assert contrast.age_range == "3mo vs 24mo"
    assert contrast.sex == "unspecified"


def test_only_control_levels_of_other_attributes_are_kept():
    meta = _meta(
        _arm("GSE1", "3 months", 3, 0, ",treatment: vehicle")
        + _arm("GSE1", "24 months", 3, 10, ",treatment: vehicle")
        + _arm("GSE1", "3 months", 3, 20, ",treatment: LPS")
        + _arm("GSE1", "24 months", 3, 30, ",treatment: LPS")
        # Not in any list, so it must exclude rather than pass as harmless.
        + _arm("GSE2", "3 months", 3, 40, ",perturbation: treadmill")
        + _arm("GSE2", "24 months", 3, 50, ",perturbation: treadmill")
    )
    (contrast,) = archs4.select_contrasts(meta, "mouse")
    assert contrast.stratum == "tissue=liver / treatment=vehicle"
    assert set(contrast.young_samples) == {"GSM0", "GSM1", "GSM2"}


def test_single_cell_small_arms_excluded_series_and_cancer_series_are_dropped():
    young = _arm("GSE1", "3 months", 3, 0)
    single_cell = [(*r, 0.9) for r in _arm("GSE1", "24 months", 3, 10)]
    assert archs4.select_contrasts(_meta(young + single_cell), "mouse") == []

    too_small = _meta(_arm("GSE1", "3 months", 3, 0) + _arm("GSE1", "24 months", 2, 10))
    assert archs4.select_contrasts(too_small, "mouse") == []

    clean = _meta(_arm("GSE1", "3 months", 3, 0) + _arm("GSE1", "24 months", 3, 10))
    assert archs4.select_contrasts(clean, "mouse", exclude_series={"GSE1"}) == []

    # Tumour-adjacent "normal" tissue is not healthy tissue: one tumour sample
    # removes the whole series.
    tumour = [("GSM99", "GSE1", "tissue: liver,age: 24 months", "hepatocellular carcinoma")]
    with_tumour = _meta(young + _arm("GSE1", "24 months", 3, 10) + tumour)
    assert archs4.select_contrasts(with_tumour, "mouse") == []


def test_source_name_splits_only_when_both_arms_share_it():
    cells = [("SST", 0), ("PV", 100)]
    rows = []
    for source, base in cells:
        rows += _arm("GSE1", "2 months", 3, base, source=source)
        rows += _arm("GSE1", "22 months", 3, base + 50, source=source)
    contrasts = archs4.select_contrasts(_meta(rows), "mouse")
    assert [c.stratum for c in contrasts] == [
        "tissue=liver / source=pv",
        "tissue=liver / source=sst",
    ]

    # A source name that encodes the age itself is not a stratum.
    by_age = _meta(
        _arm("GSE1", "3 months", 3, 0, source="young liver")
        + _arm("GSE1", "24 months", 3, 10, source="old liver")
    )
    (contrast,) = archs4.select_contrasts(by_age, "mouse")
    assert contrast.stratum == "tissue=liver"


def test_strata_whose_slugs_collide_keep_distinct_study_ids():
    rows = []
    for i, cell in enumerate(("il27ra+ hsc", "il27ra- hsc")):
        for age, base in (("3 months", 0), ("24 months", 10)):
            rows += [
                (f"GSM{i}{base + j}", "GSE1", f"cell type: {cell},age: {age}", "bone marrow")
                for j in range(3)
            ]
    ids = [c.study_id for c in archs4.select_contrasts(_meta(rows), "mouse")]
    assert len(ids) == len(set(ids)) == 2


def test_log_cpm_and_expression_filter():
    counts = pd.DataFrame({"a": [0, 10, 990_000], "b": [0, 0, 1_000_000]}, index=["x", "y", "z"])
    np.testing.assert_allclose(archs4.log_cpm(counts).loc["z", "b"], np.log2(1e6 + 1))
    assert list(archs4.expressed(counts, 1)) == [False, True, True]
    assert list(archs4.expressed(counts, 2)) == [False, False, True]


def test_the_committed_selection_matches_the_shipped_studies():
    """Every ARCHS4 study shipped must be listed, with its samples, in the audit file."""
    from geroquery.config import SOURCES_DATA

    selection_path = SOURCES_DATA / "archs4_selection.csv"
    if not selection_path.exists():
        pytest.skip("no ARCHS4 selection committed")
    selection = pd.read_csv(selection_path)
    studies = pd.read_csv(SOURCES_DATA / "studies.csv")
    shipped = set(studies.loc[studies["source"] == "ARCHS4", "study_id"])
    assert shipped and shipped <= set(selection["study_id"])
    assert selection["study_id"].is_unique
    sizes = sum(selection[arm].str.split().str.len() for arm in ("young_samples", "old_samples"))
    merged = studies.merge(selection.assign(n=sizes), on="study_id")
    assert (merged["sample_size"] == merged["n"]).all()


class _FakeRelease:
    """Just enough of an h5py.File for the ETL: path -> array."""

    def __init__(self, rows, genes, symbols, counts):
        meta = _meta(rows)
        self._data = {f"meta/samples/{c}": meta[c].to_numpy() for c in meta.columns}
        self._data["meta/genes/ensembl_gene"] = np.array([g.encode() for g in genes])
        self._data["meta/genes/symbol"] = np.array([s.encode() for s in symbols])
        self._data["data/expression"] = np.asarray(counts)

    def __getitem__(self, path):
        return self._data[path]


def test_archs4_frames_estimates_each_contrast_and_records_its_samples(monkeypatch, tmp_path):
    from geroquery import config
    from geroquery.etl import build_signatures
    from geroquery.exceptions import NetworkDisabledError

    monkeypatch.setattr(config, "CACHE_HOME", tmp_path)

    rows = _arm("GSE1", "3 months", 3, 0) + _arm("GSE1", "24 months", 3, 10)
    # Two Ensembl ids share a symbol (collapsed by mean expression, not by
    # significance); XYZ rises with age; NIL is never expressed and is dropped.
    counts = [
        [500, 520, 480, 510, 490, 505],
        [50, 60, 40, 55, 45, 52],
        [100, 110, 90, 400, 420, 390],
        [0, 0, 0, 0, 0, 0],
    ]
    genes, symbols = ["E1", "E2", "E3", "E4"], ["ABC", "ABC", "XYZ", "NIL"]
    release = _FakeRelease(rows, genes, symbols, counts)
    empty = _FakeRelease([], genes, symbols, np.zeros((4, 0)))
    monkeypatch.setattr(
        archs4, "open_release", lambda species: release if species == "mouse" else empty
    )

    signatures, studies, selection = build_signatures.archs4_frames(
        exclude_series=set(), allow_network=True, quiet=True, workers=1
    )
    assert sorted(signatures["gene_key"]) == ["SYM:ABC", "SYM:XYZ"]
    xyz = signatures.set_index("gene_key").loc["SYM:XYZ"]
    assert xyz["effect_size"] > 2 and xyz["source"] == "ARCHS4"
    (study,) = studies.to_dict("records")
    assert study["sample_size"] == 6 and study["series_id"] == "GSE1"
    assert selection["young_samples"].tolist() == ["GSM0 GSM1 GSM2"]

    # The second run is served from the per-contrast cache, identically.
    assert len(list((tmp_path / "archs4").glob("*.parquet"))) == 1
    again, _, _ = build_signatures.archs4_frames(
        exclude_series=set(), allow_network=True, quiet=True, workers=1
    )
    pd.testing.assert_frame_equal(again, signatures)

    with pytest.raises(NetworkDisabledError):
        build_signatures.archs4_frames(exclude_series=set(), allow_network=False)


def test_a_series_with_both_species_yields_two_distinct_study_ids(monkeypatch, tmp_path):
    from geroquery import config
    from geroquery.etl import build_signatures

    monkeypatch.setattr(config, "CACHE_HOME", tmp_path)
    counts = [[500, 520, 480, 510, 490, 505], [100, 110, 90, 400, 420, 390]]
    genes, symbols = ["E1", "E2"], ["ABC", "XYZ"]
    mouse = _FakeRelease(
        _arm("GSE1", "3 months", 3, 0) + _arm("GSE1", "24 months", 3, 10), genes, symbols, counts
    )
    human = _FakeRelease(
        _arm("GSE1", "25 years", 3, 20) + _arm("GSE1", "70 years", 3, 30), genes, symbols, counts
    )
    monkeypatch.setattr(archs4, "open_release", lambda s: human if s == "human" else mouse)

    _, studies, selection = build_signatures.archs4_frames(
        exclude_series=set(), allow_network=True, quiet=True, workers=1
    )
    assert studies["study_id"].is_unique and selection["study_id"].is_unique
    assert len(studies) == 2
