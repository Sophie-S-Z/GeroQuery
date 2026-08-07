"""The signature ETL and the bulk identifier resolver.

Offline throughout: the GEO fetch is redirected at a gzipped fixture and
mygene.info is stubbed. What is under test is the ETL's *rules* — which probes
survive resolution, which get an ``ENTREZ:`` fallback, how duplicate gene ids
collapse, and what the committed slice is allowed to contain.
"""

from __future__ import annotations

import gzip
import json
import math

import pandas as pd
import pytest

from geroquery.etl import build_signatures
from geroquery.exceptions import NetworkDisabledError
from geroquery.idmap import bulk
from geroquery.sources.manifest import GEO_AGING_PANEL

SOFT = """^DATASET = GDS707
!dataset_title = Aging brain: frontal cortex expression profiles at various ages
!dataset_platform = GPL8300
!dataset_sample_organism = Homo sapiens
!dataset_value_type = count
!dataset_sample_count = 6
!dataset_reference_series = GSE1572
!dataset_pubmed_id = 15190254
!dataset_update_date = Oct 10 2012
^SUBSET = GDS707_1
!subset_dataset_id = GDS707
!subset_description = 20 - 39 y
!subset_sample_id = GSM1,GSM2,GSM3
!subset_type = age
^SUBSET = GDS707_2
!subset_dataset_id = GDS707
!subset_description = 60 - 79 y
!subset_sample_id = GSM4,GSM5,GSM6
!subset_type = age
^DATASET = GDS707
#GSM1 = Value for GSM1: 26 year old male; src: human frontal cortex
!dataset_table_begin
ID_REF\tIDENTIFIER\tGSM1\tGSM2\tGSM3\tGSM4\tGSM5\tGSM6\tGene ID
p1\tAAA\t100\t110\t105\t400\t410\t405\t111
p2\tBBB\t50\t51\t52\t53\t54\t55\t222
p3\tCCC\t10\t12\t11\t9\t8\t10\t333
!dataset_table_end
"""


@pytest.fixture
def soft_path(tmp_path):
    path = tmp_path / "GDS707__GDS707_full.soft.gz"
    path.write_bytes(gzip.compress(SOFT.encode("utf-8")))
    return path


@pytest.fixture
def offline_geo(monkeypatch, soft_path):
    monkeypatch.setattr(build_signatures, "fetch_artifact", lambda *a, **k: soft_path)
    return soft_path


# ---- contrast estimation --------------------------------------------------


def test_contrast_frames_emits_one_signature_row_per_gene(offline_geo):
    signatures, studies, skipped = build_signatures.contrast_frames(
        accessions=["GDS707"], quiet=True
    )
    assert skipped == []
    assert sorted(signatures["gene_key"]) == ["111", "222", "333"]
    assert set(signatures["study_id"]) == {"GEO:GDS707"}
    assert signatures["species"].unique().tolist() == ["human"]
    assert signatures["tissue"].unique().tolist() == ["brain"]
    # p1 quadruples in the old group; p3 barely moves.
    by_gene = signatures.set_index("gene_key")
    assert by_gene.loc["111", "direction"] == "up"
    assert by_gene.loc["111", "effect_size"] > by_gene.loc["333", "effect_size"]


def test_contrast_frames_records_study_provenance(offline_geo):
    _signatures, studies, _skipped = build_signatures.contrast_frames(
        accessions=["GDS707"], quiet=True
    )
    (study,) = studies.to_dict("records")
    assert study["series_id"] == "GSE1572"
    assert study["pubmed_id"] == "15190254"
    assert study["sample_size"] == 6
    assert study["n_young"] == 3 and study["n_old"] == 3
    assert study["source"] == "GEO"
    assert "Hedges' g" in study["processing_method"]


# ---- identifier resolution ------------------------------------------------


def _raw(gene_keys, species="human"):
    return pd.DataFrame(
        {
            "gene_key": gene_keys,
            "species": species,
            "study_id": "GEO:GDSX",
            "omic_layer": "transcriptome",
            "tissue": "brain",
            "sex": "both",
            "age_range": "30y vs 70y",
            "effect_size": [0.5] * len(gene_keys),
            "direction": "up",
            "p_value": [0.01, 0.5, 0.9][: len(gene_keys)],
            "q_value": 0.5,
            "standard_error": 0.2,
            "source": "GEO",
        }
    )


def test_unresolved_entrez_keys_keep_an_entrez_id_but_unknown_symbols_are_dropped(monkeypatch):
    """An ``ENTREZ:`` id is a real, stable identifier for a gene mygene could not
    place in Ensembl. An unrecognised symbol is not an identifier at all, and a
    row keyed on one could never join to curated knowledge."""
    monkeypatch.setattr(
        build_signatures, "resolve_bulk", lambda queries, *a, **k: dict.fromkeys(queries)
    )
    out, stats, symbols = build_signatures.resolve_gene_ids(
        _raw(["111", "SYM:MYSTERY"]), quiet=True
    )
    assert out["gene_id"].tolist() == ["ENTREZ:111"]
    assert stats["rows_dropped_unresolved"] == 1
    assert symbols == {}


def test_two_probe_keys_for_one_gene_collapse_to_the_most_significant_row(monkeypatch):
    """Deprecated Entrez ids can both now point at one gene. Emitting the gene
    twice per study with different effect sizes would double-weight it in the
    pool."""
    record = {"canonical_id": "ENSG00000000001", "symbol": "AAA"}
    monkeypatch.setattr(
        build_signatures, "resolve_bulk", lambda queries, *a, **k: dict.fromkeys(queries, record)
    )
    out, stats, symbols = build_signatures.resolve_gene_ids(_raw(["111", "222"]), quiet=True)
    assert out["gene_id"].tolist() == ["ENSG00000000001"]
    assert out["p_value"].tolist() == [0.01]  # the more significant of the two
    assert stats["rows_dropped_duplicate_gene"] == 1
    assert symbols == {"ENSG00000000001": "AAA"}


def test_resolved_rows_are_rounded_and_ordered(monkeypatch):
    def _resolve(queries, species, **kwargs):
        return {q: {"canonical_id": f"ENSG{q}", "symbol": f"S{q}"} for q in queries}

    monkeypatch.setattr(build_signatures, "resolve_bulk", _resolve)
    out, _stats, _symbols = build_signatures.resolve_gene_ids(_raw(["222", "111"]), quiet=True)
    assert out["gene_id"].tolist() == ["ENSG111", "ENSG222"]
    assert list(out.columns) == build_signatures.SIGNATURE_COLUMNS


# ---- bulk identifier resolver ---------------------------------------------


def test_resolve_bulk_refuses_to_silently_drop_genes_when_offline(tmp_path):
    """Unlike the request-path resolver, a missing id here is fatal: a build that
    quietly dropped every unresolved gene would write a plausible-looking but
    incomplete signature table."""
    with pytest.raises(NetworkDisabledError, match="not in the local id map"):
        bulk.resolve_bulk(["1029"], "human", cache_dir=tmp_path, allow_network=False)


def test_resolve_bulk_caches_negatives_so_they_are_not_re_queried(tmp_path, monkeypatch):
    calls = []

    def _query(queries, species, timeout):
        calls.append(list(queries))
        return [{"query": "1029", "symbol": "CDKN2A", "taxid": 9606, "entrezgene": "1029"}]

    monkeypatch.setattr(bulk, "_query", _query)
    first = bulk.resolve_bulk(["1029", "999999"], "human", cache_dir=tmp_path, allow_network=True)
    assert first["1029"]["symbol"] == "CDKN2A"
    assert first["999999"] is None

    # Second call, still allowing network: nothing new to ask about.
    bulk.resolve_bulk(["1029", "999999"], "human", cache_dir=tmp_path, allow_network=True)
    assert len(calls) == 1

    stored = json.loads(bulk.map_path("human", tmp_path).read_text(encoding="utf-8"))
    assert stored["999999"] is None


def test_resolve_bulk_discards_records_from_the_wrong_species(tmp_path, monkeypatch):
    """Defence in depth behind mygene's own species filter."""
    monkeypatch.setattr(
        bulk,
        "_query",
        lambda queries, species, timeout: [
            {"query": "1029", "symbol": "Cdkn2a", "taxid": 10090, "entrezgene": "12578"}
        ],
    )
    out = bulk.resolve_bulk(["1029"], "human", cache_dir=tmp_path, allow_network=True)
    assert out["1029"] is None


def test_resolve_bulk_deduplicates_and_preserves_the_requested_order(tmp_path, monkeypatch):
    monkeypatch.setattr(
        bulk,
        "_query",
        lambda queries, species, timeout: [
            {"query": q, "symbol": f"S{q}", "taxid": 9606, "entrezgene": q} for q in queries
        ],
    )
    out = bulk.resolve_bulk(
        ["7157", "1029", "7157"], "human", cache_dir=tmp_path, allow_network=True
    )
    assert list(out) == ["7157", "1029"]


# ---- the whole ETL --------------------------------------------------------


def test_build_all_writes_every_table_and_a_slice_that_is_not_the_full_panel(
    tmp_path, monkeypatch, offline_geo
):
    monkeypatch.setattr(build_signatures, "GEO_AGING_PANEL", {"GDS707": GEO_AGING_PANEL["GDS707"]})

    def _resolve(queries, species, **kwargs):
        return {q: {"canonical_id": f"ENSG{q}", "symbol": f"S{q}"} for q in queries}

    monkeypatch.setattr(build_signatures, "resolve_bulk", _resolve)
    monkeypatch.setattr(
        build_signatures,
        "build_curated",
        lambda **kwargs: (
            pd.DataFrame(
                [{"gene_id": "ENSG111", "database": "GenAge", "assertion": "a", "url": "u"}]
            ),
            pd.DataFrame([{"intervention_id": "IV1", "name": "x"}]),
            {"assertions_loaded": 1},
        ),
    )

    summary = build_signatures.build_all(tmp_path, allow_network=False, quiet=True)
    assert summary["contrasts"] == 1
    assert summary["signature_rows"] == 3

    full = pd.read_csv(tmp_path / build_signatures.FULL_SIGNATURES)
    slice_ = pd.read_csv(tmp_path / build_signatures.CURATED_SIGNATURES)
    assert len(full) == 3
    assert slice_["gene_id"].tolist() == ["ENSG111"]  # only the curated gene
    for name in (
        build_signatures.STUDIES,
        build_signatures.CURATED_KNOWLEDGE,
        build_signatures.INTERVENTIONS,
    ):
        assert (tmp_path / name).exists()


def test_build_all_refuses_to_write_empty_tables(tmp_path, monkeypatch):
    monkeypatch.setattr(build_signatures, "GEO_AGING_PANEL", {})
    with pytest.raises(RuntimeError, match="refusing to write empty tables"):
        build_signatures.build_all(tmp_path, allow_network=False, quiet=True)


def test_the_exported_verdict_agrees_with_the_exported_interval():
    """A reader recomputing the verdict from the numbers on the page must get
    the printed one.

    The verdict used to be derived at full precision while the interval shipped
    rounded to four decimals, so 78 genes carried an interval of
    ``[-1.195, -0.000]`` labelled "decreases". Nothing errored and nothing looked
    wrong; the page simply asserted something its own numbers did not support.

    The negative-zero half is asserted over ``BOUND_COLUMNS`` rather than over
    two hand-named columns. The first version of this test named ``ci_low`` and
    ``ci_high`` literally and so did not notice that the prediction interval
    added in the same commit had skipped the fix — two genes shipped
    ``pi_high = -0.0``. A test that enumerates its own columns only ever covers
    the bug that was already known.
    """
    import pandas as pd

    from geroquery.etl.build_frontend_data import BOUND_COLUMNS, _pool

    # Effects engineered to land an interval bound within rounding distance of
    # zero, which is the only regime where the two definitions can differ.
    rows = []
    for i, offset in enumerate((0.0, 1e-6, -1e-6, 1e-5)):
        for study in range(4):
            rows.append(
                {
                    "gene_id": f"G{i}",
                    "species": "human",
                    "effect_size": -0.6 + offset + 0.0001 * study,
                    "standard_error": 0.3061,
                    "tissue": "blood",
                }
            )
    # And one gene whose *prediction* bound lands on zero instead. Four equal
    # effects give tau2 = 0 and se_pool = se / 2, so pi_high is
    # eff + t(0.975, k-2) * se / 2 exactly; solving that for a hair below zero
    # is what puts a -0.0 in pi_high. Derived from the quantile rather than
    # hard-coded so it stays on target if the interval's dof convention moves.
    from scipy import stats as _stats

    se_pred = 0.3
    eff_pred = -float(_stats.t.ppf(0.975, 2)) * (se_pred / 2) - 1e-7
    for _study in range(4):
        rows.append(
            {
                "gene_id": "PI",
                "species": "human",
                "effect_size": eff_pred,
                "standard_error": se_pred,
                "tissue": "blood",
            }
        )

    pooled = _pool(pd.DataFrame(rows))
    assert (pooled["gene_id"] == "PI").any(), "the prediction-bound fixture did not pool"

    for _, row in pooled.iterrows():
        expected = (
            "increases"
            if row["ci_low"] > 0
            else "decreases" if row["ci_high"] < 0 else "no_evidence"
        )
        assert row["verdict"] == expected, row.to_dict()
        # Negative zero compares as "not less than zero" in JavaScript but
        # prints as "-0.000", so it must not survive the export — on *any*
        # bound, not just the two the reported verdict is read from.
        for column in BOUND_COLUMNS:
            value = row[column]
            if value is None:
                continue
            assert not (
                value == 0 and math.copysign(1, value) < 0
            ), f"{column} shipped as negative zero for {row['gene_id']}"
