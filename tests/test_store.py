"""M4 store — curated flags, interventions, datasets, reproducible versioning."""

import pytest

from geroquery.store import DatasetNotFoundError, GeroStore


def test_curated_flags(store):
    flags = store.curated_flags("ENSG00000147889")  # CDKN2A human
    assert {f.database for f in flags} >= {"GenAge", "CellAge"}
    assert all(f.url for f in flags)


def test_interventions(store):
    rapa = store.interventions(name="rapamycin")
    assert rapa and rapa[0].lifespan_effect_pct == 14.0
    assert rapa[0].source == "NIA ITP"
    assert rapa[0].url.startswith("https://pubmed")


def test_interventions_linked_to_gene(store):
    # MTOR (human) should be linked to rapamycin.
    linked = store.interventions(gene_id="ENSG00000198793")
    assert any(iv.name == "rapamycin" for iv in linked)


def test_example_cohort_dataset(store):
    ds = {d["dataset_id"] for d in store.list_datasets()}
    assert "example_cohort_simulated" in ds
    df = store.get_dataset("example_cohort_simulated")
    assert df.shape[0] == 720
    # Carries the nine PhenoAge markers.
    assert {
        "albumin_gdl",
        "creatinine_mgdl",
        "glucose_mgdl",
        "crp_mgl",
        "lymphocyte_pct",
        "mcv_fl",
        "rdw_pct",
        "alp_ul",
        "wbc_1000ul",
        "age",
    } <= set(df.columns)


def test_missing_dataset_raises(store):
    with pytest.raises(DatasetNotFoundError):
        store.get_dataset("no_such_dataset")


def test_no_fabricated_signature_api(store):
    # The fabricated per-study signature query is gone for good.
    assert not hasattr(store, "query_signatures")
    assert not hasattr(store, "list_studies")


def test_version_is_reproducible(tmp_path):
    """Two independent builds from the same fixtures produce the same version
    and identical query results."""
    s1 = GeroStore(data_home=tmp_path / "a").build()
    s2 = GeroStore(data_home=tmp_path / "b").build()
    assert s1.version() == s2.version()
    f1 = s1.curated_flags("ENSG00000113368")
    f2 = s2.curated_flags("ENSG00000113368")
    assert [x.to_dict() for x in f1] == [x.to_dict() for x in f2]
