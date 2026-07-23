"""M4 store — round-trip, filtered query, versioning, datasets."""

import pytest

from geroquery.store import DatasetNotFoundError, GeroStore


def test_query_signatures_filters(store):
    sigs = store.query_signatures(
        gene_id="ENSG00000147889", species="human", omic_layer="transcriptome"
    )
    assert len(sigs) == 9  # 3 tissues x 3 studies
    assert all(s.species == "human" and s.omic_layer == "transcriptome" for s in sigs)
    assert all(s.direction == "up" for s in sigs)  # CDKN2A up with age


def test_partition_pruning_by_tissue(store):
    blood = store.query_signatures(gene_id="ENSG00000147889", tissue="blood")
    assert blood and all(s.tissue == "blood" for s in blood)


def test_studies_and_provenance(store):
    studies = store.list_studies()
    assert len(studies) == 228
    one = store.get_study(studies[0].study_id)
    assert one is not None and one.version == "2024.1"


def test_curated_and_interventions(store):
    flags = store.curated_flags("ENSG00000147889")
    assert {f.database for f in flags} >= {"GenAge", "CellAge"}
    rapa = store.interventions(name="rapamycin")
    assert rapa and rapa[0].lifespan_effect_pct == 14.0


def test_dataset_load_and_missing(store):
    df = store.get_dataset("clinical_nhanes_slice")
    assert df.shape[0] == 720
    with pytest.raises(DatasetNotFoundError):
        store.get_dataset("no_such_dataset")


def test_version_is_reproducible(tmp_path):
    """Two independent builds from the same fixtures produce the same version
    and identical query results."""
    s1 = GeroStore(data_home=tmp_path / "a").build()
    s2 = GeroStore(data_home=tmp_path / "b").build()
    assert s1.version() == s2.version()
    q1 = s1.query_signatures(gene_id="ENSG00000113368")
    q2 = s2.query_signatures(gene_id="ENSG00000113368")
    assert [x.to_dict() for x in q1] == [x.to_dict() for x in q2]
