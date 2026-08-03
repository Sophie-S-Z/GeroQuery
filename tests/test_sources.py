"""M2 sources — adapter contract and registry."""

from geroquery.sources import (
    CuratedKnowledgeSource,
    InterventionSource,
    LocalSignatureSource,
    all_adapters,
)


def test_local_signature_adapter_reads_the_committed_slice():
    """``prefer_full=False`` on purpose: a developer who has run `make data` has
    signatures_full.csv on disk, and without this the assertion would depend on
    their download cache rather than on the repository."""
    source = LocalSignatureSource(prefer_full=False)
    assert source.mode() == "curated_slice"
    sigs = source.signatures()
    assert len(sigs) == 40585
    assert {s.species for s in sigs} == {"human", "mouse"}
    assert all(s.effect_size is not None and s.standard_error is not None for s in sigs)
    assert all(s.source == "GEO" for s in sigs)


ONE_ROW = (
    "gene_id,study_id,omic_layer,species,effect_size,direction\nG,S,transcriptome,human,1,up\n"
)


def test_signature_adapter_prefers_the_full_panel_when_it_has_been_built(tmp_path):
    (tmp_path / "signatures_curated.csv").write_text(ONE_ROW)
    assert LocalSignatureSource(tmp_path).mode() == "curated_slice"
    (tmp_path / "signatures_full.csv").write_text(ONE_ROW)
    assert LocalSignatureSource(tmp_path).mode() == "full"


def test_curated_and_intervention_adapters():
    flags = CuratedKnowledgeSource().flags()
    assert len(flags) == 2720
    assert {f.database for f in flags} == {
        "GenAge",
        "GenAge (models)",
        "CellAge",
        "LongevityMap",
        "GenDR",
    }
    ivs = InterventionSource().interventions()
    names = {i.name.lower() for i in ivs}
    assert {"rapamycin", "metformin", "dietary restriction"} <= names


def test_registry_exposes_capabilities_and_license():
    for a in all_adapters():
        cap, lic = a.capabilities(), a.license()
        assert cap.source_name
        assert isinstance(cap.federated, bool)
        assert isinstance(lic.redistributable, bool)


def test_registry_has_cached_and_federated_adapters():
    adapters = all_adapters()
    assert any(not a.capabilities().federated for a in adapters)  # cached tier
    assert any(a.capabilities().federated for a in adapters)  # federated tier
