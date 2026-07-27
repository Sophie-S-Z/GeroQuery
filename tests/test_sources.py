"""M2 sources — adapter contract and registry."""

from geroquery.sources import (
    CuratedKnowledgeSource,
    InterventionSource,
    LocalEvidenceSource,
    all_adapters,
)


def test_curated_and_intervention_adapters():
    assert len(CuratedKnowledgeSource().flags()) == 50
    ivs = InterventionSource().interventions()
    assert {"rapamycin", "metformin"} <= {i.name for i in ivs}
    # Interventions carry a real primary-citation URL.
    rapa = next(i for i in ivs if i.name == "rapamycin")
    assert rapa.url.startswith("https://pubmed")


def test_local_evidence_source_is_cacheable():
    LocalEvidenceSource().assert_cacheable()  # redistributable, should not raise


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
