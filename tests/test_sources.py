"""M2 sources — adapter contract and registry."""

from geroquery.sources import (
    CuratedKnowledgeSource,
    InterventionSource,
    LocalSignatureSource,
    all_adapters,
)


def test_local_signature_adapter_reads_slice():
    sigs = LocalSignatureSource().signatures()
    assert len(sigs) == 228
    assert all(s.effect_size is not None and s.standard_error is not None for s in sigs)


def test_curated_and_intervention_adapters():
    assert len(CuratedKnowledgeSource().flags()) == 50
    ivs = InterventionSource().interventions()
    assert {"rapamycin", "metformin"} <= {i.name for i in ivs}


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
