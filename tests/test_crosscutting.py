"""Cross-cutting guarantees: reproducibility, provenance, licence enforcement."""

import pytest

from geroquery.exceptions import LicenseViolationError
from geroquery.sources import all_adapters
from geroquery.sources.federated import UK_BIOBANK
from geroquery.sources.local_fixture import LocalSignatureSource


def test_reproducibility_same_query_same_result(service):
    a = service.gene_signature("CDKN2A", species="human", omic_layer="transcriptome")
    b = service.gene_signature("CDKN2A", species="human", omic_layer="transcriptome")
    assert a == b


def test_every_signature_carries_provenance(service):
    body = service.gene_signature("GDF15")
    assert body["signatures"]
    for s in body["signatures"]:
        assert s["study_id"] and s["source"]


def test_every_study_is_versioned(service):
    assert all(s["version"] for s in service.studies())


def test_controlled_source_is_federate_only_and_uncacheable():
    cap = UK_BIOBANK.capabilities()
    assert cap.federated and not cap.cacheable
    assert UK_BIOBANK.license().redistributable is False
    with pytest.raises(LicenseViolationError):
        UK_BIOBANK.assert_cacheable()  # must never be cached/re-hosted


def test_redistributable_source_passes_cache_gate():
    LocalSignatureSource().assert_cacheable()  # should not raise


def test_all_controlled_adapters_refuse_caching():
    for a in all_adapters():
        if not a.license().redistributable:
            with pytest.raises(LicenseViolationError):
                a.assert_cacheable()
