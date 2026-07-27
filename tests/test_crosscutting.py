"""Cross-cutting guarantees: reproducibility, provenance, licence enforcement."""

import pytest

from geroquery.exceptions import LicenseViolationError
from geroquery.sources import all_adapters
from geroquery.sources.federated import UK_BIOBANK
from geroquery.sources.local_fixture import LocalEvidenceSource


def test_reproducibility_same_query_same_result(service):
    a = service.gene_report("CDKN2A", species="human")
    b = service.gene_report("CDKN2A", species="human")
    assert a == b


def test_every_evidence_line_carries_a_citation(service):
    report = service.gene_report("GDF15")
    assert report["knowledge"]["evidence"]
    refmap = {r["key"] for r in report["references"]}
    for e in report["knowledge"]["evidence"]:
        assert e["reference_keys"]
        assert all(k in refmap for k in e["reference_keys"])


def test_references_are_real_and_linkable(service):
    for r in service.references():
        assert r["url"].startswith("https://pubmed")
        assert r["year"] > 1990


def test_controlled_source_is_federate_only_and_uncacheable():
    cap = UK_BIOBANK.capabilities()
    assert cap.federated and not cap.cacheable
    assert UK_BIOBANK.license().redistributable is False
    with pytest.raises(LicenseViolationError):
        UK_BIOBANK.assert_cacheable()  # must never be cached/re-hosted


def test_redistributable_source_passes_cache_gate():
    LocalEvidenceSource().assert_cacheable()  # should not raise


def test_all_controlled_adapters_refuse_caching():
    for a in all_adapters():
        if not a.license().redistributable:
            with pytest.raises(LicenseViolationError):
                a.assert_cacheable()
