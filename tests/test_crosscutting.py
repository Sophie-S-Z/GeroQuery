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


def test_references_are_real_and_linkable(service):
    """From the citation layer carried over in the PR #1 merge. A reference that
    does not resolve is worse than no reference: it looks like sourcing."""
    references = service.references()
    assert references
    for reference in references:
        assert reference["url"].startswith("https://pubmed")
        assert reference["year"] > 1990
        assert reference["citation"]


def test_curated_genes_are_browsable_without_a_network(service):
    """The bundled identifier table used to name 24 genes while the curated set
    named 1,880, so browsing offline silently returned a tiny subset."""
    genes = service.list_curated_genes()
    assert len(genes) > 1000
    assert all(g["symbol"] and g["databases"] for g in genes)


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
