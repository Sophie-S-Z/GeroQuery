"""mygene.info batch resolution and the GTEx Portal v2 adapter.

Offline by default: HTTP is stubbed so behaviour is pinned without touching a
real upstream. The `live` tests actually call the APIs and are excluded from the
default run.
"""

from __future__ import annotations

import json

import pytest

from geroquery.exceptions import SourceError
from geroquery.idmap.mygene import MyGeneClient, normalize_hit
from geroquery.idmap.resolver import GeneResolver
from geroquery.sources.gtex import DEFAULT_DATASET_ID, GtexOpenSource

# ---- mygene: response normalization ---------------------------------------

_CDKN2A_HUMAN = {
    "query": "CDKN2A",
    "symbol": "CDKN2A",
    "name": "cyclin dependent kinase inhibitor 2A",
    "taxid": 9606,
    "entrezgene": "1029",
    "ensembl": {"gene": "ENSG00000147889"},
    "uniprot": {"Swiss-Prot": ["P42771", "Q8N726"]},
    "alias": ["P16", "INK4A"],
}
_CDKN2A_MOUSE = {
    "query": "CDKN2A",
    "symbol": "Cdkn2a",
    "taxid": 10090,
    "entrezgene": "12578",
    "ensembl": {"gene": "ENSMUSG00000044303"},
    "uniprot": {"Swiss-Prot": ["P51480"]},
}
_GAPDH = {
    "query": "ENSG00000111640",
    "symbol": "GAPDH",
    "taxid": 9606,
    "entrezgene": "2597",
    "ensembl": {"gene": "ENSG00000111640"},
    "uniprot": {"Swiss-Prot": "P04406"},  # scalar, not a list
}


def test_swiss_prot_list_is_collapsed_to_one_accession():
    """mygene returns Swiss-Prot as a scalar for some genes and a list for others.
    Storing the list puts an unhashable value into the uniprot index."""
    assert normalize_hit(_CDKN2A_HUMAN)["uniprot"] == "P42771"
    assert normalize_hit(_GAPDH)["uniprot"] == "P04406"


def test_notfound_hits_normalize_to_none():
    assert normalize_hit({"query": "___nope___", "notfound": True}) is None


def test_hit_without_a_usable_id_is_rejected():
    assert normalize_hit({"query": "x", "symbol": "X", "taxid": 9606}) is None


def test_unknown_taxid_is_labelled_unknown_not_guessed():
    hit = dict(_CDKN2A_HUMAN, taxid=7227)  # fly
    assert normalize_hit(hit)["species"] == "unknown"


def test_ensembl_cross_reference_list_does_not_leak_another_species_id():
    """mygene puts homology cross-references in the ``ensembl`` list, so a mouse
    gene's list can lead with a flycatcher or ferret id. Taking element zero put
    those into the canonical id space, where they silently became the join key
    for that mouse gene's signatures."""
    hit = {
        "query": "100017",
        "symbol": "Ldlrap1",
        "taxid": 10090,
        "entrezgene": "100017",
        "ensembl": [
            {"gene": "ENSFALG00000025362"},  # collared flycatcher
            {"gene": "ENSMPUG00000015874"},  # ferret
            {"gene": "ENSMUSG00000037295"},  # the gene's own id
        ],
    }
    record = normalize_hit(hit)
    assert record["canonical_id"] == "ENSMUSG00000037295"
    assert record["species"] == "mouse"


def test_ensembl_list_of_only_cross_references_falls_back_to_entrez():
    """An ``ENTREZ:`` id is at least the right gene; another species' Ensembl id
    is not."""
    hit = {
        "query": "12345",
        "symbol": "Fake1",
        "taxid": 10090,
        "entrezgene": "12345",
        "ensembl": [{"gene": "ENSFALG00000025362"}],
    }
    assert normalize_hit(hit)["canonical_id"] == "ENTREZ:12345"


# ---- mygene: client behaviour ---------------------------------------------


@pytest.fixture
def client(tmp_path):
    return MyGeneClient(cache_dir=tmp_path / "mygene", allow_network=True, timeout=5)


def _stub_post(monkeypatch, hits, counter):
    def _post(self, queries, species):
        counter.append(list(queries))
        return [h for h in hits if h["query"] in set(queries)]

    monkeypatch.setattr(MyGeneClient, "_post", _post)


def test_resolve_batches_every_miss_into_one_request(client, monkeypatch):
    """The whole point of the batch path: N unknown genes must cost one request,
    not N. Looping resolve_gene issues one round trip per gene."""
    calls: list[list[str]] = []
    _stub_post(monkeypatch, [_CDKN2A_HUMAN, _CDKN2A_MOUSE, _GAPDH], calls)

    out = client.resolve(["CDKN2A", "ENSG00000111640"])
    assert len(calls) == 1
    assert sorted(calls[0]) == ["CDKN2A", "ENSG00000111640"]
    assert len(out["CDKN2A"]) == 2
    assert len(out["ENSG00000111640"]) == 1


def test_multi_species_hits_are_ordered_deterministically(client, monkeypatch):
    """Upstream orders by relevance score. Repeated runs must pick the same
    primary regardless, so ordering is imposed here: human first."""
    _stub_post(monkeypatch, [_CDKN2A_MOUSE, _CDKN2A_HUMAN], [])
    records = client.resolve(["CDKN2A"])["CDKN2A"]
    assert [r["species"] for r in records] == ["human", "mouse"]


def test_second_call_is_served_from_cache(client, monkeypatch):
    calls: list[list[str]] = []
    _stub_post(monkeypatch, [_CDKN2A_HUMAN], calls)

    client.resolve(["CDKN2A"])
    client.resolve(["CDKN2A"])
    assert len(calls) == 1  # second call hit the cache


def test_negative_results_are_cached_too(client, monkeypatch):
    """Without this, every lookup of a typo re-hits the API forever."""
    calls: list[list[str]] = []
    _stub_post(monkeypatch, [], calls)

    assert client.resolve(["___nope___"])["___nope___"] == []
    assert client.resolve(["___nope___"])["___nope___"] == []
    assert len(calls) == 1


def test_corrupt_cache_entry_is_treated_as_a_miss(client, monkeypatch):
    calls: list[list[str]] = []
    _stub_post(monkeypatch, [_CDKN2A_HUMAN], calls)
    client.resolve(["CDKN2A"])

    client._cache_path("CDKN2A", None).write_text("{not json", encoding="utf-8")
    assert client.resolve(["CDKN2A"])["CDKN2A"]
    assert len(calls) == 2


def test_offline_client_degrades_to_not_found_rather_than_raising(tmp_path):
    """This is a fallback behind a local table. A disabled network must not turn
    a mostly-resolvable query into an exception."""
    offline = MyGeneClient(cache_dir=tmp_path / "c", allow_network=False)
    assert offline.resolve(["CDKN2A"]) == {"CDKN2A": []}


def test_cache_key_separates_species_scopes(client):
    assert client._cache_path("CDKN2A", None) != client._cache_path("CDKN2A", "human")


# ---- resolver integration --------------------------------------------------


def test_resolve_batch_only_asks_upstream_about_local_misses(tmp_path, monkeypatch):
    """Genes already in the bundled table must never generate a request."""
    asked: list[list[str]] = []

    class _Recorder:
        def resolve(self, queries, species=None):
            asked.append(list(queries))
            return {q: [] for q in queries}

    resolver = GeneResolver(mygene=_Recorder())
    out = resolver.resolve_batch(["CDKN2A", "___nope___"], species="human")

    assert out["CDKN2A"] is not None and out["CDKN2A"].symbol == "CDKN2A"
    assert out["___nope___"] is None
    assert asked == [["___nope___"]]


def test_resolve_batch_survives_an_upstream_failure(tmp_path):
    """A dead mygene must not fail the whole batch — locally resolvable genes
    still resolve, and the rest come back as None."""

    class _Broken:
        def resolve(self, queries, species=None):
            raise SourceError("upstream down")

    resolver = GeneResolver(mygene=_Broken())
    out = resolver.resolve_batch(["CDKN2A", "___nope___"], species="human")
    assert out["CDKN2A"] is not None
    assert out["___nope___"] is None


def test_resolve_batch_preserves_query_order_and_duplicates():
    resolver = GeneResolver()
    queries = ["CDKN2A", "LMNB1", "CDKN2A"]
    assert list(resolver.resolve_batch(queries, species="human")) == ["CDKN2A", "LMNB1"]


# ---- GTEx ------------------------------------------------------------------


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def _stub_gtex(monkeypatch, routes, seen):
    import httpx

    def _get(url, params=None, headers=None, timeout=None, follow_redirects=None):
        seen.append((url, dict(params or {})))
        for fragment, payload in routes.items():
            if fragment in url:
                return _FakeResponse(payload)
        return _FakeResponse({"data": []})

    monkeypatch.setattr(httpx, "get", _get)


_GENE_ROW = {
    "gencodeId": "ENSG00000147889.17",
    "geneSymbol": "CDKN2A",
    "geneSymbolUpper": "CDKN2A",
    "chromosome": "chr9",
}


def test_every_expression_request_carries_a_dataset_id(monkeypatch):
    """Omitting datasetId returns HTTP 200 with an empty list — a silent empty
    result that reads as 'this gene is not expressed anywhere'."""
    seen: list[tuple[str, dict]] = []
    _stub_gtex(
        monkeypatch,
        {
            "/reference/gene": {"data": [_GENE_ROW]},
            "/medianGeneExpression": {
                "data": [
                    {
                        "median": 0.53,
                        "tissueSiteDetailId": "Adipose_Subcutaneous",
                        "ontologyId": "UBERON:0002190",
                        "gencodeId": "ENSG00000147889.17",
                        "geneSymbol": "CDKN2A",
                        "unit": "TPM",
                        "datasetId": DEFAULT_DATASET_ID,
                    }
                ]
            },
        },
        seen,
    )
    GtexOpenSource(allow_network=True).median_expression("CDKN2A")
    expression_calls = [p for url, p in seen if "Expression" in url]
    assert expression_calls and all(p.get("datasetId") for p in expression_calls)


def test_median_expression_is_sorted_and_carries_uberon(monkeypatch):
    _stub_gtex(
        monkeypatch,
        {
            "/reference/gene": {"data": [_GENE_ROW]},
            "/medianGeneExpression": {
                "data": [
                    {
                        "median": 0.5,
                        "tissueSiteDetailId": "Adipose_Subcutaneous",
                        "ontologyId": "UBERON:0002190",
                        "gencodeId": "ENSG00000147889.17",
                        "geneSymbol": "CDKN2A",
                        "unit": "TPM",
                        "datasetId": DEFAULT_DATASET_ID,
                    },
                    {
                        "median": 9.1,
                        "tissueSiteDetailId": "Whole_Blood",
                        "ontologyId": "UBERON:0013756",
                        "gencodeId": "ENSG00000147889.17",
                        "geneSymbol": "CDKN2A",
                        "unit": "TPM",
                        "datasetId": DEFAULT_DATASET_ID,
                    },
                ]
            },
        },
        [],
    )
    out = GtexOpenSource(allow_network=True).median_expression("CDKN2A")
    assert [t.tissue_id for t in out] == ["Whole_Blood", "Adipose_Subcutaneous"]
    assert out[0].uberon_id == "UBERON:0013756"


def test_sample_expression_states_it_is_not_age_stratified(monkeypatch):
    """The open API accepts attributeSubset=ageBracket but does not actually
    split by it. The response must say so rather than let a caller assume."""
    _stub_gtex(
        monkeypatch,
        {
            "/reference/gene": {"data": [_GENE_ROW]},
            "/geneExpression": {
                "data": [
                    {
                        "data": [0.16, 0.24, 0.31],
                        "tissueSiteDetailId": "Whole_Blood",
                        "ontologyId": "UBERON:0013756",
                        "unit": "TPM",
                        "subsetGroup": None,
                    }
                ]
            },
        },
        [],
    )
    out = GtexOpenSource(allow_network=True).sample_expression("CDKN2A", "Whole_Blood")
    assert out["n_samples"] == 3
    assert out["age_stratified"] is False
    assert "dbGaP" in out["note"]


def test_unknown_gene_raises_rather_than_returning_empty(monkeypatch):
    _stub_gtex(monkeypatch, {"/reference/gene": {"data": []}}, [])
    with pytest.raises(SourceError, match="no gene matching"):
        GtexOpenSource(allow_network=True).lookup_gene("___nope___")


def test_gtex_refuses_to_query_when_network_is_disabled():
    with pytest.raises(SourceError, match="network access is disabled"):
        GtexOpenSource(allow_network=False).lookup_gene("CDKN2A")


def test_gtex_open_is_registered_and_declares_no_age_stratification():
    from geroquery.sources import all_adapters

    adapters = {a.name: a for a in all_adapters()}
    assert "gtex-open" in adapters
    notes = adapters["gtex-open"].capabilities().notes
    assert "NOT available" in notes
    # The controlled tier is where donor age lives; it must stay uncacheable.
    protected = adapters["gtex-protected"]
    assert protected.capabilities().cacheable is False


# ---- live ------------------------------------------------------------------


@pytest.mark.live
def test_live_mygene_batch_resolves_mixed_identifier_spaces(tmp_path):
    out = MyGeneClient(cache_dir=tmp_path / "c", allow_network=True).resolve(
        ["CDKN2A", "ENSG00000111640", "7157", "___nope___"]
    )
    assert {r["symbol"] for r in out["CDKN2A"]} == {"CDKN2A", "Cdkn2a"}
    assert out["ENSG00000111640"][0]["symbol"] == "GAPDH"
    assert out["7157"][0]["symbol"] == "TP53"
    assert out["___nope___"] == []


@pytest.mark.live
def test_live_gtex_median_expression():
    out = GtexOpenSource(allow_network=True).median_expression("CDKN2A")
    assert len(out) > 40  # GTEx v8 has ~54 tissues
    # Every tissue carries an ontology id, but not all are UBERON: cell-line
    # tissues (e.g. Cells_Cultured_fibroblasts) come back as EFO terms.
    assert all(t.uberon_id for t in out)
    assert sum(t.uberon_id.startswith("UBERON:") for t in out) > 40


@pytest.mark.live
def test_live_gtex_age_bracket_subsetting_still_does_not_work():
    """Pins the limitation the adapter documents. If GTEx ever starts honouring
    attributeSubset=ageBracket on the open tier, this fails and we revisit."""
    import httpx

    resp = httpx.get(
        "https://gtexportal.org/api/v2/expression/geneExpression",
        params={
            "gencodeId": "ENSG00000147889.17",
            "tissueSiteDetailId": "Whole_Blood",
            "attributeSubset": "ageBracket",
            "datasetId": DEFAULT_DATASET_ID,
        },
        headers={"Accept": "application/json"},
        timeout=60,
        follow_redirects=True,
    )
    resp.raise_for_status()
    groups = [row.get("subsetGroup") for row in resp.json().get("data", [])]
    assert groups == [None], f"GTEx now returns age subsets: {json.dumps(groups)}"


@pytest.mark.live
def test_live_geo_soft_urls_still_serve_the_pinned_artifacts():
    """The one thing a checksum cannot tell you is whether the URL still exists.
    A HEAD over a sample of the panel catches an FTP layout change before a full
    `make data` does."""
    import httpx

    from geroquery.sources.manifest import GEO_AGING_PANEL

    sample = [GEO_AGING_PANEL[k] for k in ("GDS707", "GDS5226", "GDS156")]
    with httpx.Client(timeout=60, follow_redirects=True) as http:
        for artifact in sample:
            resp = http.head(artifact.url)
            resp.raise_for_status()
            length = int(resp.headers["content-length"])
            assert length == artifact.n_bytes, (
                f"{artifact.key} changed size upstream: pinned {artifact.n_bytes}, "
                f"now {length}. Re-verify and bump the manifest."
            )


@pytest.mark.live
def test_live_hagr_zips_still_contain_the_member_we_parse():
    """HAGR renames files between releases. A checksum failure would say the
    bytes changed; this says which structural assumption broke."""
    import io
    import zipfile

    import httpx

    from geroquery.sources.manifest import HAGR, HAGR_MEMBERS

    with httpx.Client(timeout=90, follow_redirects=True) as http:
        for key, artifact in HAGR.items():
            resp = http.get(artifact.url)
            resp.raise_for_status()
            names = zipfile.ZipFile(io.BytesIO(resp.content)).namelist()
            assert HAGR_MEMBERS[key] in names, f"{key}: expected member missing, got {names}"


@pytest.mark.live
def test_live_geo_datasets_age_query_still_returns_the_panel_population():
    """The panel is the result of a GEO DataSets query. If that query stops
    returning a comparable number of records, the panel is stale rather than
    complete."""
    import httpx

    resp = httpx.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={
            "db": "gds",
            "term": (
                "age[Subset Variable Type] AND gds[Filter] AND "
                "(Homo sapiens[Organism] OR Mus musculus[Organism])"
            ),
            "retmode": "json",
            "retmax": 1,
        },
        timeout=60,
    )
    resp.raise_for_status()
    count = int(resp.json()["esearchresult"]["count"])
    assert count >= 180, f"GEO now returns {count} age-subset DataSets; the panel query changed."


@pytest.mark.live
def test_live_every_pmid_resolves_to_the_paper_we_claim():
    """A wrong PMID is worse than no PMID: it survives the skim it invites.

    Seven of the twenty-eight references in this file once pointed at unrelated
    papers — a lamin B1 study resolved to macrophages in breast cancer, the
    PhenoAge derivation to contraceptive implant bleeding — while the module
    docstring asserted that identifiers were never invented. This compares each
    recorded title against the one PubMed returns.
    """
    import difflib

    import httpx

    from geroquery.knowledge import REFERENCES

    pmids = [r.pmid for r in REFERENCES.values() if r.pmid]
    assert pmids, "the citation layer should carry PMIDs"

    resp = httpx.post(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        data={"db": "pubmed", "id": ",".join(pmids), "retmode": "json"},
        timeout=90,
    )
    resp.raise_for_status()
    summaries = resp.json()["result"]

    mismatched = []
    for reference in REFERENCES.values():
        if not reference.pmid:
            continue
        record = summaries.get(reference.pmid) or {}
        upstream = record.get("title", "")
        similarity = difflib.SequenceMatcher(
            None, reference.title.lower()[:70], upstream.lower()[:70]
        ).ratio()
        if similarity < 0.6:
            mismatched.append(f"{reference.key} (PMID {reference.pmid}) -> {upstream[:70]!r}")

    assert not mismatched, "PMIDs pointing at the wrong paper:\n  " + "\n  ".join(mismatched)
