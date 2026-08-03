"""mygene.info client — batch identifier resolution with an on-disk cache.

Used as the fallback when an identifier is not in the bundled gene table. Two
properties that matter:

**Batch, not per-gene.** ``POST /v3/query`` resolves a whole list in one request.
Resolving a 50-gene set one HTTP call at a time is ~50x the latency and is exactly
the access pattern that gets a client rate-limited.

**Cached on disk.** Identifier mappings change on the scale of Ensembl releases,
not minutes. Re-asking mygene.info for the same symbol on every process start is
wasted latency for the caller and wasted load on a free public service.

Responses are normalized into the same record shape the bundled table uses, so
:class:`~geroquery.idmap.resolver.GeneResolver` never learns where a record came
from. That is deliberate: swapping this for a local Ensembl dump should not touch
the resolver.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from pathlib import Path

from .. import config
from ..exceptions import SourceError

MYGENE_QUERY_URL = "https://mygene.info/v3/query"

# Identifier spaces we accept as input. Order does not imply priority; mygene
# scores hits and we sort deterministically ourselves.
DEFAULT_SCOPES = "symbol,alias,ensembl.gene,entrezgene,uniprot"

DEFAULT_FIELDS = "symbol,name,taxid,entrezgene,ensembl.gene,uniprot.Swiss-Prot,alias"

# The two species the rest of GeroQuery models. Anything else is returned with an
# explicit "unknown" rather than silently coerced into one of these.
TAXID_TO_SPECIES = {9606: "human", 10090: "mouse"}
SPECIES_TO_TAXID = {v: k for k, v in TAXID_TO_SPECIES.items()}

# mygene.info's documented ceiling for a POST query batch.
MAX_BATCH = 1000


def _first(value):
    """mygene returns some fields as either a scalar or a list. Normalize.

    ``uniprot.Swiss-Prot`` is the usual offender: "P04406" for GAPDH but
    ``["P42771", "Q8N726"]`` for CDKN2A. Storing the list where a string is
    expected puts an unhashable value into an identifier index.
    """
    if isinstance(value, list):
        return value[0] if value else None
    return value


# Ensembl stable-id prefixes for the species we model. Used to tell a gene's own
# id apart from a cross-species cross-reference.
ENSEMBL_PREFIX = {9606: "ENSG", 10090: "ENSMUSG"}


def _ensembl_gene(hit: dict) -> str | None:
    """The gene's own Ensembl id, or ``None``.

    ``ensembl`` is sometimes a list. Taking its first element is wrong: mygene
    puts homology cross-references in there, so a mouse gene's list can lead with
    ``ENSFALG...`` (collared flycatcher) or ``ENSMPUG...`` (ferret). That put
    other species' identifiers into the canonical id space, where they silently
    became the join key for a mouse gene's signatures.

    So candidates are filtered to the prefix that matches the hit's taxid. If the
    list contains no id from the expected assembly, the whole list is
    cross-references and ``None`` is returned — the caller then falls back to an
    ``ENTREZ:`` id, which is at least the right gene.
    """
    ens = hit.get("ensembl")
    if isinstance(ens, dict):
        candidates = [_first(ens.get("gene"))]
    elif isinstance(ens, list):
        candidates = [_first(e.get("gene")) for e in ens if isinstance(e, dict)]
    else:
        return None

    genes = sorted(c for c in candidates if isinstance(c, str) and c)
    if not genes:
        return None
    taxid = hit.get("taxid")
    prefix = ENSEMBL_PREFIX.get(taxid) if isinstance(taxid, int) else None
    if prefix is None:
        return genes[0]
    matching = [g for g in genes if g.startswith(prefix)]
    return matching[0] if matching else None


def normalize_hit(hit: dict) -> dict | None:
    """Convert one mygene hit into a resolver record. ``None`` if unusable."""
    if hit.get("notfound"):
        return None
    symbol = hit.get("symbol")
    entrez = hit.get("entrezgene")
    ensembl = _ensembl_gene(hit)
    if not symbol or not (ensembl or entrez):
        return None

    taxid = hit.get("taxid")
    species = TAXID_TO_SPECIES.get(taxid, "unknown") if isinstance(taxid, int) else "unknown"
    uniprot = hit.get("uniprot")
    swiss = _first(uniprot.get("Swiss-Prot")) if isinstance(uniprot, dict) else None

    return {
        "canonical_id": (ensembl or f"ENTREZ:{entrez}").upper(),
        "symbol": symbol,
        "species": species,
        "entrez": str(entrez) if entrez else None,
        "ensembl": ensembl,
        "uniprot": swiss,
        "name": hit.get("name"),
        "aliases": [a for a in (hit.get("alias") or []) if isinstance(a, str)],
        # Ortholog grouping upstream is by upper-cased symbol, which is how the
        # bundled table groups human/mouse pairs (CDKN2A / Cdkn2a -> "CDKN2A").
        "ortholog_group": symbol.upper(),
        "source": "mygene.info",
    }


class MyGeneClient:
    """Batch resolver against mygene.info, with a JSON file cache."""

    def __init__(
        self,
        cache_dir: Path | None = None,
        allow_network: bool | None = None,
        timeout: float | None = None,
    ):
        self.cache_dir = Path(cache_dir) if cache_dir else Path(config.CACHE_HOME) / "mygene"
        self._allow_network = allow_network
        self._timeout = timeout

    @property
    def allow_network(self) -> bool:
        return config.ALLOW_NETWORK if self._allow_network is None else self._allow_network

    @property
    def timeout(self) -> float:
        return config.HTTP_TIMEOUT if self._timeout is None else self._timeout

    # ---- cache -----------------------------------------------------------

    def _cache_path(self, query: str, species: str | None) -> Path:
        key = f"{query.strip().lower()}|{(species or '*').lower()}"
        digest = hashlib.sha256(key.encode()).hexdigest()[:20]
        return self.cache_dir / f"{digest}.json"

    def _read_cache(self, query: str, species: str | None) -> list[dict] | None:
        path = self._cache_path(query, species)
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            # A corrupt cache entry is a cache miss, not a failure.
            return None

    def _write_cache(self, query: str, species: str | None, records: list[dict]) -> None:
        path = self._cache_path(query, species)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(records), encoding="utf-8")
        tmp.replace(path)

    # ---- query -----------------------------------------------------------

    def _post(self, queries: Sequence[str], species: str | None) -> list[dict]:
        import httpx

        data = {
            "q": ",".join(queries),
            "scopes": DEFAULT_SCOPES,
            "fields": DEFAULT_FIELDS,
            "species": species if species else ",".join(SPECIES_TO_TAXID),
        }
        try:
            resp = httpx.post(MYGENE_QUERY_URL, data=data, timeout=self.timeout)
            resp.raise_for_status()
            payload = resp.json()
        except Exception as exc:  # httpx errors, JSON errors
            raise SourceError(
                f"mygene.info batch query failed: {exc}",
                detail={"n_queries": len(queries), "species": species},
            ) from exc
        # A single-item batch can come back as a bare object.
        if isinstance(payload, dict):
            payload = [payload]
        return payload

    def resolve(self, queries: Sequence[str], species: str | None = None) -> dict[str, list[dict]]:
        """Resolve identifiers to lists of records, keyed by the original query.

        A query that resolves to nothing maps to an empty list. Cached queries are
        served without network access; if network is disabled, uncached queries
        simply return empty rather than raising — this is a *fallback* path, and a
        disabled network should degrade to "not found", not to an error.
        """
        wanted = [q for q in dict.fromkeys(q.strip() for q in queries) if q]
        out: dict[str, list[dict]] = {}
        misses: list[str] = []

        for q in wanted:
            cached = self._read_cache(q, species)
            if cached is None:
                misses.append(q)
            else:
                out[q] = cached

        if misses and self.allow_network:
            for start in range(0, len(misses), MAX_BATCH):
                chunk = misses[start : start + MAX_BATCH]
                grouped: dict[str, list[dict]] = {q: [] for q in chunk}
                for hit in self._post(chunk, species):
                    key = hit.get("query")
                    record = normalize_hit(hit)
                    if isinstance(key, str) and key in grouped and record is not None:
                        grouped[key].append(record)
                for q, records in grouped.items():
                    # Deterministic order: human first, then by canonical id, so
                    # repeated runs pick the same primary regardless of scoring.
                    records.sort(key=lambda r: (r["species"] != "human", r["canonical_id"]))
                    self._write_cache(q, species, records)
                    out[q] = records
        else:
            for q in misses:
                out[q] = []

        return {q: out.get(q, []) for q in wanted}
