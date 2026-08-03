"""Bulk identifier resolution for the offline ETL.

:class:`~geroquery.idmap.mygene.MyGeneClient` caches one JSON file per query,
which is the right shape for the API path: a request resolves a handful of
identifiers and wants each of them individually re-usable. The ETL resolves the
entire probe annotation of thirty microarray platforms at once — tens of
thousands of identifiers — and one file each would put ~50,000 sub-kilobyte
files in the cache directory.

So this module keeps one consolidated map per species instead. It reuses
:func:`~geroquery.idmap.mygene.normalize_hit`, so the record shape stays
identical to the request-path resolver and a gene resolved here joins to a gene
resolved there. Only the caching strategy differs.

Unresolvable identifiers are stored as ``None`` rather than omitted. Without
that, every re-run would re-ask mygene.info for the same identifiers that failed
last time — which for a fixed panel of arrays is thousands of pointless queries
against a free public service on every build.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from .. import config
from ..exceptions import NetworkDisabledError
from .mygene import MAX_BATCH, MYGENE_QUERY_URL, normalize_hit

# Entrez first: it is what the GEO platform annotations supply and what HAGR
# supplies for four of its five gene tables. Symbol and alias catch LongevityMap,
# which indexes by symbol only.
BULK_SCOPES = "entrezgene,symbol,alias"
BULK_FIELDS = "symbol,name,taxid,entrezgene,ensembl.gene,uniprot.Swiss-Prot,alias"


def map_path(species: str, cache_dir: Path | None = None) -> Path:
    directory = Path(cache_dir) if cache_dir else Path(config.CACHE_HOME) / "idmap"
    return directory / f"bulk_{species}.json"


def _load(path: Path) -> dict[str, dict | None]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _save(path: Path, mapping: dict[str, dict | None]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(mapping, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def _query(queries: list[str], species: str, timeout: float) -> list[dict]:
    import httpx

    response = httpx.post(
        MYGENE_QUERY_URL,
        data={
            "q": ",".join(queries),
            "scopes": BULK_SCOPES,
            "fields": BULK_FIELDS,
            "species": species,
        },
        timeout=timeout,
    )
    response.raise_for_status()
    payload = response.json()
    return [payload] if isinstance(payload, dict) else payload


def resolve_bulk(
    queries: Iterable[str],
    species: str,
    *,
    cache_dir: Path | None = None,
    allow_network: bool | None = None,
    timeout: float | None = None,
    progress: bool = False,
) -> dict[str, dict | None]:
    """Resolve many identifiers of one species to canonical gene records.

    Args:
        queries: Entrez ids or gene symbols.
        species: ``"human"`` or ``"mouse"``.
        cache_dir: override for the consolidated map location.
        allow_network: override ``GEROQUERY_ALLOW_NETWORK`` for this call.
        timeout: per-request timeout in seconds.
        progress: print one line per batch. The ETL runs for minutes; silence
            makes a stalled network indistinguishable from slow progress.

    Returns:
        ``{query: record or None}`` for every input, deduplicated.

    Raises:
        NetworkDisabledError: identifiers are missing from the map and network
            access is off. Unlike the request-path resolver this is fatal: a
            build that silently dropped every unresolved gene would produce a
            plausible-looking but incomplete signature table.
    """
    allow = config.ALLOW_NETWORK if allow_network is None else allow_network
    timeout = config.HTTP_TIMEOUT if timeout is None else timeout

    path = map_path(species, cache_dir)
    mapping = _load(path)
    wanted = list(dict.fromkeys(q.strip() for q in queries if q and q.strip()))
    missing = [q for q in wanted if q not in mapping]

    if missing and not allow:
        raise NetworkDisabledError(
            f"{len(missing)} {species} identifiers are not in the local id map and network "
            f"access is disabled. Run `make data` (or set GEROQUERY_ALLOW_NETWORK=1).",
            detail={"species": species, "missing": len(missing), "map_path": str(path)},
        )

    for start in range(0, len(missing), MAX_BATCH):
        chunk = missing[start : start + MAX_BATCH]
        resolved: dict[str, dict | None] = dict.fromkeys(chunk, None)
        for hit in _query(chunk, species, timeout):
            key = hit.get("query")
            record = normalize_hit(hit)
            if not (isinstance(key, str) and key in resolved and record is not None):
                continue
            # Defence in depth behind the species filter in the query itself: a
            # record whose taxid is not the species we asked for is discarded
            # rather than written into that species' map.
            if record.get("species") != species:
                continue
            # mygene can return several hits per query; the first accepted one
            # wins, and the batch order is stable, so a rebuild picks the same
            # primary rather than whichever scored marginally higher.
            if resolved[key] is None:
                resolved[key] = record
        mapping.update(resolved)
        if progress:
            done = min(start + MAX_BATCH, len(missing))
            print(f"  idmap {species}: {done}/{len(missing)} resolved", flush=True)

    if missing:
        _save(path, mapping)
    return {q: mapping.get(q) for q in wanted}
