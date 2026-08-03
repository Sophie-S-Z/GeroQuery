"""M1 idmap — identifier & ontology harmonization.

Resolves any gene identifier (symbol, alias, Ensembl, Entrez, UniProt) to one
canonical internal gene, maps tissue labels to UBERON, and normalizes ages to
fractional lifespan via AnAge.

Design: the interface is tiny and stable (``resolve_gene``, ``map_tissue``,
``fractional_age``); the internals — which indices exist, whether we fall back
to mygene.info, how we cache — can change freely. Everything downstream depends
on canonical IDs, so this contract must not churn.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Sequence
from functools import lru_cache

from .. import config
from ..exceptions import (
    GeneNotFoundError,
    SourceError,
    TissueNotFoundError,
    UnknownSpeciesError,
)
from ..models import CanonicalGene, UberonTerm
from .mygene import MyGeneClient

_ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d{6,}$", re.IGNORECASE)
_ENTREZ_RE = re.compile(r"^\d+$")


def _norm(s: str) -> str:
    return s.strip()


class GeneResolver:
    """Loads the bundled gene table once and answers resolution queries.

    A single instance is cheap to keep around (built via :func:`get_resolver`).
    Network fallback to mygene.info is attempted only when the identifier is not
    in the bundled table *and* ``config.ALLOW_NETWORK`` is true; it never runs in
    tests/CI by default, keeping resolution hermetic.
    """

    def __init__(self, genes_path=None, tissues_path=None, anage_path=None, mygene=None):
        self._mygene = mygene if mygene is not None else MyGeneClient()
        genes_path = genes_path or (config.IDMAP_DATA / "genes.json")
        tissues_path = tissues_path or (config.IDMAP_DATA / "tissues.json")
        anage_path = anage_path or (config.IDMAP_DATA / "anage.json")

        with open(genes_path) as fh:
            gene_doc = json.load(fh)
        self._records: list[dict] = gene_doc["genes"]
        self._deprecated: dict[str, str] = {
            k.upper(): v for k, v in gene_doc.get("deprecated_ids", {}).items()
        }

        # Build lookup indices. Symbol/alias are one-to-many (multiple species),
        # so those map to *lists*; unique IDs map one-to-one.
        self._by_canonical: dict[str, dict] = {}
        self._by_entrez: dict[str, dict] = {}
        self._by_uniprot: dict[str, dict] = {}
        self._by_symbol: dict[str, list[dict]] = {}
        self._by_alias: dict[str, list[dict]] = {}
        self._by_group: dict[str, list[dict]] = {}

        for rec in self._records:
            self._by_canonical[rec["canonical_id"].upper()] = rec
            if rec.get("ensembl"):
                self._by_canonical[rec["ensembl"].upper()] = rec
            if rec.get("entrez"):
                self._by_entrez[str(rec["entrez"])] = rec
            if rec.get("uniprot"):
                self._by_uniprot[rec["uniprot"].upper()] = rec
            self._by_symbol.setdefault(rec["symbol"].lower(), []).append(rec)
            for alias in rec.get("aliases", []):
                self._by_alias.setdefault(alias.lower(), []).append(rec)
            self._by_group.setdefault(rec["ortholog_group"], []).append(rec)

        with open(tissues_path) as fh:
            self._tissue_terms: list[dict] = json.load(fh)["terms"]
        self._tissue_index: dict[str, dict] = {}
        for term in self._tissue_terms:
            for syn in term["synonyms"]:
                self._tissue_index[syn.lower()] = term

        with open(anage_path) as fh:
            anage = json.load(fh)["species"]
        self._species: dict[str, float] = {}
        for sp in anage:
            names = {sp["key"], sp["latin"], *sp.get("common", [])}
            for n in names:
                self._species[n.lower()] = sp["max_lifespan_years"]

    # ---- gene resolution -------------------------------------------------

    def _candidates(self, query: str) -> list[dict]:
        """Return every record the raw query could refer to, most-specific first."""
        q = _norm(query)
        qu = q.upper()

        # Deprecated identifier redirect happens before anything else.
        if qu in self._deprecated:
            target = self._deprecated[qu].upper()
            rec = self._by_canonical.get(target)
            return [rec] if rec else []

        # Unique-ID spaces: exact, unambiguous.
        if _ENSEMBL_RE.match(qu) and qu in self._by_canonical:
            return [self._by_canonical[qu]]
        if _ENTREZ_RE.match(q) and q in self._by_entrez:
            return [self._by_entrez[q]]
        if qu in self._by_uniprot:
            return [self._by_uniprot[qu]]

        # Symbol beats alias (a symbol is more authoritative than someone else's
        # alias of the same string).
        if q.lower() in self._by_symbol:
            return list(self._by_symbol[q.lower()])
        if q.lower() in self._by_alias:
            return list(self._by_alias[q.lower()])
        return []

    def resolve_gene(self, query: str, species: str | None = None) -> CanonicalGene:
        """Resolve ``query`` to a canonical gene.

        If ``species`` is given, the result is scoped to it. If omitted and the
        query matches more than one species (e.g. a symbol shared by human and
        mouse orthologs), a deterministic primary is chosen (human preferred)
        and the alternatives are reported in ``ambiguous_candidates`` rather than
        silently dropped.
        """
        if not query or not query.strip():
            raise GeneNotFoundError("Empty gene query.")

        candidates = self._local_candidates(query, species)

        if not candidates:
            fetched = self._fetch_remote(query, species)
            if fetched is not None:
                candidates = [fetched]

        if not candidates:
            raise GeneNotFoundError(
                f"Could not resolve gene identifier {query!r}"
                + (f" for species {species!r}" if species else "")
                + ".",
                detail={"query": query, "species": species},
            )

        return self._from_candidates(candidates, query)

    def resolve_batch(
        self, queries: Iterable[str], species: str | None = None
    ) -> dict[str, CanonicalGene | None]:
        """Resolve many identifiers; unresolved ones map to ``None``.

        Everything resolvable from the bundled table is answered locally, and the
        remainder goes upstream in a *single* batched request. Calling
        ``resolve_gene`` in a loop instead would issue one HTTP round trip per
        missing gene — for a 50-gene set that is 50 sequential requests to a free
        public API, which is both slow and the fastest way to get rate-limited.
        """
        queries = list(queries)
        out: dict[str, CanonicalGene | None] = {}
        pending: list[str] = []

        for q in queries:
            local = self._local_candidates(q, species)
            if local:
                out[q] = self._from_candidates(local, q)
            else:
                pending.append(q)

        remote = self._remote_records(pending, species)
        for q in pending:
            records = remote.get(q) or []
            if species is not None:
                sp = species.strip().lower()
                records = [r for r in records if r["species"].lower() == sp]
            out[q] = self._from_candidates(records, q) if records else None

        return {q: out.get(q) for q in queries}

    def _local_candidates(self, query: str, species: str | None) -> list[dict]:
        """Candidates from the bundled table only, scoped to species if given."""
        if not query or not query.strip():
            return []
        candidates = self._candidates(query)
        if species is not None:
            sp = species.strip().lower()
            candidates = [c for c in candidates if c["species"].lower() == sp]
        return candidates

    def _from_candidates(self, candidates: list[dict], query: str) -> CanonicalGene:
        primary = self._pick_primary(candidates)
        others = tuple(
            c["canonical_id"] for c in candidates if c["canonical_id"] != primary["canonical_id"]
        )
        return self._to_gene(primary, matched_from=query, ambiguous=others)

    def orthologs(self, canonical_id: str) -> list[CanonicalGene]:
        """All canonical genes sharing an ortholog group (cross-species view)."""
        rec = self._by_canonical.get(canonical_id.upper())
        if rec is None:
            raise GeneNotFoundError(f"Unknown canonical id {canonical_id!r}.")
        return [
            self._to_gene(r, matched_from=canonical_id)
            for r in self._by_group[rec["ortholog_group"]]
        ]

    def genes_in_group(self, group: str) -> list[CanonicalGene]:
        """All canonical genes in a named ortholog group (empty list if unknown)."""
        return [self._to_gene(r, matched_from=group) for r in self._by_group.get(group, [])]

    @staticmethod
    def _pick_primary(candidates: list[dict]) -> dict:
        # Deterministic: human first, then alphabetical by canonical_id.
        return sorted(candidates, key=lambda c: (c["species"] != "human", c["canonical_id"]))[0]

    @staticmethod
    def _to_gene(
        rec: dict, matched_from: str | None, ambiguous: tuple[str, ...] = ()
    ) -> CanonicalGene:
        return CanonicalGene(
            canonical_id=rec["canonical_id"],
            symbol=rec["symbol"],
            species=rec["species"],
            entrez=rec.get("entrez"),
            ensembl=rec.get("ensembl"),
            uniprot=rec.get("uniprot"),
            name=rec.get("name"),
            aliases=tuple(rec.get("aliases", [])),
            matched_from=matched_from,
            ambiguous_candidates=ambiguous,
        )

    def _remote_records(self, queries: Sequence[str], species: str | None) -> dict[str, list[dict]]:
        """Batch mygene.info fallback for identifiers not in the bundled table.

        Returns ``{query: [record, ...]}``. A disabled network or an upstream
        failure yields empty lists rather than an exception: this is a fallback
        behind a local table, so it must degrade to "not found" instead of taking
        down a query that was mostly resolvable locally.
        """
        if not queries:
            return {}
        try:
            return self._mygene.resolve(list(queries), species)
        except SourceError:
            return {q: [] for q in queries}

    def _fetch_remote(self, query: str, species: str | None) -> dict | None:
        """Single-identifier convenience wrapper over the batch path."""
        records = self._remote_records([query], species).get(query) or []
        return records[0] if records else None

    # ---- tissue mapping --------------------------------------------------

    def map_tissue(self, label: str) -> UberonTerm:
        if not label or not label.strip():
            raise TissueNotFoundError("Empty tissue label.")
        key = label.strip().lower()
        if key in self._tissue_index:
            term = self._tissue_index[key]
            conf = 1.0 if key == term["label"].lower() else 0.9
            return UberonTerm(term["term_id"], term["label"], conf, matched_from=label)
        # Loose contains-match as a lower-confidence fallback.
        for syn, term in self._tissue_index.items():
            if syn in key or key in syn:
                return UberonTerm(term["term_id"], term["label"], 0.6, matched_from=label)
        raise TissueNotFoundError(f"No UBERON mapping for tissue label {label!r}.", detail=label)

    # ---- age normalization ----------------------------------------------

    def max_lifespan(self, species: str) -> float:
        key = species.strip().lower()
        if key not in self._species:
            raise UnknownSpeciesError(
                f"No AnAge maximum-lifespan value for species {species!r}.", detail=species
            )
        return self._species[key]

    def fractional_age(self, age_years: float, species: str) -> float:
        """Chronological age as a fraction of the species maximum lifespan."""
        if age_years < 0:
            raise UnknownSpeciesError("Age must be non-negative.", detail=age_years)
        return float(age_years) / self.max_lifespan(species)


@lru_cache(maxsize=1)
def get_resolver() -> GeneResolver:
    """Process-wide singleton over the bundled reference data."""
    return GeneResolver()
