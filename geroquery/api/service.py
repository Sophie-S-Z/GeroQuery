"""Orchestration layer: composes the lower modules into product operations.

This is the single place that knows how a gene card is assembled, how meta
signatures are pooled from raw signatures, and how a gene maps to interventions.
The HTTP layer (app.py) stays thin and just translates request/response.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .. import __version__
from ..clocks import ClockService
from ..exceptions import (
    GeneNotFoundError,
    InterventionNotFoundError,
    ResilienceInputError,
    UnknownSpeciesError,
)
from ..harmonize import random_effects
from ..idmap import get_resolver
from ..knowledge import HALLMARKS, REFERENCES
from ..models import GeneCard, MetaSignature
from ..resilience import ResilienceService
from ..sources import all_adapters
from ..store import GeroStore
from .cache import VersionedLRU

# Identifier and study-design columns, never biomarkers. `survey_weight` matters
# here: NHANES carries WTMEC2YR alongside the markers, and letting it fall into
# an auto-inferred biomarker list would put the sampling design into the health
# state itself — a variance trend that is pure survey artefact.
NON_BIOMARKER_COLUMNS = frozenset({"subject_id", "age", "sex", "survey_weight"})

# Ordering for an intervention's per-organism records. Not a claim that a mouse
# result is better evidence than a worm one — only that when a caller asks about
# a compound in a human-aging context, the mammalian result should not be buried
# behind three invertebrates by alphabetical accident.
_ORGANISM_RANK = {
    "Homo sapiens": 0,
    "Mus musculus": 1,
    "Rattus norvegicus": 2,
    "Danio rerio": 3,
    "Drosophila melanogaster": 4,
    "Caenorhabditis elegans": 5,
    "Saccharomyces cerevisiae": 6,
}


class GeroService:
    def __init__(self, store: GeroStore | None = None):
        self.store = (store or GeroStore()).ensure_built()
        self.resolver = get_resolver()
        self.clocks = ClockService()
        self.resilience = ResilienceService()
        self._cache = VersionedLRU()

    # ---- version / provenance -------------------------------------------

    def version(self) -> dict:
        return {"code_version": __version__, "data_version": self.store.version()}

    # ---- meta-analysis ---------------------------------------------------

    def _meta_from_signatures(self, signatures) -> list[MetaSignature]:
        groups: dict[tuple, list] = {}
        for s in signatures:
            if s.standard_error is None:
                continue
            groups.setdefault((s.gene_id, s.omic_layer, s.species), []).append(s)
        metas = []
        for (gene_id, omic, species), sigs in sorted(groups.items()):
            pooled = random_effects([s.effect_size for s in sigs], [s.standard_error for s in sigs])
            metas.append(
                MetaSignature(
                    gene_id=gene_id,
                    omic_layer=omic,
                    species=species,
                    pooled_effect=round(pooled.pooled_effect, 4),
                    standard_error=round(pooled.standard_error, 4),
                    ci_low=round(pooled.ci_low, 4),
                    ci_high=round(pooled.ci_high, 4),
                    p_value=pooled.p_value,
                    heterogeneity_i2=round(pooled.i2, 2),
                    tau2=round(pooled.tau2, 4),
                    n_studies=pooled.n_studies,
                    direction=pooled.direction,
                    ci_low_dl=round(pooled.ci_low_dl, 4),
                    ci_high_dl=round(pooled.ci_high_dl, 4),
                    pi_low=round(pooled.pi_low, 4) if pooled.pi_low is not None else None,
                    pi_high=round(pooled.pi_high, 4) if pooled.pi_high is not None else None,
                    verdict=pooled.verdict,
                )
            )
        return metas

    # ---- gene queries ----------------------------------------------------

    def _group_canonical_ids(self, gene, species: str | None) -> list[str]:
        orthologs = self.resolver.orthologs(gene.canonical_id)
        if species is not None:
            return [g.canonical_id for g in orthologs if g.species.lower() == species.lower()]
        return [g.canonical_id for g in orthologs]

    def gene_signature(
        self,
        gene_query: str,
        species: str | None = None,
        tissue: str | None = None,
        omic_layer: str | None = None,
        sex: str | None = None,
    ) -> dict:
        gene = self.resolver.resolve_gene(gene_query, species)
        ids = self._group_canonical_ids(gene, species)
        signatures = []
        for cid in ids:
            signatures.extend(
                self.store.query_signatures(
                    gene_id=cid, species=species, tissue=tissue, omic_layer=omic_layer, sex=sex
                )
            )
        metas = self._meta_from_signatures(signatures)
        return {
            "gene": gene.to_dict(),
            "meta_signatures": [m.to_dict() for m in metas],
            "signatures": [s.to_dict() for s in signatures],
            "n_signatures": len(signatures),
        }

    def gene_card(self, gene_query: str, species: str | None = None) -> dict:
        gene = self.resolver.resolve_gene(gene_query, species)
        cache_key = self._cache.key(self.store.version(), "card", gene.canonical_id, species)
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        ids = self._group_canonical_ids(gene, species)
        signatures = []
        for cid in ids:
            signatures.extend(self.store.query_signatures(gene_id=cid, species=species))
        metas = self._meta_from_signatures(signatures)

        curated, interventions, seen_iv = [], [], set()
        for cid in ids:
            curated.extend(self.store.curated_flags(cid))
            for iv in self.store.interventions(gene_id=cid):
                if iv.intervention_id not in seen_iv:
                    seen_iv.add(iv.intervention_id)
                    interventions.append(iv)

        card = GeneCard(
            gene=gene,
            meta_signatures=metas,
            signatures=signatures,
            curated_flags=curated,
            interventions=interventions,
        )
        out = card.to_dict()
        self._cache.set(cache_key, out)
        return out

    def gene_report(self, gene_query: str, species: str | None = None) -> dict:
        """The assembled aging profile of one gene, for the dashboard.

        A superset of :meth:`gene_card`: same measured evidence, plus the
        ortholog list, the citation layer, and a ``hallmarks`` slot.

        ``hallmarks`` is currently always empty. The hand-written per-gene
        hallmark annotation that used to fill it was dropped along with the rest
        of the curated `KNOWLEDGE` table, because the ingested GEO panel now
        measures what that table asserted — and disagrees with it for CDKN2A.
        The key is kept rather than removed so the shape is stable when a
        sourced gene-to-hallmark mapping replaces it; see
        ``docs/HANDOFF_2026-08-03.md``.
        """
        card = self.gene_card(gene_query, species)
        gene = self.resolver.resolve_gene(gene_query, species)
        orthologs = self.resolver.orthologs(gene.canonical_id)
        return {
            **card,
            "group": gene.ortholog_group or gene.symbol.upper(),
            "orthologs": [
                {"species": o.species, "symbol": o.symbol, "canonical_id": o.canonical_id}
                for o in orthologs
            ],
            "knowledge": None,
            "hallmarks": [],
            "hallmark_vocabulary": [
                {"key": key, "name": key.replace("_", " ").title(), "description": text}
                for key, text in HALLMARKS.items()
            ],
            "references": self.references(),
            "n_signatures": len(card.get("signatures", [])),
        }

    def list_curated_genes(self, limit: int | None = None) -> list[dict]:
        """Genes carrying curated evidence, for browsing and autocomplete.

        Sourced from the ingested HAGR tables rather than a hand-written list,
        so the set grows when HAGR does. Ordered by how many independent
        databases assert the gene, which is the closest thing the curated layer
        has to a confidence ranking.
        """
        counts: dict[str, dict] = {}
        for flag in self.store.all_curated_flags():
            entry = counts.setdefault(
                flag.gene_id,
                {
                    "gene_id": flag.gene_id,
                    "symbol": flag.symbol,
                    "species": flag.species,
                    "databases": set(),
                    "n_assertions": 0,
                },
            )
            entry["databases"].add(flag.database)
            entry["n_assertions"] += 1
            entry["symbol"] = entry["symbol"] or flag.symbol
            entry["species"] = entry["species"] or flag.species

        out = [
            {**entry, "databases": sorted(entry["databases"]), "symbol": entry["symbol"] or "?"}
            for entry in counts.values()
        ]
        out.sort(key=lambda r: (-len(r["databases"]), -r["n_assertions"], r["symbol"]))
        return out[:limit] if limit else out

    def references(self) -> list[dict]:
        """The literature citation layer. Real papers; never an invented PMID."""
        return [r.to_dict() for r in REFERENCES.values()]

    def geneset_signature(
        self,
        gene_queries: Sequence[str],
        species: str | None = None,
        omic_layer: str | None = None,
    ) -> dict:
        resolved, unresolved, per_gene = [], [], []
        unresolved_detail: dict[str, str] = {}
        for q in gene_queries:
            try:
                gene = self.resolver.resolve_gene(q, species)
            except (GeneNotFoundError, UnknownSpeciesError) as exc:
                # Narrow on purpose. A bare `except Exception` here also swallowed
                # network failures, bad species arguments, and genuine bugs in the
                # resolver, reporting all of them as "gene not found" — so a broken
                # id-mapping backend looked exactly like a typo in the gene list.
                unresolved.append(q)
                unresolved_detail[q] = exc.message
                continue
            resolved.append(gene.symbol)
            sigs = self.store.query_signatures(
                gene_id=gene.canonical_id, species=species, omic_layer=omic_layer
            )
            metas = self._meta_from_signatures(sigs)
            per_gene.append(
                {"gene": gene.to_dict(), "meta_signatures": [m.to_dict() for m in metas]}
            )

        pooled_effects = [m["pooled_effect"] for g in per_gene for m in g["meta_signatures"]]
        aggregate = (sum(pooled_effects) / len(pooled_effects)) if pooled_effects else None
        return {
            "resolved": resolved,
            "unresolved": unresolved,
            # Per-query reason, so a caller can tell a typo from a bad species
            # filter without re-running each query on its own.
            "unresolved_detail": unresolved_detail,
            "aggregate_pooled_effect": round(aggregate, 4) if aggregate is not None else None,
            "n_genes": len(per_gene),
            "per_gene": per_gene,
        }

    # ---- clocks ----------------------------------------------------------

    def list_clocks(self) -> list[dict]:
        return [c.to_dict() for c in self.clocks.list_clocks()]

    def clock_tiers(self) -> dict:
        """Which clock tiers loaded, and why any of them did not.

        Surfaced on /v1/clocks because "three reference clocks and nothing else"
        has two very different causes — the optional libraries are not installed,
        or they are installed and broke on import. Without this the caller cannot
        tell those apart.
        """
        registry = self.clocks.registry
        return {
            "reference": {"n_clocks": 3, "always_available": True},
            "biolearn": registry.library_status.to_dict(),
            "pyaging": registry.pyaging_status.to_dict(),
        }

    def apply_clock(
        self,
        clock_id: str,
        matrix: pd.DataFrame,
        chronological_age: Sequence[float] | None = None,
    ) -> dict:
        return self.clocks.apply_clock(clock_id, matrix, chronological_age).to_dict()

    def compare_clocks(
        self,
        clock_ids: Sequence[str],
        matrix: pd.DataFrame,
        chronological_age: Sequence[float] | None = None,
    ) -> list[dict]:
        return [
            r.to_dict() for r in self.clocks.compare_clocks(clock_ids, matrix, chronological_age)
        ]

    # ---- interventions ---------------------------------------------------

    def intervention(self, name: str) -> dict:
        """Every organism a named intervention was tested in, not just one.

        DrugAge is a table of experiments, so one compound has a record per
        organism: rapamycin extends median lifespan ~13% in mouse and ~20% in
        *C. elegans*. Returning a single record meant returning whichever
        organism happened to sort first, which for rapamycin was the worm — a
        number about nematodes presented as the answer to a question about a
        drug. All matches are returned, ordered so mammals lead.
        """
        matches = self.store.interventions(name=name)
        if not matches:
            raise InterventionNotFoundError(f"Unknown intervention {name!r}.", detail=name)
        matches = sorted(
            matches, key=lambda iv: (_ORGANISM_RANK.get(iv.organism or "", 9), iv.organism or "")
        )

        linked = []
        for cid in dict.fromkeys(c for iv in matches for c in iv.linked_gene_ids):
            metas = self._meta_from_signatures(self.store.query_signatures(gene_id=cid))
            if metas:
                linked.append({"gene_id": cid, "meta_signatures": [m.to_dict() for m in metas]})
        return {
            "name": matches[0].name,
            "interventions": [iv.to_dict() for iv in matches],
            "n": len(matches),
            "linked_signatures": linked,
        }

    # ---- resilience ------------------------------------------------------

    def resilience_csd(
        self,
        dataset_id: str | None = None,
        data: pd.DataFrame | None = None,
        biomarker_cols: Sequence[str] | None = None,
        age_col: str = "age",
        n_strata: int = 6,
        longitudinal: bool = False,
    ) -> dict:
        if data is None:
            if dataset_id is None:
                raise ResilienceInputError("Provide either dataset_id or inline data.")
            data = self.store.get_dataset(dataset_id)
            if biomarker_cols is None:
                biomarker_cols = [c for c in data.columns if c not in NON_BIOMARKER_COLUMNS]
        if biomarker_cols is None:
            raise ResilienceInputError("biomarker_cols is required for inline data.")
        result = self.resilience.csd(
            data, biomarker_cols, age_col=age_col, n_strata=n_strata, longitudinal=longitudinal
        )
        out = result.to_dict()
        out["dataset_id"] = dataset_id
        out["biomarker_cols"] = list(biomarker_cols)
        return out

    def resilience_recovery(self, series: Sequence[float]) -> dict:
        return self.resilience.recovery(series).to_dict()

    # ---- provenance ------------------------------------------------------

    def studies(self) -> list[dict]:
        return [s.to_dict() for s in self.store.list_studies()]

    def sources(self) -> list[dict]:
        out = []
        for a in all_adapters():
            cap, lic = a.capabilities(), a.license()
            out.append(
                {
                    "name": a.name,
                    "omics": list(cap.omics),
                    "species": list(cap.species),
                    "federated": cap.federated,
                    "cacheable": cap.cacheable,
                    "notes": cap.notes,
                    "license": {
                        "name": lic.name,
                        "redistributable": lic.redistributable,
                        "attribution": lic.attribution,
                    },
                }
            )
        return out

    def datasets(self) -> list[dict]:
        return self.store.list_datasets()


_SERVICE: GeroService | None = None


def get_service() -> GeroService:
    global _SERVICE
    if _SERVICE is None:
        _SERVICE = GeroService()
    return _SERVICE
