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
from ..harmonize import random_effects
from ..idmap import get_resolver
from ..models import GeneCard, MetaSignature
from ..resilience import ResilienceService
from ..sources import all_adapters
from ..store import GeroStore
from .cache import VersionedLRU


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

    def geneset_signature(
        self,
        gene_queries: Sequence[str],
        species: str | None = None,
        omic_layer: str | None = None,
    ) -> dict:
        resolved, unresolved, per_gene = [], [], []
        for q in gene_queries:
            try:
                gene = self.resolver.resolve_gene(q, species)
            except Exception:
                unresolved.append(q)
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
            "aggregate_pooled_effect": round(aggregate, 4) if aggregate is not None else None,
            "n_genes": len(per_gene),
            "per_gene": per_gene,
        }

    # ---- clocks ----------------------------------------------------------

    def list_clocks(self) -> list[dict]:
        return [c.to_dict() for c in self.clocks.list_clocks()]

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
        matches = self.store.interventions(name=name)
        if not matches:
            from ..exceptions import GeroQueryError

            class InterventionNotFound(GeroQueryError):
                code = "intervention_not_found"
                http_status = 404

            raise InterventionNotFound(f"Unknown intervention {name!r}.", detail=name)
        iv = matches[0]
        # Linked signatures: pooled effect per linked gene.
        linked = []
        for cid in iv.linked_gene_ids:
            sigs = self.store.query_signatures(gene_id=cid)
            metas = self._meta_from_signatures(sigs)
            if metas:
                linked.append({"gene_id": cid, "meta_signatures": [m.to_dict() for m in metas]})
        return {"intervention": iv.to_dict(), "linked_signatures": linked}

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
                from ..exceptions import ResilienceInputError

                raise ResilienceInputError("Provide either dataset_id or inline data.")
            data = self.store.get_dataset(dataset_id)
            if biomarker_cols is None:
                biomarker_cols = [c for c in data.columns if c not in ("subject_id", "age", "sex")]
        if biomarker_cols is None:
            from ..exceptions import ResilienceInputError

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
