"""Orchestration layer: composes the lower modules into product operations.

This is the single place that knows how a gene report is assembled from the
curated knowledge base, gene-ID resolution, curated database flags, and linked
interventions. The HTTP layer (app.py) and the dashboard both call this.

Everything a gene report asserts is real, curated, and cited — there is no
meta-analysis of fabricated per-study statistics anymore, because those data no
longer exist in the project.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .. import __version__
from ..clocks import ClockService
from ..idmap import get_resolver
from ..knowledge import HALLMARKS, REFERENCES, gene_knowledge, interventions_for_group
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

    # ---- helpers ---------------------------------------------------------

    def _references_for(self, keys: Sequence[str]) -> list[dict]:
        out, seen = [], set()
        for k in keys:
            if k in REFERENCES and k not in seen:
                seen.add(k)
                out.append(REFERENCES[k].to_dict())
        return out

    def _curated_flags(self, ortholog_ids: Sequence[str]) -> list[dict]:
        flags, seen = [], set()
        for cid in ortholog_ids:
            for f in self.store.curated_flags(cid):
                key = (f.database, f.assertion)
                if key not in seen:
                    seen.add(key)
                    flags.append(f.to_dict())
        return flags

    # ---- gene report -----------------------------------------------------

    def gene_report(self, gene_query: str, species: str | None = None) -> dict:
        """The assembled, real-evidence aging profile of one gene."""
        gene = self.resolver.resolve_gene(gene_query, species)
        cache_key = self._cache.key(self.store.version(), "report", gene.canonical_id, species)
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        group = gene.ortholog_group
        orthologs = self.resolver.orthologs(gene.canonical_id)
        ortholog_ids = [o.canonical_id for o in orthologs]

        know = gene_knowledge(group) if group else None

        # Interventions linked to this gene, enriched with their citations.
        interventions = []
        ref_keys: list[str] = []
        if know:
            ref_keys.extend(know.reference_keys())
        if group:
            for iv in interventions_for_group(group):
                d = iv.to_dict()
                d["references"] = self._references_for(iv.reference_keys)
                interventions.append(d)
                ref_keys.extend(iv.reference_keys)

        hallmarks = []
        if know:
            for h in know.hallmarks:
                hallmarks.append(
                    {
                        "key": h,
                        "name": h.replace("_", " ").title(),
                        "description": HALLMARKS.get(h, ""),
                    }
                )

        report = {
            "gene": gene.to_dict(),
            "group": group,
            "has_knowledge": know is not None,
            "orthologs": [
                {"species": o.species, "symbol": o.symbol, "canonical_id": o.canonical_id}
                for o in orthologs
            ],
            "knowledge": know.to_dict() if know else None,
            "hallmarks": hallmarks,
            "curated_flags": self._curated_flags(ortholog_ids),
            "interventions": interventions,
            "references": self._references_for(ref_keys),
        }
        self._cache.set(cache_key, report)
        return report

    # Back-compat alias for the assembled view.
    def gene_card(self, gene_query: str, species: str | None = None) -> dict:
        return self.gene_report(gene_query, species)

    def list_curated_genes(self) -> list[dict]:
        """The curated gene set, for browsing / autocomplete."""
        from ..knowledge import KNOWLEDGE

        out = []
        for group, know in KNOWLEDGE.items():
            genes = self.resolver.genes_in_group(group)
            human = next((g for g in genes if g.species == "human"), genes[0] if genes else None)
            if human is None:
                continue
            out.append(
                {
                    "group": group,
                    "symbol": human.symbol,
                    "name": human.name,
                    "direction_with_age": know.direction_with_age,
                    "confidence": know.confidence,
                    "one_liner": know.one_liner,
                }
            )
        return sorted(out, key=lambda d: d["symbol"])

    def geneset_summary(self, gene_queries: Sequence[str], species: str | None = None) -> dict:
        """Direction-of-change summary across a set of genes."""
        resolved, unresolved, per_gene = [], [], []
        counts = {"up": 0, "down": 0, "context-dependent": 0}
        for q in gene_queries:
            try:
                gene = self.resolver.resolve_gene(q, species)
            except Exception:  # noqa: BLE001
                unresolved.append(q)
                continue
            resolved.append(gene.symbol)
            know = gene_knowledge(gene.ortholog_group) if gene.ortholog_group else None
            direction = know.direction_with_age if know else "unknown"
            if direction in counts:
                counts[direction] += 1
            per_gene.append(
                {
                    "symbol": gene.symbol,
                    "group": gene.ortholog_group,
                    "direction_with_age": direction,
                    "confidence": know.confidence if know else None,
                }
            )
        return {
            "resolved": resolved,
            "unresolved": unresolved,
            "n_genes": len(per_gene),
            "direction_counts": counts,
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
        from ..knowledge import INTERVENTIONS

        know = INTERVENTIONS.get(iv.name)
        linked_genes = []
        seen = set()
        for cid in iv.linked_gene_ids:
            try:
                orths = self.resolver.orthologs(cid)
            except Exception:  # noqa: BLE001
                continue
            for o in orths:
                if o.ortholog_group not in seen:
                    seen.add(o.ortholog_group)
                    linked_genes.append({"group": o.ortholog_group})
        out = iv.to_dict()
        out["display_name"] = know.display_name if know else iv.name
        if know:
            out["summary"] = know.summary
            out["references"] = self._references_for(know.reference_keys)
        out["linked_groups"] = [g["group"] for g in linked_genes]
        return {"intervention": out}

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

    def references(self) -> list[dict]:
        return [r.to_dict() for r in REFERENCES.values()]

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
