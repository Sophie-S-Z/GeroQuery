"""Cached, redistributable local sources backed by the bundled data slice.

These stand in for the cache tier of the open aging databases (HAGR/GenAge,
CellAge, LongevityMap, OpenGenes, DrugAge, NIA ITP). They expose only real,
redistributable curated facts — database memberships and lifespan interventions
— with the same adapter contract the federated sources use, so the store and API
never learn where a row came from.

The fabricated per-study "signature" adapter that earlier versions shipped has
been removed: GeroQuery no longer synthesises effect sizes or GEO accessions.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..config import SOURCES_DATA
from ..models import CuratedFlag, Intervention
from .base import Capabilities, License, SourceAdapter


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


class LocalEvidenceSource(SourceAdapter):
    """The cached, redistributable tier: curated aging knowledge + interventions.

    Kept as a first-class adapter so the licence/cache gate and the ``/sources``
    provenance view still describe a real cached source.
    """

    name = "local-curated"

    def __init__(self, data_dir: Path | None = None):
        self.data_dir = data_dir or SOURCES_DATA

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("annotation", "intervention"),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes="Bundled curated aging knowledge (GenAge/CellAge/LongevityMap/OpenGenes) "
            "and lifespan interventions (DrugAge/NIA ITP). Real, cited, redistributable.",
        )

    def license(self) -> License:
        return License(
            "HAGR/OpenGenes (free, attribute)",
            redistributable=True,
            attribution="HAGR (GenAge, CellAge, LongevityMap, DrugAge), OpenGenes, NIA ITP",
        )


class CuratedKnowledgeSource(SourceAdapter):
    name = "curated-knowledge"

    def __init__(self, data_dir: Path | None = None):
        self.data_dir = data_dir or SOURCES_DATA

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("annotation",),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes="GenAge/OpenGenes/CellAge/LongevityMap curated flags.",
        )

    def license(self) -> License:
        return License(
            "HAGR/OpenGenes (free, attribute)", redistributable=True, attribution="HAGR, OpenGenes"
        )

    def flags(self) -> list[CuratedFlag]:
        return [
            CuratedFlag(
                gene_id=r["gene_id"],
                database=r["database"],
                assertion=r["assertion"],
                url=r.get("url"),
            )
            for r in _read_csv(self.data_dir / "curated_knowledge.csv")
        ]


class InterventionSource(SourceAdapter):
    name = "interventions"

    def __init__(self, data_dir: Path | None = None):
        self.data_dir = data_dir or SOURCES_DATA

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("intervention",),
            species=("mouse",),
            federated=False,
            cacheable=True,
            notes="DrugAge / NIA ITP / GenDR lifespan-intervention slice.",
        )

    def license(self) -> License:
        return License(
            "DrugAge/ITP/GenDR (free, attribute)",
            redistributable=True,
            attribution="HAGR DrugAge, NIA ITP",
        )

    def interventions(self) -> list[Intervention]:
        out = []
        for r in _read_csv(self.data_dir / "interventions.csv"):
            linked = tuple(x for x in (r.get("linked_gene_ids") or "").split("|") if x)
            out.append(
                Intervention(
                    intervention_id=r["intervention_id"],
                    name=r["name"],
                    itype=r["itype"],
                    source=r["source"],
                    organism=r.get("organism"),
                    lifespan_effect_pct=(
                        float(r["lifespan_effect_pct"]) if r.get("lifespan_effect_pct") else None
                    ),
                    linked_gene_ids=linked,
                    url=r.get("url"),
                )
            )
        return out
