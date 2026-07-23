"""Cached, redistributable local sources backed by the bundled CSV slice.

These read the harmonized demonstration data produced by the ETL. They stand in
for the cache tier of GEO/GTEx/HAGR/etc.: same shape the federated adapters will
return, so the store and API never learn where a row came from.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..config import SOURCES_DATA
from ..models import AgingSignature, CuratedFlag, Intervention, Study
from .base import Capabilities, License, SourceAdapter


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


class LocalSignatureSource(SourceAdapter):
    name = "local-harmonized"

    def __init__(self, data_dir: Path | None = None):
        self.data_dir = data_dir or SOURCES_DATA

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("transcriptome", "methylome", "proteome"),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes="Bundled harmonized aging-signature slice (GEO-derived demo).",
        )

    def license(self) -> License:
        return License(
            "public-attribute", redistributable=True, attribution="NCBI GEO / harmonized"
        )

    def signatures(self) -> list[AgingSignature]:
        out = []
        for r in _read_csv(self.data_dir / "signatures.csv"):
            out.append(
                AgingSignature(
                    gene_id=r["gene_id"],
                    study_id=r["study_id"],
                    omic_layer=r["omic_layer"],
                    species=r["species"],
                    effect_size=float(r["effect_size"]),
                    direction=r["direction"],
                    tissue=r.get("tissue") or None,
                    sex=r.get("sex") or None,
                    age_range=r.get("age_range") or None,
                    p_value=float(r["p_value"]) if r.get("p_value") else None,
                    q_value=float(r["q_value"]) if r.get("q_value") else None,
                    standard_error=float(r["standard_error"]) if r.get("standard_error") else None,
                    source=r.get("source"),
                )
            )
        return out

    def studies(self) -> list[Study]:
        out = []
        for r in _read_csv(self.data_dir / "studies.csv"):
            out.append(
                Study(
                    study_id=r["study_id"],
                    source=r["source"],
                    omic_layer=r["omic_layer"],
                    species=r["species"],
                    tissue=r.get("tissue") or None,
                    sample_size=int(r["sample_size"]) if r.get("sample_size") else None,
                    processing_method=r.get("processing_method"),
                    license=r.get("license"),
                    url=r.get("url"),
                    version=r.get("version"),
                )
            )
        return out


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
