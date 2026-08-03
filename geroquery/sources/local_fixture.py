"""Cached, redistributable local sources backed by the tables the ETL writes.

Every table these read is now real: signatures come from the checksum-pinned GEO
DataSets panel via :mod:`geroquery.etl.build_signatures`, curated knowledge and
interventions from the HAGR releases. Nothing here is generated.

Same shape a federated adapter would return, so the store and API never learn
where a row came from.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..config import SOURCES_DATA
from ..exceptions import SourceError
from ..models import AgingSignature, CuratedFlag, Intervention, Study
from .base import Capabilities, License, SourceAdapter
from .manifest import GEO_LICENSE

# Written by `make signatures`, git-ignored: reproducible from the manifest, so a
# committed copy would be an unverifiable second one.
FULL_SIGNATURES = "signatures_full.csv"

# Committed. The same real estimates, restricted to genes carrying a HAGR
# assertion, so tests and CI run without a 313 MB download.
CURATED_SIGNATURES = "signatures_curated.csv"


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


class LocalSignatureSource(SourceAdapter):
    name = "local-harmonized"

    def __init__(self, data_dir: Path | None = None, *, prefer_full: bool = True):
        self.data_dir = data_dir or SOURCES_DATA
        self.prefer_full = prefer_full

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("transcriptome",),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes=(
                "GEO DataSets aging panel: young-vs-old Hedges' g per gene, "
                f"served from {self.signature_path().name}."
            ),
        )

    def license(self) -> License:
        return License(GEO_LICENSE, redistributable=True, attribution="NCBI GEO")

    def signature_path(self) -> Path:
        """The table in use: the full panel when built, else the committed slice."""
        full = self.data_dir / FULL_SIGNATURES
        if self.prefer_full and full.exists():
            return full
        return self.data_dir / CURATED_SIGNATURES

    def mode(self) -> str:
        """``"full"`` or ``"curated_slice"``.

        Returned rather than only logged, for the same reason the NHANES adapter
        returns its mode: a number computed on the curated slice covers ~2,700
        genes, not ~20,000, and must not be reported as the full-panel result.
        """
        return "full" if self.signature_path().name == FULL_SIGNATURES else "curated_slice"

    def signatures(self) -> list[AgingSignature]:
        path = self.signature_path()
        rows = _read_csv(path)
        if not rows:
            raise SourceError(
                f"No signature table at {path}. Run `make signatures` to build the full "
                f"panel, or restore the committed {CURATED_SIGNATURES}.",
                detail={"path": str(path)},
            )
        out = []
        for r in rows:
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
                    series_id=r.get("series_id") or None,
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
