"""Internal domain types shared across modules.

These are plain dataclasses, not Pydantic models: the API layer (M7) owns the
HTTP-facing Pydantic schemas and translates to/from these. Keeping the domain
free of the web framework is what lets `store`, `harmonize`, etc. be tested with
no FastAPI import.

The fields mirror the conceptual schema in the PRD (§6 "Schema").
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field


def _clean(d: dict) -> dict:
    """asdict() but dropping None values, for compact JSON envelopes."""
    return {k: v for k, v in d.items() if v is not None}


@dataclass(frozen=True)
class CanonicalGene:
    """A gene resolved to one canonical internal identity."""

    canonical_id: str  # Ensembl gene ID, our canonical space
    symbol: str
    species: str
    entrez: str | None = None
    ensembl: str | None = None
    uniprot: str | None = None
    name: str | None = None
    aliases: tuple[str, ...] = ()
    matched_from: str | None = None  # what the caller queried with
    ambiguous_candidates: tuple[str, ...] = ()  # other canonical_ids the query could mean

    def to_dict(self) -> dict:
        d = asdict(self)
        d["aliases"] = list(self.aliases)
        d["ambiguous_candidates"] = list(self.ambiguous_candidates)
        return _clean(d)


@dataclass(frozen=True)
class UberonTerm:
    term_id: str  # e.g. "UBERON:0000178"
    label: str
    confidence: float
    matched_from: str | None = None

    def to_dict(self) -> dict:
        return _clean(asdict(self))


@dataclass(frozen=True)
class Study:
    study_id: str
    source: str
    omic_layer: str  # transcriptome | methylome | proteome | clinical
    species: str
    tissue: str | None = None
    sample_size: int | None = None
    processing_method: str | None = None
    license: str | None = None
    url: str | None = None
    version: str | None = None

    def to_dict(self) -> dict:
        return _clean(asdict(self))


@dataclass(frozen=True)
class AgingSignature:
    """One gene's age-associated change within a single study."""

    gene_id: str
    study_id: str
    omic_layer: str
    species: str
    effect_size: float  # standardized (Hedges g): +ve = up with age
    direction: str  # "up" | "down"
    tissue: str | None = None
    sex: str | None = None
    age_range: str | None = None
    p_value: float | None = None
    q_value: float | None = None
    standard_error: float | None = None
    source: str | None = None

    def to_dict(self) -> dict:
        return _clean(asdict(self))


@dataclass(frozen=True)
class MetaSignature:
    """Random-effects pool of one gene's signatures within an omic layer."""

    gene_id: str
    omic_layer: str
    species: str
    pooled_effect: float
    standard_error: float
    ci_low: float
    ci_high: float
    p_value: float
    heterogeneity_i2: float
    tau2: float
    n_studies: int
    direction: str

    def to_dict(self) -> dict:
        return _clean(asdict(self))


@dataclass(frozen=True)
class ClockInfo:
    clock_id: str
    name: str
    library: str
    predicted_outcome: str  # chronological_age | mortality | pace_of_aging
    training_population: str
    input_type: str  # methylation | expression | clinical
    units: str
    required_features: tuple[str, ...]
    notes: str | None = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["required_features"] = list(self.required_features)
        return _clean(d)


@dataclass(frozen=True)
class Intervention:
    intervention_id: str
    name: str
    itype: str  # drug | dietary | genetic
    source: str  # DrugAge | ITP | GenDR
    organism: str | None = None
    lifespan_effect_pct: float | None = None  # median lifespan change, %
    linked_gene_ids: tuple[str, ...] = ()
    url: str | None = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["linked_gene_ids"] = list(self.linked_gene_ids)
        return _clean(d)


@dataclass(frozen=True)
class CuratedFlag:
    """A curated-knowledge assertion that a gene is aging/longevity related."""

    gene_id: str
    database: str  # GenAge | OpenGenes | CellAge | LongevityMap
    assertion: str
    url: str | None = None

    def to_dict(self) -> dict:
        return _clean(asdict(self))


@dataclass
class GeneCard:
    """Assembled multi-omic + curated + intervention view of one gene."""

    gene: CanonicalGene
    meta_signatures: list[MetaSignature] = field(default_factory=list)
    signatures: list[AgingSignature] = field(default_factory=list)
    curated_flags: list[CuratedFlag] = field(default_factory=list)
    interventions: list[Intervention] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "gene": self.gene.to_dict(),
            "meta_signatures": [m.to_dict() for m in self.meta_signatures],
            "signatures": [s.to_dict() for s in self.signatures],
            "curated_flags": [c.to_dict() for c in self.curated_flags],
            "interventions": [i.to_dict() for i in self.interventions],
        }
