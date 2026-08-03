"""NCBI GEO DataSets (GDS) adapter — real gene aging signatures.

This is the source behind GeroQuery's multi-omic signature path. It replaces the
generated demonstration slice that used to live in ``etl/build_fixtures.py``.

**Why GEO DataSets and not GEO Series.** A GSE series matrix carries sample
metadata as free text: age arrives as ``age: 67``, ``Age (yrs): 67``, ``67y``, or
buried in a sample title, differently in every series. A GDS is the curator-built
view of a series: age is a declared *subset variable* with named levels and
explicit sample lists, and the value matrix ships with the platform annotation
already joined, so probes carry a gene symbol and an Entrez id. Both of the
fragile steps in GEO ingestion — deciding who is old, and deciding what a probe
measures — are answered by GEO's curators rather than by a regex of ours.

The cost is coverage: GEO stopped producing new GDS records, so this panel is
microarray-era and the per-study sample sizes are small. That shows up honestly
as wide confidence intervals and high heterogeneity in the pooled estimates.

**Confounding is handled by rule, not by hand.** Most aging datasets vary
something else too — a disease, a drug, a knockout, a diet. Each non-age subset
variable is dispatched by :data:`STRATIFIED_VARIABLES` (produce one contrast per
level, because these are genuinely different tissues), :data:`CONTROLLED_VARIABLES`
(restrict to a control level, or drop the dataset if none is recognisable), or
pooled. Every restriction applied is recorded on the resulting contrast, so a
signature row can be traced back to the exact sample set that produced it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import pandas as pd

from ..exceptions import SourceError
from .base import Capabilities, License, SourceAdapter

GDS_FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/geo/datasets"
GDS_BROWSER = "https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc="

LICENSE = "NCBI GEO — US public domain, attribute"
ATTRIBUTION = (
    "National Center for Biotechnology Information (NCBI), Gene Expression Omnibus (GEO). "
    "Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and "
    "hybridization array data repository. Nucleic Acids Res. 2002;30(1):207-10."
)

# --- age vocabulary ---------------------------------------------------------

# Multiplier from each unit to years. GDS age levels are free-text but drawn
# from a small vocabulary; every form below appears in the real records.
_YEARS_PER: dict[str, float] = {
    "y": 1.0,
    "yr": 1.0,
    "yrs": 1.0,
    "year": 1.0,
    "years": 1.0,
    "m": 1 / 12,
    "mo": 1 / 12,
    "mon": 1 / 12,
    "month": 1 / 12,
    "months": 1 / 12,
    "w": 1 / 52.18,
    "wk": 1 / 52.18,
    "wks": 1 / 52.18,
    "week": 1 / 52.18,
    "weeks": 1 / 52.18,
    "d": 1 / 365.25,
    "day": 1 / 365.25,
    "days": 1 / 365.25,
}
_UNIT_RE = "|".join(sorted(_YEARS_PER, key=len, reverse=True))

QUALITATIVE_YOUNG = frozenset({"young", "young adult", "younger", "juvenile"})
QUALITATIVE_OLD = frozenset({"old", "older", "aged", "elderly", "very old", "senescent"})

# (young_lower, young_upper, old_lower) in years, per organism.
#
# Human: the lower bound excludes children for the same reason the NHANES adapter
# does — paediatric physiology differs enough that including it would make part
# of any age trend a developmental artefact rather than aging. Mouse: 2-8 months
# is young adult and 18 months is the conventional start of old age for the
# inbred strains these studies use.
AGE_BANDS: dict[str, tuple[float, float, float]] = {
    "human": (18.0, 40.0, 60.0),
    "mouse": (2 / 12, 8 / 12, 18 / 12),
}

ORGANISM_TO_SPECIES = {"Homo sapiens": "human", "Mus musculus": "mouse"}

# --- subset-variable dispatch ----------------------------------------------

# Variables that describe *where* the measurement was taken. One contrast per
# level: an age effect in bone-marrow adipocytes and one in epididymal adipocytes
# are two findings, and merging them would put between-tissue variance into the
# within-group SD and shrink both.
STRATIFIED_VARIABLES = frozenset({"tissue", "cell type", "specimen"})

# Variables that introduce a non-age difference between samples. The contrast is
# restricted to the control level; a dataset with no recognisable control level
# is dropped rather than guessed at.
CONTROLLED_VARIABLES = frozenset(
    {
        "disease state",
        "agent",
        "genotype/variation",
        "protocol",
        "stress",
        "dose",
        "infection",
        "time",
        "development stage",
    }
)

# Level descriptions that identify the untouched arm of a controlled variable.
# Matched after lowercasing and collapsing whitespace; a level also matches if it
# *starts with* one of these followed by a space or bracket, which is what
# catches "control (healthy, age-matched)" and "wild type control".
CONTROL_LEVELS = frozenset(
    {
        "control",
        "controls",
        "control fed",
        "control diet",
        "normal diet",
        "ad libitum",
        "ad lib",
        "normal",
        "normal karyotype",
        "healthy",
        "non-diseased",
        "untreated",
        "unaffected",
        "none",
        "wild type",
        "wildtype",
        "wt",
        "vehicle",
        "sham",
        "mock",
        "uninfected",
        "baseline",
        "basal",
        "sedentary",
        "non-smoker",
        "nonsmoker",
        "never smoker",
    }
)

# Pooled rather than split: sex is carried onto the signature row instead, and
# splitting every contrast by it would halve already small groups.
SEX_VARIABLES = frozenset({"gender", "sex"})

# Recorded but neither split nor restricted. "individual" is sample identity;
# "strain" and "other" cover background variation that pooling absorbs into the
# within-group SD, which biases the effect size towards zero — a conservative
# direction, and preferable to dropping the study.
POOLED_VARIABLES = frozenset({"individual", "strain", "other", "age"})


@dataclass(frozen=True)
class Subset:
    """One level of one GDS subset variable."""

    description: str
    variable: str
    sample_ids: tuple[str, ...]


@dataclass(frozen=True)
class GdsHeader:
    """Everything in a GDS SOFT file except the value matrix."""

    accession: str
    title: str
    organism: str
    platform: str
    series_id: str
    pubmed_id: str | None
    value_type: str
    sample_count: int
    update_date: str | None
    subsets: tuple[Subset, ...]
    sample_sources: dict[str, str] = field(default_factory=dict)

    @property
    def species(self) -> str | None:
        return ORGANISM_TO_SPECIES.get(self.organism)

    @property
    def url(self) -> str:
        return f"{GDS_BROWSER}{self.accession}"

    def by_variable(self, variable: str) -> list[Subset]:
        return [s for s in self.subsets if s.variable == variable]

    @property
    def variables(self) -> list[str]:
        seen: dict[str, None] = {}
        for s in self.subsets:
            seen.setdefault(s.variable, None)
        return list(seen)


@dataclass(frozen=True)
class AgeContrast:
    """A young-vs-old sample split within one dataset, ready to estimate."""

    accession: str
    species: str
    tissue: str | None
    sex: str
    age_range: str
    young_samples: tuple[str, ...]
    old_samples: tuple[str, ...]
    restrictions: tuple[str, ...]
    # The subset level this contrast was stratified to, verbatim. Kept apart from
    # ``tissue`` because tissue is normalized to a coarse vocabulary: two strata
    # of one dataset can share a tissue ("bone marrow adipocyte" and "epididymal
    # adipocyte" are both adipose tissue) and would then collide on study_id.
    stratum: str | None = None

    @property
    def study_id(self) -> str:
        """``GEO:GDS5216`` for a whole dataset, suffixed when it was stratified."""
        base = f"GEO:{self.accession}"
        return f"{base}:{_slug(self.stratum)}" if self.stratum else base


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def gds_directory(accession: str) -> str:
    """FTP shard for an accession: ``GDS707`` -> ``GDSnnn``, ``GDS5216`` -> ``GDS5nnn``."""
    digits = accession[3:]
    if not accession.startswith("GDS") or not digits.isdigit():
        raise ValueError(f"Not a GDS accession: {accession!r}")
    return f"GDS{digits[:-3]}nnn" if len(digits) > 3 else "GDSnnn"


def soft_url(accession: str) -> str:
    """Canonical URL of a GDS full SOFT file (matrix plus platform annotation)."""
    return f"{GDS_FTP_ROOT}/{gds_directory(accession)}/{accession}/soft/{accession}_full.soft.gz"


# --- parsing ----------------------------------------------------------------

_TABLE_BEGIN = "!dataset_table_begin"
_TABLE_END = "!dataset_table_end"
_SAMPLE_SRC = re.compile(r"^#(GSM\d+) = .*?;\s*src:\s*(.+?)\s*$")


def parse_header(text: str) -> GdsHeader:
    """Parse the metadata block of a GDS SOFT file.

    Accepts the whole file or just the part before ``!dataset_table_begin``, so
    the panel can be triaged from range-requested prefixes without downloading
    every matrix.
    """
    head = text.split(_TABLE_BEGIN, 1)[0]
    meta: dict[str, str] = {}
    subsets: list[Subset] = []
    sources: dict[str, str] = {}
    current: dict[str, str] = {}

    for raw in head.splitlines():
        line = raw.rstrip()
        if line.startswith("^SUBSET"):
            current = {}
        elif line.startswith("!subset_"):
            key, _, value = line[1:].partition(" = ")
            current[key[len("subset_") :]] = value
            if {"description", "type", "sample_id"} <= current.keys():
                subsets.append(
                    Subset(
                        description=current["description"],
                        variable=current["type"],
                        sample_ids=tuple(s.strip() for s in current["sample_id"].split(",")),
                    )
                )
                current = {}
        elif line.startswith("!dataset_"):
            key, _, value = line[1:].partition(" = ")
            meta.setdefault(key[len("dataset_") :], value)
        elif line.startswith("#GSM"):
            match = _SAMPLE_SRC.match(line)
            if match:
                sources[match.group(1)] = match.group(2)

    accession = _first_dataset_accession(head)
    if not accession:
        raise SourceError("SOFT text has no ^DATASET accession.", detail={"head": head[:200]})

    return GdsHeader(
        accession=accession,
        title=meta.get("title", ""),
        organism=meta.get("sample_organism", ""),
        platform=meta.get("platform", ""),
        series_id=meta.get("reference_series", ""),
        pubmed_id=meta.get("pubmed_id") or None,
        value_type=meta.get("value_type", ""),
        sample_count=int(meta.get("sample_count", "0") or 0),
        update_date=meta.get("update_date") or None,
        subsets=tuple(subsets),
        sample_sources=sources,
    )


def _first_dataset_accession(head: str) -> str:
    for line in head.splitlines():
        if line.startswith("^DATASET"):
            return line.partition(" = ")[2].strip()
    return ""


def parse_table(text: str) -> pd.DataFrame:
    """Parse the value matrix, indexed by probe id.

    Keeps ``IDENTIFIER`` (gene symbol) and ``Gene ID`` (Entrez) where present and
    drops the remaining annotation columns — UniGene, GenBank, and chromosome
    strings are megabytes of text this pipeline never reads.
    """
    if _TABLE_BEGIN not in text:
        raise SourceError("SOFT text contains no dataset table.", detail={"len": len(text)})
    body = text.split(_TABLE_BEGIN, 1)[1].split(_TABLE_END, 1)[0].strip("\n")
    frame = pd.read_csv(StringIO(body), sep="\t", dtype=str, na_values=["null", "NULL", ""])
    if "ID_REF" not in frame.columns:
        raise SourceError("Dataset table has no ID_REF column.", detail=list(frame.columns)[:10])

    keep_meta = [c for c in ("IDENTIFIER", "Gene ID", "Gene symbol") if c in frame.columns]
    sample_cols = [c for c in frame.columns if c.startswith("GSM")]
    frame = frame.set_index("ID_REF")[keep_meta + sample_cols]
    for col in sample_cols:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def read_soft(path: Path) -> tuple[GdsHeader, pd.DataFrame]:
    """Read a (optionally gzipped) GDS SOFT file into a header and a table."""
    import gzip

    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    else:
        text = path.read_text(encoding="utf-8", errors="replace")
    return parse_header(text), parse_table(text)


# --- age classification -----------------------------------------------------


def age_in_years(description: str) -> float | None:
    """Convert a GDS age-level description to years, or ``None`` if it is not a
    postnatal numeric age.

    Prenatal levels (``E16.5``, "embryonic day 14") return ``None`` rather than a
    negative age: they are not part of any aging contrast.
    """
    text = " ".join(description.strip().lower().split())
    if re.match(r"^e\s?\d", text) or any(w in text for w in ("embryo", "fetal", "gestation")):
        return None

    match = re.match(r"^p\s?(\d+(?:\.\d+)?)$", text)  # P7 = postnatal day 7
    if match:
        return float(match.group(1)) * _YEARS_PER["d"]

    match = re.search(rf"(\d+(?:\.\d+)?)\s*(?:-|to|~)\s*(\d+(?:\.\d+)?)\s*({_UNIT_RE})\b", text)
    if match:  # "20 - 39 y" -> midpoint
        return (float(match.group(1)) + float(match.group(2))) / 2 * _YEARS_PER[match.group(3)]

    match = re.search(r"(?:postnatal\s+)?days?\s+(\d+(?:\.\d+)?)", text)
    if match:
        return float(match.group(1)) * _YEARS_PER["d"]

    match = re.search(rf"(\d+(?:\.\d+)?)\s*({_UNIT_RE})\b", text)
    if match:
        return float(match.group(1)) * _YEARS_PER[match.group(2)]
    return None


def age_group(description: str, species: str) -> str | None:
    """Classify one age level as ``"young"``, ``"old"``, or neither.

    Qualitative levels are honoured only as themselves: "adult" and
    "middle-aged" fall in the gap between the bands on purpose, because a
    contrast built out of them would not be a young-vs-old contrast.
    """
    text = " ".join(description.strip().lower().split())
    if text in QUALITATIVE_YOUNG:
        return "young"
    if text in QUALITATIVE_OLD:
        return "old"

    years = age_in_years(description)
    if years is None or species not in AGE_BANDS:
        return None
    young_lo, young_hi, old_lo = AGE_BANDS[species]
    if young_lo <= years <= young_hi:
        return "young"
    if years >= old_lo:
        return "old"
    return None


def is_control_level(description: str) -> bool:
    """True if a subset level names the untouched arm of its variable."""
    text = " ".join(description.strip().lower().split())
    if text in CONTROL_LEVELS:
        return True
    return any(
        text.startswith(f"{level} ") or text.startswith(f"{level}(") for level in CONTROL_LEVELS
    )


# --- contrast construction --------------------------------------------------

MIN_GROUP_SIZE = 3

# Coarse tissue vocabulary, aligned with ``idmap/data/tissues.json`` so a
# signature's tissue label resolves to an UBERON term through the existing
# ``GeneResolver.map_tissue``.
#
# GDS carries tissue only as per-sample free text (``src: human frontal cortex``,
# but also ``src: vastus lateralis biopsy, 69 yr old male``), so the label has to
# be recognised rather than read. Order matters: the first entry whose keywords
# appear wins, so specific anatomy precedes the organ it sits in.
TISSUE_KEYWORDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("adipose tissue", ("adipocyte", "adipose", "epididymal fat", "white fat")),
    (
        "skeletal muscle",
        (
            "vastus lateralis",
            "gastrocnemius",
            "quadriceps",
            "biceps",
            "skeletal muscle",
            "muscle precursor",
            "myoblast",
            "satellite cell",
            "muscle fiber",
            "muscle fibre",
            "mhc i",
            "mhc iia",
            "muscle",
        ),
    ),
    (
        "brain",
        (
            "frontal cortex",
            "prefrontal",
            "neocortex",
            "hippocamp",
            "cerebell",
            "cerebral cortex",
            "ba10",
            "ba22",
            "temporal cortex",
            "brain",
            "cortex",
        ),
    ),
    (
        "bone marrow",
        ("bone marrow", "hematopoietic stem cell", "hsc", "lymphoid progenitor", "clp"),
    ),
    (
        "blood",
        (
            "whole blood",
            "peripheral blood",
            "pbmc",
            "lymphoblastoid",
            "lymphocyte",
            "monocyte",
            "dendritic cell",
            "leukocyte",
            "blood",
        ),
    ),
    ("thymus", ("thymus", "thymi", "thymocyte")),
    ("liver", ("liver", "hepat")),
    ("lung", ("lung", "pulmonary", "upper lobe", "lower lobe")),
    ("kidney", ("kidney", "renal")),
    ("heart", ("heart", "cardiac", "myocardium", "ventricle")),
    ("cochlea", ("cochlea",)),
    ("spleen", ("spleen", "splenic")),
    ("skin", ("skin", "dermis", "epidermis")),
)


def normalize_tissue(text: str) -> str | None:
    """Map free text to the coarse tissue vocabulary, or ``None`` if unrecognised."""
    haystack = " ".join(text.lower().split())
    for tissue, keywords in TISSUE_KEYWORDS:
        if any(word in haystack for word in keywords):
            return tissue
    return None


def _contrast_tissue(
    header: GdsHeader, samples: tuple[str, ...], stratum: str | None
) -> str | None:
    """Best available tissue label: the stratum, then sample sources, then the title.

    Falling back to the dataset title is safe here because it is only consulted
    when neither the subset level nor any sample source named a tissue.
    """
    if stratum:
        return normalize_tissue(stratum) or stratum
    sources = " ; ".join(header.sample_sources.get(s, "") for s in samples)
    return normalize_tissue(sources) or normalize_tissue(header.title)


def _format_age(years: float, species: str) -> str:
    return f"{years:.0f}y" if species == "human" else f"{years * 12:.0f}mo"


def _age_range_label(descriptions: dict[str, list[str]], species: str) -> str:
    parts = []
    for group in ("young", "old"):
        values = [v for v in (age_in_years(d) for d in descriptions[group]) if v is not None]
        if not values:
            parts.append(group)
            continue
        bounds = dict.fromkeys(_format_age(v, species) for v in (min(values), max(values)))
        parts.append("-".join(bounds))
    return f"{parts[0]} vs {parts[1]}"


# Emitted when the dataset declares no sex variable at all. Distinct from
# "both": one says the contrast mixes sexes, the other says GEO never recorded
# which sexes are in it, and a caller filtering on sex must not conflate them.
SEX_UNSPECIFIED = "unspecified"

_SEX_ALIASES = {
    "male": "male",
    "men": "male",
    "man": "male",
    "m": "male",
    "female": "female",
    "women": "female",
    "woman": "female",
    "f": "female",
}


def _sex_label(header: GdsHeader, samples: set[str]) -> str:
    levels = [s for variable in SEX_VARIABLES for s in header.by_variable(variable)]
    if not levels:
        return SEX_UNSPECIFIED
    matched = {
        _SEX_ALIASES.get(s.description.strip().lower())
        for s in levels
        if samples <= set(s.sample_ids)
    }
    labels = {m for m in matched if m is not None}
    return labels.pop() if len(labels) == 1 else "both"


def build_contrasts(header: GdsHeader) -> tuple[list[AgeContrast], list[str]]:
    """Derive every usable young-vs-old contrast from a GDS header.

    Returns ``(contrasts, skipped)``. ``skipped`` holds one human-readable reason
    per rejection, so a dataset that produces nothing says why rather than
    silently vanishing from the panel.
    """
    species = header.species
    if species is None:
        return [], [f"{header.accession}: organism {header.organism!r} is not modelled"]

    buckets: dict[str, set[str]] = {"young": set(), "old": set()}
    labels: dict[str, list[str]] = {"young": [], "old": []}
    for subset in header.by_variable("age"):
        group = age_group(subset.description, species)
        if group:
            buckets[group].update(subset.sample_ids)
            labels[group].append(subset.description)
    if not buckets["young"] or not buckets["old"]:
        return [], [f"{header.accession}: no young and old age levels for {species}"]
    # A sample assigned to both bands (overlapping curator levels) is unusable.
    ambiguous = buckets["young"] & buckets["old"]
    buckets["young"] -= ambiguous
    buckets["old"] -= ambiguous

    allowed: set[str] | None = None
    restrictions: list[str] = []
    for variable in header.variables:
        if variable not in CONTROLLED_VARIABLES:
            continue
        controls = [s for s in header.by_variable(variable) if is_control_level(s.description)]
        if not controls:
            levels = ", ".join(s.description for s in header.by_variable(variable))
            return [], [
                f"{header.accession}: {variable!r} has no recognisable control level ({levels})"
            ]
        control_ids = {sid for s in controls for sid in s.sample_ids}
        allowed = control_ids if allowed is None else (allowed & control_ids)
        restrictions.append(f"{variable} = {' | '.join(s.description for s in controls)}")

    strata: list[tuple[str | None, set[str] | None, str | None]] = [(None, None, None)]
    for variable in header.variables:
        if variable not in STRATIFIED_VARIABLES:
            continue
        stratum_levels = header.by_variable(variable)
        strata = [
            (
                level.description if base_label is None else f"{base_label} / {level.description}",
                set(level.sample_ids) if base_ids is None else base_ids & set(level.sample_ids),
                variable if base_var is None else f"{base_var}, {variable}",
            )
            for base_label, base_ids, base_var in strata
            for level in stratum_levels
        ]

    contrasts: list[AgeContrast] = []
    skipped: list[str] = []
    for label, ids, stratified_by in strata:
        young, old = buckets["young"], buckets["old"]
        for restrict in (allowed, ids):
            if restrict is not None:
                young = young & restrict
                old = old & restrict
        if len(young) < MIN_GROUP_SIZE or len(old) < MIN_GROUP_SIZE:
            where = f" [{label}]" if label else ""
            skipped.append(
                f"{header.accession}{where}: young={len(young)}, old={len(old)} "
                f"(need {MIN_GROUP_SIZE} each after restrictions)"
            )
            continue

        applied = list(restrictions)
        if stratified_by:
            applied.append(f"stratified by {stratified_by} = {label}")
        samples = young | old
        contrasts.append(
            AgeContrast(
                accession=header.accession,
                species=species,
                tissue=_contrast_tissue(header, tuple(sorted(samples)), label),
                sex=_sex_label(header, samples),
                age_range=_age_range_label(labels, species),
                young_samples=tuple(sorted(young)),
                old_samples=tuple(sorted(old)),
                restrictions=tuple(applied),
                stratum=label,
            )
        )
    return contrasts, skipped


# --- probe annotation -------------------------------------------------------

# GEO writes multi-gene probes as "1029///1030". A probe that measures more than
# one gene cannot be attributed to either, so it is dropped rather than assigned
# to the first id, which would fabricate evidence for whichever gene sorts first.
_MULTI = "///"

_UNINFORMATIVE = frozenset({"", "none", "null", "---", "nan"})


def probe_annotation(table: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    """Split a parsed GDS table into ``(probe -> gene key, key -> identifiers)``.

    The key is the Entrez id where the platform annotation supplies one and the
    gene symbol otherwise. Entrez is preferred because symbols get renamed
    between annotation releases while Entrez ids are stable, and the store joins
    signatures to curated knowledge on the resolved identifier.
    """
    symbols = table.get("IDENTIFIER")
    if symbols is None:
        symbols = table.get("Gene symbol")
    if symbols is None:
        raise SourceError(
            "GDS table has neither IDENTIFIER nor 'Gene symbol'.", detail=list(table.columns)[:10]
        )
    symbols = symbols.fillna("").astype(str).str.strip()
    entrez = (
        table["Gene ID"].fillna("").astype(str).str.strip()
        if "Gene ID" in table.columns
        else pd.Series("", index=table.index)
    )

    usable = (
        ~symbols.str.lower().isin(_UNINFORMATIVE)
        & ~symbols.str.contains(_MULTI, regex=False)
        & ~entrez.str.contains(_MULTI, regex=False)
    )
    symbols, entrez = symbols[usable], entrez[usable]
    keys = entrez.where(entrez.str.isdigit(), "SYM:" + symbols)

    identifiers = (
        pd.DataFrame(
            {"gene_key": keys, "symbol": symbols, "entrez": entrez.where(entrez.str.isdigit())}
        )
        .drop_duplicates("gene_key")
        .set_index("gene_key")
    )
    return keys, identifiers


class GeoDataSetsSource(SourceAdapter):
    """Adapter over the pinned GEO DataSets aging panel.

    The ETL (:mod:`geroquery.etl.build_signatures`) does the fetching, so this
    class exists for the provenance surface: ``/v1/sources`` should name every
    upstream the store's contents actually came from, and a source that is
    invisible there is a source a caller cannot audit.
    """

    name = "geo-datasets"

    def capabilities(self) -> Capabilities:
        from .manifest import GEO_AGING_PANEL, VERIFIED_ON

        return Capabilities(
            source_name=self.name,
            omics=("transcriptome",),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes=(
                f"{len(GEO_AGING_PANEL)} curated GEO DataSets with an age subset "
                f"variable, harmonized to young-vs-old Hedges' g per gene. "
                f"Checksums last verified {VERIFIED_ON}."
            ),
        )

    def license(self) -> License:
        return License(LICENSE, redistributable=True, attribution=ATTRIBUTION)
