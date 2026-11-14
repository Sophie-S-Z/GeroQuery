"""ARCHS4 — rule-selected young-vs-old RNA-seq contrasts from GEO.

The GEO DataSets panel is exhausted: NCBI stopped curating DataSets, all 189
that declare an ``age`` subset were audited, and the 158 the contrast rules
reject are developmental time courses, "old" arms that are middle-aged, disease
designs with no control arm, or groups under three. None hides a usable aging
contrast, so the panel cannot grow from that source without loosening a rule.

ARCHS4 (Lachmann et al., *Nat Commun* 2018) re-aligned every public human and
mouse RNA-seq sample in GEO through one pipeline. That removes the two limits of
the DataSets panel at once: it reaches the RNA-seq era, and every sample went
through the same processing, so cross-study heterogeneity is biology and design
rather than pipeline.

**Selection is a rule, not a list.** A sample enters when its GEO
characteristics carry a parseable age in the young or old band of
:mod:`geroquery.sources.geo` (the same bands, the same parser) and the profile
is bulk. Within a series, samples are grouped by stratum (tissue, cell type,
region) and by every other recorded attribute, so a young and an old arm differ
in age and in nothing else that was recorded. Attributes fall into three tiers:

* **neutral** — sex, strain, subject and batch identifiers, RNA quality, race,
  BMI: descriptors that vary between people or animals without being a
  manipulation. Pooled, as the DataSets rules pool ``individual`` and ``strain``.
* **strata** — tissue, organ, cell type, region: one contrast per level.
* **everything else** — treatment, genotype, diet, disease, time point, and any
  key not recognised — must sit at a control level (``wild type``,
  ``healthy``, ``vehicle``, ``baseline`` ...), or the sample is excluded. An
  unknown key excludes rather than passes: ``perturbation: treadmill`` is not
  in any list, and that is exactly why it must not be assumed harmless.

Groups need :data:`geo.MIN_GROUP_SIZE` per arm, as in the DataSets panel.

**Pinning.** ARCHS4 publishes only a moving ``latest`` file (64 GB human, 51 GB
mouse), so a whole-file checksum is neither computable here nor stable. The pin
is the release string inside the file plus its exact byte size; a moved upstream
fails the build rather than silently changing published numbers. Because the
file will be replaced, the selected sample lists are committed
(``archs4_selection.csv``): they are the durable record of what was measured.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from ..exceptions import SourceError
from . import geo
from .base import Capabilities, License, SourceAdapter

ROOT = "https://s3.dev.maayanlab.cloud/archs4/files"
VERSION = "2.7"
RELEASED = "2026-11-06"
# species -> (file name, exact byte size of the pinned release)
FILES: dict[str, tuple[str, int]] = {
    "human": ("human_gene_v2.latest.h5", 64_432_723_113),
    "mouse": ("mouse_gene_v2.latest.h5", 51_458_003_701),
}

LICENSE = (
    "ARCHS4 - uniformly processed counts of public GEO/SRA samples; the licence of "
    "the processed files is not stated (confirm with the Ma'ayan Laboratory before "
    "commercial use)"
)
ATTRIBUTION = (
    "Lachmann A, Torre D, Keenan AB, et al. Massive mining of publicly available RNA-seq "
    "data from human and mouse. Nat Commun. 2018;9:1366. https://archs4.org/"
)

# A profile ARCHS4 scores as more likely single-cell than bulk is excluded: its
# count distribution is a different measurement from the bulk arrays and
# RNA-seq the rest of the panel pools.
MAX_SINGLE_CELL_PROBABILITY = 0.5

# A gene is analysed in a contrast when at least this many CPM are seen in as
# many samples as the smaller arm has (edgeR's filterByExpr idea). Below that a
# standardized difference measures counting noise, not expression.
MIN_CPM = 1.0

_AGE_KEY = re.compile(r"^age(?:\s*[(_]?\s*(?:in\s+)?(years?|months?|weeks?|days?)\s*\)?)?$")
_NEUTRAL_KEY = re.compile(
    r"^(?:sex|gender|genotyped sex|race|ethnicity|bmi|rin|pmi"
    r"|strain|mouse strain|strain background|background strain|strain/background"
    r"|background|genetic background|dev stage|developmental stage"
    r"|individual|subject|donor|patient|participant|mouse|animal"
    r"|batch|enrollment batch|lane|plate|plate position|plate_location|well address"
    r"|flowcell|barcode|block|sequencing pool|library construction|spike-in"
    r"|replicate|biological replicate|lab|site|provider|link|biomaterial_type"
    r"|biosample encode accession)$"
    r"|(?:^|[\s_])id(?: number)?$"
)
_STRATUM_KEYS = frozenset(
    {
        "tissue",
        "tissue type",
        "organ",
        "cell type",
        "celltype",
        "cell_type",
        "region",
        "brain region",
    }
)
_SEX = {"male": "male", "m": "male", "female": "female", "f": "female"}

# Tumour material is sometimes recorded as a cell type or a source name rather
# than as a disease state. A series that profiles any tumour is a cancer study,
# and its "normal" samples are tumour-adjacent tissue from patients, so the
# whole series is excluded, not just the tumour samples.
_TUMOUR = re.compile(
    r"tumou?r|cancer|carcinoma|adenoma|sarcoma|lymphoma|leuka?emia|melanoma|glioma|metasta"
)


def characteristics(text: str) -> dict[str, str]:
    """``"age: 3 months,tissue: liver"`` -> ``{"age": "3 months", "tissue": "liver"}``."""
    out: dict[str, str] = {}
    for part in text.split(","):
        key, sep, value = part.partition(":")
        if sep and key.strip():
            out[key.strip().lower()] = value.strip()
    return out


def sample_age(fields: dict[str, str], species: str) -> str | None:
    """The age description, with a unit carried over from the key when the value lacks one.

    A bare number is read as years for human samples only, where it is the
    universal convention ("age: 67"). For mouse it stays unreadable: weeks and
    months are both common, and guessing would put a 24-week animal in the old
    band.
    """
    for key, value in fields.items():
        match = _AGE_KEY.match(key)
        if match:
            unit = match.group(1) or ("years" if species == "human" else None)
            return f"{value} {unit}" if unit and not re.search(r"[a-z]", value.lower()) else value
    return None


@dataclass(frozen=True)
class SeriesContrast:
    """A young-vs-old split within one GEO series, ready to estimate."""

    series: str
    species: str
    tissue: str | None
    stratum: str | None
    sex: str
    age_range: str
    young_samples: tuple[str, ...]
    old_samples: tuple[str, ...]

    @property
    def study_id(self) -> str:
        base = f"ARCHS4:{self.series}"
        return f"{base}:{geo.slug(self.stratum)}" if self.stratum else base


def _split_by_source(young: list, old: list) -> list[tuple[str | None, list, list]]:
    """Split an arm pair by GEO source name when it is a covariate both arms share.

    Cell type is often recorded only in the source name ("SST Neurons", "PV
    Neurons"); pooling those would put between-cell-type variance into the
    within-group SD and bury the age effect. A source name is a stratum only
    when the young and old arms carry exactly the same set of them: one that
    differs by arm encodes age itself ("old_liver"), and one unique per sample is
    an identifier.
    """
    young_sources = {s for *_, s in young}
    if len(young_sources) > 1 and young_sources == {s for *_, s in old}:
        return [
            (src, [x for x in young if x[3] == src], [x for x in old if x[3] == src])
            for src in sorted(young_sources)
        ]
    return [(None, young, old)]


def select_contrasts(
    meta: pd.DataFrame, species: str, *, exclude_series: Iterable[str] = ()
) -> list[SeriesContrast]:
    """Apply the selection rule to ARCHS4 sample metadata for one species.

    ``meta`` needs ``geo_accession``, ``series_id``, ``characteristics_ch1``,
    ``source_name_ch1`` and ``singlecellprobability``. Series in
    ``exclude_series`` (the DataSets panel's) are skipped, so one experiment is
    never pooled twice.
    """
    excluded = set(exclude_series)
    for series, text, source in zip(
        meta["series_id"], meta["characteristics_ch1"], meta["source_name_ch1"], strict=True
    ):
        if _TUMOUR.search(f"{text} {source}".lower()):
            excluded.update(series.split(","))
    groups: dict[tuple[str, tuple, tuple], dict[str, list]] = {}
    for accession, series, text, source, single_cell in zip(
        meta["geo_accession"],
        meta["series_id"],
        meta["characteristics_ch1"],
        meta["source_name_ch1"],
        meta["singlecellprobability"],
        strict=True,
    ):
        if single_cell >= MAX_SINGLE_CELL_PROBABILITY or excluded & set(series.split(",")):
            continue
        fields = characteristics(text)
        age = sample_age(fields, species)
        group = geo.age_group(age, species) if age else None
        if group is None:
            continue
        strata, held = [], []
        for key, value in sorted(fields.items()):
            if key.startswith("age") or _NEUTRAL_KEY.search(key):
                continue
            if key in _STRATUM_KEYS:
                strata.append((key, value.lower()))
            elif geo.is_control_level(value):
                held.append((key, value.lower()))
            else:
                break  # a manipulated or unrecognised attribute: not a control sample
        else:
            # Keyed on the primary series: a sample also listed under a
            # SuperSeries belongs to the same experiment as its siblings.
            primary = series.split(",")[0]
            entry = groups.setdefault(
                (primary, tuple(strata), tuple(held)), {"young": [], "old": []}
            )
            entry[group].append((accession, age, fields, source))

    contrasts = []
    for (series, stratum_fields, held_fields), arms in sorted(groups.items()):
        for source, young, old in _split_by_source(arms["young"], arms["old"]):
            if len(young) < geo.MIN_GROUP_SIZE or len(old) < geo.MIN_GROUP_SIZE:
                continue
            samples = young + old
            sexes = [
                _SEX.get((f.get("sex") or f.get("gender") or "").lower()) for _, _, f, _ in samples
            ]
            recorded = sorted({s for s in sexes if s})
            parts = [
                *stratum_fields,
                *held_fields,
                *([("source", source.lower())] if source else []),
            ]
            label = " / ".join(f"{k}={v}" for k, v in parts) or None
            contrasts.append(
                SeriesContrast(
                    series=series,
                    species=species,
                    tissue=geo.normalize_tissue(
                        " ".join([v for _, v in stratum_fields] + [samples[0][3]])
                    ),
                    stratum=label,
                    # "both" only when every sample records its sex; a partly recorded
                    # contrast cannot claim to be single-sex or mixed.
                    sex=(
                        geo.SEX_UNSPECIFIED
                        if None in sexes
                        else recorded[0] if len(recorded) == 1 else "both"
                    ),
                    age_range=geo.age_range_label(
                        {"young": [a for _, a, _, _ in young], "old": [a for _, a, _, _ in old]},
                        species,
                    ),
                    young_samples=tuple(sorted(a for a, _, _, _ in young)),
                    old_samples=tuple(sorted(a for a, _, _, _ in old)),
                )
            )
    # Distinct strata can share a slug ("il27ra+ hsc", "il27ra- hsc"); number
    # the repeats so every contrast keeps its own study id.
    seen: dict[str, int] = {}
    for i, contrast in enumerate(contrasts):
        n = seen[contrast.study_id] = seen.get(contrast.study_id, 0) + 1
        if n > 1:
            contrasts[i] = replace(contrast, stratum=f"{contrast.stratum} #{n}")
    return contrasts


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts -> log2(CPM + 1), the scale the Hedges' g is computed on."""
    library = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts / library * 1e6 + 1.0)


def expressed(counts: pd.DataFrame, min_samples: int) -> pd.Series:
    """Genes with at least :data:`MIN_CPM` in ``min_samples`` samples."""
    cpm = counts / counts.sum(axis=0).replace(0, np.nan) * 1e6
    return (cpm >= MIN_CPM).sum(axis=1) >= min_samples


# --- remote access ----------------------------------------------------------


def open_release(species: str):
    """Open the pinned ARCHS4 file over HTTP range requests, checking the pin.

    Needs the optional ``h5py`` and ``fsspec[http]`` packages. Range reads keep
    a 64 GB file out of local disk: only metadata columns and the selected
    samples' count columns are transferred.
    """
    try:
        import fsspec
        import h5py
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SourceError(
            "ARCHS4 needs the optional 'h5py' and 'fsspec[http]' packages: "
            "pip install -e '.[archs4]'"
        ) from exc

    name, size = FILES[species]
    url = f"{ROOT}/{name}"
    fs = fsspec.filesystem("http", skip_instance_cache=True)  # one session per process
    remote_size = fs.size(url)
    handle = h5py.File(fs.open(url, "rb", block_size=4 * 2**20, cache_type="blockcache"), "r")
    version = handle["meta/info/version"][()].decode()
    if (version, remote_size) != (VERSION, size):
        raise SourceError(
            f"ARCHS4 {species} moved from release {VERSION} ({size} bytes) to {version} "
            f"({remote_size} bytes). Re-run the selection and update the pin in "
            "geroquery/sources/archs4.py; published numbers must not change silently.",
            detail={"url": url, "pinned": [VERSION, size], "found": [version, remote_size]},
        )
    return handle


_META_COLUMNS = (
    "geo_accession",
    "series_id",
    "characteristics_ch1",
    "source_name_ch1",
    "singlecellprobability",
)


def read_metadata(handle) -> pd.DataFrame:
    """The sample columns the selection rule reads, in file order."""
    out = {}
    for column in _META_COLUMNS:
        values = handle[f"meta/samples/{column}"][:]
        out[column] = [v.decode("utf-8", "replace") if isinstance(v, bytes) else v for v in values]
    return pd.DataFrame(out)


def read_counts(handle, columns: list[int], accessions: list[str]) -> pd.DataFrame:
    """Gene-by-sample counts for the given file columns, indexed by Ensembl id.

    Columns are read in file order: neighbouring samples share storage blocks,
    so one series costs roughly one remote read instead of one per sample.
    """
    order = np.argsort(columns)
    data = handle["data/expression"][:, np.asarray(columns)[order]]
    genes = [g.decode() for g in handle["meta/genes/ensembl_gene"][:]]
    frame = pd.DataFrame(data, index=pd.Index(genes, name="ensembl"))
    frame.columns = np.asarray(accessions)[order]
    return frame[accessions]


def gene_keys(handle) -> pd.Series:
    """Ensembl id -> ``SYM:<symbol>``, the key the shared resolver canonicalizes.

    Several Ensembl ids can share a symbol; passing this map to
    :func:`~geroquery.harmonize.differential.contrast_effects` collapses them by
    highest mean expression (MaxMean, as for array probes), never by
    significance, which would select on the outcome.
    """
    genes = [g.decode() for g in handle["meta/genes/ensembl_gene"][:]]
    symbols = [s.decode().strip() for s in handle["meta/genes/symbol"][:]]
    return pd.Series([f"SYM:{s}" if s else None for s in symbols], index=genes)


class Archs4Source(SourceAdapter):
    """Provenance surface for the ARCHS4 contrasts (the ETL does the fetching)."""

    name = "archs4"

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("transcriptome",),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes=(
                f"ARCHS4 release {VERSION} ({RELEASED}): bulk RNA-seq samples from GEO series "
                "with a recorded age, rule-selected into young-vs-old contrasts in which "
                "every other recorded attribute is held at a control level."
            ),
        )

    def license(self) -> License:
        return License(LICENSE, redistributable=True, attribution=ATTRIBUTION)
