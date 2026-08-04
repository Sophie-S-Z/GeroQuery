"""GEO DNA-methylation series adapter — real beta matrices for the clock layer.

GeroQuery wires 236 published aging clocks. Until this module existed, none of
them had run on real data: the repository contained no methylation matrix, so
"236 clocks" meant 236 objects that had been shown not to crash. This closes
that loop.

**Why GEO Series here, when the expression panel deliberately uses GEO DataSets.**
The expression side avoided free-text age parsing by taking the curator-built GDS
view, where age is a declared subset variable. No GDS exists for methylation —
GEO stopped curating them before 450K arrays became common — so the age metadata
has to be read out of ``!Sample_characteristics_ch1``. That is tractable here and
was not there, for one reason: this is **two hand-verified series**, not thirty.
Each one's characteristic keys are checked against
:data:`SERIES` at parse time, so a change in upstream formatting fails loudly
rather than yielding a plausible frame full of NaN ages.

**Why these two series.**

``GSE64495``
    450K, whole blood, n=113, ages 0–94. The primary validation set: a wide,
    continuous age range is what makes a mean absolute error meaningful. It also
    ships the authors' **own Horvath-clock output** per sample
    (``dna methylation age``), which is a far stronger check than MAE against
    chronological age — it compares our wrapper against a published number
    computed from the same input by the clock's own conventions.

``GSE30870``
    450K, newborns vs nonagenarians, n=40. An extreme-contrast set. A clock that
    cannot separate cord blood from 90-year-old blood is broken in a way no
    correlation coefficient would hide.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ..exceptions import SourceError

GSE_FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GSE_BROWSER = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="

LICENSE = "NCBI GEO — US public domain, attribute"

_TABLE_BEGIN = "!series_matrix_table_begin"
_TABLE_END = "!series_matrix_table_end"

# Illumina 450K. Both series use it; a different platform would mean different
# probe identifiers and is refused rather than silently accepted.
PLATFORM_450K = "GPL13534"


@dataclass(frozen=True)
class MethylationSeries:
    """A hand-verified GEO methylation series and how to read its metadata."""

    accession: str
    title: str
    platform: str
    n_samples: int
    tissue: str
    # Characteristic key -> the column name it becomes. Verified against the
    # actual file; a missing key raises rather than producing a NaN column.
    characteristics: dict[str, str]
    # Restrict to samples whose named column equals one of these values. Same
    # control-arm discipline the expression panel applies.
    control_filter: tuple[str, str] | None = None
    notes: str = ""
    reference_age_column: str | None = None

    @property
    def url(self) -> str:
        digits = self.accession[len("GSE") :]
        shard = f"GSE{digits[:-3]}nnn" if len(digits) > 3 else "GSEnnn"
        return (
            f"{GSE_FTP_ROOT}/{shard}/{self.accession}/matrix/"
            f"{self.accession}_series_matrix.txt.gz"
        )

    @property
    def browser_url(self) -> str:
        return f"{GSE_BROWSER}{self.accession}"


SERIES: dict[str, MethylationSeries] = {
    "GSE64495": MethylationSeries(
        accession="GSE64495",
        title=(
            "DNA methylation profiles of human blood samples "
            "(developmental disorder and controls)"
        ),
        platform=PLATFORM_450K,
        n_samples=113,
        tissue="whole blood",
        characteristics={
            "age": "age",
            "Sex": "sex",
            "tissue": "tissue",
            "disease status": "disease_status",
            "dna methylation age": "published_horvath_age",
        },
        control_filter=("disease_status", "Control"),
        reference_age_column="published_horvath_age",
        notes=(
            "Ages carry one decimal place. The authors' own Horvath-clock output "
            "ships in the metadata, so the wrapper can be checked against a "
            "published per-sample number rather than only against chronology."
        ),
    ),
    "GSE30870": MethylationSeries(
        accession="GSE30870",
        title="DNA methylomes of newborns and nonagenarians",
        platform=PLATFORM_450K,
        n_samples=40,
        tissue="whole blood / cord blood",
        characteristics={
            "age": "age_label",
            "disease status": "disease_status",
            "tissue": "tissue",
            "cell type": "cell_type",
        },
        control_filter=("disease_status", "healthy sample"),
        notes=(
            "Ages are qualitative for the newborn arm ('Newborn') and in years "
            "for the nonagenarians ('97 years'), so age is parsed from a label "
            "rather than read as a float."
        ),
    ),
}


@dataclass
class MethylationDataset:
    """A parsed series: beta matrix plus per-sample metadata."""

    series: MethylationSeries
    betas: pd.DataFrame  # CpG id x sample
    metadata: pd.DataFrame  # sample x (age, sex, ...)
    dropped: dict[str, float] = field(default_factory=dict)

    @property
    def accession(self) -> str:
        return self.series.accession

    def matrix(self) -> pd.DataFrame:
        """Samples as rows, CpGs as columns — the orientation the clocks take.

        Stored the other way round because that is how GEO ships it and how the
        data is naturally indexed, but every clock wrapper in this project
        expects one row per sample. Converting here rather than at each call
        site is what stops the two conventions being mixed up, which is exactly
        the mistake that made the first run of this panel fail 236 times.
        """
        return self.betas.T

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"<MethylationDataset {self.accession} "
            f"{self.betas.shape[0]} CpGs x {self.betas.shape[1]} samples>"
        )


# "Newborn" is age 0. Cord blood at birth is the intended meaning, and treating
# it as missing would discard half of GSE30870.
_NEWBORN = re.compile(r"^\s*(newborn|new born|cord blood|birth)\s*$", re.IGNORECASE)
_YEARS = re.compile(r"(-?\d+(?:\.\d+)?)\s*(?:years?|yrs?|y)?\s*$", re.IGNORECASE)


def parse_age(label: str) -> float | None:
    """Parse a GEO age characteristic into years.

    Returns ``None`` rather than guessing when the label is not an age — an
    unparsed age must not become a silent zero in a clock's error statistics.
    """
    if label is None:
        return None
    text = str(label).strip()
    if not text or text.upper() in {"NA", "N/A", "NONE", "--"}:
        return None
    if _NEWBORN.match(text):
        return 0.0
    match = _YEARS.match(text)
    return float(match.group(1)) if match else None


def _header_and_offset(path: Path) -> tuple[list[str], int]:
    """Return the metadata lines and the line index where the table starts."""
    header: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            if line.startswith(_TABLE_BEGIN):
                return header, index + 1
            header.append(line.rstrip("\n"))
    raise SourceError(
        f"{path.name} contains no {_TABLE_BEGIN!r} marker.", detail={"path": str(path)}
    )


def _split_values(line: str) -> list[str]:
    """Split a `!Key\\t"v1"\\t"v2"` line into its unquoted values."""
    return [part.strip().strip('"') for part in line.split("\t")[1:]]


def parse_metadata(header_lines: list[str], series: MethylationSeries) -> pd.DataFrame:
    """Build the per-sample metadata frame from a series-matrix header.

    Raises if a characteristic key declared in :data:`SERIES` is absent, so an
    upstream formatting change is a loud failure rather than a column of NaN.
    """
    accessions: list[str] = []
    platform: str | None = None
    characteristics: dict[str, list[str]] = {}

    for line in header_lines:
        if line.startswith("!Sample_geo_accession"):
            accessions = _split_values(line)
        elif line.startswith("!Series_platform_id"):
            platform = _split_values(line)[0] if _split_values(line) else None
        elif line.startswith("!Sample_characteristics_ch1"):
            values = _split_values(line)
            keys = {v.split(":", 1)[0].strip() for v in values if ":" in v}
            if len(keys) != 1:
                continue  # a mixed row is not a single characteristic
            key = keys.pop()
            characteristics[key] = [v.split(":", 1)[1].strip() if ":" in v else "" for v in values]

    if not accessions:
        raise SourceError(
            f"{series.accession}: no !Sample_geo_accession row.", detail=series.accession
        )
    if platform and platform != series.platform:
        raise SourceError(
            f"{series.accession}: expected platform {series.platform}, found {platform}. "
            f"Probe identifiers would not match the pinned expectation.",
            detail={"expected": series.platform, "found": platform},
        )

    frame = pd.DataFrame(index=pd.Index(accessions, name="sample_id"))
    missing = [k for k in series.characteristics if k not in characteristics]
    if missing:
        raise SourceError(
            f"{series.accession}: characteristics {missing} are no longer in the series "
            f"matrix. Available: {sorted(characteristics)}",
            detail={"missing": missing, "available": sorted(characteristics)},
        )
    for key, column in series.characteristics.items():
        frame[column] = characteristics[key]

    age_source = "age" if "age" in frame.columns else "age_label"
    frame["age"] = [parse_age(v) for v in frame[age_source]]
    if series.reference_age_column and series.reference_age_column in frame.columns:
        frame[series.reference_age_column] = pd.to_numeric(
            frame[series.reference_age_column], errors="coerce"
        )
    return frame


def read_series(
    path: Path,
    series: MethylationSeries,
    *,
    cpgs: list[str] | None = None,
    impute: bool = True,
) -> MethylationDataset:
    """Parse a pinned series matrix into betas + metadata.

    Args:
        path: the verified ``*_series_matrix.txt.gz``.
        series: its :class:`MethylationSeries` declaration.
        cpgs: optional CpG allow-list. Passing the union of the clocks' required
            probes turns a 485,000-row read into a few-thousand-row one, which is
            the difference between ~450 MB of resident beta matrix and ~5 MB.

    Returns:
        A :class:`MethylationDataset` restricted to the control arm and to
        samples with a parseable age.
    """
    header_lines, skiprows = _header_and_offset(path)
    metadata = parse_metadata(header_lines, series)

    wanted = set(cpgs) if cpgs else None
    frame = pd.read_csv(
        path,
        sep="\t",
        skiprows=skiprows,
        index_col=0,
        na_values=["null", "NA", ""],
        comment="!",  # drops the trailing !series_matrix_table_end row
        low_memory=False,
    )
    frame.index = frame.index.astype(str).str.strip().str.strip('"')
    frame.columns = [str(c).strip().strip('"') for c in frame.columns]
    if wanted:
        frame = frame.loc[frame.index.isin(wanted)]

    dropped: dict[str, float] = {"samples_total": len(metadata)}

    keep = metadata.index
    if series.control_filter:
        column, value = series.control_filter
        keep = metadata.index[metadata[column].str.strip().str.casefold() == value.casefold()]
        dropped["samples_not_control"] = len(metadata) - len(keep)
    metadata = metadata.loc[keep]

    unparsed = metadata["age"].isna()
    dropped["samples_no_age"] = int(unparsed.sum())
    metadata = metadata.loc[~unparsed]

    shared = [s for s in metadata.index if s in frame.columns]
    dropped["samples_absent_from_matrix"] = len(metadata) - len(shared)
    metadata = metadata.loc[shared]
    betas = frame[shared].astype("float32")

    before = len(betas)
    all_missing = betas.isna().all(axis=1)
    dropped["cpgs_all_missing"] = int(all_missing.sum())
    betas = betas.loc[~all_missing]

    if impute:
        # Mean-impute the remaining gaps, per CpG, across samples.
        #
        # Dropping every CpG with any missing value instead — which this did
        # first — silently deletes probes the clocks need: 7,510 of GSE30870's
        # CpGs have at least one NA, and three of them are Horvath coefficients.
        # The result was a clock that "requires features not present", i.e. a
        # data-preparation failure wearing the costume of an incompatible input.
        #
        # Imputing a handful of probes from the cohort mean is the standard
        # handling and is what the clock libraries do internally anyway. The
        # fraction is reported so it can never be large without being visible.
        gaps = betas.isna()
        n_gaps = int(gaps.to_numpy().sum())
        dropped["cpgs_with_any_gap"] = int(gaps.any(axis=1).sum())
        dropped["cells_imputed"] = n_gaps
        dropped["fraction_imputed"] = round(n_gaps / float(betas.size), 6) if betas.size else 0.0
        if n_gaps:
            # numpy, not `betas.T.fillna(means).T`: pandas aligns fillna per
            # column, and with ~480,000 CpG columns after the transpose that is
            # minutes rather than milliseconds.
            values = betas.to_numpy(dtype="float32", copy=True)
            rows = np.where(gaps.to_numpy())[0]
            means = np.nanmean(values, axis=1)
            values[gaps.to_numpy()] = means[rows]
            betas = pd.DataFrame(values, index=betas.index, columns=betas.columns)
    dropped["cpgs_kept"] = len(betas)
    dropped["cpgs_total"] = before

    return MethylationDataset(series=series, betas=betas, metadata=metadata, dropped=dropped)


def get_series(accession: str) -> MethylationSeries:
    try:
        return SERIES[accession]
    except KeyError:
        raise KeyError(
            f"Unknown methylation series {accession!r}. Known: {sorted(SERIES)}"
        ) from None
