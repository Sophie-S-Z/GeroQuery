"""NHANES 2017-2018 clinical adapter — real data, not a generator.

This is the source behind GeroQuery's clinical/resilience path. It downloads the
four pinned NHANES XPORT files (see :mod:`geroquery.sources.manifest`), joins them
on ``SEQN``, and emits the six-marker clinical frame the clocks and resilience
modules already expect.

Two operating modes, both real data:

``full``
    All four XPT files, fetched and checksum-verified. n = 4,895 complete cases
    aged 20+. This is what ``make data`` produces and what any reported result
    should be computed on.

``sample``
    A committed, seeded subset of the same real rows (``data/nhanes_2017_2018_sample.csv``)
    so tests and CI run offline without shipping 7 MB of XPORT into git. It is a
    subset of real measurements, not synthetic data — but it is a *sample*, so
    estimates computed on it are not the reported result.

Known limitations, carried in :data:`CAVEATS` rather than left implicit:

* ``RIDAGEYR`` is topcoded at 80, which compresses exactly the oldest stratum —
  the one where any critical-slowing-down signal should be strongest.
* NHANES is cross-sectional. Age strata are a *proxy* for within-individual
  trajectories; no relaxation time is actually observed.
* Survey weights are carried through as ``survey_weight`` but **not applied**.
  Estimates are unweighted and therefore not nationally representative.
* No fasting-status, medication, or comorbidity adjustment.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..config import SOURCES_DATA
from ..exceptions import SourceError
from .base import Capabilities, License, SourceAdapter
from .fetch import fetch_artifact
from .manifest import NHANES_2017_2018, VERIFIED_ON

RELEASE = "NHANES 2017-2018"

# Only adults: paediatric reference ranges differ enough that pooling them would
# make an age trend partly an artefact of changing physiology, not of aging.
MIN_AGE = 20.0

# CDC topcodes RIDAGEYR at 80 for disclosure control; everyone older is coded 80.
AGE_TOPCODE = 80.0

# geroquery column -> (manifest key, NHANES variable). These six markers are
# exactly clocks.registry.CLINICAL_FEATURES, so the clock path needs no change.
MARKER_MAP: dict[str, tuple[str, str]] = {
    "albumin": ("BIOPRO_J", "LBXSAL"),  # g/dL
    "creatinine": ("BIOPRO_J", "LBXSCR"),  # mg/dL
    "glucose": ("BIOPRO_J", "LBXSGL"),  # mg/dL
    "crp": ("HSCRP_J", "LBXHSCRP"),  # mg/L, high-sensitivity
    "lymphocyte_pct": ("CBC_J", "LBXLYPCT"),  # %
    "rdw": ("CBC_J", "LBXRDW"),  # %
}

MARKERS: tuple[str, ...] = tuple(MARKER_MAP)

DEMOGRAPHICS: dict[str, tuple[str, str]] = {
    "age": ("DEMO_J", "RIDAGEYR"),
    "_sex_code": ("DEMO_J", "RIAGENDR"),  # 1 = male, 2 = female
    "survey_weight": ("DEMO_J", "WTMEC2YR"),
}

SEX_LABELS = {1: "male", 2: "female"}

SAMPLE_FILENAME = "nhanes_2017_2018_sample.csv"

# Full harmonized table written by `make data`. Git-ignored: it is reproducible
# from the pinned manifest, so a committed copy would just be an unverifiable one.
FULL_FILENAME = "clinical_nhanes_full.csv"

# Seeded so the committed sample is reproducible from the full download.
SAMPLE_SEED = 20240117
SAMPLE_SIZE = 600

CAVEATS: tuple[str, ...] = (
    f"Age is topcoded at {AGE_TOPCODE:.0f}, compressing the oldest stratum.",
    "Cross-sectional: age strata proxy for within-individual trajectories; "
    "no relaxation time is observed.",
    "Survey weights (WTMEC2YR) are carried but NOT applied; estimates are unweighted "
    "and not nationally representative.",
    "No fasting-status, medication, or comorbidity adjustment.",
)


def _read_xpt(path: Path, columns: list[str]) -> pd.DataFrame:
    """Read the requested columns out of an NHANES XPORT file."""
    try:
        df = pd.read_sas(path, format="xport")
    except Exception as exc:  # pandas raises assorted ValueError/struct errors
        raise SourceError(
            f"Could not parse {path.name} as XPORT. The upstream URL may now serve "
            f"HTML instead of a data file.",
            detail={"path": str(path), "error": str(exc)},
        ) from exc
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise SourceError(
            f"{path.name} is missing expected NHANES variables {missing}.",
            detail={"path": str(path), "missing": missing, "available": list(df.columns)[:40]},
        )
    return df[["SEQN", *columns]]


def _columns_by_file() -> dict[str, list[str]]:
    """Invert the variable maps into {manifest key: [NHANES vars]}."""
    wanted: dict[str, list[str]] = {}
    for _, (file_key, var) in {**DEMOGRAPHICS, **MARKER_MAP}.items():
        wanted.setdefault(file_key, []).append(var)
    return wanted


def load_full(
    *,
    allow_network: bool | None = None,
    cache_directory: Path | None = None,
    min_age: float = MIN_AGE,
) -> pd.DataFrame:
    """Fetch, verify, join, and harmonize the full NHANES 2017-2018 clinical frame.

    Returns complete cases only — a subject missing any of the six markers is
    dropped rather than imputed, because imputing the covariance structure is
    precisely what the resilience metrics are trying to measure.
    """
    frames: dict[str, pd.DataFrame] = {}
    for file_key, variables in _columns_by_file().items():
        path = fetch_artifact(
            NHANES_2017_2018[file_key],
            directory=cache_directory,
            allow_network=allow_network,
        )
        frames[file_key] = _read_xpt(path, variables)

    merged = frames["DEMO_J"]
    for file_key, frame in frames.items():
        if file_key != "DEMO_J":
            merged = merged.merge(frame, on="SEQN", how="inner")

    rename = {var: col for col, (_, var) in {**DEMOGRAPHICS, **MARKER_MAP}.items()}
    out = merged.rename(columns=rename)

    out["subject_id"] = "NHANES:" + out["SEQN"].astype("int64").astype(str)
    out["sex"] = out["_sex_code"].map(SEX_LABELS)

    columns = ["subject_id", "age", "sex", *MARKERS, "survey_weight"]
    out = out[columns]
    out = out[out["age"] >= min_age]
    out = out.dropna(subset=["age", "sex", *MARKERS])
    return out.sort_values("subject_id").reset_index(drop=True)


def full_path(data_dir: Path | None = None) -> Path:
    """Path to the full harmonized real table written by ``make data``."""
    return (data_dir or SOURCES_DATA) / FULL_FILENAME


def sample_path(data_dir: Path | None = None) -> Path:
    """Path to the committed offline sample of real NHANES rows."""
    return (data_dir or SOURCES_DATA) / SAMPLE_FILENAME


def load_sample(data_dir: Path | None = None) -> pd.DataFrame:
    """Load the committed sample of real NHANES rows (no network required)."""
    path = sample_path(data_dir)
    if not path.exists():
        raise SourceError(
            f"NHANES offline sample missing at {path}. Regenerate it with "
            f"`python -m geroquery.etl.build_data --write-sample`.",
            detail={"path": str(path)},
        )
    return pd.read_csv(path)


def write_sample(
    frame: pd.DataFrame,
    data_dir: Path | None = None,
    *,
    n: int = SAMPLE_SIZE,
    seed: int = SAMPLE_SEED,
) -> Path:
    """Write a seeded subset of ``frame`` as the committed offline sample."""
    take = frame if len(frame) <= n else frame.sample(n=n, random_state=seed)
    path = sample_path(data_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    take.sort_values("subject_id").to_csv(path, index=False)
    return path


class NhanesClinicalSource(SourceAdapter):
    """Adapter over NHANES 2017-2018 clinical chemistry and CBC."""

    name = "nhanes-2017-2018"

    def __init__(self, data_dir: Path | None = None, cache_directory: Path | None = None):
        self.data_dir = data_dir or SOURCES_DATA
        self.cache_directory = cache_directory

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("clinical",),
            species=("human",),
            # Not federated: fetched once by ETL, then served from local storage.
            federated=False,
            cacheable=True,
            notes=(
                f"{RELEASE} clinical markers ({', '.join(MARKERS)}), joined on SEQN. "
                f"Checksums last verified {VERIFIED_ON}."
            ),
        )

    def license(self) -> License:
        artifact = NHANES_2017_2018["DEMO_J"]
        return License(artifact.license, redistributable=True, attribution=artifact.attribution)

    def clinical_frame(
        self, *, prefer_full: bool = True, allow_network: bool | None = None
    ) -> tuple[pd.DataFrame, str]:
        """Return ``(frame, mode)`` where mode is ``"full"`` or ``"sample"``.

        Resolution order: the table already built by ``make data``, then a fresh
        verified fetch, then the committed offline sample. The mode is *returned*
        rather than just logged so callers can record which one produced a number
        — an estimate from the 600-row sample must never be reported as the
        full-cohort result.
        """
        if prefer_full:
            built = full_path(self.data_dir)
            if built.exists():
                return pd.read_csv(built), "full"
            try:
                frame = load_full(allow_network=allow_network, cache_directory=self.cache_directory)
                return frame, "full"
            except SourceError:
                pass
        return load_sample(self.data_dir), "sample"
