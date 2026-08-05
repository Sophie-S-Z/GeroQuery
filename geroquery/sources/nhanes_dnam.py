"""NHANES 1999-2002 cross-layer adapter: DNAm clocks + health state + mortality.

This is the only place in GeroQuery where all three layers land on the same
person. Everything else in the repo validates a measurement against another
measurement; this joins an epigenetic clock, a six-marker health state, and a
death certificate on ``SEQN``, and is therefore the first path here that can be
wrong in a way an outcome would notice.

Three upstreams, all pinned:

``dnmepi.sas7bdat``
    NCHS-computed DNAm biomarkers, released 2024-07-31. 4,449 rows of which
    **2,532 carry measurements**: adults aged 50+ from the 1999-2000 and
    2001-2002 cycles, whole blood. Twelve clocks plus cell-type fractions.
    NCHS ran the clocks; we read the outputs. The CpG betas are RDC-only.

1999-2002 demographics, biochemistry, CBC and CRP
    The same six markers the resilience module measures, plus the three extra
    ones PhenoAge needs, from the cycles the methylation subsample was drawn
    from.

NCHS 2019 public-use linked mortality files
    Follow-up through 2019-12-31, fixed-width, joined on ``SEQN``.

Why this cohort and not NHANES 2017-2018: 2017-2018 has the better clinical
panel and a wider age range, and no methylation at all. These two cycles are the
only NHANES with both.

What it costs, carried in :data:`CAVEATS` rather than left implicit:

* **Age floor of 50.** The methylation subsample is adults 50+, so the age range
  is 50-85 rather than 20-80. The critical-slowing-down design measures variance
  across age strata, and roughly half its dynamic range is gone. Any CSD result
  computed here is *less* powered than the published 2017-2018 one, not more.
* **CRP is a different assay.** 1999-2002 measured standard CRP in mg/dL;
  2017-2018 measured high-sensitivity CRP in mg/L. The values are converted to
  mg/L here so the column means the same thing, but the assays have different
  detection limits at the low end, and a standard-CRP cohort is not
  interchangeable with an hs-CRP one.
* **Public-use mortality records are perturbed.** NCHS substitutes synthetic
  follow-up time or cause of death for a subset of records to prevent
  re-identification. Aggregate estimates survive this; individual rows do not.
* **Variable names differ between the two cycles.** Creatinine is ``LBXSCR`` in
  1999-2000 and ``LBDSCR`` in 2001-2002; alkaline phosphatase ``LBXSAPSI`` then
  ``LBDSAPSI``. Getting this wrong silently halves the cohort instead of
  raising, so :data:`CYCLES` is pinned by a test.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..config import SOURCES_DATA
from ..exceptions import SourceError
from .base import Capabilities, License, SourceAdapter
from .fetch import fetch_artifact
from .manifest import NHANES_1999_2002, NHANES_DNAM, NHANES_MORTALITY, VERIFIED_ON

RELEASE = "NHANES 1999-2002 (DNAm subsample)"

# The methylation subsample is adults 50+. Not a choice we made; stating it so
# nobody reads a truncated age range as a cohort definition.
MIN_AGE = 50.0

# RIDAGEYR is topcoded at 85 in these cycles (80 in 2017-2018).
AGE_TOPCODE = 85.0

# mg/dL -> mg/L. 1999-2002 reports standard CRP in mg/dL; the rest of GeroQuery
# carries CRP in mg/L because that is what hs-CRP is reported in.
CRP_MG_DL_TO_MG_L = 10.0


@dataclass(frozen=True)
class Cycle:
    """One NHANES two-year cycle's file keys and variable names."""

    name: str
    demo: str
    bio: str
    cbc: str
    crp: str
    # geroquery column -> (file attribute above, NHANES variable)
    markers: dict[str, tuple[str, str]]


# Shared between cycles. Only creatinine and ALP move, and they are overridden
# per cycle below.
_COMMON_MARKERS: dict[str, tuple[str, str]] = {
    "albumin": ("bio", "LBXSAL"),  # g/dL
    "creatinine": ("bio", "LBXSCR"),  # mg/dL
    "glucose": ("bio", "LBXSGL"),  # mg/dL
    "crp": ("crp", "LBXCRP"),  # mg/dL -> converted to mg/L
    "lymphocyte_pct": ("cbc", "LBXLYPCT"),  # %
    "rdw": ("cbc", "LBXRDW"),  # %
    "mcv": ("cbc", "LBXMCVSI"),  # fL
    "wbc": ("cbc", "LBXWBCSI"),  # 10^3/uL
    "alp": ("bio", "LBXSAPSI"),  # IU/L
}

CYCLES: tuple[Cycle, ...] = (
    Cycle(
        name="1999-2000",
        demo="DEMO_1999",
        bio="LAB18_1999",
        cbc="LAB25_1999",
        crp="LAB11_1999",
        markers=dict(_COMMON_MARKERS),
    ),
    Cycle(
        name="2001-2002",
        demo="DEMO_2001",
        bio="L40_2001",
        cbc="L25_2001",
        crp="L11_2001",
        # The two renames. LBXSCR/LBXSAPSI simply do not exist in L40_B, so
        # using the 1999-2000 names here drops every 2001-2002 subject.
        markers={
            **_COMMON_MARKERS,
            "creatinine": ("bio", "LBDSCR"),
            "alp": ("bio", "LBDSAPSI"),
        },
    ),
)

# The six-marker health state, identical to nhanes.MARKERS so the resilience
# module sees the same state variable it was validated on. Do not widen: see the
# warning on nhanes.MARKER_MAP.
MARKERS: tuple[str, ...] = (
    "albumin",
    "creatinine",
    "glucose",
    "crp",
    "lymphocyte_pct",
    "rdw",
)

PHENOAGE_ONLY_MARKERS: tuple[str, ...] = ("mcv", "wbc", "alp")
PHENOAGE_MARKERS: tuple[str, ...] = MARKERS + PHENOAGE_ONLY_MARKERS

DEMOGRAPHIC_VARS: dict[str, str] = {
    "age": "RIDAGEYR",
    "_sex_code": "RIAGENDR",
    # Four-year MEC weight: correct for a 1999-2002 pooled analysis. The DNAm
    # subsample has its own weight (WTDN4YR) carried off the methylation file.
    "survey_weight_mec": "WTMEC4YR",
    "sdmvpsu": "SDMVPSU",
    "sdmvstra": "SDMVSTRA",
}

SEX_LABELS = {1: "male", 2: "female"}


@dataclass(frozen=True)
class ClockColumn:
    """One DNAm biomarker column, and what it actually predicts.

    ``predicts`` matters: a clock trained on chronological age and a clock
    trained on mortality are not interchangeable, and calling both "age" is how
    an age-acceleration residual gets computed against the wrong target.
    """

    column: str
    clock_id: str
    predicts: str
    units: str
    note: str


# The twelve clocks NCHS released. `predicts` drives whether age acceleration is
# even a meaningful quantity for a given column — for DunedinPoAm (a rate) and
# HorvathTelo (a length) it is not, so they are excluded from the residual path.
CLOCKS: tuple[ClockColumn, ...] = (
    ClockColumn(
        "HorvathAge",
        "nhanes-horvath",
        "chronological_age",
        "years",
        "Horvath 2013 multi-tissue clock.",
    ),
    ClockColumn(
        "HannumAge", "nhanes-hannum", "chronological_age", "years", "Hannum 2013 whole-blood clock."
    ),
    ClockColumn(
        "SkinBloodAge",
        "nhanes-skinblood",
        "chronological_age",
        "years",
        "Horvath 2018 skin-and-blood clock.",
    ),
    ClockColumn(
        "ZhangAge",
        "nhanes-zhang",
        "chronological_age",
        "years",
        "Zhang 2019 elastic-net age predictor.",
    ),
    ClockColumn("LinAge", "nhanes-lin", "chronological_age", "years", "Lin 2016 age predictor."),
    ClockColumn(
        "WeidnerAge",
        "nhanes-weidner",
        "chronological_age",
        "years",
        "Weidner 2014 three-CpG age predictor.",
    ),
    ClockColumn(
        "VidalBraloAge",
        "nhanes-vidalbralo",
        "chronological_age",
        "years",
        "Vidal-Bralo 2016 age predictor.",
    ),
    ClockColumn(
        "PhenoAge",
        "nhanes-phenoage",
        "phenotypic_age",
        "years",
        "Levine 2018 DNAm PhenoAge. Trained on a clinical phenotype, not on age.",
    ),
    ClockColumn(
        "GrimAgeMort",
        "nhanes-grimage",
        "mortality",
        "years",
        "Lu 2019 DNAm GrimAge. Trained on time-to-death.",
    ),
    ClockColumn(
        "GrimAge2Mort", "nhanes-grimage2", "mortality", "years", "Lu 2022 DNAm GrimAge version 2."
    ),
    ClockColumn(
        "DunedinPoAm",
        "nhanes-dunedinpoam",
        "pace_of_aging",
        "ratio",
        "Belsky 2020 pace of aging. A rate, not an age: 1.0 is one year per year.",
    ),
    ClockColumn(
        "HorvathTelo",
        "nhanes-dnamtl",
        "telomere_length",
        "kb",
        "Lu 2019 DNAm telomere length. Shortens with age; the sign is the check.",
    ),
)

CLOCK_COLUMNS: tuple[str, ...] = tuple(c.column for c in CLOCKS)

# Clocks for which "age acceleration" (clock minus chronological age) is a
# meaningful residual. Excludes the pace-of-aging rate and the telomere length,
# whose units are not years.
AGE_LIKE_CLOCKS: tuple[str, ...] = tuple(c.column for c in CLOCKS if c.units == "years")

# Estimated leukocyte fractions, released alongside the clocks. Carried because
# cell composition confounds every blood clock, and a result that does not
# adjust for it should at least be able to show that it could have.
CELL_FRACTIONS: tuple[str, ...] = (
    "CD8TPP",
    "CD4TPP",
    "Nkcell",
    "Bcell",
    "MonoPP",
    "NeuPP",
)

# DNAm-predicted sex. A free QC check: it should agree with reported sex on
# essentially every sample, and where it does not, the sample was probably
# swapped or mislabelled upstream.
#
# **The coding is the inverse of NHANES's.** RIAGENDR is 1 = male, 2 = female;
# XY_Estimation is 1 = female, 2 = male. Reusing SEX_LABELS here produced 2.4%
# agreement instead of 97.6% — that is, it was wrong for every subject, and it
# failed silently because both codings are valid-looking 1/2 integers. Verified
# by crosstab against RIAGENDR on all 2,532 assayed rows:
#
#     XY_Estimation=1 -> 1217 female, 30 male
#     XY_Estimation=2 -> 1255 male,   30 female
#
# The 60 discordant samples are what the check is for; they are counted by
# :func:`sex_discordance` rather than dropped, because a 2.4% mismatch rate is
# a property of the upstream assay that a reader should see.
DNAM_SEX_COLUMN = "XY_Estimation"

DNAM_SEX_LABELS = {1: "female", 2: "male"}

# Four-year survey weight specific to the DNAm subsample.
DNAM_WEIGHT_COLUMN = "WTDN4YR"

# --- mortality file layout ---------------------------------------------------
#
# Fixed-width, no header, one record per survey participant. Positions are
# 0-indexed half-open slices of a line padded to 48 characters — lines in the
# published files are 46-48 characters because trailing blanks are not padded,
# so slicing an unpadded line silently truncates PERMTH_EXM.
#
# Columns 21-42 hold the NHIS-only weight fields and are blank for NHANES.
MORTALITY_LINE_WIDTH = 48

MORTALITY_LAYOUT: tuple[tuple[str, int, int], ...] = (
    ("SEQN", 0, 6),
    ("eligstat", 14, 15),
    ("mortstat", 15, 16),
    ("ucod_leading", 16, 19),
    ("diabetes_flag", 19, 20),
    ("hypertension_flag", 20, 21),
    ("permth_int", 42, 45),
    ("permth_exm", 45, 48),
)

# ELIGSTAT: 1 = eligible for linkage, 2 = under 18 / not available, 3 = ineligible.
ELIGIBLE = 1

# MORTSTAT: 0 = assumed alive at end of follow-up, 1 = assumed deceased.
DECEASED = 1

# NCHS leading-cause groupings, for cause-specific work later. 001-010 only.
UCOD_LEADING_LABELS: dict[str, str] = {
    "001": "Diseases of heart",
    "002": "Malignant neoplasms",
    "003": "Chronic lower respiratory diseases",
    "004": "Accidents (unintentional injuries)",
    "005": "Cerebrovascular diseases",
    "006": "Alzheimer's disease",
    "007": "Diabetes mellitus",
    "008": "Influenza and pneumonia",
    "009": "Nephritis, nephrotic syndrome and nephrosis",
    "010": "All other causes",
}

SAMPLE_FILENAME = "nhanes_dnam_1999_2002_sample.csv"
FULL_FILENAME = "crosslayer_nhanes_dnam_full.csv"

SAMPLE_SEED = 20260805
SAMPLE_SIZE = 600

CAVEATS: tuple[str, ...] = (
    f"Age floor of {MIN_AGE:.0f}: the DNAm subsample is adults 50+, so the age range is "
    f"50-{AGE_TOPCODE:.0f} rather than the 20-80 of the 2017-2018 cohort. The "
    "critical-slowing-down design loses roughly half its dynamic range; any CSD "
    "estimate here is less powered than the published one, not more.",
    f"Age is topcoded at {AGE_TOPCODE:.0f}.",
    "CRP is the standard assay in mg/dL (converted to mg/L here), not the "
    "high-sensitivity assay used in 2017-2018. Detection limits differ at the low end.",
    "Public-use mortality records are perturbed by NCHS: a subset carry synthetic "
    "follow-up time or cause of death. Aggregate estimates are designed to survive "
    "this; individual records are not trustworthy on their own.",
    "Survey weights (WTMEC4YR, WTDN4YR) and the design variables (SDMVPSU, SDMVSTRA) "
    "are carried but NOT applied by default; unweighted estimates are not nationally "
    "representative.",
    "Clocks are NCHS's own computed values, not recomputed here. That removes our "
    "normalization from the picture and also means we cannot audit it.",
    "Cross-sectional exposure: markers and methylation are measured once, at baseline. "
    "Only the outcome is longitudinal.",
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
            f"{path.name} is missing expected NHANES variables {missing}. "
            f"Variable names are not stable across NHANES cycles - check the "
            f"per-cycle map in sources.nhanes_dnam.CYCLES.",
            detail={"path": str(path), "missing": missing, "available": list(df.columns)[:40]},
        )
    return df[["SEQN", *columns]]


def parse_mortality(text: str) -> pd.DataFrame:
    """Parse one NCHS public-use linked mortality file.

    Takes text rather than a path so the layout is testable against a literal
    without a fixture file. Lines are padded before slicing because the
    published files do not pad trailing blanks, and an unpadded 46-character
    line would yield a one-character ``permth_exm``.
    """
    rows = [line.ljust(MORTALITY_LINE_WIDTH) for line in text.splitlines() if line.strip()]
    frame = pd.DataFrame(
        {name: [row[start:stop].strip() for row in rows] for name, start, stop in MORTALITY_LAYOUT}
    )
    # NCHS writes missing as "." (SAS) or blank depending on the field.
    frame = frame.replace({"": None, ".": None})
    for column in frame.columns:
        if column != "ucod_leading":
            frame[column] = pd.to_numeric(frame[column])
    frame["SEQN"] = frame["SEQN"].astype("int64")
    return frame


def load_mortality(
    *, allow_network: bool | None = None, cache_directory: Path | None = None
) -> pd.DataFrame:
    """Both cycles' mortality follow-up, concatenated."""
    frames = []
    for key in NHANES_MORTALITY:
        path = fetch_artifact(
            NHANES_MORTALITY[key], directory=cache_directory, allow_network=allow_network
        )
        frames.append(parse_mortality(path.read_text(encoding="latin-1")))
    return pd.concat(frames, ignore_index=True)


def load_dnam(
    *, allow_network: bool | None = None, cache_directory: Path | None = None
) -> pd.DataFrame:
    """NCHS-computed DNAm biomarkers, measured rows only.

    The published file carries a row for every participant in the sampling
    frame, most of them empty. Rows without a Horvath age were never assayed;
    keeping them would turn "not measured" into "missing at random" three joins
    downstream.
    """
    path = fetch_artifact(NHANES_DNAM, directory=cache_directory, allow_network=allow_network)
    try:
        frame = pd.read_sas(path, format="sas7bdat")
    except Exception as exc:
        raise SourceError(
            f"Could not parse {path.name} as SAS7BDAT.",
            detail={"path": str(path), "error": str(exc)},
        ) from exc

    wanted = [
        "SEQN",
        DNAM_SEX_COLUMN,
        DNAM_WEIGHT_COLUMN,
        *CLOCK_COLUMNS,
        *CELL_FRACTIONS,
    ]
    missing = [c for c in wanted if c not in frame.columns]
    if missing:
        raise SourceError(
            f"dnmepi is missing expected columns {missing}.",
            detail={"missing": missing, "available": list(frame.columns)},
        )
    frame = frame[wanted].dropna(subset=["HorvathAge"]).copy()
    frame["SEQN"] = frame["SEQN"].astype("int64")
    return frame.reset_index(drop=True)


def load_clinical(
    *, allow_network: bool | None = None, cache_directory: Path | None = None
) -> pd.DataFrame:
    """Demographics and the nine markers for both cycles, harmonized.

    Each cycle is assembled with its own variable map and only then concatenated,
    because the names differ and a single map would drop one cycle entirely.
    """
    per_cycle = []
    for cycle in CYCLES:
        files: dict[str, pd.DataFrame] = {}
        needed: dict[str, list[str]] = {}
        for _column, (attr, variable) in cycle.markers.items():
            needed.setdefault(attr, []).append(variable)
        needed.setdefault("demo", []).extend(DEMOGRAPHIC_VARS.values())

        for attr, variables in needed.items():
            key = getattr(cycle, attr)
            path = fetch_artifact(
                NHANES_1999_2002[key], directory=cache_directory, allow_network=allow_network
            )
            files[attr] = _read_xpt(path, sorted(set(variables)))

        frame = files["demo"].rename(columns={v: k for k, v in DEMOGRAPHIC_VARS.items()})
        frame = frame[["SEQN", *DEMOGRAPHIC_VARS]]
        for column, (attr, variable) in cycle.markers.items():
            frame = frame.merge(
                files[attr][["SEQN", variable]].rename(columns={variable: column}),
                on="SEQN",
                how="left",
            )
        frame["cycle"] = cycle.name
        per_cycle.append(frame)

    out = pd.concat(per_cycle, ignore_index=True)
    out["crp"] = out["crp"] * CRP_MG_DL_TO_MG_L
    out["SEQN"] = out["SEQN"].astype("int64")
    return out


def load_full(
    *,
    allow_network: bool | None = None,
    cache_directory: Path | None = None,
    min_age: float = MIN_AGE,
) -> pd.DataFrame:
    """The joined cross-layer cohort: clocks + health state + mortality.

    Complete cases on the six resilience markers, as in
    :mod:`geroquery.sources.nhanes` — a subject missing a marker is dropped
    rather than imputed, because imputing the covariance structure is precisely
    what the resilience metrics measure.

    Restricted to ``eligstat == 1``. A participant not eligible for mortality
    linkage has no observable outcome; carrying them as "alive" would code
    missing follow-up as survival.
    """
    dnam = load_dnam(allow_network=allow_network, cache_directory=cache_directory)
    clinical = load_clinical(allow_network=allow_network, cache_directory=cache_directory)
    mortality = load_mortality(allow_network=allow_network, cache_directory=cache_directory)

    frame = dnam.merge(clinical, on="SEQN", how="inner").merge(mortality, on="SEQN", how="left")

    frame["subject_id"] = "NHANES:" + frame["SEQN"].astype(str)
    frame["sex"] = frame["_sex_code"].map(SEX_LABELS)
    frame["dnam_sex"] = frame[DNAM_SEX_COLUMN].map(DNAM_SEX_LABELS)
    frame = frame.rename(columns={DNAM_WEIGHT_COLUMN: "survey_weight"})

    # Follow-up in years, from the MEC exam — the visit the blood was drawn at,
    # so it is the correct origin for a biomarker measured on that blood.
    frame["followup_years"] = frame["permth_exm"] / 12.0
    frame["died"] = (frame["mortstat"] == DECEASED).astype("int64")

    columns = [
        "subject_id",
        "SEQN",
        "cycle",
        "age",
        "sex",
        "dnam_sex",
        *PHENOAGE_MARKERS,
        *CLOCK_COLUMNS,
        *CELL_FRACTIONS,
        "eligstat",
        "mortstat",
        "died",
        "followup_years",
        "permth_exm",
        "permth_int",
        "ucod_leading",
        "survey_weight",
        "survey_weight_mec",
        "sdmvpsu",
        "sdmvstra",
    ]
    frame = frame[columns]
    frame = frame[frame["age"] >= min_age]
    frame = frame[frame["eligstat"] == ELIGIBLE]
    frame = frame.dropna(subset=["age", "sex", "followup_years", *MARKERS])
    return frame.sort_values("SEQN").reset_index(drop=True)


def age_acceleration(frame: pd.DataFrame, clock: str, age_col: str = "age") -> pd.Series:
    """Age acceleration for one clock: the residual of clock age on chronological age.

    The residual, not the raw difference. A raw ``clock - age`` difference
    inherits the clock's calibration error, which is age-dependent for every
    published clock — so the difference correlates with age by construction and
    any association with an age-related outcome is partly that artefact. The
    OLS residual is the standard definition and is orthogonal to age by
    construction.

    Only defined for clocks measured in years (:data:`AGE_LIKE_CLOCKS`).
    DunedinPoAm is a rate and HorvathTelo is a length; neither has a meaningful
    difference from a chronological age.
    """
    if clock not in AGE_LIKE_CLOCKS:
        raise SourceError(
            f"{clock!r} is not measured in years, so age acceleration is undefined for it. "
            f"Age-like clocks: {list(AGE_LIKE_CLOCKS)}.",
            detail={"clock": clock, "age_like": list(AGE_LIKE_CLOCKS)},
        )
    import numpy as np

    ok = frame[[clock, age_col]].notna().all(axis=1)
    x = frame.loc[ok, age_col].to_numpy(dtype=float)
    y = frame.loc[ok, clock].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    out = pd.Series(np.nan, index=frame.index, dtype=float)
    out.loc[ok] = y - (slope * x + intercept)
    return out


def sex_discordance(frame: pd.DataFrame) -> dict[str, float | int]:
    """How often DNAm-predicted sex disagrees with reported sex.

    Reported as a QC number rather than used to filter. A rising discordance
    rate on a future release is the cheapest available signal that samples were
    mixed up upstream; on this release it is 60 of 2,517.
    """
    ok = frame[["sex", "dnam_sex"]].notna().all(axis=1)
    n = int(ok.sum())
    discordant = int((frame.loc[ok, "sex"] != frame.loc[ok, "dnam_sex"]).sum())
    return {
        "n_compared": n,
        "n_discordant": discordant,
        "rate": round(discordant / n, 4) if n else 0.0,
    }


def resilience_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """The six-marker health state, in the shape ``ResilienceService`` expects.

    Narrowed here for the same reason as in :mod:`geroquery.sources.nhanes`: the
    service infers its biomarker list by exclusion, so any extra numeric column
    left in the frame joins the health state without anyone choosing it. This
    frame has twelve clock columns in it, which would be a spectacular version
    of that bug.
    """
    return frame[["subject_id", "age", "sex", *MARKERS, "survey_weight"]].copy()


def full_path(data_dir: Path | None = None) -> Path:
    """Path to the joined table written by ``make crosslayer``."""
    return (data_dir or SOURCES_DATA) / FULL_FILENAME


def sample_path(data_dir: Path | None = None) -> Path:
    """Path to the committed offline sample of real joined rows."""
    return (data_dir or SOURCES_DATA) / SAMPLE_FILENAME


def load_sample(data_dir: Path | None = None) -> pd.DataFrame:
    """Load the committed sample of real cross-layer rows (no network required)."""
    path = sample_path(data_dir)
    if not path.exists():
        raise SourceError(
            f"Cross-layer offline sample missing at {path}. Regenerate it with "
            f"`python -m geroquery.etl.build_crosslayer --write-sample`.",
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
    take.sort_values("SEQN").to_csv(path, index=False)
    return path


class NhanesDnamSource(SourceAdapter):
    """Adapter over the NHANES 1999-2002 methylation subsample."""

    name = "nhanes-dnam-1999-2002"

    def __init__(self, data_dir: Path | None = None, cache_directory: Path | None = None):
        self.data_dir = data_dir or SOURCES_DATA
        self.cache_directory = cache_directory

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("clinical", "methylation"),
            species=("human",),
            federated=False,
            cacheable=True,
            notes=(
                f"{RELEASE}: {len(CLOCKS)} NCHS-computed DNAm clocks joined to the "
                f"six-marker health state and to NCHS 2019 linked mortality on SEQN. "
                f"Checksums last verified {VERIFIED_ON}."
            ),
        )

    def license(self) -> License:
        return License(
            NHANES_DNAM.license, redistributable=True, attribution=NHANES_DNAM.attribution
        )

    def crosslayer_frame(
        self, *, prefer_full: bool = True, allow_network: bool | None = None
    ) -> tuple[pd.DataFrame, str]:
        """Return ``(frame, mode)`` where mode is ``"full"`` or ``"sample"``.

        Same contract as :meth:`NhanesClinicalSource.clinical_frame`: the mode is
        returned rather than logged, because a hazard ratio computed on the
        600-row sample must never be reported as the cohort result.
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
