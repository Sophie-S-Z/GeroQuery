"""The cross-layer analysis: does an epigenetic clock, a dysregulated health
state, or both predict death — and does either add anything over the other.

This module answers the question `docs/HANDOFF.md` called "the reason this
exists" and could not answer, because clocks and health state used to live in
different cohorts. They no longer do; see
:mod:`geroquery.sources.nhanes_dnam`.

Three questions, three model comparisons, all nested:

======  ==========================================  ================================
Model   Covariates                                  What it is for
======  ==========================================  ================================
A       age, sex                                    Baseline. Everything must beat it
B       age, sex, clock age acceleration            Q1: does the clock predict death?
C       age, sex, dysregulation                     Q2: does the health state?
D       age, sex, both                              Q3: does either add over the other?
======  ==========================================  ================================

Q3 is the one nothing else reports, and it is the only one that needs the models
to be nested rather than merely compared: two hazard ratios side by side cannot
tell you whether the second predictor carries information the first does not.

**On the resilience measure.** Critical slowing down as implemented in
:mod:`geroquery.resilience` is a *population* statistic — variance and
cross-correlation across age strata. It has no per-subject value, so it cannot
enter a Cox model. The individual-level analogue used here is Mahalanobis
distance from a young-reference centroid over the same six markers, the
"physiological dysregulation" statistic of Cohen et al. (2013). It is computed
from the same covariance structure that CSD measures the drift of, which is why
it is the right per-subject stand-in and also why the two must not be described
as the same quantity. The population CSD is computed alongside and reported
separately.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .cox import CoxResult, SurvivalInputError, cox_regression, likelihood_ratio_test

# Markers whose absolute skew in the reference sample exceeds this are
# log-transformed before the covariance is estimated. Mahalanobis distance
# assumes approximate multivariate normality; CRP in particular is so
# right-skewed that without this the statistic is mostly "how high is this
# person's CRP" rather than a multi-system measure.
SKEW_LOG_THRESHOLD = 1.0

# Reference population for the dysregulation centroid: the youngest stratum of
# this cohort. Cohen et al. use a healthy reference; we do not have a health
# screen here, so the reference is defined by age alone and that is stated
# rather than dressed up. Everyone in this cohort is 50+, so "young" is
# relative — which weakens the statistic and is a limitation, not a choice.
REFERENCE_AGE_MAX = 60.0

# Minimum reference-sample size for a stable 6x6 covariance estimate.
MIN_REFERENCE_N = 100


@dataclass(frozen=True)
class DysregulationResult:
    """Per-subject Mahalanobis dysregulation, plus how it was defined."""

    values: pd.Series
    markers: list[str]
    log_transformed: list[str]
    reference_n: int
    reference_age_max: float
    skew: dict[str, float]

    def to_dict(self) -> dict:
        return {
            "markers": self.markers,
            "log_transformed": self.log_transformed,
            "reference_n": self.reference_n,
            "reference_age_max": self.reference_age_max,
            "skew": {k: round(v, 3) for k, v in self.skew.items()},
            "n": int(self.values.notna().sum()),
            "median": round(float(self.values.median()), 4),
        }


def mahalanobis_dysregulation(
    frame: pd.DataFrame,
    markers: list[str],
    *,
    age_col: str = "age",
    reference_age_max: float = REFERENCE_AGE_MAX,
    skew_threshold: float = SKEW_LOG_THRESHOLD,
) -> DysregulationResult:
    """Per-subject distance from the young-reference marker centroid.

    Skewed markers are log-transformed first, chosen by measured skew in the
    reference sample rather than by a hard-coded list, so the choice is
    auditable and moves with the data instead of with an author's memory.

    The covariance is estimated on the reference sample only. Estimating it on
    the whole cohort would let the oldest, most dysregulated subjects define
    what "normal covariance" means, which shrinks exactly the distances the
    statistic exists to detect.
    """
    missing = [m for m in markers if m not in frame.columns]
    if missing:
        raise SurvivalInputError(
            f"Frame is missing markers {missing}.", detail={"missing": missing}
        )

    values = frame[markers].astype(float).copy()
    reference_mask = frame[age_col] <= reference_age_max
    if int(reference_mask.sum()) < MIN_REFERENCE_N:
        raise SurvivalInputError(
            f"Reference sample (age <= {reference_age_max}) has "
            f"{int(reference_mask.sum())} subjects; need at least {MIN_REFERENCE_N} "
            f"for a stable {len(markers)}x{len(markers)} covariance.",
            detail={"reference_n": int(reference_mask.sum())},
        )

    skew = {m: float(values.loc[reference_mask, m].skew()) for m in markers}
    logged = [m for m in markers if abs(skew[m]) > skew_threshold]
    for marker in logged:
        shift = values[marker].min()
        # log1p of a shifted value, so a marker with legitimate zeros (CRP below
        # the detection limit is reported as 0) does not become -inf.
        values[marker] = np.log1p(values[marker] - min(shift, 0.0))

    reference = values.loc[reference_mask]
    centroid = reference.mean().to_numpy()
    covariance = np.cov(reference.to_numpy(), rowvar=False)
    try:
        precision = np.linalg.inv(covariance)
    except np.linalg.LinAlgError as exc:
        raise SurvivalInputError(
            "Reference covariance is singular; markers are collinear.",
            detail={"markers": markers},
        ) from exc

    centred = values.to_numpy() - centroid
    # einsum rather than a loop: n x p times p x p times p x n, keeping only the
    # diagonal, which is all the distance needs.
    squared = np.einsum("ij,jk,ik->i", centred, precision, centred)
    distance = pd.Series(np.sqrt(np.clip(squared, 0.0, None)), index=frame.index)
    distance[values.isna().any(axis=1)] = np.nan

    return DysregulationResult(
        values=distance,
        markers=list(markers),
        log_transformed=logged,
        reference_n=int(reference_mask.sum()),
        reference_age_max=reference_age_max,
        skew=skew,
    )


@dataclass(frozen=True)
class CrossLayerResult:
    """One clock's worth of the three-question analysis."""

    clock: str
    predicts: str
    n_subjects: int
    n_events: int
    baseline: CoxResult
    clock_model: CoxResult
    dysregulation_model: CoxResult
    joint_model: CoxResult
    clock_adds_over_baseline: dict
    dysregulation_adds_over_baseline: dict
    dysregulation_adds_over_clock: dict
    clock_adds_over_dysregulation: dict

    def to_dict(self) -> dict:
        return {
            "clock": self.clock,
            "predicts": self.predicts,
            "n_subjects": self.n_subjects,
            "n_events": self.n_events,
            "models": {
                "A_baseline": self.baseline.to_dict(),
                "B_clock": self.clock_model.to_dict(),
                "C_dysregulation": self.dysregulation_model.to_dict(),
                "D_joint": self.joint_model.to_dict(),
            },
            "tests": {
                "clock_vs_baseline": self.clock_adds_over_baseline,
                "dysregulation_vs_baseline": self.dysregulation_adds_over_baseline,
                "dysregulation_adds_over_clock": self.dysregulation_adds_over_clock,
                "clock_adds_over_dysregulation": self.clock_adds_over_dysregulation,
            },
        }

    def headline(self) -> str:
        """One line a human can read without opening the nested dictionaries."""
        joint = self.joint_model.summary_rows()
        by_name = {row["covariate"]: row for row in joint}
        clock_row = by_name.get("clock_acceleration", {})
        dys_row = by_name.get("dysregulation", {})
        return (
            f"{self.clock}: clock HR/SD {clock_row.get('hazard_ratio')} "
            f"[{clock_row.get('ci_low')}, {clock_row.get('ci_high')}]; "
            f"dysregulation HR/SD {dys_row.get('hazard_ratio')} "
            f"[{dys_row.get('ci_low')}, {dys_row.get('ci_high')}]; "
            f"dysregulation adds over clock p={self.dysregulation_adds_over_clock['p_value']:.3g}"
        )


def crosslayer_analysis(
    frame: pd.DataFrame,
    clock_acceleration: pd.Series,
    dysregulation: pd.Series,
    *,
    clock_name: str,
    predicts: str = "",
    age_col: str = "age",
    sex_col: str = "sex",
    time_col: str = "followup_years",
    event_col: str = "died",
) -> CrossLayerResult:
    """Fit the four nested models for one clock.

    One design matrix is built and then sliced, so every model is fitted on
    exactly the same rows. Building four matrices independently would let each
    model drop a different set of incomplete cases, and a likelihood ratio
    between models fitted on different subsets is not a test of anything.
    """
    data = pd.DataFrame(
        {
            "age": frame[age_col].astype(float),
            "is_female": (frame[sex_col] == "female").astype(float),
            "clock_acceleration": pd.to_numeric(clock_acceleration, errors="coerce"),
            "dysregulation": pd.to_numeric(dysregulation, errors="coerce"),
            "time": frame[time_col].astype(float),
            "event": frame[event_col].astype(int),
        }
    ).dropna()

    if data.empty:
        raise SurvivalInputError("No complete cases for the cross-layer analysis.")

    time = data["time"].to_numpy()
    event = data["event"].to_numpy()

    def fit(columns: list[str]) -> CoxResult:
        return cox_regression(data[columns].to_numpy(), time, event, columns)

    base_cols = ["age", "is_female"]
    baseline = fit(base_cols)
    clock_model = fit([*base_cols, "clock_acceleration"])
    dysregulation_model = fit([*base_cols, "dysregulation"])
    joint = fit([*base_cols, "clock_acceleration", "dysregulation"])

    return CrossLayerResult(
        clock=clock_name,
        predicts=predicts,
        n_subjects=int(len(data)),
        n_events=int(event.sum()),
        baseline=baseline,
        clock_model=clock_model,
        dysregulation_model=dysregulation_model,
        joint_model=joint,
        clock_adds_over_baseline=likelihood_ratio_test(clock_model, baseline),
        dysregulation_adds_over_baseline=likelihood_ratio_test(dysregulation_model, baseline),
        # The two that matter: each predictor tested against a model that
        # already contains the other.
        dysregulation_adds_over_clock=likelihood_ratio_test(joint, clock_model),
        clock_adds_over_dysregulation=likelihood_ratio_test(joint, dysregulation_model),
    )
