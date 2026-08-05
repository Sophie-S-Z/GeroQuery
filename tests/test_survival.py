"""Tests for the Cox estimator and the cross-layer analysis.

The Cox model is hand-rolled (no lifelines, no statsmodels), so it has to earn
trust rather than inherit it. Three independent checks, none of which can pass
by accident together:

* a closed form small enough to differentiate by hand,
* finite-difference agreement between the analytic derivatives and the
  likelihood they claim to differentiate,
* recovery of a planted log hazard ratio from simulated data.

The third is the same argument the synthetic CSD fixture exists to make: an
estimator that cannot find an effect known to be there has not earned the right
to report a null.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from geroquery.survival import (
    concordance_index,
    cox_regression,
    crosslayer_analysis,
    likelihood_ratio_test,
    mahalanobis_dysregulation,
)
from geroquery.survival.cox import SurvivalInputError, _breslow_terms

# --- the estimator itself ----------------------------------------------------


def test_cox_matches_a_hand_computable_closed_form():
    """Three subjects, one binary covariate, no censoring.

    With x = [0, 1, 0] and event times 1 < 2 < 3, the Breslow log partial
    likelihood is

        L(b) = -log(e^b + 2) + b - log(e^b + 1)

    whose derivative is (2 - e^{2b}) / ((e^b + 2)(e^b + 1)). That is zero at
    e^b = sqrt(2), so the MLE is exactly b = ln(2)/2. Any tie handling, risk-set
    ordering, or sign error moves this.
    """
    result = cox_regression(
        np.array([[0.0], [1.0], [0.0]]), [1, 2, 3], [1, 1, 1], standardize=False
    )
    assert result.converged
    assert result.coefficients[0] == pytest.approx(math.log(2) / 2, abs=1e-9)
    assert result.hazard_ratios[0] == pytest.approx(math.sqrt(2), abs=1e-9)


def test_analytic_derivatives_match_finite_differences():
    """The gradient and Hessian must differentiate the likelihood they ship with.

    This is the check that catches a risk-set built in the wrong direction: such
    a bug still produces a smooth, convergent, entirely wrong fit, and only
    disagrees with the numerical derivative of its own objective.
    """
    rng = np.random.default_rng(7)
    x = rng.normal(size=(60, 3))
    time = rng.exponential(size=60) * 10
    event = (rng.random(60) < 0.6).astype(int)
    beta = rng.normal(scale=0.3, size=3)

    _, gradient, neg_hessian = _breslow_terms(x, time, event, beta)

    eye = np.eye(3)
    numeric_gradient = np.array(
        [
            (
                _breslow_terms(x, time, event, beta + eye[i] * 1e-6)[0]
                - _breslow_terms(x, time, event, beta - eye[i] * 1e-6)[0]
            )
            / 2e-6
            for i in range(3)
        ]
    )
    numeric_hessian = np.array(
        [
            (
                _breslow_terms(x, time, event, beta + eye[i] * 1e-5)[1]
                - _breslow_terms(x, time, event, beta - eye[i] * 1e-5)[1]
            )
            / 2e-5
            for i in range(3)
        ]
    )
    assert np.abs(gradient - numeric_gradient).max() < 1e-6
    # _breslow_terms returns the NEGATIVE Hessian (the observed information),
    # which is what Newton-Raphson solves against.
    assert np.abs(neg_hessian + numeric_hessian).max() < 1e-6


def test_cox_recovers_a_planted_hazard_ratio():
    """Simulated exponential survival with known coefficients."""
    rng = np.random.default_rng(42)
    n = 6000
    x = rng.normal(size=(n, 2))
    true = np.array([0.7, -0.4])
    latent = rng.exponential(1.0 / np.exp(x @ true))
    censor = rng.exponential(2.0, n)
    time = np.minimum(latent, censor)
    event = (latent <= censor).astype(int)

    result = cox_regression(x, time, event, ["a", "b"], standardize=False)
    assert result.converged
    for estimated, expected in zip(result.coefficients, true, strict=True):
        assert estimated == pytest.approx(expected, abs=0.06)
    # A predictor that genuinely drives the hazard must beat a coin flip.
    assert result.concordance > 0.6


def test_standardized_hazard_ratio_is_per_standard_deviation():
    """Scaling a covariate must scale its coefficient, not change the fit."""
    rng = np.random.default_rng(3)
    n = 2000
    x = rng.normal(size=(n, 1)) * 7.0
    latent = rng.exponential(1.0 / np.exp(x[:, 0] * 0.1))
    event = np.ones(n, dtype=int)

    raw = cox_regression(x, latent, event, standardize=False)
    scaled = cox_regression(x, latent, event, standardize=True)
    assert scaled.standardized and not raw.standardized
    assert scaled.coefficients[0] == pytest.approx(
        raw.coefficients[0] * x[:, 0].std(ddof=1), rel=1e-6
    )


def test_concordance_is_one_for_a_perfect_ranker_and_half_for_a_constant():
    time = np.arange(1.0, 11.0)
    event = np.ones(10, dtype=int)
    assert concordance_index(time, event, -time) == pytest.approx(1.0)
    assert concordance_index(time, event, time) == pytest.approx(0.0)
    # A constant predictor orders nothing; ties count as half, so it is exactly
    # a coin flip rather than accidentally perfect.
    assert concordance_index(time, event, np.zeros(10)) == pytest.approx(0.5)


def test_cox_rejects_input_it_cannot_fit():
    x = np.zeros((5, 1))
    with pytest.raises(SurvivalInputError, match="No events"):
        cox_regression(x, [1, 2, 3, 4, 5], [0, 0, 0, 0, 0])
    with pytest.raises(SurvivalInputError, match="non-negative"):
        cox_regression(x, [-1, 2, 3, 4, 5], [1, 0, 0, 0, 0])
    with pytest.raises(SurvivalInputError, match="same number of rows"):
        cox_regression(x, [1, 2, 3], [1, 0, 1])


def test_likelihood_ratio_refuses_models_fitted_on_different_rows():
    """A likelihood ratio between different complete-case subsets is meaningless.

    Nothing downstream can detect it, so it is refused here rather than
    reported as a p-value.
    """
    rng = np.random.default_rng(11)
    x = rng.normal(size=(300, 2))
    time = rng.exponential(size=300)
    event = (rng.random(300) < 0.7).astype(int)

    full = cox_regression(x, time, event, ["a", "b"])
    reduced = cox_regression(x[:200, :1], time[:200], event[:200], ["a"])
    with pytest.raises(SurvivalInputError, match="identical rows"):
        likelihood_ratio_test(full, reduced)


# --- the cross-layer composition --------------------------------------------


def _toy_cohort(n: int = 900, seed: int = 5) -> tuple[pd.DataFrame, pd.Series]:
    """A cohort where dysregulation drives death and the clock does not.

    Built so the analysis has a known right answer: the health state carries the
    signal, the clock is noise. If the composition ever wires the two together,
    this flips.
    """
    rng = np.random.default_rng(seed)
    age = rng.uniform(50, 85, n)
    markers = {
        "albumin": rng.normal(4.2, 0.3, n),
        "creatinine": np.exp(rng.normal(0.0, 0.25, n)),
        "glucose": np.exp(rng.normal(4.6, 0.2, n)),
        "crp": np.exp(rng.normal(0.5, 1.0, n)),
        "lymphocyte_pct": rng.normal(30, 7, n),
        "rdw": np.exp(rng.normal(2.6, 0.1, n)),
    }
    frame = pd.DataFrame({"age": age, "sex": rng.choice(["male", "female"], n), **markers})
    dysregulation = mahalanobis_dysregulation(frame, list(markers))
    signal = (dysregulation.values - dysregulation.values.mean()) / dysregulation.values.std()

    latent = rng.exponential(1.0 / np.exp(0.05 * (age - 50) + 0.8 * signal))
    censor = rng.exponential(25.0, n)
    frame["followup_years"] = np.minimum(latent, censor)
    frame["died"] = (latent <= censor).astype(int)
    return frame, dysregulation.values


def test_dysregulation_is_orthogonal_to_nothing_but_is_finite_and_positive():
    frame, distance = _toy_cohort()
    assert distance.notna().all()
    assert (distance >= 0).all()
    # A Mahalanobis distance over six markers has a mode well away from zero;
    # a near-zero median would mean the covariance was estimated on the wrong
    # rows and everyone landed on the centroid.
    assert distance.median() > 1.0


def test_crosslayer_finds_the_planted_signal_and_not_the_planted_noise():
    frame, distance = _toy_cohort()
    rng = np.random.default_rng(99)
    noise_clock = pd.Series(rng.normal(size=len(frame)), index=frame.index)

    result = crosslayer_analysis(
        frame, noise_clock, distance, clock_name="noise", predicts="nothing"
    )
    rows = {row["covariate"]: row for row in result.joint_model.summary_rows()}

    # The planted health-state effect is found...
    assert rows["dysregulation"]["excludes_null"]
    assert rows["dysregulation"]["hazard_ratio"] > 1.0
    assert result.dysregulation_adds_over_baseline["p_value"] < 1e-6
    # ...and the pure-noise clock is not claimed.
    assert not rows["clock_acceleration"]["excludes_null"]
    assert result.clock_adds_over_dysregulation["p_value"] > 0.01


def test_crosslayer_models_are_fitted_on_identical_rows():
    """All four models must share a complete-case set, or the LR tests are void."""
    frame, distance = _toy_cohort()
    clock = pd.Series(np.linspace(-3, 3, len(frame)), index=frame.index)
    # Punch a hole in each predictor in a different place.
    clock.iloc[:20] = np.nan
    distance = distance.copy()
    distance.iloc[-30:] = np.nan

    result = crosslayer_analysis(frame, clock, distance, clock_name="c")
    counts = {
        result.baseline.n_subjects,
        result.clock_model.n_subjects,
        result.dysregulation_model.n_subjects,
        result.joint_model.n_subjects,
    }
    assert len(counts) == 1
    assert result.n_subjects == len(frame) - 50


def test_dysregulation_reference_must_be_large_enough():
    frame, _ = _toy_cohort(n=900)
    with pytest.raises(SurvivalInputError, match="Reference sample"):
        mahalanobis_dysregulation(
            frame,
            ["albumin", "creatinine", "glucose", "crp", "lymphocyte_pct", "rdw"],
            reference_age_max=50.5,
        )
