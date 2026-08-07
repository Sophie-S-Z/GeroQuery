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
    wald_test,
)
from geroquery.survival.cox import SurvivalInputError, _breslow_terms, _score_residuals

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


# --- survey weights and the design-based variance ----------------------------
#
# NHANES is not a simple random sample. WTDN4YR, SDMVSTRA and SDMVPSU were
# carried into the cross-layer table and never used, which made every hazard
# ratio a claim about 2,517 volunteers rather than about US adults aged 50+.
# Applying them changes both the point estimate (unequal selection) and the
# interval (clustering), and each of those needs its own check.


def test_unit_weights_reproduce_the_unweighted_fit_exactly():
    """Adding a weights argument must not perturb the default path."""
    rng = np.random.default_rng(13)
    x = rng.normal(size=(400, 2))
    time = rng.exponential(size=400)
    event = (rng.random(400) < 0.7).astype(int)

    plain = cox_regression(x, time, event, standardize=False)
    unit = cox_regression(x, time, event, weights=np.ones(400), standardize=False)
    assert unit.coefficients == pytest.approx(plain.coefficients, abs=1e-12)
    assert unit.log_likelihood == pytest.approx(plain.log_likelihood, abs=1e-9)


def test_a_weight_of_two_is_the_same_as_the_row_appearing_twice():
    """The pseudo-likelihood's defining property, including at tied times.

    A weight enters both the event term and the risk-set sum. Getting it into
    only one of the two produces a fit that converges and is wrong, and this is
    the check that separates them: duplicating a row is the one case where the
    right answer is known without deriving anything.
    """
    rng = np.random.default_rng(21)
    x = rng.normal(size=(150, 2))
    time = rng.exponential(size=150)
    event = (rng.random(150) < 0.6).astype(int)

    weights = np.ones(150)
    weights[:40] = 2.0
    weighted = cox_regression(x, time, event, weights=weights, standardize=False)

    duplicated = cox_regression(
        np.vstack([x, x[:40]]),
        np.concatenate([time, time[:40]]),
        np.concatenate([event, event[:40]]),
        standardize=False,
    )
    assert weighted.coefficients == pytest.approx(duplicated.coefficients, abs=1e-8)
    assert weighted.population_size == pytest.approx(190.0)


def test_score_residuals_sum_to_the_score_and_so_vanish_at_the_fit():
    """The sandwich's meat is built from per-subject scores; they must be scores.

    A mis-derived residual still yields a plausible-looking variance and is only
    caught by the identity it is supposed to satisfy: summed over subjects it is
    the gradient of the log pseudo-likelihood, which is zero at the MLE. Same
    argument as the finite-difference check above, one derivative further out.
    """
    rng = np.random.default_rng(5)
    x = rng.normal(size=(300, 2))
    time = rng.exponential(size=300)
    event = (rng.random(300) < 0.65).astype(int)
    weights = rng.uniform(0.5, 4.0, 300)

    fit = cox_regression(x, time, event, weights=weights, standardize=False)
    beta = np.asarray(fit.coefficients)
    centred = x - x.mean(axis=0)

    residuals = _score_residuals(centred, time, event, beta, weights)
    assert residuals.shape == (300, 2)

    # The identity, stated exactly rather than as "close to zero": the residuals
    # sum to the gradient. Asserting they sum to zero instead would be asserting
    # the optimizer's tolerance, which is a different claim.
    _, gradient_at_fit, _ = _breslow_terms(centred, time, event, beta, weights)
    assert residuals.sum(axis=0) == pytest.approx(gradient_at_fit, abs=1e-9)
    assert np.abs(gradient_at_fit).max() < 1e-4  # and the fit really is stationary

    # Away from the MLE they must still equal the gradient, where it is large.
    off = beta + 0.3
    _, gradient, _ = _breslow_terms(centred, time, event, off, weights)
    off_residuals = _score_residuals(centred, time, event, off, weights)
    assert np.abs(gradient).max() > 1.0
    assert off_residuals.sum(axis=0) == pytest.approx(gradient, abs=1e-8)


def test_weights_recover_the_population_effect_from_a_deliberately_biased_sample():
    """The planted-effect argument, applied to the survey design itself.

    A population made of two strata whose hazard ratios differ, sampled at
    different rates. The unweighted sample fit answers a question about the
    sample; the weighted one has to answer the question about the population,
    and the population answer is known because the whole population is in hand.
    """
    rng = np.random.default_rng(2024)
    n = 24_000
    stratum = (rng.random(n) < 0.2).astype(int)  # 20% in the oversampled stratum
    x = rng.normal(size=n)
    beta_true = np.where(stratum == 1, 1.1, 0.15)
    latent = rng.exponential(1.0 / np.exp(beta_true * x))
    censor = rng.exponential(3.0, n)
    time = np.minimum(latent, censor)
    event = (latent <= censor).astype(int)

    population = cox_regression(x[:, None], time, event, standardize=False)

    # Sample stratum 1 five times as heavily as stratum 0.
    probability = np.where(stratum == 1, 0.5, 0.1)
    picked = rng.random(n) < probability
    sx, st, se = x[picked, None], time[picked], event[picked]
    weights = 1.0 / probability[picked]

    naive = cox_regression(sx, st, se, standardize=False)
    weighted = cox_regression(sx, st, se, weights=weights, standardize=False)

    target = population.coefficients[0]
    assert weighted.coefficients[0] == pytest.approx(target, abs=0.05)
    # And the unweighted fit is visibly the wrong answer, or the test proves
    # nothing about the weighting.
    assert abs(naive.coefficients[0] - target) > 0.10
    assert weighted.weighted and weighted.variance == "robust"


def test_clustering_widens_the_interval_and_reports_the_design_effect():
    """A PSU-level effect the model does not carry makes residuals correlate
    inside a cluster, and that is exactly when the design variance must exceed
    the independence one.

    Worth stating why the obvious construction does not work: making the
    *covariate* cluster-level is not enough. If the model is correctly specified
    the score residuals inside a PSU stay uncorrelated and the design effect
    comes out below 1 — which is legitimate, and was the first version of this
    test failing for the right reason. What inflates a clustered variance is
    shared variation the model has not accounted for, so that is what is planted.
    """
    rng = np.random.default_rng(77)
    n_strata, psus_per_stratum, per_psu = 24, 2, 30
    strata, psu, x, frailty = [], [], [], []
    for h in range(n_strata):
        for a in range(psus_per_stratum):
            shared, unmodelled = rng.normal(), rng.normal()
            strata += [h] * per_psu
            psu += [a] * per_psu
            x += list(shared + rng.normal(scale=0.3, size=per_psu))
            frailty += [unmodelled] * per_psu
    strata, psu, x, frailty = (np.array(v) for v in (strata, psu, x, frailty))
    n = len(x)
    latent = rng.exponential(1.0 / np.exp(0.6 * x + 1.2 * frailty))
    censor = rng.exponential(3.0, n)
    time = np.minimum(latent, censor)
    event = (latent <= censor).astype(int)
    weights = np.full(n, 400.0)

    independent = cox_regression(x[:, None], time, event, weights=weights, standardize=False)
    design = cox_regression(
        x[:, None], time, event, weights=weights, strata=strata, psu=psu, standardize=False
    )

    assert independent.variance == "robust"
    assert design.variance == "design"
    assert design.n_psu == n_strata * psus_per_stratum
    assert design.n_strata == n_strata
    assert design.lonely_psu_strata == 0
    # Same point estimate — the design changes the interval, not the answer.
    assert design.coefficients == pytest.approx(independent.coefficients, abs=1e-12)
    assert design.standard_errors[0] > independent.standard_errors[0] * 1.5
    assert design.design_effect[0] > 2.0
    assert independent.design_effect[0] == 1.0
    assert design.population_size == pytest.approx(400.0 * n)
    assert design.population_events == pytest.approx(400.0 * event.sum())


def test_a_stratum_with_one_psu_is_reported_rather_than_silently_contributing_zero():
    rng = np.random.default_rng(4)
    n = 400
    x = rng.normal(size=(n, 1))
    time = rng.exponential(size=n)
    event = (rng.random(n) < 0.7).astype(int)
    strata = np.repeat(np.arange(10), 40)
    psu = np.tile(np.repeat([0, 1], 20), 10)
    psu[:40] = 0  # stratum 0 now has a single PSU

    fit = cox_regression(
        x, time, event, weights=np.ones(n), strata=strata, psu=psu, standardize=False
    )
    assert fit.lonely_psu_strata == 1


def test_a_weighted_fit_refuses_a_likelihood_ratio_test():
    """2*delta log pseudo-likelihood is not chi-squared under survey weighting.

    Reporting it anyway is the kind of error that produces a confident p-value
    from an invalid reference distribution, so it is refused and the design-based
    Wald test is what the cross-layer analysis uses instead.
    """
    rng = np.random.default_rng(9)
    x = rng.normal(size=(300, 2))
    time = rng.exponential(size=300)
    event = (rng.random(300) < 0.7).astype(int)
    weights = rng.uniform(0.5, 3.0, 300)

    full = cox_regression(x, time, event, ["a", "b"], weights=weights)
    reduced = cox_regression(x[:, :1], time, event, ["a"], weights=weights)
    with pytest.raises(SurvivalInputError, match="pseudo-likelihood"):
        likelihood_ratio_test(full, reduced)

    wald = wald_test(full, ["b"])
    assert wald["df"] == 1
    assert 0.0 <= wald["p_value"] <= 1.0
    assert wald["tested"] == ["b"]


def test_the_wald_test_agrees_with_the_likelihood_ratio_when_unweighted():
    """Two asymptotically equivalent tests on the same well-behaved fit.

    Not a tautology: they are computed from different quantities — one from the
    curvature at the fit, one from the height of the likelihood at two points —
    so agreement is evidence the covariance being handed to the Wald test is the
    covariance of the coefficients it is testing.
    """
    rng = np.random.default_rng(31)
    n = 4000
    x = rng.normal(size=(n, 2))
    latent = rng.exponential(1.0 / np.exp(x @ np.array([0.5, 0.3])))
    censor = rng.exponential(2.0, n)
    time = np.minimum(latent, censor)
    event = (latent <= censor).astype(int)

    full = cox_regression(x, time, event, ["a", "b"], standardize=False)
    reduced = cox_regression(x[:, :1], time, event, ["a"], standardize=False)
    lr = likelihood_ratio_test(full, reduced)
    wald = wald_test(full, ["b"])
    assert wald["statistic"] == pytest.approx(lr["statistic"], rel=0.15)


def test_weighted_concordance_uses_the_population_not_the_sample():
    """A C-index over an unequal-probability sample describes the sample."""
    time = np.array([1.0, 2.0, 3.0, 4.0])
    event = np.array([1, 1, 1, 0])
    risk = np.array([4.0, 3.0, 1.0, 2.0])
    assert concordance_index(time, event, risk) == pytest.approx(5.0 / 6.0)
    # Weighting the one discordant pair's members up must lower the value.
    heavy = np.array([1.0, 1.0, 8.0, 8.0])
    assert concordance_index(time, event, risk, weights=heavy) < 0.7


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


def test_crosslayer_carries_the_design_through_the_complete_case_filter():
    """Weights must be filtered with their rows, not aligned afterwards.

    A weight vector one row out of step produces a fit that converges, reports a
    plausible hazard ratio, and is wrong for every subject — the same failure
    shape as bug #12 and #13, and equally invisible. So the design columns go
    through the same ``dropna`` as the covariates and the test punches holes to
    prove it.
    """
    frame, distance = _toy_cohort()
    rng = np.random.default_rng(17)
    frame = frame.copy()
    frame["w"] = rng.uniform(1_000, 40_000, len(frame))
    frame["stratum"] = np.repeat(np.arange(15), len(frame) // 15)[: len(frame)]
    frame["unit"] = np.tile([1, 2], len(frame))[: len(frame)]
    clock = pd.Series(rng.normal(size=len(frame)), index=frame.index)
    clock.iloc[:25] = np.nan  # holes that shift every later row

    result = crosslayer_analysis(
        frame,
        clock,
        distance,
        clock_name="c",
        weight_col="w",
        strata_col="stratum",
        psu_col="unit",
    )
    assert result.weighted
    assert result.n_subjects == len(frame) - 25
    assert result.joint_model.variance == "design"
    assert result.joint_model.n_strata == 15
    # The population the weights describe, not the number of people measured.
    kept = frame["w"].iloc[25:]
    assert result.joint_model.population_size == pytest.approx(kept.sum())

    # And the nested tests switched to the one that survives pseudo-likelihood.
    for test in (
        result.dysregulation_adds_over_clock,
        result.clock_adds_over_dysregulation,
        result.clock_adds_over_baseline,
    ):
        assert test["test"] == "design_wald"
        assert test["variance"] == "design"


def test_an_unweighted_crosslayer_still_uses_the_likelihood_ratio():
    """The dispatch must not quietly change the unweighted result."""
    frame, distance = _toy_cohort()
    clock = pd.Series(np.linspace(-2, 2, len(frame)), index=frame.index)
    result = crosslayer_analysis(frame, clock, distance, clock_name="c")
    assert not result.weighted
    assert result.dysregulation_adds_over_clock["test"] == "likelihood_ratio"
    assert result.joint_model.variance == "model"


def test_dysregulation_reference_must_be_large_enough():
    frame, _ = _toy_cohort(n=900)
    with pytest.raises(SurvivalInputError, match="Reference sample"):
        mahalanobis_dysregulation(
            frame,
            ["albumin", "creatinine", "glucose", "crp", "lymphocyte_pct", "rdw"],
            reference_age_max=50.5,
        )
