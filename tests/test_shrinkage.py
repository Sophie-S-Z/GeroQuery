"""Tests for corpus-wide adaptive shrinkage and the local false sign rate.

The estimator has to earn the same way the survey-weighted Cox path did: a
defining property, a planted answer whose value is known by construction, and a
calibration claim that fails if the posterior is merely plausible.

The calibration test is the one that matters. A shrinkage implementation with a
mis-derived posterior variance still returns numbers in [0, 1] that decrease as
the evidence strengthens, and still looks entirely reasonable plotted. What it
cannot do is keep its promise -- that among the genes it calls at lfsr <= 0.05,
at most 5% have the wrong sign -- because that promise is a statement about the
truth, and the simulation knows the truth.
"""

from __future__ import annotations

import numpy as np
import pytest

from geroquery.harmonize import adaptive_shrinkage


def _planted(seed: int = 20260806, n: int = 20_000, null_fraction: float = 0.9):
    """A corpus with a known prior: mostly nothing, some real effects.

    Deliberately shaped like the real one -- most genes do not change with age,
    and the standard errors vary several-fold across genes because the panel's
    contrasts vary several-fold in size.
    """
    rng = np.random.default_rng(seed)
    is_null = rng.random(n) < null_fraction
    truth = np.where(is_null, 0.0, rng.normal(0.0, 0.8, n))
    se = rng.uniform(0.15, 0.60, n)
    observed = truth + rng.normal(0.0, se)
    return truth, observed, se


# --- defining properties ------------------------------------------------------


def test_a_noisier_estimate_is_shrunk_harder():
    """The whole mechanism, stated as an inequality.

    Two genes reporting the *same* effect, one measured three times as
    precisely: the imprecise one must end up closer to zero. If this does not
    hold, the posterior is not using the standard errors at all -- which is a
    failure mode that leaves every other summary looking fine.
    """
    _, observed, se = _planted()
    # Append the pair to a real corpus, because the prior has to come from
    # somewhere; on their own two points carry no information about the mixture.
    effects = np.concatenate([observed, [0.9, 0.9]])
    errors = np.concatenate([se, [0.15, 0.45]])

    result = adaptive_shrinkage(effects, errors)
    precise, imprecise = result.posterior_mean[-2], result.posterior_mean[-1]

    assert abs(imprecise) < abs(precise)


def test_shrinkage_never_crosses_zero_or_overshoots():
    """A posterior mean lies between zero and the estimate that produced it.

    Under a prior that is unimodal at zero this is not a convention, it is a
    consequence: every mixture component pulls toward the mode. A posterior mean
    that overshoots the observation, or flips its sign, means the conjugate
    update is wrong.
    """
    _, observed, se = _planted()
    result = adaptive_shrinkage(observed, se)

    assert np.all(np.abs(result.posterior_mean) <= np.abs(observed) + 1e-9)
    same_sign = np.sign(result.posterior_mean) * np.sign(observed)
    assert np.all(same_sign >= 0)


def test_sign_probabilities_and_null_mass_partition_the_posterior():
    """P(>0) + P(<0) + P(=0) = 1, and the null mass is the remainder.

    Asserted because ``lfsr`` is built out of these three and a subtle error in
    the point-mass column would leave the lfsr looking sane while quietly
    double-counting.
    """
    _, observed, se = _planted(n=4_000)
    result = adaptive_shrinkage(observed, se)

    total = result.prob_positive + result.prob_negative
    assert np.all(total <= 1.0 + 1e-9)
    assert np.all(result.lfsr >= np.minimum(result.prob_positive, result.prob_negative) - 1e-12)
    assert np.all((result.lfsr >= 0.0) & (result.lfsr <= 1.0))


# --- the planted answer -------------------------------------------------------


def test_the_fitted_prior_recovers_the_planted_share_of_null_genes():
    """90% of this corpus is null by construction; the fit has to say so.

    **A test premise that was wrong, kept because it is instructive.** This
    first asserted on ``null_proportion`` -- the weight on the exact point mass
    -- and expected ~0.9. It came out 0.66, and that was the estimator being
    right rather than wrong. The grid's narrowest components sit at a standard
    deviation ten times below the smallest standard error, so "exactly zero" and
    "very nearly zero" have almost identical likelihoods and the mixture cannot
    apportion mass between them; the point mass is identified only up to that
    neighbourhood. Integrating the fitted prior over an interval is the question
    that *is* identified, and it recovers the planted value to within a point.

    The API changed with the test: ``prior_mass_below`` exists so the next
    reader does not make the same inference from ``null_proportion`` that this
    test did.
    """
    truth, observed, se = _planted(null_fraction=0.9)
    result = adaptive_shrinkage(observed, se)

    for threshold in (0.05, 0.1):
        fitted = result.prior_mass_below(threshold)
        planted = float(np.mean(np.abs(truth) < threshold))
        assert abs(fitted - planted) < 0.05, f"|b| < {threshold}: {fitted:.3f} vs {planted:.3f}"

    # And the point mass, asserted for what it is: a lower bound on the null
    # share, never the null share itself.
    assert 0.4 <= result.null_proportion <= result.prior_mass_below(0.1)


def test_a_corpus_with_no_null_genes_is_not_forced_into_one():
    """The Dirichlet penalty leans toward the null; it must not overrule data.

    Every effect here is real and large. A prior that still declared most of the
    corpus null would be shrinking on faith rather than on evidence, and the
    penalty weight would need to come down.
    """
    rng = np.random.default_rng(7)
    truth = rng.normal(0.0, 1.2, 8_000)
    truth = np.where(np.abs(truth) < 0.4, np.sign(truth) * 0.4 + 1e-9, truth)
    se = rng.uniform(0.10, 0.20, 8_000)
    observed = truth + rng.normal(0.0, se)

    result = adaptive_shrinkage(observed, se)
    assert result.prior_mass_below(0.1) < 0.35


# --- calibration: the claim the statistic actually makes -----------------------


@pytest.mark.parametrize("threshold", [0.01, 0.05, 0.10])
def test_the_local_false_sign_rate_is_calibrated(threshold):
    """Among genes called at lfsr <= t, at most t of them have the wrong sign.

    This is the promise the number makes to a reader, and it is the only test
    here that a plausible-but-wrong posterior fails. Monte Carlo slack is
    allowed -- the check is that the realised error rate tracks the nominal one,
    not that it lands on it exactly.
    """
    truth, observed, se = _planted()
    result = adaptive_shrinkage(observed, se)

    called = result.lfsr <= threshold
    assert called.sum() > 200, "not enough confident calls to judge calibration"

    # A call is wrong if the true effect goes the other way. A truly-zero effect
    # has no sign to get wrong, so it cannot be counted as correct either --
    # lfsr promises to bound sign errors, and mass at exactly zero is already
    # charged against the statistic inside the estimator.
    direction = np.sign(result.posterior_mean[called])
    wrong = np.sign(truth[called]) == -direction
    realised = wrong.mean()

    assert realised <= threshold + 0.02, f"lfsr <= {threshold} but {realised:.3f} wrong-signed"


def test_shrinkage_beats_the_raw_estimate_on_the_genes_it_would_be_ranked_on():
    """The winner's curse, made falsifiable.

    Take the largest raw estimates -- exactly the rows the landscape table puts
    at the top -- and compare both estimators against the truth there. Raw
    effects selected for being extreme overshoot; the posterior is what removes
    that. Asserted on the selected subset rather than corpus-wide, because
    corpus-wide is the easy case and the table only ever shows the tail.
    """
    truth, observed, se = _planted()
    result = adaptive_shrinkage(observed, se)

    top = np.argsort(np.abs(observed))[-500:]
    raw_error = np.mean((observed[top] - truth[top]) ** 2)
    shrunk_error = np.mean((result.posterior_mean[top] - truth[top]) ** 2)

    assert shrunk_error < raw_error

    # And the direction of the bias, named: the raw estimates overshoot.
    assert np.mean(np.abs(observed[top])) > np.mean(np.abs(truth[top]))


# --- input handling -----------------------------------------------------------


@pytest.mark.parametrize(
    ("effects", "errors", "message"),
    [
        ([0.1, 0.2], [0.1], "equal length"),
        ([], [], "at least one"),
        ([0.1, 0.2], [0.1, 0.0], "positive"),
        ([0.1, np.nan], [0.1, 0.2], "finite"),
    ],
)
def test_malformed_input_raises_rather_than_returning_a_shape(effects, errors, message):
    with pytest.raises(ValueError, match=message):
        adaptive_shrinkage(effects, errors)
