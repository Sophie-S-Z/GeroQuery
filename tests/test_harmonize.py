"""M3 harmonize — scientific correctness on planted-truth synthetic data."""

import numpy as np
import pytest
from scipy import stats

from geroquery.harmonize import (
    batch_variance,
    hedges_g,
    random_effects,
    remove_batch_effect,
)


def test_hedges_g_recovers_planted_effect():
    rng = np.random.default_rng(0)
    young = rng.normal(0.0, 1.0, 400)
    old = rng.normal(2.0, 1.0, 400)  # true standardized diff ~= 2
    es = hedges_g(young, old)
    assert es.effect == pytest.approx(2.0, abs=0.2)
    assert es.direction == "up"
    assert es.standard_error > 0


def test_hedges_g_sign_flips_for_down_regulation():
    rng = np.random.default_rng(1)
    young = rng.normal(2.0, 1.0, 200)
    old = rng.normal(0.0, 1.0, 200)
    es = hedges_g(young, old)
    assert es.effect < 0
    assert es.direction == "down"


def test_random_effects_homogeneous_studies():
    pooled = random_effects([0.5, 0.5, 0.5, 0.5], [0.1, 0.1, 0.1, 0.1])
    assert pooled.pooled_effect == pytest.approx(0.5, abs=1e-6)
    assert pooled.i2 == pytest.approx(0.0, abs=1e-6)
    assert pooled.tau2 == pytest.approx(0.0, abs=1e-9)
    assert pooled.ci_low < 0.5 < pooled.ci_high
    assert pooled.n_studies == 4


def test_random_effects_detects_heterogeneity():
    pooled = random_effects([0.1, 0.9, 0.5, 0.0], [0.1, 0.1, 0.1, 0.1])
    assert pooled.tau2 > 0
    assert pooled.i2 > 50  # substantial heterogeneity


def test_random_effects_single_study_passthrough():
    pooled = random_effects([0.7], [0.2])
    assert pooled.pooled_effect == pytest.approx(0.7)
    assert pooled.standard_error == pytest.approx(0.2)
    assert pooled.i2 == 0.0


def test_random_effects_rejects_nonpositive_se():
    with pytest.raises(ValueError):
        random_effects([0.5, 0.5], [0.1, 0.0])


# --- Hartung-Knapp and the prediction interval -------------------------------
#
# The DerSimonian-Laird interval treats tau-squared as if it were known rather
# than estimated, which makes it too narrow whenever k is small. A quarter of
# this corpus pools fewer than ten contrasts, and the product's entire claim is
# that its intervals are honest, so an interval that is systematically too
# narrow is the worst available failure mode.


def test_hartung_knapp_collapses_to_a_one_sample_t_interval():
    """The property that makes the implementation checkable without restating it.

    With equal within-study variances and no between-study variance, the
    Hartung-Knapp standard error reduces exactly to the ordinary standard error
    of the mean of the study estimates, and the interval to Student's t on k-1
    degrees of freedom. Nothing about the meta-analytic machinery survives that
    reduction, so matching it is evidence the machinery is right.
    """
    effects = [0.20, 0.55, 0.35, 0.62, 0.41]
    pooled = random_effects(effects, [0.30] * 5)

    y = np.array(effects)
    expected_se = y.std(ddof=1) / np.sqrt(len(y))
    t_crit = float(stats.t.ppf(0.975, len(y) - 1))

    assert pooled.se_hk == pytest.approx(expected_se, rel=1e-9)
    assert pooled.ci_low_hk == pytest.approx(y.mean() - t_crit * expected_se, rel=1e-9)
    assert pooled.ci_high_hk == pytest.approx(y.mean() + t_crit * expected_se, rel=1e-9)


def test_the_reported_interval_is_never_narrower_than_the_dersimonian_laird_one():
    """Hartung-Knapp can be anti-conservative when the studies agree unusually
    well — the correction divides by an observed dispersion that happens to be
    small. The modified rule (Roever et al. 2015) takes the wider of the two, and
    this asserts it over a range that includes both regimes.
    """
    cases = [
        ([0.5, 0.5, 0.5, 0.5], [0.1] * 4),  # perfectly homogeneous: HK is narrower
        ([0.1, 0.9, 0.5, 0.0], [0.1] * 4),  # heterogeneous: HK is wider
        ([0.3, 0.31, 0.29, 0.30, 0.30], [0.2] * 5),
        ([-1.2, 0.4, 0.9], [0.3, 0.5, 0.2]),
    ]
    narrower_somewhere = False
    for effects, ses in cases:
        pooled = random_effects(effects, ses)
        dl_width = pooled.ci_high_dl - pooled.ci_low_dl
        reported = pooled.ci_high - pooled.ci_low
        hk_width = pooled.ci_high_hk - pooled.ci_low_hk
        assert reported >= dl_width - 1e-12
        assert reported == pytest.approx(max(dl_width, hk_width), rel=1e-12)
        narrower_somewhere |= hk_width < dl_width
    # If Hartung-Knapp were never narrower the guard would be dead code.
    assert narrower_somewhere


def test_the_prediction_interval_covers_a_new_study_and_the_confidence_interval_does_not():
    """A confidence interval is about the mean effect; a prediction interval is
    about the next study. Under real heterogeneity they answer different
    questions, and conflating them is the most common misreading of a forest
    plot — so this is a coverage simulation, not an algebraic identity.
    """
    rng = np.random.default_rng(4)
    k, mu, tau = 8, 0.4, 0.35
    inside_ci = inside_pi = 0
    trials = 1500
    for _ in range(trials):
        se = rng.uniform(0.15, 0.35, k)
        truths = rng.normal(mu, tau, k)
        observed = rng.normal(truths, se)
        pooled = random_effects(observed.tolist(), se.tolist())
        # One more study from the same distribution, observed without error, so
        # the target is the study's true effect rather than its estimate.
        new_study = rng.normal(mu, tau)
        inside_ci += pooled.ci_low <= new_study <= pooled.ci_high
        inside_pi += pooled.pi_low <= new_study <= pooled.pi_high

    assert 0.90 < inside_pi / trials < 0.99
    # The confidence interval is not a prediction interval and must visibly fail
    # at the job, or the distinction this repo draws is decorative.
    assert inside_ci / trials < 0.75


def test_the_prediction_interval_needs_three_studies():
    """With two studies the between-study variance has one degree of freedom and
    a prediction interval is not estimable. Returning a number anyway would be
    the same error as pooling two contrasts and calling the result a summary."""
    two = random_effects([0.4, 0.8], [0.2, 0.2])
    assert two.pi_low is None and two.pi_high is None
    assert two.ci_low_hk is not None  # Hartung-Knapp is defined from k=2

    one = random_effects([0.7], [0.2])
    assert one.pi_low is None and one.se_hk is None
    assert one.ci_low == pytest.approx(one.ci_low_dl)


def test_widening_the_interval_can_only_retract_a_verdict_never_create_one():
    """The corpus-level consequence, stated as an invariant.

    Swapping in a wider interval must never turn a null into a claim. If it ever
    did, the wider interval would not be the more conservative one and the reason
    for adopting it would be gone.
    """
    rng = np.random.default_rng(11)
    for _ in range(400):
        k = int(rng.integers(3, 16))
        se = rng.uniform(0.1, 0.6, k)
        effects = rng.normal(rng.uniform(-0.6, 0.6), 0.4, k)
        pooled = random_effects(effects.tolist(), se.tolist())
        dl_claims = pooled.ci_low_dl > 0 or pooled.ci_high_dl < 0
        reported_claims = pooled.ci_low > 0 or pooled.ci_high < 0
        assert not (reported_claims and not dl_claims)


def test_batch_correction_reduces_batch_signal_and_preserves_biology():
    rng = np.random.default_rng(3)
    n_features, per = 30, 10
    # 40 samples: 2 biological groups x 2 batches.
    group = np.array(["young"] * (2 * per) + ["old"] * (2 * per))
    batch = np.array((["A"] * per + ["B"] * per) * 2)
    base = rng.normal(0, 1, (n_features, 40))
    group_effect = np.where(group == "old", 3.0, 0.0)  # real biology
    batch_effect = np.where(batch == "B", 5.0, 0.0)  # nuisance
    X = base + group_effect + batch_effect

    before = batch_variance(X, batch)
    corrected = remove_batch_effect(X, batch, biological_group=group)
    after = batch_variance(corrected, batch)
    assert after < before * 0.2  # batch signal largely removed

    # Biological difference (old - young) is preserved.
    diff_before = X[:, group == "old"].mean() - X[:, group == "young"].mean()
    diff_after = corrected[:, group == "old"].mean() - corrected[:, group == "young"].mean()
    assert diff_after == pytest.approx(diff_before, abs=0.3)
