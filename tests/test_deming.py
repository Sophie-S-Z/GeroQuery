"""Deming (orthogonal) regression for the cross-species scatter (§4.2).

The cross-species view asks whether a gene's age effect agrees between human and
mouse. Both axes are pooled Hedges' g, so both carry measurement error, and
ordinary least squares is biased toward a flat slope when the x variable is
noisy — it would systematically *understate* conservation and the understatement
would look like a result. Deming regression accounts for error in both
variables; with delta = 1 it is orthogonal/total-least-squares regression.

The tests are the §3 rule-3 contract: defining identities with known answers,
not "it runs".
"""

from __future__ import annotations

import numpy as np

from geroquery.harmonize.deming import deming_regression


def test_exact_recovery_on_a_perfectly_collinear_line():
    """Points on y = a + b·x recover (a, b) exactly, and the jackknife SE is ~0."""
    x = np.linspace(-2, 2, 40)
    y = 0.5 + 1.7 * x
    fit = deming_regression(x, y)
    assert np.isclose(fit.slope, 1.7, atol=1e-9)
    assert np.isclose(fit.intercept, 0.5, atol=1e-9)
    # No scatter about the line, so the leave-one-out slope never moves.
    assert fit.slope_se < 1e-6


def test_orthogonal_regression_is_symmetric():
    """With delta = 1, slope(y~x) · slope(x~y) = 1 — the property OLS lacks.

    This is exactly what makes Deming the right choice here: it does not privilege
    one axis. OLS of y on x and of x on y give slopes whose product is R² < 1;
    orthogonal regression gives reciprocals.
    """
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, 200)
    y = 0.8 * x + rng.normal(0, 0.4, 200)

    forward = deming_regression(x, y).slope
    backward = deming_regression(y, x).slope
    assert np.isclose(forward * backward, 1.0, atol=1e-9)


def test_slope_is_steeper_than_ols_when_x_is_noisy():
    """OLS attenuates toward zero under x-error; Deming corrects it.

    The whole reason §4.2 forbids an OLS line: with noise on the x-axis the OLS
    slope is pulled below the true 1.0, and drawing that would understate
    cross-species conservation. Deming must sit above OLS and nearer the truth.
    """
    rng = np.random.default_rng(1)
    truth = rng.normal(0, 1, 500)
    x = truth + rng.normal(0, 0.5, 500)  # both axes measure `truth` with error
    y = truth + rng.normal(0, 0.5, 500)

    ols = float(np.polyfit(x, y, 1)[0])
    deming = deming_regression(x, y).slope
    assert ols < deming  # OLS attenuated
    assert abs(deming - 1.0) < abs(ols - 1.0)  # Deming closer to the true slope 1


def test_confidence_interval_brackets_the_slope():
    x = np.linspace(-1, 1, 60)
    y = 0.3 + 0.9 * x + np.sin(np.arange(60))  # some structured scatter
    fit = deming_regression(x, y)
    assert fit.slope_ci_low < fit.slope < fit.slope_ci_high
