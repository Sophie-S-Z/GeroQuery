"""Deming (errors-in-both-variables) regression.

The cross-species scatter (§4.2 of ``VISUAL_AND_STATISTICAL_PLAN``) plots human
pooled effect against mouse pooled effect, one point per ortholog group. Both
axes are Hedges' g and both carry measurement error, so the fitted line must be
Deming/orthogonal, not ordinary least squares: OLS is biased toward a flat slope
when the predictor is noisy, which here would understate cross-species
conservation and make the understatement look like a result.

``delta`` is the ratio of the error variance in y to the error variance in x.
``delta = 1`` is orthogonal / total-least-squares regression, the default and
the right choice when both axes are on the same scale (both are Hedges' g).

Only the estimator lives here; the exporter computes the fit once over the
corpus and ships slope/intercept/CI in ``crossspecies.json`` so the browser
draws the line rather than fitting one (the standing "JS renders, never derives"
rule).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DemingFit:
    slope: float
    intercept: float
    slope_se: float
    slope_ci_low: float
    slope_ci_high: float
    n: int


def _slope_intercept(x: np.ndarray, y: np.ndarray, delta: float) -> tuple[float, float]:
    xbar = float(x.mean())
    ybar = float(y.mean())
    sxx = float(np.sum((x - xbar) ** 2))
    syy = float(np.sum((y - ybar) ** 2))
    sxy = float(np.sum((x - xbar) * (y - ybar)))
    if sxy == 0.0:
        # No linear relationship to orient the line; slope is undefined. Report
        # a flat line through the means rather than dividing by zero.
        return 0.0, ybar
    # Standard Deming closed form (Deming 1943); delta = Var(e_y)/Var(e_x).
    term = syy - delta * sxx
    slope = (term + np.sqrt(term**2 + 4.0 * delta * sxy**2)) / (2.0 * sxy)
    intercept = ybar - slope * xbar
    return float(slope), float(intercept)


def deming_regression(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    *,
    delta: float = 1.0,
    ci: float = 0.95,
) -> DemingFit:
    """Fit a Deming line with a jackknife interval for the slope.

    The jackknife (leave-one-out) SE is used rather than a closed-form variance
    because it needs no distributional assumption on the two error terms and is
    exact-zero on collinear data, which is the property the test pins.
    """
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    if xa.shape != ya.shape or xa.ndim != 1:
        raise ValueError("x and y must be 1-D and equal length.")
    n = len(xa)
    if n < 3:
        raise ValueError("Need at least three points to fit a line with an interval.")
    if delta <= 0:
        raise ValueError("delta (error-variance ratio) must be positive.")

    slope, intercept = _slope_intercept(xa, ya, delta)

    # Leave-one-out jackknife on the slope.
    loo = np.empty(n, dtype=float)
    for i in range(n):
        keep = np.arange(n) != i
        loo[i], _ = _slope_intercept(xa[keep], ya[keep], delta)
    mean_loo = float(loo.mean())
    # Jackknife variance of the estimator.
    var = (n - 1) / n * float(np.sum((loo - mean_loo) ** 2))
    se = float(np.sqrt(max(var, 0.0)))

    from scipy import stats

    t_crit = float(stats.t.ppf(1.0 - (1.0 - ci) / 2.0, n - 2))
    return DemingFit(
        slope=slope,
        intercept=intercept,
        slope_se=se,
        slope_ci_low=slope - t_crit * se,
        slope_ci_high=slope + t_crit * se,
        n=n,
    )
