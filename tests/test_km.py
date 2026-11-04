"""Weighted Kaplan-Meier and restricted mean survival time (§4.3, §3.4).

The mortality view reports hazard ratios and a concordance, but no survival
curve and no answer in years. RMST — the area under the survival curve to a
horizon — is model-free and says how many event-free years separate two groups,
which a hazard ratio cannot. These tests are the §3.4 contract: a closed form
the estimator must recover, and the defining property that unit weights collapse
the weighted estimator onto the unweighted one.
"""

from __future__ import annotations

import numpy as np

from geroquery.survival.km import rmst, weighted_kaplan_meier


def test_km_without_censoring_is_one_minus_the_ecdf():
    """With every subject observed, S(t) is just the survival ECDF."""
    rng = np.random.default_rng(0)
    time = rng.exponential(5.0, 2000)
    event = np.ones_like(time)

    grid = np.linspace(0, 10, 21)
    curve = weighted_kaplan_meier(time, event, grid=grid)
    empirical = np.array([np.mean(time > t) for t in grid])
    assert np.allclose(curve.survival, empirical, atol=0.02)


def test_rmst_recovers_the_exponential_closed_form():
    """Under S(t) = exp(-λt), RMST to τ is (1 − exp(−λτ)) / λ."""
    rng = np.random.default_rng(1)
    lam = 0.2
    time = rng.exponential(1.0 / lam, 50_000)
    event = np.ones_like(time)

    for tau in (5.0, 10.0):
        closed = (1.0 - np.exp(-lam * tau)) / lam
        got = rmst(time, event, tau=tau).value
        assert abs(got - closed) < 0.1, (tau, got, closed)


def test_unit_weights_match_the_unweighted_estimator_exactly():
    """The defining property the weighted Cox path is already tested on."""
    rng = np.random.default_rng(2)
    time = rng.exponential(4.0, 500)
    event = (rng.random(500) < 0.7).astype(float)

    a = rmst(time, event, tau=8.0)
    b = rmst(time, event, tau=8.0, weight=np.ones_like(time))
    assert np.isclose(a.value, b.value)

    grid = np.linspace(0, 8, 17)
    ca = weighted_kaplan_meier(time, event, grid=grid)
    cb = weighted_kaplan_meier(time, event, weight=np.ones_like(time), grid=grid)
    assert np.allclose(ca.survival, cb.survival)


def test_rmst_is_the_area_under_the_km_curve():
    """RMST and the KM curve are the same object seen two ways."""
    rng = np.random.default_rng(3)
    time = rng.exponential(3.0, 3000)
    event = (rng.random(3000) < 0.8).astype(float)
    tau = 6.0

    # Integrate the curve on a fine grid by the trapezoid of the step function.
    grid = np.linspace(0, tau, 2001)
    curve = weighted_kaplan_meier(time, event, grid=grid)
    area = np.trapezoid(curve.survival, grid)
    assert abs(area - rmst(time, event, tau=tau).value) < 0.05


def test_weighting_shifts_the_curve_toward_the_upweighted_subjects():
    """Upweighting the long survivors must raise the curve and the RMST."""
    time = np.array([1.0, 2.0, 3.0, 10.0, 11.0, 12.0])
    event = np.ones_like(time)
    heavy_long = np.array([1.0, 1.0, 1.0, 5.0, 5.0, 5.0])

    base = rmst(time, event, tau=12.0).value
    up = rmst(time, event, tau=12.0, weight=heavy_long).value
    assert up > base
