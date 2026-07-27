"""M6 resilience — CSD, recovery, control energy on engineered ground truth."""

import numpy as np
import pandas as pd
import pytest

from geroquery.exceptions import ResilienceInputError
from geroquery.resilience import (
    ResilienceService,
    control_energy,
    csd_indicators,
    recovery_rate,
)


def _tipping_point_dataset(seed=0, n=600):
    """Shared aging factor whose variance grows with age -> rising variance and
    rising cross-correlation among markers (critical slowing down)."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        age = rng.integers(20, 86)
        aging = (age - 20) / 65.0
        common = rng.normal(0, 0.4 + 1.8 * aging)
        markers = [0.5 * (i + 1) * common + rng.normal(0, 0.5) for i in range(5)]
        rows.append([age, *markers])
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1:]


def test_csd_indicators_rise_toward_tipping_point():
    ages, values = _tipping_point_dataset()
    res = csd_indicators(values, ages, n_strata=6)
    assert res.variance_trend_slope > 0
    assert res.crosscorr_trend_slope > 0
    assert res.resilience_declines
    assert res.variance_trend_tau > 0  # monotone increase


def test_csd_labels_fallback_when_not_longitudinal():
    ages, values = _tipping_point_dataset()
    res = csd_indicators(values, ages, n_strata=5, longitudinal=False)
    assert res.fallback_used
    assert "FALLBACK" in " ".join(res.assumptions)
    assert res.method == "age_stratified_csd_fallback"


def test_csd_needs_enough_samples():
    with pytest.raises(ValueError):
        csd_indicators(np.zeros((5, 3)), np.arange(5), n_strata=6)


def test_recovery_rate_higher_for_resilient_system():
    def ar1(a, n=500, seed=1):
        rng = np.random.default_rng(seed)
        x = np.zeros(n)
        for t in range(1, n):
            x[t] = a * x[t - 1] + rng.normal(0, 1)
        return x

    young = recovery_rate(ar1(0.3))  # resilient
    old = recovery_rate(ar1(0.9))  # frail
    assert young.recovery_rate > old.recovery_rate
    assert young.relaxation_time < old.relaxation_time


def test_recovery_rate_input_guards():
    with pytest.raises(ResilienceInputError):
        recovery_rate([1.0, 2.0])  # too short
    with pytest.raises(ResilienceInputError):
        recovery_rate([2.0, 2.0, 2.0, 2.0, 2.0])  # constant


def test_control_energy_increases_with_distance():
    A = np.array([[-1.0, 0.2], [0.0, -1.5]])
    B = np.eye(2)
    near = control_energy(A, B, np.zeros(2), np.array([1.0, 0.0]))
    far = control_energy(A, B, np.zeros(2), np.array([3.0, 0.0]))
    assert 0 < near < far


def test_service_csd_missing_columns_raises():
    df = pd.DataFrame({"age": range(30), "x": range(30)})
    with pytest.raises(ResilienceInputError):
        ResilienceService().csd(df, ["x", "not_here"], age_col="age", n_strata=3)


def test_csd_reports_dnb_and_bootstrap_cis():
    from geroquery.resilience import csd_indicators

    ages, values = _tipping_point_dataset()
    res = csd_indicators(values, ages, n_strata=6, n_bootstrap=100, seed=0)
    # DNB series + trend are present.
    assert len(res.dnb) == len(res.strata_midpoints)
    assert np.isfinite(res.dnb_trend_slope)
    # Bootstrap CIs are 2-element [low, high] on each trend slope.
    for c in (res.variance_trend_ci, res.crosscorr_trend_ci, res.dnb_trend_ci):
        assert c is None or (len(c) == 2 and c[0] <= c[1])
    assert res.n_bootstrap == 100
    # Reframing keeps the contested-EWS caveat in the assumptions.
    joined = " ".join(res.assumptions).lower()
    assert "false positive" in joined and "critical slowing down" in joined


def test_dnb_index_rises_toward_tipping_point():
    from geroquery.resilience import auto_module, dnb_index

    ages, values = _tipping_point_dataset()
    order = np.argsort(ages)
    young = values[order[:150]]
    old = values[order[-150:]]
    part = auto_module(values)
    assert dnb_index(old, part) > dnb_index(young, part)


def test_csd_bootstrap_can_be_disabled():
    from geroquery.resilience import csd_indicators

    ages, values = _tipping_point_dataset()
    res = csd_indicators(values, ages, n_strata=6, n_bootstrap=0)
    assert res.variance_trend_ci is None and res.n_bootstrap == 0
