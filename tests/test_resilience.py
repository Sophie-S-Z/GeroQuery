"""M6 resilience — CSD and recovery metrics on engineered ground truth."""

import numpy as np
import pandas as pd
import pytest

from geroquery.exceptions import ResilienceInputError
from geroquery.resilience import (
    ResilienceService,
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


def test_service_csd_missing_columns_raises():
    df = pd.DataFrame({"age": range(30), "x": range(30)})
    with pytest.raises(ResilienceInputError):
        ResilienceService().csd(df, ["x", "not_here"], age_col="age", n_strata=3)


# ---- service: log transform -------------------------------------------------


def _skewed_frame(seed=11):
    rng = np.random.default_rng(seed)
    n = 300
    age = rng.uniform(20, 80, n)
    return pd.DataFrame(
        {
            "age": age,
            "crp": rng.lognormal(0.0, 1.0, n),
            "albumin": rng.normal(4.2, 0.3, n),
        }
    )


def test_service_log_transform_changes_the_estimate():
    """Right-skewed markers let a few extreme values dominate a stratum's
    variance, making the trend a story about outliers rather than the cohort."""
    df = _skewed_frame()
    raw = ResilienceService().csd(df, ["crp", "albumin"], n_strata=6, n_bootstrap=100)
    logged = ResilienceService().csd(
        df, ["crp", "albumin"], n_strata=6, n_bootstrap=100, log_columns=["crp"]
    )
    assert raw.variance_evidence.slope != logged.variance_evidence.slope


def test_service_log_transform_rejects_non_positive_values():
    df = _skewed_frame()
    df.loc[0, "crp"] = 0.0
    with pytest.raises(ResilienceInputError, match="non-positive"):
        ResilienceService().csd(df, ["crp", "albumin"], n_strata=6, log_columns=["crp"])


def test_service_log_transform_rejects_unknown_columns():
    with pytest.raises(ResilienceInputError, match="log_columns"):
        ResilienceService().csd(
            _skewed_frame(), ["crp", "albumin"], n_strata=6, log_columns=["nope"]
        )


def test_service_log_transform_does_not_mutate_the_callers_frame():
    df = _skewed_frame()
    before = df["crp"].copy()
    ResilienceService().csd(df, ["crp", "albumin"], n_strata=6, n_bootstrap=50, log_columns=["crp"])
    pd.testing.assert_series_equal(df["crp"], before)


def test_service_detrend_flag_is_plumbed_through():
    df = _skewed_frame()
    on = ResilienceService().csd(df, ["crp", "albumin"], n_strata=6, n_bootstrap=50)
    off = ResilienceService().csd(df, ["crp", "albumin"], n_strata=6, n_bootstrap=50, detrend=False)
    assert on.detrended is True and off.detrended is False
