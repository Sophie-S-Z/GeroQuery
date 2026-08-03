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


# ---- service: control energy ------------------------------------------------


def test_service_control_energy_returns_diagnostics_not_a_bare_number():
    """Regression: this method referenced an undefined name and raised NameError
    on every call, because nothing exercised it. On a near-uncontrollable system
    the energy is numerically huge and meaningless, so the caller needs the
    conditioning diagnostics alongside it, not just the number."""
    A = np.array([[-1.0, 0.0], [0.0, -2.0]])
    B = np.eye(2)
    out = ResilienceService().control_energy(A, B, np.zeros(2), np.array([1.0, 0.0]))

    assert out["control_energy"] > 0
    assert out["horizon"] == 1.0
    for key in ("gramian_condition_number", "gramian_rank", "dimension", "well_conditioned"):
        assert key in out
    assert len(out["assumptions"]) >= 4


def test_service_control_energy_warns_when_the_gramian_is_ill_conditioned():
    """An input direction with no actuation is unreachable; the response must
    say the energy is a lower bound rather than present it as calibrated."""
    A = np.array([[-1.0, 0.0], [0.0, -1.0]])
    B = np.array([[1.0], [0.0]])  # second state has no input at all
    service = ResilienceService()
    # Default is strict: an ill-conditioned system is refused outright.
    with pytest.raises(ResilienceInputError, match="ill-conditioned"):
        service.control_energy(A, B, np.zeros(2), np.array([0.0, 1.0]))

    # strict=False returns the estimate, flagged.
    out = service.control_energy(A, B, np.zeros(2), np.array([0.0, 1.0]), strict=False)
    assert out["unreachable_fraction"] > 0.5
    assert out["well_conditioned"] is False
    assert any("WARNING" in a for a in out["assumptions"])


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
