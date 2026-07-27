"""Age-stratified dispersion & co-fluctuation indicators of resilience.

As a body loses resilience its biomarker state tends to become **more variable**
and its markers **more coupled**. GeroQuery measures three such indicators across
age strata:

  * **dispersion** — the variance of a scalar health state,
  * **co-fluctuation** — the mean absolute cross-correlation among markers, and
  * a **Dynamic Network Biomarker (DNB)** composite (:mod:`dnb`).

These are *inspired by* critical-slowing-down (CSD) early-warning-signal theory
(Scheffer et al., 2009) but this is a **cross-sectional proxy computed across age
strata, not validated critical slowing down.** Early-warning signals are
contested: they produce systematic false positives and are sensitive to
detrending, window, and bandwidth choices (Boettiger & Hastings, 2012; Dakos et
al., 2012). The sharper, defensible signal — lag-1 temporal autocorrelation /
Ornstein-Uhlenbeck relaxation — requires within-individual longitudinal data and
lives in :mod:`recovery`. Everything here is labelled as a proxy and shipped with
bootstrap confidence intervals so trends are not over-read.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from .dnb import auto_module, dnb_per_block


@dataclass
class CSDResult:
    strata_midpoints: list[float]
    n_per_stratum: list[int]
    variance: list[float]
    cross_correlation: list[float]
    dnb: list[float]
    variance_trend_slope: float
    variance_trend_tau: float
    variance_trend_p: float
    crosscorr_trend_slope: float
    crosscorr_trend_tau: float
    crosscorr_trend_p: float
    dnb_trend_slope: float
    dnb_trend_tau: float
    dnb_trend_p: float
    method: str
    fallback_used: bool
    # Bootstrap 95% CIs on the trend slopes ([low, high]); None if not computed.
    variance_trend_ci: list[float] | None = None
    crosscorr_trend_ci: list[float] | None = None
    dnb_trend_ci: list[float] | None = None
    n_bootstrap: int = 0
    assumptions: list[str] = field(default_factory=list)

    @property
    def resilience_declines(self) -> bool:
        """Both classic early-warning signals point toward loss of resilience."""
        return self.variance_trend_slope > 0 and self.crosscorr_trend_slope > 0

    def to_dict(self) -> dict:
        d = self.__dict__.copy()
        d["resilience_declines"] = self.resilience_declines
        return d


def _zscore_columns(X: np.ndarray) -> np.ndarray:
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def _mean_abs_offdiag_corr(X: np.ndarray) -> float:
    if X.shape[0] < 3 or X.shape[1] < 2:
        return float("nan")
    c = np.corrcoef(X, rowvar=False)
    iu = np.triu_indices_from(c, k=1)
    vals = np.abs(c[iu])
    vals = vals[np.isfinite(vals)]
    return float(vals.mean()) if vals.size else float("nan")


def _trend(x, y) -> tuple[float, float, float]:
    """OLS slope plus Kendall tau (robust monotonic-trend test)."""
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    mask = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[mask], ys[mask]
    if len(xs) < 3:
        return float("nan"), float("nan"), float("nan")
    slope = float(np.polyfit(xs, ys, 1)[0])
    tau, p = stats.kendalltau(xs, ys)
    return slope, float(tau), float(p)


def _series(X: np.ndarray, ages: np.ndarray, n_strata: int, partition):
    """Per-stratum midpoints, counts, and the three indicators."""
    Z = _zscore_columns(X)
    state = Z.mean(axis=1)
    order = np.argsort(ages, kind="mergesort")
    splits = [idx for idx in np.array_split(order, n_strata) if len(idx) > 0]

    mids = [float(np.mean(ages[idx])) for idx in splits]
    ns = [int(len(idx)) for idx in splits]
    variances = [
        float(np.var(state[idx], ddof=1)) if len(idx) > 1 else float("nan") for idx in splits
    ]
    crosscorrs = [_mean_abs_offdiag_corr(X[idx]) for idx in splits]
    dnbs = dnb_per_block(X, splits, partition)
    return mids, ns, variances, crosscorrs, dnbs


def _slopes(X: np.ndarray, ages: np.ndarray, n_strata: int, partition):
    mids, _, v, c, d = _series(X, ages, n_strata, partition)
    return _trend(mids, v)[0], _trend(mids, c)[0], _trend(mids, d)[0]


def _bootstrap_cis(X, ages, n_strata, partition, n_boot, seed):
    """Percentile 95% CIs on the three trend slopes via subject resampling."""
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    vs, cs, ds = [], [], []
    for _ in range(n_boot):
        samp = rng.integers(0, n, n)
        try:
            sv, sc, sd = _slopes(X[samp], ages[samp], n_strata, partition)
        except Exception:  # noqa: BLE001 — a degenerate resample just gets skipped
            continue
        if np.isfinite(sv):
            vs.append(sv)
        if np.isfinite(sc):
            cs.append(sc)
        if np.isfinite(sd):
            ds.append(sd)

    def ci(a):
        return (
            [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))] if len(a) >= 10 else None
        )

    return ci(vs), ci(cs), ci(ds)


def csd_indicators(
    values: np.ndarray,
    ages: np.ndarray,
    n_strata: int = 6,
    longitudinal: bool = False,
    n_bootstrap: int = 200,
    seed: int = 0,
) -> CSDResult:
    """Compute age-stratified dispersion, co-fluctuation, and DNB indicators.

    Parameters
    ----------
    values: (n_samples, n_markers) biomarker matrix.
    ages: (n_samples,) age per sample.
    n_strata: number of equal-count age bins.
    longitudinal: set True only if `values` are repeated measures over time for
        one system; controls the method label (temporal autocorrelation via
        :mod:`recovery` is then the recommended companion metric).
    n_bootstrap: subject-resampling replicates for trend-slope CIs (0 to skip).
    seed: RNG seed for the bootstrap.
    """
    X = np.asarray(values, dtype=float)
    ages = np.asarray(ages, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    if X.shape[0] != ages.shape[0]:
        raise ValueError("values and ages must have the same number of samples.")
    if X.shape[0] < n_strata * 3:
        raise ValueError(f"Need >= {n_strata * 3} samples for {n_strata} strata.")

    # Fix the DNB module/background partition once, on the full cohort, so the
    # index is comparable across strata and bootstrap resamples.
    partition = auto_module(X)

    mids, ns, variances, crosscorrs, dnbs = _series(X, ages, n_strata, partition)
    v_slope, v_tau, v_p = _trend(mids, variances)
    c_slope, c_tau, c_p = _trend(mids, crosscorrs)
    d_slope, d_tau, d_p = _trend(mids, dnbs)

    v_ci = c_ci = d_ci = None
    if n_bootstrap and n_bootstrap > 0:
        v_ci, c_ci, d_ci = _bootstrap_cis(X, ages, n_strata, partition, n_bootstrap, seed)

    fallback = not longitudinal
    method = "age_stratified_csd_fallback" if fallback else "longitudinal_csd"
    dnb_mode = (
        "full three-term index"
        if partition is not None
        else "no-background reduction (SD_in x PCC_in)"
    )
    assumptions = [
        "Cross-sectional age strata are treated as a proxy for a within-individual "
        "temporal trajectory; individual longitudinal recovery is not observed.",
        "Rising dispersion (variance) and rising co-fluctuation (cross-correlation) "
        "are interpreted as loss-of-resilience indicators, inspired by — but not a "
        "validated measurement of — critical slowing down.",
        "Early-warning signals are contested: they produce systematic false positives "
        "and are sensitive to detrending/window choices (Boettiger & Hastings 2012; "
        "Dakos et al. 2012). Read trends alongside the bootstrap confidence intervals.",
        f"DNB composite computed as: {dnb_mode}.",
        "Markers are z-scored per column; the health 'state' is their mean. Absolute "
        "magnitudes are not comparable across datasets.",
    ]
    if fallback:
        assumptions.append(
            "FALLBACK: no dense longitudinal data supplied — this is the age-stratified "
            "proxy, not a true temporal-autocorrelation measurement (use the recovery-rate "
            "metric on within-individual time series for that)."
        )
    return CSDResult(
        strata_midpoints=mids,
        n_per_stratum=ns,
        variance=variances,
        cross_correlation=crosscorrs,
        dnb=dnbs,
        variance_trend_slope=v_slope,
        variance_trend_tau=v_tau,
        variance_trend_p=v_p,
        crosscorr_trend_slope=c_slope,
        crosscorr_trend_tau=c_tau,
        crosscorr_trend_p=c_p,
        dnb_trend_slope=d_slope,
        dnb_trend_tau=d_tau,
        dnb_trend_p=d_p,
        method=method,
        fallback_used=fallback,
        variance_trend_ci=v_ci,
        crosscorr_trend_ci=c_ci,
        dnb_trend_ci=d_ci,
        n_bootstrap=int(n_bootstrap or 0),
        assumptions=assumptions,
    )
