"""Critical-slowing-down (CSD) indicators of resilience loss.

Two classic early-warning signals of an approaching tipping point (Scheffer et
al.), applied across age strata as the cross-sectional proxy used in aging:

  * **rising variance** of the health state, and
  * **rising cross-correlation** among fluctuating markers (the system's
    fluctuations become lower-dimensional / more coupled as recovery slows).

When true longitudinal data is present, lag-1 temporal autocorrelation
(:mod:`recovery`) is the sharper signal; here we compute the age-stratified
proxy and label it as such so nothing is over-claimed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats


@dataclass
class CSDResult:
    strata_midpoints: list[float]
    n_per_stratum: list[int]
    variance: list[float]
    cross_correlation: list[float]
    variance_trend_slope: float
    variance_trend_tau: float
    variance_trend_p: float
    crosscorr_trend_slope: float
    crosscorr_trend_tau: float
    crosscorr_trend_p: float
    method: str
    fallback_used: bool
    assumptions: list[str] = field(default_factory=list)

    @property
    def resilience_declines(self) -> bool:
        """Both early-warning signals point toward loss of resilience with age."""
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


def _trend(x: list[float], y: list[float]) -> tuple[float, float, float]:
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


def csd_indicators(
    values: np.ndarray,
    ages: np.ndarray,
    n_strata: int = 6,
    longitudinal: bool = False,
) -> CSDResult:
    """Compute age-stratified CSD indicators.

    Parameters
    ----------
    values: (n_samples, n_markers) biomarker matrix.
    ages: (n_samples,) age per sample.
    n_strata: number of equal-count age bins.
    longitudinal: set True only if `values` are repeated measures over time for
        one system; controls the method label (temporal autocorrelation is then
        the recommended companion metric).
    """
    X = np.asarray(values, dtype=float)
    ages = np.asarray(ages, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    if X.shape[0] != ages.shape[0]:
        raise ValueError("values and ages must have the same number of samples.")
    if X.shape[0] < n_strata * 3:
        raise ValueError(f"Need >= {n_strata * 3} samples for {n_strata} strata.")

    # Health "state" = mean z-scored marker per sample (a scalar order parameter).
    Z = _zscore_columns(X)
    state = Z.mean(axis=1)

    # Equal-count strata by age quantiles.
    order = np.argsort(ages, kind="mergesort")
    splits = np.array_split(order, n_strata)

    mids, ns, variances, crosscorrs = [], [], [], []
    for idx in splits:
        if len(idx) == 0:
            continue
        mids.append(float(np.mean(ages[idx])))
        ns.append(int(len(idx)))
        variances.append(float(np.var(state[idx], ddof=1)) if len(idx) > 1 else float("nan"))
        crosscorrs.append(_mean_abs_offdiag_corr(X[idx]))

    v_slope, v_tau, v_p = _trend(mids, variances)
    c_slope, c_tau, c_p = _trend(mids, crosscorrs)

    fallback = not longitudinal
    method = "age_stratified_csd_fallback" if fallback else "longitudinal_csd"
    assumptions = [
        "Cross-sectional age strata are treated as a proxy for a within-individual "
        "temporal trajectory; individual longitudinal recovery is not observed.",
        "Rising variance and rising cross-correlation are interpreted as loss-of-"
        "resilience early-warning signals (critical slowing down).",
        "Markers are z-scored per column; the health 'state' is their mean. Absolute "
        "magnitudes are not comparable across datasets.",
    ]
    if fallback:
        assumptions.append(
            "FALLBACK: no dense longitudinal data supplied — this is the age-stratified "
            "proxy, not a true temporal-autocorrelation CSD measurement."
        )
    return CSDResult(
        strata_midpoints=mids,
        n_per_stratum=ns,
        variance=variances,
        cross_correlation=crosscorrs,
        variance_trend_slope=v_slope,
        variance_trend_tau=v_tau,
        variance_trend_p=v_p,
        crosscorr_trend_slope=c_slope,
        crosscorr_trend_tau=c_tau,
        crosscorr_trend_p=c_p,
        method=method,
        fallback_used=fallback,
        assumptions=assumptions,
    )
