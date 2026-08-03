"""Critical-slowing-down (CSD) indicators of resilience loss.

Two classic early-warning signals of an approaching tipping point (Scheffer et
al.), applied across age strata as the cross-sectional proxy used in aging:

  * **rising variance** of the health state, and
  * **rising cross-correlation** among fluctuating markers (the system's
    fluctuations become lower-dimensional / more coupled as recovery slows).

When true longitudinal data is present, lag-1 temporal autocorrelation
(:mod:`recovery`) is the sharper signal; here we compute the age-stratified
proxy and label it as such so nothing is over-claimed.

Two properties matter for this to be a measurement rather than an artifact:

**Within-stratum age detrending.** Markers are z-scored globally, then the
variance of the mean-z state is computed *within* each age bin. If markers drift
monotonically with age, the age heterogeneity inside a bin inflates that bin's
variance — manufacturing a rising-variance trend with no critical slowing down
present at all. We therefore regress age out within each stratum before
computing variance and cross-correlation, and report the undetrended series
alongside so the difference is visible rather than assumed away.

**Bootstrap confidence intervals over subjects.** A Kendall tau over ~6 strata
has a p-value floor near 0.003 and no power; a bare ``slope > 0`` test fires on
pure noise roughly a quarter of the time. Trend significance is therefore
established by resampling *subjects* (the actual sampling unit) and taking a
percentile interval on the slope. ``resilience_declines`` requires both
intervals to exclude zero — it is a claim about evidence, not about sign.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import stats

DEFAULT_N_BOOTSTRAP = 500
DEFAULT_RANDOM_STATE = 20240117


@dataclass
class TrendEvidence:
    """Slope of one early-warning indicator across age strata, with evidence."""

    indicator: str
    values: list[float]  # per-stratum indicator values (detrended)
    values_undetrended: list[float]  # same, without within-stratum age regression
    slope: float
    slope_undetrended: float
    tau: float
    tau_p: float
    ci_low: float | None
    ci_high: float | None
    n_bootstrap: int
    #: True only when the bootstrap CI lies entirely above zero.
    supported: bool
    verdict: str

    def to_dict(self) -> dict:
        return asdict(self)


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
    variance_evidence: TrendEvidence
    crosscorr_evidence: TrendEvidence
    detrended: bool
    n_bootstrap: int
    n_samples: int
    verdict: str
    assumptions: list[str] = field(default_factory=list)

    @property
    def resilience_declines(self) -> bool:
        """Both early-warning signals show *evidence* of rising with age.

        Requires each bootstrap CI to exclude zero. A positive slope alone is not
        sufficient — with a handful of strata that happens on noise routinely.
        """
        return self.variance_evidence.supported and self.crosscorr_evidence.supported

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["variance_evidence"] = self.variance_evidence.to_dict()
        d["crosscorr_evidence"] = self.crosscorr_evidence.to_dict()
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


def _residualize_on_age(y: np.ndarray, ages: np.ndarray) -> np.ndarray:
    """Remove the linear age trend from ``y`` within a stratum.

    Falls back to mean-centering when age has no spread inside the bin (a single
    age value, or a topcoded bin), because a slope is not identifiable there.
    """
    squeeze = y.ndim == 1
    if squeeze:
        y = y[:, None]
    if len(ages) < 3 or float(np.ptp(ages)) == 0.0:
        out = y - y.mean(axis=0, keepdims=True)
    else:
        design = np.column_stack([np.ones_like(ages), ages])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        out = y - design @ coef
    return out[:, 0] if squeeze else out


def _trend(x: Sequence[float], y: Sequence[float]) -> tuple[float, float, float]:
    """OLS slope plus Kendall tau (robust monotonic-trend test).

    The tau p-value is reported for continuity but is not what gates any claim —
    with a small number of strata it has essentially no power. See the bootstrap.
    """
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    mask = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[mask], ys[mask]
    if len(xs) < 3:
        return float("nan"), float("nan"), float("nan")
    slope = float(np.polyfit(xs, ys, 1)[0])
    tau, p = stats.kendalltau(xs, ys)
    return slope, float(tau), float(p)


def _stratum_indicators(
    X: np.ndarray, ages: np.ndarray, n_strata: int, detrend: bool
) -> tuple[list[float], list[int], list[float], list[float]]:
    """Per-stratum midpoint, count, state variance, and mean |cross-correlation|."""
    Z = _zscore_columns(X)
    state = Z.mean(axis=1)
    order = np.argsort(ages, kind="mergesort")

    mids: list[float] = []
    ns: list[int] = []
    variances: list[float] = []
    crosscorrs: list[float] = []
    for idx in np.array_split(order, n_strata):
        if len(idx) == 0:
            continue
        a = ages[idx]
        s = state[idx]
        Xi = X[idx]
        if detrend:
            s = _residualize_on_age(s, a)
            Xi = _residualize_on_age(Xi, a)
        mids.append(float(np.mean(a)))
        ns.append(int(len(idx)))
        variances.append(float(np.var(s, ddof=1)) if len(idx) > 1 else float("nan"))
        crosscorrs.append(_mean_abs_offdiag_corr(Xi))
    return mids, ns, variances, crosscorrs


def _bootstrap_slopes(
    X: np.ndarray,
    ages: np.ndarray,
    n_strata: int,
    detrend: bool,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Resample SUBJECTS with replacement and recompute both trend slopes."""
    n = len(ages)
    var_slopes, cc_slopes = [], []
    for _ in range(n_bootstrap):
        pick = rng.integers(0, n, n)
        try:
            mids, _ns, var, cc = _stratum_indicators(X[pick], ages[pick], n_strata, detrend)
        except (np.linalg.LinAlgError, ValueError):
            continue
        v_slope, _t, _p = _trend(mids, var)
        c_slope, _t2, _p2 = _trend(mids, cc)
        if np.isfinite(v_slope):
            var_slopes.append(v_slope)
        if np.isfinite(c_slope):
            cc_slopes.append(c_slope)
    return np.asarray(var_slopes), np.asarray(cc_slopes)


def _evidence(
    indicator: str,
    values: list[float],
    values_undetrended: list[float],
    mids: list[float],
    boot: np.ndarray,
) -> TrendEvidence:
    slope, tau, tau_p = _trend(mids, values)
    slope_undet, _t, _p = _trend(mids, values_undetrended)

    ci_low: float | None
    ci_high: float | None
    if boot.size >= 20:
        ci_low = float(np.percentile(boot, 2.5))
        ci_high = float(np.percentile(boot, 97.5))
        supported = ci_low > 0.0
        ci_text = f"bootstrap 95% CI {ci_low:+.5g} to {ci_high:+.5g}"
        if supported:
            verdict = f"rises with age ({ci_text}, excludes 0)"
        elif ci_high < 0.0:
            verdict = f"falls with age ({ci_text}, excludes 0)"
        else:
            verdict = f"no evidence of a trend ({ci_text} includes 0; point slope {slope:+.5g})"
    else:
        ci_low = ci_high = None
        supported = False
        verdict = "not assessed (bootstrap disabled or too few usable resamples)"

    return TrendEvidence(
        indicator=indicator,
        values=values,
        values_undetrended=values_undetrended,
        slope=slope,
        slope_undetrended=slope_undet,
        tau=tau,
        tau_p=tau_p,
        ci_low=ci_low,
        ci_high=ci_high,
        n_bootstrap=int(boot.size),
        supported=supported,
        verdict=verdict,
    )


def csd_indicators(
    values: np.ndarray,
    ages: np.ndarray,
    n_strata: int = 6,
    longitudinal: bool = False,
    detrend: bool = True,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    random_state: int = DEFAULT_RANDOM_STATE,
) -> CSDResult:
    """Compute age-stratified CSD indicators with trend evidence.

    Parameters
    ----------
    values: (n_samples, n_markers) biomarker matrix.
    ages: (n_samples,) age per sample.
    n_strata: number of equal-count age bins.
    longitudinal: set True only if `values` are repeated measures over time for
        one system; controls the method label (temporal autocorrelation is then
        the recommended companion metric).
    detrend: regress age out within each stratum before computing indicators.
        Leave on unless you specifically want the confounded comparison; the
        undetrended series is returned either way.
    n_bootstrap: subject-level resamples used for the trend CIs. Set 0 to skip,
        in which case no trend is reported as supported.
    random_state: seed, so a reported CI is reproducible.
    """
    X = np.asarray(values, dtype=float)
    ages = np.asarray(ages, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    if X.shape[0] != ages.shape[0]:
        raise ValueError("values and ages must have the same number of samples.")
    if X.shape[0] < n_strata * 3:
        raise ValueError(f"Need >= {n_strata * 3} samples for {n_strata} strata.")

    mids, ns, variances, crosscorrs = _stratum_indicators(X, ages, n_strata, detrend)
    _m, _n, var_undet, cc_undet = _stratum_indicators(X, ages, n_strata, detrend=False)

    rng = np.random.default_rng(random_state)
    if n_bootstrap > 0:
        var_boot, cc_boot = _bootstrap_slopes(X, ages, n_strata, detrend, n_bootstrap, rng)
    else:
        var_boot = cc_boot = np.asarray([])

    var_ev = _evidence("variance", variances, var_undet, mids, var_boot)
    cc_ev = _evidence("cross_correlation", crosscorrs, cc_undet, mids, cc_boot)

    fallback = not longitudinal
    method = "age_stratified_csd_fallback" if fallback else "longitudinal_csd"

    if var_ev.supported and cc_ev.supported:
        verdict = "Both early-warning indicators rise with age: evidence of resilience loss."
    elif var_ev.supported or cc_ev.supported:
        rising, other = (var_ev, cc_ev) if var_ev.supported else (cc_ev, var_ev)
        verdict = (
            f"Mixed: {rising.indicator} rises with age, but {other.indicator} shows "
            f"{other.verdict.split('(')[0].strip()}. A single indicator is not sufficient "
            "evidence of critical slowing down."
        )
    else:
        verdict = "No evidence of rising early-warning indicators with age."

    assumptions = [
        "Cross-sectional age strata are treated as a proxy for a within-individual "
        "temporal trajectory; individual longitudinal recovery is not observed.",
        "Rising variance and rising cross-correlation are interpreted as loss-of-"
        "resilience early-warning signals (critical slowing down).",
        "Markers are z-scored per column; the health 'state' is their mean. Absolute "
        "magnitudes are not comparable across datasets.",
        "Trend significance comes from resampling subjects, not strata. The Kendall tau "
        "p-value is reported but has almost no power at this number of strata and gates "
        "no claim.",
    ]
    if detrend:
        assumptions.append(
            "Age is regressed out linearly within each stratum before computing variance "
            "and cross-correlation, so within-bin mean drift cannot masquerade as critical "
            "slowing down. Undetrended values are reported alongside for comparison."
        )
    else:
        assumptions.append(
            "WARNING: detrending is disabled. Within-stratum age drift can inflate variance "
            "and mimic critical slowing down. Results are not confound-controlled."
        )
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
        variance_trend_slope=var_ev.slope,
        variance_trend_tau=var_ev.tau,
        variance_trend_p=var_ev.tau_p,
        crosscorr_trend_slope=cc_ev.slope,
        crosscorr_trend_tau=cc_ev.tau,
        crosscorr_trend_p=cc_ev.tau_p,
        method=method,
        fallback_used=fallback,
        variance_evidence=var_ev,
        crosscorr_evidence=cc_ev,
        detrended=detrend,
        n_bootstrap=n_bootstrap,
        n_samples=int(X.shape[0]),
        verdict=verdict,
        assumptions=assumptions,
    )
