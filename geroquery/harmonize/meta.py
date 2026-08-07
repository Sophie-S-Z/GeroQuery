"""Random-effects meta-analysis (DerSimonian-Laird, Hartung-Knapp interval).

Pools per-study effect sizes into one signature *without* ever concatenating raw
expression matrices (which would confound platform and batch with biology). The
random-effects model additionally estimates between-study heterogeneity (tau^2,
I^2) so callers can see how consistent the studies actually are.

**Three intervals, because they answer three questions.**

``ci_low``/``ci_high``
    The reported interval for the mean effect. Hartung-Knapp-Sidik-Jonkman,
    widened to the DerSimonian-Laird interval whenever HK comes out narrower.

``ci_low_dl``/``ci_high_dl``
    The plain DerSimonian-Laird interval, kept because it is what almost every
    other tool reports and a reader comparing against one deserves to see it.

``pi_low``/``pi_high``
    The 95% prediction interval: where the true effect of *the next study* is
    expected to fall. Under real heterogeneity this is much wider than the
    confidence interval, and treating the two as interchangeable is the most
    common misreading of a forest plot.

Why the reported interval is not plain DerSimonian-Laird. DL substitutes the
estimated tau^2 into the weights and then proceeds as though it were known,
which makes the interval too narrow — badly so when k is small. Hartung-Knapp
replaces the normal quantile with Student's t on k-1 degrees of freedom and
rescales by the observed dispersion of the studies around the pooled mean;
simulation work (IntHout et al. 2014) finds it recovers nominal coverage where
DL does not. A quarter of this corpus pools fewer than ten contrasts, so this is
not a marginal correction here.

The one hazard is that HK can come out *narrower* than DL when the studies agree
unusually closely, purely because the observed dispersion happened to be small.
The modified rule (Roever, Knapp & Friede 2015) takes the wider of the two, and
that is what ships: **the reported interval is never narrower than the interval
every other tool would give.**
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class PooledEffect:
    pooled_effect: float
    standard_error: float
    ci_low: float
    ci_high: float
    p_value: float
    tau2: float  # between-study variance
    i2: float  # % of variance due to heterogeneity (0..100)
    q: float  # Cochran's Q
    n_studies: int
    # The DerSimonian-Laird interval, for comparison with other tools.
    ci_low_dl: float = 0.0
    ci_high_dl: float = 0.0
    # Hartung-Knapp. None when k == 1, where there is no dispersion to observe.
    se_hk: float | None = None
    ci_low_hk: float | None = None
    ci_high_hk: float | None = None
    # Prediction interval. None below k == 3, where it is not estimable.
    pi_low: float | None = None
    pi_high: float | None = None

    @property
    def direction(self) -> str:
        return "up" if self.pooled_effect >= 0 else "down"

    @property
    def verdict(self) -> str:
        """increases / decreases / no_evidence, from the reported interval.

        Computed here rather than at each call site so the API, the static
        export and the MCP payload cannot disagree about what a gene's answer is.
        """
        if self.ci_low > 0:
            return "increases"
        if self.ci_high < 0:
            return "decreases"
        return "no_evidence"


def random_effects(
    effects: Sequence[float], standard_errors: Sequence[float], ci: float = 0.95
) -> PooledEffect:
    """DerSimonian-Laird pooling with a Hartung-Knapp interval.

    Parameters
    ----------
    effects, standard_errors:
        Per-study effect sizes and their SEs (same length, >=1).
    ci:
        Confidence level for every interval returned.
    """
    y = np.asarray(effects, dtype=float)
    se = np.asarray(standard_errors, dtype=float)
    if y.shape != se.shape or y.ndim != 1:
        raise ValueError("effects and standard_errors must be 1-D and equal length.")
    k = len(y)
    if k == 0:
        raise ValueError("Need at least one study to pool.")
    if np.any(se <= 0):
        raise ValueError("All standard errors must be positive.")

    v = se**2  # within-study variances
    w_fixed = 1.0 / v
    alpha = 1.0 - ci

    if k == 1:
        eff, se_pool = float(y[0]), float(se[0])
        tau2 = q = i2 = 0.0
        w_random = w_fixed
    else:
        # Fixed-effect pooled mean, used to compute Q.
        mean_fixed = np.sum(w_fixed * y) / np.sum(w_fixed)
        q = float(np.sum(w_fixed * (y - mean_fixed) ** 2))
        c = np.sum(w_fixed) - np.sum(w_fixed**2) / np.sum(w_fixed)
        tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0

        w_random = 1.0 / (v + tau2)
        eff = float(np.sum(w_random * y) / np.sum(w_random))
        se_pool = float(np.sqrt(1.0 / np.sum(w_random)))
        i2 = float(max(0.0, (q - (k - 1)) / q) * 100.0) if q > 0 else 0.0

    z = eff / se_pool if se_pool > 0 else 0.0
    p = float(2.0 * stats.norm.sf(abs(z)))
    z_crit = float(stats.norm.ppf(1.0 - alpha / 2.0))
    dl_low, dl_high = eff - z_crit * se_pool, eff + z_crit * se_pool

    se_hk: float | None = None
    hk_low: float | None = None
    hk_high: float | None = None
    low, high = dl_low, dl_high
    if k > 1:
        # Hartung-Knapp: rescale by how far the studies actually sit from the
        # pooled mean, then use t rather than z. With equal within-study
        # variances and tau2 = 0 this collapses to the ordinary standard error
        # of the mean, which is how the implementation is checked.
        dispersion = float(np.sum(w_random * (y - eff) ** 2) / ((k - 1) * np.sum(w_random)))
        se_hk = float(np.sqrt(max(dispersion, 0.0)))
        t_crit = float(stats.t.ppf(1.0 - alpha / 2.0, k - 1))
        hk_low, hk_high = eff - t_crit * se_hk, eff + t_crit * se_hk
        # The modified rule: take the wider. HK is anti-conservative when the
        # studies happen to agree closely, and this repo cannot afford an
        # interval narrower than the one every other tool reports.
        if (hk_high - hk_low) > (dl_high - dl_low):
            low, high = hk_low, hk_high

    pi_low: float | None = None
    pi_high: float | None = None
    if k > 2:
        # Higgins-Thompson-Spiegelhalter: the spread of the mean *plus* the
        # spread between studies, on k-2 degrees of freedom.
        t_pred = float(stats.t.ppf(1.0 - alpha / 2.0, k - 2))
        spread = float(np.sqrt(tau2 + se_pool**2))
        pi_low, pi_high = eff - t_pred * spread, eff + t_pred * spread

    return PooledEffect(
        pooled_effect=eff,
        standard_error=se_pool,
        ci_low=low,
        ci_high=high,
        p_value=p,
        tau2=float(tau2),
        i2=i2,
        q=float(q),
        n_studies=k,
        ci_low_dl=dl_low,
        ci_high_dl=dl_high,
        se_hk=se_hk,
        ci_low_hk=hk_low,
        ci_high_hk=hk_high,
        pi_low=pi_low,
        pi_high=pi_high,
    )
