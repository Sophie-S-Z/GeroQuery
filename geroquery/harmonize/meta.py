"""Random-effects meta-analysis (DerSimonian-Laird).

Pools per-study effect sizes into one signature *without* ever concatenating raw
expression matrices (which would confound platform and batch with biology). The
random-effects model additionally estimates between-study heterogeneity (tau^2,
I^2) so callers can see how consistent the studies actually are.
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

    @property
    def direction(self) -> str:
        return "up" if self.pooled_effect >= 0 else "down"


def random_effects(
    effects: Sequence[float], standard_errors: Sequence[float], ci: float = 0.95
) -> PooledEffect:
    """DerSimonian-Laird pooling of standardized effects.

    Parameters
    ----------
    effects, standard_errors:
        Per-study effect sizes and their SEs (same length, >=1).
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

    if k == 1:
        eff = float(y[0])
        se_pool = float(se[0])
    else:
        # Fixed-effect pooled mean, used to compute Q.
        mean_fixed = np.sum(w_fixed * y) / np.sum(w_fixed)
        q = float(np.sum(w_fixed * (y - mean_fixed) ** 2))
        c = np.sum(w_fixed) - np.sum(w_fixed**2) / np.sum(w_fixed)
        tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0

        w_random = 1.0 / (v + tau2)
        eff = float(np.sum(w_random * y) / np.sum(w_random))
        se_pool = float(np.sqrt(1.0 / np.sum(w_random)))

    # For k==1 there is no heterogeneity to estimate.
    if k == 1:
        tau2, q, i2 = 0.0, 0.0, 0.0
    else:
        i2 = float(max(0.0, (q - (k - 1)) / q) * 100.0) if q > 0 else 0.0

    z = eff / se_pool if se_pool > 0 else 0.0
    p = float(2.0 * stats.norm.sf(abs(z)))
    z_crit = float(stats.norm.ppf(1.0 - (1.0 - ci) / 2.0))
    return PooledEffect(
        pooled_effect=eff,
        standard_error=se_pool,
        ci_low=eff - z_crit * se_pool,
        ci_high=eff + z_crit * se_pool,
        p_value=p,
        tau2=float(tau2),
        i2=i2,
        q=float(q),
        n_studies=k,
    )
