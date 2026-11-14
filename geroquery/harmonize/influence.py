"""Leave-one-out influence for a random-effects pool.

Two charts stand on this estimator, both from §3.3 of ``VISUAL_AND_STATISTICAL_PLAN``:
the leave-one-out strip beside the forest plot (§4.4), which shows whether one
study is carrying a gene's result, and the corpus-influence chart (§4.5), which
the exporter (:func:`geroquery.etl.build_frontend_data.add_influence`) builds by
aggregating these per-gene fits for every contrast across the genes it measures.

Nothing here re-derives the pool. Every leave-one-out estimate is
:func:`geroquery.harmonize.random_effects` on the k-1 subset — the same
estimator the API and the static export already serve — so the strip cannot
disagree with the forest plot it sits beside. The browser only renders these
numbers (the standing repo rule); they are computed once in the ETL.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from .meta import random_effects

# A gene needs at least this many contrasts for a leave-one-out to make sense:
# with k contrasts each subset has k-1, and below three the remaining pool is
# not a random-effects fit. Matches the exporter's MIN_CONTRASTS pooling gate.
MIN_CONTRASTS = 3


@dataclass(frozen=True)
class LeaveOneOut:
    """Per-contrast influence on one gene's pooled effect.

    Every array is aligned to the input order, so the §4.4 strip can share the
    forest plot's y-axis row for row.
    """

    full_effect: float
    estimates: np.ndarray  # pooled effect with contrast i removed
    shifts: np.ndarray  # full_effect - estimates[i] (signed; how far it moves)
    cooks: np.ndarray  # shifts[i]**2 / se_full**2 (Cook's distance, intercept-only)
    baujat_q: np.ndarray  # x: contrast i's contribution to Cochran's Q
    baujat_influence: np.ndarray  # y: standardised influence on the pooled estimate


def leave_one_out(effects: Sequence[float], standard_errors: Sequence[float]) -> LeaveOneOut:
    """Recompute the pool with each study removed, and the influence it carries.

    Parameters
    ----------
    effects, standard_errors:
        Per-contrast effect sizes and SEs (same length, k >= 2).
    """
    y = np.asarray(effects, dtype=float)
    se = np.asarray(standard_errors, dtype=float)
    if y.shape != se.shape or y.ndim != 1:
        raise ValueError("effects and standard_errors must be 1-D and equal length.")
    k = len(y)
    if k < 2:
        raise ValueError("Leave-one-out needs at least two studies.")

    full = random_effects(y.tolist(), se.tolist())
    full_effect = full.pooled_effect
    se_full = full.standard_error

    estimates = np.empty(k, dtype=float)
    se_without = np.empty(k, dtype=float)
    for i in range(k):
        keep = np.arange(k) != i
        sub = random_effects(y[keep].tolist(), se[keep].tolist())
        estimates[i] = sub.pooled_effect
        se_without[i] = sub.standard_error

    shifts = full_effect - estimates
    cooks = shifts**2 / se_full**2 if se_full > 0 else np.zeros(k)

    # Baujat (2002): x is the fixed-effect contribution to heterogeneity Q, y is
    # the standardised influence on the pooled estimate (Viechtbauer & Cheung
    # 2010). Both read off quantities already computed above.
    v = se**2
    w_fixed = 1.0 / v
    mean_fixed = float(np.sum(w_fixed * y) / np.sum(w_fixed))
    baujat_q = w_fixed * (y - mean_fixed) ** 2
    baujat_influence = np.where(se_without > 0, shifts**2 / se_without**2, 0.0)

    return LeaveOneOut(
        full_effect=full_effect,
        estimates=estimates,
        shifts=shifts,
        cooks=cooks,
        baujat_q=baujat_q,
        baujat_influence=baujat_influence,
    )
