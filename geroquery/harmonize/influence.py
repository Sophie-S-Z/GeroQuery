"""Leave-one-out influence for a random-effects pool.

Two charts stand on this module, both from §3.3 of ``VISUAL_AND_STATISTICAL_PLAN``:

- **§4.4, the leave-one-out strip** beside the forest plot — for one gene, the
  pooled effect recomputed with each contrast removed, so a reader can see at a
  glance whether one study is carrying the result.
- **§4.5, the corpus-influence bar chart** — because this panel is near-balanced
  (§2.1), every gene is pooled from very nearly the same contrasts, so summing a
  contrast's leave-one-out shift across all genes answers a question about the
  *panel itself*: is any one contrast systematically bending the corpus? The
  aggregation is ``etl.build_frontend_data.add_influence``.

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
    for i in range(k):
        keep = np.arange(k) != i
        estimates[i] = random_effects(y[keep].tolist(), se[keep].tolist()).pooled_effect

    shifts = full_effect - estimates
    cooks = shifts**2 / se_full**2 if se_full > 0 else np.zeros(k)

    return LeaveOneOut(full_effect=full_effect, estimates=estimates, shifts=shifts, cooks=cooks)
