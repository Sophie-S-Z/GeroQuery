"""Leave-one-out influence for a random-effects pool, and its corpus aggregate.

Two charts stand on this module, both from §3.3 of ``VISUAL_AND_STATISTICAL_PLAN``:

- **§4.4, the leave-one-out strip** beside the forest plot — for one gene, the
  pooled effect recomputed with each contrast removed, so a reader can see at a
  glance whether one study is carrying the result.
- **§4.5, the corpus-influence bar chart** — because this panel is near-balanced
  (§2.1), every gene is pooled from very nearly the same contrasts, so summing a
  contrast's leave-one-out shift across all genes answers a question about the
  *panel itself*: is any one contrast systematically bending the corpus?

Nothing here re-derives the pool. Every leave-one-out estimate is
:func:`geroquery.harmonize.random_effects` on the k-1 subset — the same
estimator the API and the static export already serve — so the strip cannot
disagree with the forest plot it sits beside. The browser only renders these
numbers (the standing repo rule); they are computed once in the ETL.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

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


def corpus_influence(
    contrasts: Sequence[Mapping[str, object]],
) -> list[dict]:
    """Aggregate leave-one-out shift per contrast across every pooled gene.

    ``contrasts`` is the long table of per-study effects with at least
    ``gene_id``, ``study_id``, ``effect_size`` and ``standard_error``. Genes with
    fewer than :data:`MIN_CONTRASTS` contrasts are skipped (they carry no pooled
    row either). The result is one row per ``study_id``, ranked by the total
    absolute pooled-estimate shift it is responsible for — the §4.5 bar chart.
    """
    by_gene: dict[str, list[tuple[str, float, float]]] = defaultdict(list)
    for row in contrasts:
        eff = row.get("effect_size")
        se = row.get("standard_error")
        if eff is None or se is None:
            continue
        se_f = float(se)  # type: ignore[arg-type]
        if not np.isfinite(float(eff)) or not np.isfinite(se_f) or se_f <= 0:  # type: ignore[arg-type]
            continue
        by_gene[str(row["gene_id"])].append(
            (str(row["study_id"]), float(eff), se_f)  # type: ignore[arg-type]
        )

    total_shift: dict[str, float] = defaultdict(float)
    total_cook: dict[str, float] = defaultdict(float)
    max_shift: dict[str, float] = defaultdict(float)
    n_genes: dict[str, int] = defaultdict(int)

    for members in by_gene.values():
        if len(members) < MIN_CONTRASTS:
            continue
        studies = [m[0] for m in members]
        loo = leave_one_out([m[1] for m in members], [m[2] for m in members])
        for study, shift, cook in zip(studies, loo.shifts, loo.cooks, strict=True):
            abs_shift = float(abs(shift))
            total_shift[study] += abs_shift
            total_cook[study] += float(cook)
            max_shift[study] = max(max_shift[study], abs_shift)
            n_genes[study] += 1

    ranking = [
        {
            "study_id": study,
            "total_abs_shift": round(total_shift[study], 6),
            "mean_abs_shift": round(total_shift[study] / n_genes[study], 6),
            "max_abs_shift": round(max_shift[study], 6),
            "total_cooks": round(total_cook[study], 6),
            "n_genes": n_genes[study],
        }
        for study in total_shift
    ]
    ranking.sort(key=lambda r: -cast(float, r["total_abs_shift"]))
    return ranking
