"""Aligning a knockdown's downstream response with the pooled aging direction.

The question this answers, for one CRISPRi knockdown: *do the genes it moves,
move the way they move with age?* Both sides are vectors over the same readout
genes — the perturbation's pseudobulk response on one side, GeroQuery's pooled
Hedges' g on the other — so the comparison is a cosine between them.

Two corrections are not optional, and both were forced by the data rather than
chosen for elegance:

**The perturbed gene is excluded from its own score.** A CRISPRi knockdown
drives its target's own readout down (median z of -0.64 across the screens,
against an off-target spread of 0.29). That entry is the knockdown working, not
a downstream effect, and leaving it in lets a gene's own aging estimate answer a
question about its consequences. The effect on the ranking is small (r = 0.9999)
— it is excluded because the alternative cannot be defended, not because it
changes the answer.

**The score is adjusted for how hard the knockdown hit.** Raw alignment
correlates +0.39 to +0.45 with the size of the response, and the mean alignment
over perturbations is positive (+0.06 to +0.14), not zero. Knocking out anything
important makes a cell sick, and a sick cell's transcriptome resembles an old
tissue's, so the raw score ranks damage and would present that as biology. The
adjustment regresses alignment on log magnitude and standardizes the residual,
which leaves a score answering a narrower and honest question: *among knockdowns
that disturbed the cell about this much, was this one unusually aging-like?*

What the adjusted score is not: a significance test. The null it is measured
against is the other perturbations, so it ranks and cannot single a gene out —
no gene in either screen reaches an FDR-significant alignment, and the
cross-screen agreement of the full ranking is 52.7%, a coin flip. Only the
extremes reproduce (the 300 most extreme agree at rho = +0.36, and the top-50
lists overlap 4-5x above chance), so the exporter gates on replication across
cell lines and the interface reports agreement per gene. See
``docs/PERTURBATION.md``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

# Below this many perturbations the magnitude slope is fitted on noise.
MIN_FOR_ADJUSTMENT = 20


@dataclass(frozen=True)
class Alignment:
    """One knockdown scored against the aging direction."""

    alignment: float
    """Cosine of the centred response with the centred aging direction, in
    [-1, +1]. NaN when the response carries no direction to compare."""

    magnitude: float
    """Euclidean norm of the centred response: how hard the knockdown hit."""


def _centred(values: np.ndarray) -> tuple[np.ndarray, float]:
    centred = values - values.mean()
    return centred, float(np.linalg.norm(centred))


def align(
    response: Sequence[float] | np.ndarray,
    direction: Sequence[float] | np.ndarray,
    skip: int | None = None,
) -> Alignment:
    """Score one knockdown's response against the aging direction.

    ``response`` and ``direction`` are over the same readout genes, in the same
    order. ``skip`` drops one index from both — the perturbed gene's own
    readout, which is the knockdown rather than its consequence.

    A flat ``direction`` raises: a caller asking to project onto nothing has a
    bug. A flat ``response`` returns NaN: that is a property of one knockdown in
    a batch of thousands, not of the call.
    """
    y = np.asarray(response, dtype=float)
    a = np.asarray(direction, dtype=float)
    if y.shape != a.shape or y.ndim != 1:
        raise ValueError(f"response {y.shape} and direction {a.shape} must be one matching axis.")
    if skip is not None:
        keep = np.ones(len(a), dtype=bool)
        keep[skip] = False
        y, a = y[keep], a[keep]
    a_c, a_norm = _centred(a)
    if a_norm == 0.0:
        raise ValueError("the aging direction has no direction: every value is identical.")
    y_c, y_norm = _centred(y)
    if y_norm == 0.0:
        return Alignment(alignment=float("nan"), magnitude=0.0)
    return Alignment(alignment=float(y_c @ a_c / (y_norm * a_norm)), magnitude=y_norm)


def align_many(
    responses: np.ndarray,
    direction: Sequence[float] | np.ndarray,
    skip: Sequence[int | None] | None = None,
) -> list[Alignment]:
    """:func:`align` over a (perturbations x genes) matrix, one ``skip`` each."""
    matrix = np.asarray(responses, dtype=float)
    if matrix.ndim != 2:
        raise ValueError(f"expected a (perturbations x genes) matrix, got {matrix.shape}.")
    skips: Sequence[int | None] = [None] * len(matrix) if skip is None else skip
    if len(skips) != len(matrix):
        raise ValueError(f"{len(skips)} skip entries for {len(matrix)} perturbations.")
    return [align(row, direction, skip=s) for row, s in zip(matrix, skips, strict=True)]


def adjust_for_magnitude(
    alignment: Sequence[float] | np.ndarray,
    magnitude: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Remove the response-magnitude confound and standardize what is left.

    Regresses alignment on ``log(magnitude)``, takes the residual, and divides
    by its standard deviation. The result is a z-score against the other
    perturbations in the same screen: positive means more aging-like than
    knockdowns that hit about as hard, negative means less.

    A straight line is enough — on both screens it leaves the residual
    uncorrelated with magnitude to machine precision, where decile matching
    leaves 0.009-0.015.
    """
    y = np.asarray(alignment, dtype=float)
    m = np.asarray(magnitude, dtype=float)
    if y.shape != m.shape or y.ndim != 1:
        raise ValueError(f"alignment {y.shape} and magnitude {m.shape} must be one matching axis.")
    if len(y) < MIN_FOR_ADJUSTMENT:
        raise ValueError(
            f"need at least {MIN_FOR_ADJUSTMENT} perturbations to fit the magnitude "
            f"slope, got {len(y)}."
        )
    if not np.all(m > 0.0):
        raise ValueError("every response magnitude must be positive to take a logarithm.")
    log_m = np.log(m)
    slope, intercept = np.polyfit(log_m, y, 1)
    residual = y - (slope * log_m + intercept)
    spread = residual.std()
    if spread == 0.0:
        raise ValueError("alignment is explained exactly by magnitude: no residual to score.")
    return residual / spread
