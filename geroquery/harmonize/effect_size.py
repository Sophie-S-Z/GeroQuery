"""Standardized effect sizes for age-associated change.

We express every study's signal as Hedges' g (bias-corrected standardized mean
difference between an *old* and a *young* group), with a standard error. This is
the currency the meta-analysis pools in: it is unit-free, so studies measured on
different platforms/scales become comparable.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class EffectSize:
    effect: float  # Hedges' g: +ve => higher in old group
    standard_error: float
    n_young: int
    n_old: int

    @property
    def direction(self) -> str:
        return "up" if self.effect >= 0 else "down"


def cohens_d(young: np.ndarray, old: np.ndarray) -> float:
    """Cohen's d with a pooled standard deviation."""
    young = np.asarray(young, dtype=float)
    old = np.asarray(old, dtype=float)
    n1, n2 = len(young), len(old)
    if n1 < 2 or n2 < 2:
        raise ValueError("Each group needs at least 2 samples for an effect size.")
    v1, v2 = young.var(ddof=1), old.var(ddof=1)
    pooled_sd = math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    if pooled_sd == 0:
        return 0.0
    return (old.mean() - young.mean()) / pooled_sd


def hedges_g(young: np.ndarray, old: np.ndarray) -> EffectSize:
    """Bias-corrected standardized mean difference and its SE.

    The small-sample correction factor J removes the upward bias of Cohen's d.
    SE follows the standard large-sample formula used in meta-analysis.
    """
    young = np.asarray(young, dtype=float)
    old = np.asarray(old, dtype=float)
    n1, n2 = len(young), len(old)
    d = cohens_d(young, old)
    df = n1 + n2 - 2
    # Hedges' correction J (exact form uses gamma; the 1 - 3/(4df-1) approx is
    # standard and accurate to <1% for df>=2).
    j = 1.0 - (3.0 / (4.0 * df - 1.0))
    g = j * d
    se = math.sqrt((n1 + n2) / (n1 * n2) + (g * g) / (2.0 * (n1 + n2)))
    return EffectSize(effect=g, standard_error=se, n_young=n1, n_old=n2)
