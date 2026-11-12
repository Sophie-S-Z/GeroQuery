"""M3 harmonize — cross-study normalization & meta-analysis (public surface)."""

from .deming import DemingFit, deming_regression
from .effect_size import EffectSize, cohens_d, hedges_g
from .influence import LeaveOneOut, leave_one_out
from .meta import PooledEffect, random_effects
from .shrinkage import Shrunk, adaptive_shrinkage

__all__ = [
    "EffectSize",
    "cohens_d",
    "hedges_g",
    "PooledEffect",
    "random_effects",
    "Shrunk",
    "adaptive_shrinkage",
    "LeaveOneOut",
    "leave_one_out",
    "DemingFit",
    "deming_regression",
]
