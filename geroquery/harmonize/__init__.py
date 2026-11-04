"""M3 harmonize — cross-study normalization & meta-analysis (public surface)."""

from .batch import batch_variance, remove_batch_effect
from .effect_size import EffectSize, cohens_d, hedges_g
from .influence import LeaveOneOut, corpus_influence, leave_one_out
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
    "remove_batch_effect",
    "batch_variance",
    "LeaveOneOut",
    "leave_one_out",
    "corpus_influence",
]
