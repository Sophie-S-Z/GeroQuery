"""M3 harmonize — cross-study normalization & meta-analysis (public surface)."""

from .batch import batch_variance, remove_batch_effect
from .effect_size import EffectSize, cohens_d, hedges_g
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
]
