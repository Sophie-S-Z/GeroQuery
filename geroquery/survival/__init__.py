"""M-survival — time-to-event estimators for the mortality-linked cohort.

A peer of :mod:`geroquery.harmonize` and :mod:`geroquery.resilience`: an
estimator layer with no knowledge of sources or storage. Depends on numpy and
scipy only.
"""

from .cox import (
    CoxResult,
    SurvivalInputError,
    concordance_index,
    cox_regression,
    likelihood_ratio_test,
)
from .crosslayer import (
    CrossLayerResult,
    DysregulationResult,
    crosslayer_analysis,
    mahalanobis_dysregulation,
)

__all__ = [
    "CoxResult",
    "CrossLayerResult",
    "DysregulationResult",
    "SurvivalInputError",
    "concordance_index",
    "cox_regression",
    "crosslayer_analysis",
    "likelihood_ratio_test",
    "mahalanobis_dysregulation",
]
