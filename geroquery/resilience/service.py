"""M6 resilience — orchestration over CSD, recovery, and control energy."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from ..exceptions import ResilienceInputError
from .control import control_energy
from .csd import CSDResult, csd_indicators
from .recovery import RecoveryResult, recovery_rate


class ResilienceService:
    def csd(
        self,
        data: pd.DataFrame,
        biomarker_cols: Sequence[str],
        age_col: str = "age",
        n_strata: int = 6,
        longitudinal: bool = False,
    ) -> CSDResult:
        missing = [c for c in [*biomarker_cols, age_col] if c not in data.columns]
        if missing:
            raise ResilienceInputError(
                "Dataset is missing required columns.",
                detail={"missing": missing, "available": list(data.columns)},
            )
        values = data[list(biomarker_cols)].to_numpy(dtype=float)
        ages = data[age_col].to_numpy(dtype=float)
        return csd_indicators(values, ages, n_strata=n_strata, longitudinal=longitudinal)

    def recovery(self, series: Sequence[float]) -> RecoveryResult:
        return recovery_rate(np.asarray(series, dtype=float))

    def control_energy(self, A, B, x0, xf, T: float = 1.0) -> dict:
        energy = control_energy(A, B, x0, xf, T)
        return {
            "control_energy": energy,
            "horizon": T,
            "assumptions": [
                "Aging network is approximated as a linear time-invariant system.",
                "Energy is the minimum-norm control input driving x0 -> xf by time T.",
                "Interpretation is comparative (relative difficulty of steering), not absolute.",
            ],
        }
