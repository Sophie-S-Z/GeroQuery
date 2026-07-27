"""M5 clocks — application service.

Validates input, applies a clock, computes age acceleration, and (for clocks
that expose it, such as PhenoAge) the associated mortality risk. The interface
is fixed; new clocks arrive via the registry.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from ..exceptions import ClockInputError
from ..models import ClockInfo
from .registry import ClockRegistry, get_registry


@dataclass
class ClockResult:
    clock_id: str
    predicted_outcome: str
    units: str
    n_samples: int
    predictions: list[float]
    sample_ids: list[str]
    age_acceleration: list[float] | None = None
    mean_age_acceleration: float | None = None
    # Bootstrap 95% CI on the cohort mean age acceleration ([low, high]).
    mean_age_acceleration_ci: list[float] | None = None
    mortality_risk_10yr: list[float] | None = None

    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if v is not None}


def _bootstrap_mean_ci(values: np.ndarray, n_boot: int = 1000, seed: int = 0) -> list[float] | None:
    """Percentile 95% CI on the mean of ``values`` via resampling."""
    n = len(values)
    if n < 10:
        return None
    rng = np.random.default_rng(seed)
    means = values[rng.integers(0, n, (n_boot, n))].mean(axis=1)
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


class ClockService:
    def __init__(self, registry: ClockRegistry | None = None):
        self.registry = registry or get_registry()

    def list_clocks(self) -> list[ClockInfo]:
        return self.registry.list_clocks()

    def apply_clock(
        self,
        clock_id: str,
        matrix: pd.DataFrame,
        chronological_age: Sequence[float] | None = None,
    ) -> ClockResult:
        clock = self.registry.get(clock_id)
        if matrix is None or len(matrix) == 0:
            raise ClockInputError("Input matrix is empty.", detail={"clock_id": clock_id})

        predictions = clock.predict(matrix)
        sample_ids = [str(i) for i in matrix.index]

        result = ClockResult(
            clock_id=clock_id,
            predicted_outcome=clock.info.predicted_outcome,
            units=clock.info.units,
            n_samples=len(predictions),
            predictions=[float(x) for x in predictions],
            sample_ids=sample_ids,
        )

        # PhenoAge (and any clock exposing it) reports a companion mortality risk.
        if hasattr(clock, "mortality_risk"):
            try:
                risk = clock.mortality_risk(matrix)
                result.mortality_risk_10yr = [float(x) for x in risk]
            except Exception:  # noqa: BLE001 — mortality is auxiliary, never fatal
                result.mortality_risk_10yr = None

        # Chronological age lets us report age acceleration for age-predicting clocks.
        chrono = chronological_age
        if chrono is None and "age" in matrix.columns:
            chrono = matrix["age"].tolist()
        if chrono is not None and clock.info.predicted_outcome == "chronological_age":
            chrono_arr = np.asarray(chrono, dtype=float)
            if len(chrono_arr) != len(predictions):
                raise ClockInputError(
                    "chronological_age length does not match number of samples.",
                    detail={"n_ages": len(chrono_arr), "n_samples": len(predictions)},
                )
            accel = predictions - chrono_arr
            result.age_acceleration = [float(x) for x in accel]
            result.mean_age_acceleration = float(np.mean(accel))
            result.mean_age_acceleration_ci = _bootstrap_mean_ci(accel)
        return result

    def compare_clocks(
        self,
        clock_ids: Sequence[str],
        matrix: pd.DataFrame,
        chronological_age: Sequence[float] | None = None,
    ) -> list[ClockResult]:
        return [self.apply_clock(cid, matrix, chronological_age) for cid in clock_ids]
