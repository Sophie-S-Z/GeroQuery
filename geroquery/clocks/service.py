"""M5 clocks — application service.

Validates input, applies clocks, computes age acceleration, and compares clocks
on one input. The interface is fixed; new clocks arrive via the registry.
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

    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if v is not None}


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

        if chronological_age is not None:
            chrono = np.asarray(chronological_age, dtype=float)
            if len(chrono) != len(predictions):
                raise ClockInputError(
                    "chronological_age length does not match number of samples.",
                    detail={"n_ages": len(chrono), "n_samples": len(predictions)},
                )
            # Age acceleration is only meaningful for age-predicting clocks.
            if clock.info.predicted_outcome == "chronological_age":
                accel = predictions - chrono
                result.age_acceleration = [float(x) for x in accel]
                result.mean_age_acceleration = float(np.mean(accel))
        return result

    def compare_clocks(
        self,
        clock_ids: Sequence[str],
        matrix: pd.DataFrame,
        chronological_age: Sequence[float] | None = None,
    ) -> list[ClockResult]:
        return [self.apply_clock(cid, matrix, chronological_age) for cid in clock_ids]
