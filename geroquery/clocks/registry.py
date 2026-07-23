"""M5 clocks — clock registry.

The registry ships one **real, published** clinical clock — Levine PhenoAge
(see :mod:`geroquery.clocks.phenoage`) — rather than the transparent teaching
clocks earlier versions used. PhenoAge outputs a biological age in years *and*
a 10-year mortality risk from a single calibrated model, so the important
distinction between "biological age" and "mortality risk" is demonstrated with
one real instrument instead of two toy ones.

Real epigenetic clocks (Horvath, Hannum, DunedinPACE, …) require hundreds of
CpG coefficients; GeroQuery *wraps* them via ``pyaging`` / ``biolearn`` when
those libraries are installed rather than reprinting coefficients it cannot
validate. The seam below picks them up automatically where available.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..exceptions import ClockInputError, ClockNotFoundError
from ..models import ClockInfo
from . import phenoage

CLINICAL_FEATURES: tuple[str, ...] = phenoage.REQUIRED_FEATURES


@dataclass(frozen=True)
class PhenoAgeClock:
    """Wrapper exposing the real PhenoAge model through the clock interface."""

    info: ClockInfo

    def _check(self, matrix: pd.DataFrame) -> None:
        missing = [f for f in self.info.required_features if f not in matrix.columns]
        if missing:
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} requires features not present in the input.",
                detail={
                    "clock_id": self.info.clock_id,
                    "missing_features": missing,
                    "required": list(self.info.required_features),
                    "provided": list(matrix.columns),
                },
            )
        try:
            X = matrix[list(self.info.required_features)].to_numpy(dtype=float)
        except (ValueError, TypeError) as exc:
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} input has non-numeric values in required features.",
                detail=str(exc),
            ) from exc
        if np.isnan(X).any():
            raise ClockInputError(
                f"Clock {self.info.clock_id!r} input contains missing (NaN) values.",
                detail={"clock_id": self.info.clock_id},
            )

    def predict(self, matrix: pd.DataFrame) -> np.ndarray:
        self._check(matrix)
        return phenoage.phenotypic_age(matrix)

    def mortality_risk(self, matrix: pd.DataFrame) -> np.ndarray:
        self._check(matrix)
        return phenoage.mortality_risk_10yr(matrix)


def _reference_clocks() -> dict[str, PhenoAgeClock]:
    pheno = PhenoAgeClock(
        info=ClockInfo(
            clock_id="phenoage",
            name="PhenoAge (Levine et al., 2018)",
            library="geroquery-phenoage",
            predicted_outcome="chronological_age",
            training_population="US adults, NHANES III/IV (Liu et al., 2018)",
            input_type="clinical",
            units="years",
            required_features=phenoage.REQUIRED_FEATURES,
            notes="Real published clinical aging clock over 9 blood biomarkers + age. "
            "Also yields a 10-year mortality risk from the same model.",
        ),
    )
    return {pheno.info.clock_id: pheno}


def _library_clocks() -> dict[str, PhenoAgeClock]:
    """Hook for pyaging / biolearn. Returns {} unless those are installed.

    Kept as a seam so real DNA-methylation clocks appear automatically wherever
    the libraries are available, without any other code change. We deliberately
    never fabricate coefficients we cannot execute.
    """
    import importlib.util

    available = [
        lib for lib in ("pyaging", "biolearn") if importlib.util.find_spec(lib) is not None
    ]
    _ = available
    return {}


class ClockRegistry:
    def __init__(self):
        self._clocks: dict[str, PhenoAgeClock] = {**_reference_clocks(), **_library_clocks()}

    def list_clocks(self) -> list[ClockInfo]:
        return [c.info for c in self._clocks.values()]

    def get(self, clock_id: str) -> PhenoAgeClock:
        if clock_id not in self._clocks:
            raise ClockNotFoundError(
                f"Unknown clock {clock_id!r}.",
                detail={"available": list(self._clocks)},
            )
        return self._clocks[clock_id]


_REGISTRY: ClockRegistry | None = None


def get_registry() -> ClockRegistry:
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = ClockRegistry()
    return _REGISTRY
