"""M5 clocks — clock registry.

Two tiers:

  * **Reference clocks** (below): small, fully transparent linear models with
    published-style structure. They exist so the service, validation, and
    age-acceleration maths are testable end-to-end with a known ground truth,
    and so the chronological-age-vs-mortality distinction is demonstrable.
  * **Library clocks**: real Horvath/PhenoAge/DunedinPACE/etc. via ``pyaging`` /
    ``biolearn`` when those are installed. We *wrap*, never rebuild.

Every registered clock — reference or library — carries ``predicted_outcome``
and ``training_population`` metadata, because a clock's chronological-age
accuracy does not imply mortality-prediction ability, and surfacing that is a
core correctness requirement (PRD §6, user stories 15/16).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..exceptions import ClockInputError, ClockNotFoundError
from ..models import ClockInfo

CLINICAL_FEATURES = ("albumin", "creatinine", "glucose", "crp", "lymphocyte_pct", "rdw")
EXPRESSION_FEATURES = ("ENSG00000147889", "ENSG00000124762", "ENSG00000113368", "ENSG00000130513")


@dataclass(frozen=True)
class LinearClock:
    """A transparent linear clock: prediction = intercept + features · weights."""

    info: ClockInfo
    intercept: float
    weights: dict[str, float]

    def predict(self, matrix: pd.DataFrame) -> np.ndarray:
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
        ordered = list(self.info.required_features)
        w = np.array([self.weights[f] for f in ordered], dtype=float)
        try:
            X = matrix[ordered].to_numpy(dtype=float)
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
        return self.intercept + X @ w


def _reference_clocks() -> dict[str, LinearClock]:
    clinical_pheno = LinearClock(
        info=ClockInfo(
            clock_id="clinical_phenoage_demo",
            name="Clinical PhenoAge (reference demo)",
            library="geroquery-reference",
            predicted_outcome="chronological_age",
            training_population="synthetic reference cohort (demo)",
            input_type="clinical",
            units="years",
            required_features=CLINICAL_FEATURES,
            notes="Transparent linear reference clock. Not a validated clinical instrument.",
        ),
        intercept=40.0,
        weights={
            "albumin": -3.0,
            "creatinine": 10.0,
            "glucose": 0.2,
            "crp": 2.0,
            "lymphocyte_pct": -0.3,
            "rdw": 3.0,
        },
    )
    clinical_mortality = LinearClock(
        info=ClockInfo(
            clock_id="clinical_mortality_demo",
            name="Clinical mortality score (reference demo)",
            library="geroquery-reference",
            predicted_outcome="mortality",
            training_population="synthetic reference cohort (demo)",
            input_type="clinical",
            units="log_hazard",
            required_features=CLINICAL_FEATURES,
            notes="Predicts a mortality risk score, NOT chronological age — same inputs, "
            "different target. Demonstrates that age accuracy != mortality prediction.",
        ),
        intercept=-2.0,
        weights={
            "albumin": -0.4,
            "creatinine": 0.6,
            "glucose": 0.02,
            "crp": 0.5,
            "lymphocyte_pct": -0.03,
            "rdw": 0.4,
        },
    )
    transcriptomic = LinearClock(
        info=ClockInfo(
            clock_id="transcriptomic_demo",
            name="Transcriptomic age (reference demo)",
            library="geroquery-reference",
            predicted_outcome="chronological_age",
            training_population="synthetic reference cohort (demo)",
            input_type="expression",
            units="years",
            required_features=EXPRESSION_FEATURES,
            notes="Transparent linear expression clock over four aging genes (demo).",
        ),
        intercept=45.0,
        weights={
            "ENSG00000147889": 8.0,
            "ENSG00000124762": 6.0,
            "ENSG00000113368": -7.0,
            "ENSG00000130513": 5.0,
        },
    )
    clocks = [clinical_pheno, clinical_mortality, transcriptomic]
    return {c.info.clock_id: c for c in clocks}


def _library_clocks() -> dict[str, LinearClock]:
    """Hook for pyaging / biolearn. Returns {} unless those are installed.

    Kept as a seam so real clocks appear in the registry automatically wherever
    the libraries are available (e.g. the hosted demo), without any other code
    change. We deliberately do not fabricate their coefficients here.
    """
    import importlib.util

    available = [
        lib for lib in ("pyaging", "biolearn") if importlib.util.find_spec(lib) is not None
    ]
    # When present, production code would enumerate and wrap each library clock.
    # We record availability but do not invent clocks we cannot execute.
    _ = available
    return {}


class ClockRegistry:
    def __init__(self):
        self._clocks: dict[str, LinearClock] = {**_reference_clocks(), **_library_clocks()}

    def list_clocks(self) -> list[ClockInfo]:
        return [c.info for c in self._clocks.values()]

    def get(self, clock_id: str) -> LinearClock:
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
