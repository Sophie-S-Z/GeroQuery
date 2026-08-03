"""M6 resilience — orchestration over CSD, recovery, and control energy."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import pandas as pd

from ..exceptions import ResilienceInputError
from .control import DEFAULT_COND_LIMIT, control_energy_detailed
from .csd import DEFAULT_N_BOOTSTRAP, CSDResult, csd_indicators
from .recovery import RecoveryResult, recovery_rate


class ResilienceService:
    def csd(
        self,
        data: pd.DataFrame,
        biomarker_cols: Sequence[str],
        age_col: str = "age",
        n_strata: int = 6,
        longitudinal: bool = False,
        detrend: bool = True,
        n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
        log_columns: Sequence[str] | None = None,
    ) -> CSDResult:
        """Cross-sectional critical-slowing-down indicators over age strata.

        Args:
            log_columns: markers to natural-log transform first. Real clinical
                markers like hs-CRP are strongly right-skewed, and on the raw
                scale a handful of extreme values can dominate a stratum's
                variance — making the variance trend a story about outliers
                rather than about the population.
        """
        missing = [c for c in [*biomarker_cols, age_col] if c not in data.columns]
        if missing:
            raise ResilienceInputError(
                "Dataset is missing required columns.",
                detail={"missing": missing, "available": list(data.columns)},
            )

        frame = data[list(biomarker_cols)]
        if log_columns:
            unknown = [c for c in log_columns if c not in biomarker_cols]
            if unknown:
                raise ResilienceInputError(
                    "log_columns must name biomarker columns.",
                    detail={"unknown": unknown, "biomarker_cols": list(biomarker_cols)},
                )
            frame = frame.copy()
            for col in log_columns:
                if (frame[col] <= 0).any():
                    raise ResilienceInputError(
                        f"Cannot log-transform {col!r}: it contains non-positive values.",
                        detail={"column": col, "min": float(frame[col].min())},
                    )
                frame[col] = np.log(frame[col])

        values = frame.to_numpy(dtype=float)
        ages = data[age_col].to_numpy(dtype=float)
        return csd_indicators(
            values,
            ages,
            n_strata=n_strata,
            longitudinal=longitudinal,
            detrend=detrend,
            n_bootstrap=n_bootstrap,
        )

    def recovery(self, series: Sequence[float]) -> RecoveryResult:
        return recovery_rate(np.asarray(series, dtype=float))

    def control_energy(
        self,
        A,
        B,
        x0,
        xf,
        T: float = 1.0,
        cond_limit: float | None = None,
        strict: bool = True,
    ) -> dict:
        """Minimum control energy to steer x0 -> xf, with conditioning diagnostics.

        Returns the diagnostics rather than a bare float: on a near-uncontrollable
        system the energy is numerically enormous and meaningless, and a caller
        that only sees the number has no way to tell that apart from a genuinely
        expensive but well-posed transition.

        Args:
            strict: raise when the Gramian is worse conditioned than ``cond_limit``.
                Set ``False`` to get the estimate anyway, flagged in ``assumptions``
                — useful for comparing two equally ill-posed systems, where the
                ranking is still informative even though the magnitudes are not.
        """
        limit = DEFAULT_COND_LIMIT if cond_limit is None else cond_limit
        result = control_energy_detailed(A, B, x0, xf, T, cond_limit=limit, strict=strict)
        out = result.to_dict()  # already carries `horizon`
        out["assumptions"] = [
            "Aging network is approximated as a linear time-invariant system.",
            "Energy is the minimum-norm control input driving x0 -> xf by time T.",
            "Interpretation is comparative (relative difficulty of steering), not absolute.",
            "The Gramian is inverted on its numerically-supported subspace; directions "
            "below the rank tolerance are reported as unreachable rather than assigned "
            "an arbitrarily large finite energy.",
        ]
        if not result.well_conditioned or (
            isinstance(result.control_energy, float) and math.isinf(result.control_energy)
        ):
            out["assumptions"].append(
                "WARNING: the Gramian is ill-conditioned at this horizon; treat the "
                "energy as a lower bound on difficulty, not a calibrated quantity."
            )
        return out
