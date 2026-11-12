"""M6 resilience — orchestration over CSD and recovery metrics."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from ..exceptions import ResilienceInputError
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
        if len(biomarker_cols) < 2:
            # Cross-correlation is between markers; with one it is NaN, which
            # the API cannot even serialize.
            raise ResilienceInputError(
                "CSD needs at least two biomarkers.",
                detail={"biomarker_cols": list(biomarker_cols)},
            )
        missing = [c for c in [*biomarker_cols, age_col] if c not in data.columns]
        if missing:
            raise ResilienceInputError(
                "Dataset is missing required columns.",
                detail={"missing": missing, "available": list(data.columns)},
            )

        try:
            frame = data[list(biomarker_cols)].astype(float)
            ages = data[age_col].to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise ResilienceInputError(
                "Biomarker and age columns must be numeric.", detail={"error": str(exc)}
            ) from exc
        if log_columns:
            unknown = [c for c in log_columns if c not in biomarker_cols]
            if unknown:
                raise ResilienceInputError(
                    "log_columns must name biomarker columns.",
                    detail={"unknown": unknown, "biomarker_cols": list(biomarker_cols)},
                )
            for col in log_columns:
                if (frame[col] <= 0).any():
                    raise ResilienceInputError(
                        f"Cannot log-transform {col!r}: it contains non-positive values.",
                        detail={"column": col, "min": float(frame[col].min())},
                    )
                frame[col] = np.log(frame[col])

        values = frame.to_numpy()
        # A NaN poisons every stratum variance it falls in, and the bootstrap
        # then silently keeps only the resamples that missed it.
        incomplete = ~(np.isfinite(values).all(axis=1) & np.isfinite(ages))
        if incomplete.any():
            raise ResilienceInputError(
                f"{int(incomplete.sum())} rows have missing or non-finite values; "
                "drop or impute them first.",
                detail={"incomplete_rows": int(incomplete.sum())},
            )
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
