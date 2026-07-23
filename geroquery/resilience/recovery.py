"""DOSI-style recovery rate from longitudinal data.

Following the Pyrkov/Gero "dynamic organism state indicator" template: fit a
first-order autoregressive relaxation to a fluctuating state variable. The AR(1)
coefficient ``a`` is the fraction of a deviation that persists to the next step;
the recovery rate ``-ln(a)`` (and its inverse, the relaxation time) quantifies
how fast the system returns to baseline after perturbation. As organisms age,
``a -> 1`` (recovery slows, relaxation time diverges) — the dynamical signature
of resilience loss.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ..exceptions import ResilienceInputError


@dataclass
class RecoveryResult:
    ar1_coef: float
    recovery_rate: float  # -ln(ar1_coef); larger = more resilient
    relaxation_time: float  # 1 / recovery_rate (steps)
    n_points: int
    interpretation: str

    def to_dict(self) -> dict:
        return asdict(self)


def recovery_rate(series: np.ndarray) -> RecoveryResult:
    """Estimate AR(1) relaxation from a single evenly-sampled time series."""
    x = np.asarray(series, dtype=float).ravel()
    if x.size < 4:
        raise ResilienceInputError(
            "Need at least 4 time points to estimate a recovery rate.", detail={"n": int(x.size)}
        )
    if not np.all(np.isfinite(x)):
        raise ResilienceInputError("Time series contains non-finite values.")

    # Regress x[t+1] on x[t] about the mean: (x1-mu) = a (x0-mu).
    mu = x.mean()
    x0 = x[:-1] - mu
    x1 = x[1:] - mu
    denom = float(x0 @ x0)
    if denom == 0:
        raise ResilienceInputError("Time series is constant; AR(1) is undefined.")
    a = float((x0 @ x1) / denom)
    a_clipped = min(max(a, 1e-6), 1 - 1e-9)  # keep in (0,1) for a decaying process

    rate = float(-np.log(a_clipped))
    relax = float(1.0 / rate) if rate > 0 else float("inf")
    if a >= 0.9:
        interp = "slow recovery (high autocorrelation) — low resilience"
    elif a <= 0.5:
        interp = "fast recovery (low autocorrelation) — high resilience"
    else:
        interp = "intermediate recovery"
    return RecoveryResult(
        ar1_coef=a,
        recovery_rate=rate,
        relaxation_time=relax,
        n_points=int(x.size),
        interpretation=interp,
    )
