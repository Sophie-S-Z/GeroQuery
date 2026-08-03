"""DOSI-style recovery rate from longitudinal data.

Following the Pyrkov/Gero "dynamic organism state indicator" template: fit a
first-order autoregressive relaxation to a fluctuating state variable. The AR(1)
coefficient ``a`` is the fraction of a deviation that persists to the next step;
the recovery rate ``-ln(a)`` (and its inverse, the relaxation time) quantifies
how fast the system returns to baseline after perturbation. As organisms age,
``a -> 1`` (recovery slows, relaxation time diverges) — the dynamical signature
of resilience loss.

**The AR(1) relaxation reading is only valid for ``a`` in (0, 1).** Outside that
interval the fitted process is not a decaying one and ``-ln(a)`` is meaningless,
so we classify the dynamical regime explicitly and refuse to report a recovery
rate rather than coercing the estimate into the decaying range. An oscillatory
(anti-persistent, ``a < 0``) series is *not* a fast-recovering one, and reporting
it as "high resilience" would be confidently wrong on exactly the noisy,
artifact-prone biomarker series this is meant to handle.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np

from ..exceptions import ResilienceInputError

# Regime labels for the fitted AR(1) coefficient.
REGIME_DECAYING = "decaying"  # 0 < a < 1 — relaxation reading is valid
REGIME_OSCILLATORY = "oscillatory"  # a < 0 — anti-persistent, AR(1) relaxation undefined
REGIME_NONSTATIONARY = "nonstationary"  # a >= 1 — no return to baseline
REGIME_DEGENERATE = "degenerate"  # a == 0 — no memory; relaxation instantaneous

_Z95 = 1.959963984540054


@dataclass
class RecoveryResult:
    ar1_coef: float
    ar1_ci_low: float
    ar1_ci_high: float
    regime: str
    #: ``-ln(a)``; larger = more resilient. ``None`` outside the decaying regime.
    recovery_rate: float | None
    #: ``1 / recovery_rate`` in steps. ``None`` outside the decaying regime.
    relaxation_time: float | None
    recovery_rate_ci_low: float | None
    recovery_rate_ci_high: float | None
    n_points: int
    valid: bool
    interpretation: str
    assumptions: list[str]

    def to_dict(self) -> dict:
        return asdict(self)


def _ar1_fit(x: np.ndarray) -> tuple[float, float]:
    """OLS AR(1) coefficient about the mean, plus its standard error."""
    mu = x.mean()
    x0 = x[:-1] - mu
    x1 = x[1:] - mu
    denom = float(x0 @ x0)
    if denom == 0:
        raise ResilienceInputError("Time series is constant; AR(1) is undefined.")
    a = float((x0 @ x1) / denom)
    resid = x1 - a * x0
    dof = max(len(x0) - 1, 1)
    se = math.sqrt(float(resid @ resid) / dof / denom)
    return a, se


def _classify(a: float) -> str:
    if a < 0:
        return REGIME_OSCILLATORY
    if a == 0:
        return REGIME_DEGENERATE
    if a >= 1:
        return REGIME_NONSTATIONARY
    return REGIME_DECAYING


def recovery_rate(series: np.ndarray) -> RecoveryResult:
    """Estimate AR(1) relaxation from a single evenly-sampled time series.

    Returns a result whose ``valid`` flag and ``regime`` say whether the
    relaxation-time reading applies at all. Callers must check ``valid`` before
    interpreting ``recovery_rate``; it is ``None`` whenever the fit falls outside
    the decaying regime.
    """
    x = np.asarray(series, dtype=float).ravel()
    if x.size < 4:
        raise ResilienceInputError(
            "Need at least 4 time points to estimate a recovery rate.", detail={"n": int(x.size)}
        )
    if not np.all(np.isfinite(x)):
        raise ResilienceInputError("Time series contains non-finite values.")

    a, se = _ar1_fit(x)
    ci_low, ci_high = a - _Z95 * se, a + _Z95 * se
    regime = _classify(a)

    assumptions = [
        "Series is evenly sampled; unequal spacing invalidates the step interpretation.",
        "Fluctuations are modelled as a stationary first-order autoregressive relaxation "
        "about a constant mean; any trend should be removed before fitting.",
        "Recovery rate is reported in units of 1/step, not 1/time — convert using the "
        "sampling interval.",
    ]

    if regime != REGIME_DECAYING:
        detail = {
            REGIME_OSCILLATORY: (
                "AR(1) coefficient is negative (anti-persistent / oscillatory dynamics). "
                "The relaxation model does not apply and no recovery rate is reported — "
                "anti-persistence is not resilience. Common causes: measurement noise "
                "dominating signal, alternating acquisition artifacts, or over-differencing."
            ),
            REGIME_NONSTATIONARY: (
                "AR(1) coefficient is >= 1 (non-stationary / random-walk or divergent). "
                "The system does not return to baseline, so relaxation time is undefined."
            ),
            REGIME_DEGENERATE: (
                "AR(1) coefficient is exactly zero (no serial memory). Relaxation is "
                "instantaneous under the model, which usually means the sampling interval "
                "is long relative to the true relaxation time."
            ),
        }[regime]
        return RecoveryResult(
            ar1_coef=a,
            ar1_ci_low=ci_low,
            ar1_ci_high=ci_high,
            regime=regime,
            recovery_rate=None,
            relaxation_time=None,
            recovery_rate_ci_low=None,
            recovery_rate_ci_high=None,
            n_points=int(x.size),
            valid=False,
            interpretation=detail,
            assumptions=assumptions,
        )

    rate = float(-math.log(a))
    # -ln is monotone decreasing, so the CI endpoints swap. Clip the coefficient
    # CI into (0,1) first: outside it the transform is undefined.
    lo_a = min(max(ci_high, 1e-12), 1.0 - 1e-12)
    hi_a = min(max(ci_low, 1e-12), 1.0 - 1e-12)
    rate_lo, rate_hi = float(-math.log(lo_a)), float(-math.log(hi_a))

    if a >= 0.9:
        interp = "slow recovery (high autocorrelation) — low resilience"
    elif a <= 0.5:
        interp = "fast recovery (low autocorrelation) — high resilience"
    else:
        interp = "intermediate recovery"
    if ci_low <= 0.0 or ci_high >= 1.0:
        interp += (
            " [wide confidence interval crosses the boundary of the decaying regime; "
            "treat the point estimate as weakly identified]"
        )

    return RecoveryResult(
        ar1_coef=a,
        ar1_ci_low=ci_low,
        ar1_ci_high=ci_high,
        regime=regime,
        recovery_rate=rate,
        relaxation_time=float(1.0 / rate),
        recovery_rate_ci_low=rate_lo,
        recovery_rate_ci_high=rate_hi,
        n_points=int(x.size),
        valid=True,
        interpretation=interp,
        assumptions=assumptions,
    )
