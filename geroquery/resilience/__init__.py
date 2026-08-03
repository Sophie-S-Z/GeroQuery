"""M6 resilience — dynamical-systems / criticality metrics (public surface)."""

from .control import (
    ControlEnergyResult,
    control_energy,
    control_energy_detailed,
    controllability_gramian,
)
from .csd import CSDResult, TrendEvidence, csd_indicators
from .recovery import (
    REGIME_DECAYING,
    REGIME_DEGENERATE,
    REGIME_NONSTATIONARY,
    REGIME_OSCILLATORY,
    RecoveryResult,
    recovery_rate,
)
from .service import ResilienceService

__all__ = [
    "csd_indicators",
    "CSDResult",
    "TrendEvidence",
    "recovery_rate",
    "RecoveryResult",
    "REGIME_DECAYING",
    "REGIME_OSCILLATORY",
    "REGIME_NONSTATIONARY",
    "REGIME_DEGENERATE",
    "control_energy",
    "control_energy_detailed",
    "ControlEnergyResult",
    "controllability_gramian",
    "ResilienceService",
]
