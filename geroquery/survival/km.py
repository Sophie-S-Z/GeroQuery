"""Weighted Kaplan-Meier survival and restricted mean survival time.

The mortality view's hazard ratios describe a scale nobody has intuition for and
are only interpretable if proportional hazards holds over twenty years. RMST is
model-free — the area under the survival curve to a horizon τ — and its
difference between two groups is a difference in event-free *years* (Royston &
Parmar 2013). The curve itself is the chart every reader of a mortality result
expects and this one did not have (§4.3).

Both estimators are survey-weighted: the NHANES cohort was drawn to represent a
population, so a curve that ignores the weights describes the volunteers, not
the population. With unit weights they collapse exactly onto the ordinary
estimators — the same defining property the weighted Cox path is tested on.

Variance is left to the design-based machinery in :mod:`cox`; what ships to the
front end is the curve, the at-risk counts, and the RMST point estimate, all
computed here and drawn there.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class KMCurve:
    """A weighted survival curve sampled on a grid, with at-risk counts."""

    grid: np.ndarray  # time points the curve is sampled at
    survival: np.ndarray  # S(t) at each grid point
    at_risk: np.ndarray  # weighted number still at risk just before each grid point
    n_at_risk: np.ndarray  # unweighted subject count at risk (the table under the plot)

    def to_dict(self) -> dict:
        return {
            "grid": [round(float(t), 3) for t in self.grid],
            "survival": [round(float(s), 4) for s in self.survival],
            "at_risk": [round(float(w), 1) for w in self.at_risk],
            "n_at_risk": [int(n) for n in self.n_at_risk],
        }


@dataclass(frozen=True)
class RMST:
    value: float  # event-free years to the horizon
    tau: float


def _prepare(
    time: Sequence[float] | np.ndarray,
    event: Sequence[float] | np.ndarray,
    weight: Sequence[float] | np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.asarray(time, dtype=float)
    e = np.asarray(event, dtype=float)
    w = np.ones_like(t) if weight is None else np.asarray(weight, dtype=float)
    if not (t.shape == e.shape == w.shape) or t.ndim != 1:
        raise ValueError("time, event and weight must be 1-D and equal length.")
    if np.any(t < 0):
        raise ValueError("times must be non-negative.")
    if np.any(w < 0):
        raise ValueError("weights must be non-negative.")
    return t, e, w


def _steps(t: np.ndarray, e: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the distinct event times and the survival step value after each.

    Weighted Kaplan-Meier: at each event time, S multiplies by (1 − d/n) where
    d is the weighted number of events there and n the weighted number at risk.
    """
    order = np.argsort(t, kind="mergesort")
    t, e, w = t[order], e[order], w[order]
    event_times = np.unique(t[e > 0])
    surv = np.empty(len(event_times), dtype=float)
    running = 1.0
    for i, et in enumerate(event_times):
        at_risk = float(np.sum(w[t >= et]))
        deaths = float(np.sum(w[(t == et) & (e > 0)]))
        if at_risk > 0:
            running *= 1.0 - deaths / at_risk
        surv[i] = running
    return event_times, surv


def weighted_kaplan_meier(
    time: Sequence[float] | np.ndarray,
    event: Sequence[float] | np.ndarray,
    *,
    weight: Sequence[float] | np.ndarray | None = None,
    grid: Sequence[float] | np.ndarray | None = None,
) -> KMCurve:
    """Weighted KM survival sampled on ``grid`` (default: a 40-point even grid)."""
    t, e, w = _prepare(time, event, weight)
    if grid is None:
        grid = np.linspace(0.0, float(t.max()), 40)
    g = np.asarray(grid, dtype=float)

    event_times, surv = _steps(t, e, w)
    # S(t) is right-continuous step: the value at grid point x is the survival
    # after the last event time <= x (1.0 before the first event).
    s_at = np.ones(len(g), dtype=float)
    if len(event_times):
        idx = np.searchsorted(event_times, g, side="right") - 1
        s_at = np.where(idx >= 0, surv[np.clip(idx, 0, len(surv) - 1)], 1.0)

    at_risk = np.array([float(np.sum(w[t >= x])) for x in g])
    n_at_risk = np.array([int(np.sum(t >= x)) for x in g])
    return KMCurve(grid=g, survival=s_at, at_risk=at_risk, n_at_risk=n_at_risk)


def rmst(
    time: Sequence[float] | np.ndarray,
    event: Sequence[float] | np.ndarray,
    *,
    tau: float,
    weight: Sequence[float] | np.ndarray | None = None,
) -> RMST:
    """Restricted mean survival time: the area under the weighted KM curve to τ.

    Computed exactly from the step function (sum of S·Δt over the intervals the
    survival is constant on), not by sampling a grid, so it does not depend on a
    grid resolution.
    """
    t, e, w = _prepare(time, event, weight)
    if tau <= 0:
        raise ValueError("tau must be positive.")

    event_times, surv = _steps(t, e, w)
    # Breakpoints inside (0, tau): survival is 1 up to the first event time.
    inside = event_times[event_times < tau]
    s_inside = surv[event_times < tau]
    # Build the piecewise-constant survival with its left edges.
    edges = np.concatenate(([0.0], inside, [tau]))
    levels = np.concatenate(([1.0], s_inside))  # value on each [edge_i, edge_{i+1})
    widths = np.diff(edges)
    area = float(np.sum(levels * widths))
    return RMST(value=area, tau=float(tau))
