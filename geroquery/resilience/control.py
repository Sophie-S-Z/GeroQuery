"""Network control energy (stretch) — the explicit bridge to control theory.

Minimum control energy to drive a linear aging-network model from a current
state to a target ('youthful') state over horizon T:

    E* = (x_f - e^{A T} x_0)^T  W_c(T)^{-1}  (x_f - e^{A T} x_0)

where W_c(T) = \\int_0^T e^{A t} B B^T e^{A^T t} dt is the controllability
Gramian. Larger E* = harder to steer the network there. Deliberately small and
numerically explicit so it is testable on a toy network.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from ..exceptions import ResilienceInputError


def controllability_gramian(
    A: np.ndarray, B: np.ndarray, T: float, n_steps: int = 400
) -> np.ndarray:
    """Finite-horizon controllability Gramian via trapezoidal integration."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    n = A.shape[0]
    if A.shape != (n, n):
        raise ResilienceInputError("A must be square.")
    if B.shape[0] != n:
        raise ResilienceInputError("B must have the same number of rows as A.")
    ts = np.linspace(0.0, T, n_steps + 1)
    dt = ts[1] - ts[0]
    W = np.zeros((n, n))
    prev = None
    for t in ts:
        M = expm(A * t) @ B
        integrand = M @ M.T
        if prev is not None:
            W += 0.5 * (prev + integrand) * dt  # trapezoid
        prev = integrand
    return W


def control_energy(
    A: np.ndarray,
    B: np.ndarray,
    x0: np.ndarray,
    xf: np.ndarray,
    T: float = 1.0,
    n_steps: int = 400,
) -> float:
    """Minimum energy to steer the linear system from ``x0`` to ``xf`` by time T."""
    A = np.asarray(A, dtype=float)
    x0 = np.asarray(x0, dtype=float).ravel()
    xf = np.asarray(xf, dtype=float).ravel()
    if x0.shape[0] != A.shape[0] or xf.shape[0] != A.shape[0]:
        raise ResilienceInputError("x0 and xf must match the dimension of A.")
    W = controllability_gramian(A, B, T, n_steps)
    delta = xf - expm(A * T) @ x0
    try:
        sol = np.linalg.solve(W, delta)
    except np.linalg.LinAlgError as exc:
        raise ResilienceInputError(
            "System is not controllable to the target (singular Gramian).", detail=str(exc)
        ) from exc
    return float(delta @ sol)
