"""Network control energy (stretch) — the explicit bridge to control theory.

Minimum control energy to drive a linear aging-network model from a current
state to a target ('youthful') state over horizon T:

    E* = (x_f - e^{A T} x_0)^T  W_c(T)^{-1}  (x_f - e^{A T} x_0)

where W_c(T) = \\int_0^T e^{A t} B B^T e^{A^T t} dt is the controllability
Gramian. Larger E* = harder to steer the network there. Deliberately small and
numerically explicit so it is testable on a toy network.

**Conditioning is the whole ballgame here.** Real biological networks are
*near*-uncontrollable: the Gramian is ill-conditioned rather than exactly
singular, so a naive ``solve`` succeeds and returns energies that look like
numbers but are numerically meaningless, differing by orders of magnitude under
rounding. We therefore report ``cond(W)``, invert through a rank-truncated
eigendecomposition, and refuse above a conditioning cutoff instead of returning
a confident answer we cannot support.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.linalg import expm

from ..exceptions import ResilienceInputError

#: Above this Gramian condition number the energy is not numerically trustworthy.
DEFAULT_COND_LIMIT = 1e10


@dataclass
class ControlEnergyResult:
    control_energy: float
    horizon: float
    gramian_condition_number: float
    gramian_rank: int
    dimension: int
    #: Fraction of the target displacement lying outside the reachable subspace.
    #: Large values mean most of the requested move is simply not achievable.
    unreachable_fraction: float
    truncated_modes: int
    well_conditioned: bool
    interpretation: str

    def to_dict(self) -> dict:
        return asdict(self)


def controllability_gramian(
    A: np.ndarray, B: np.ndarray, T: float, n_steps: int = 400
) -> np.ndarray:
    """Finite-horizon controllability Gramian via trapezoidal integration."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    n = A.shape[0]
    if A.shape != (n, n):
        raise ResilienceInputError("A must be square.")
    if B.ndim == 1:
        B = B[:, None]
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
    return 0.5 * (W + W.T)  # symmetrize away integration round-off


def control_energy_detailed(
    A: np.ndarray,
    B: np.ndarray,
    x0: np.ndarray,
    xf: np.ndarray,
    T: float = 1.0,
    n_steps: int = 400,
    cond_limit: float = DEFAULT_COND_LIMIT,
    strict: bool = True,
) -> ControlEnergyResult:
    """Minimum energy to steer ``x0 -> xf`` by time T, with conditioning diagnostics.

    Parameters
    ----------
    cond_limit: refuse (or flag, if ``strict`` is False) when the Gramian's
        condition number exceeds this. The default is deliberately conservative:
        past ~1e10 in double precision the inverse has lost most of its
        significant digits.
    strict: raise :class:`ResilienceInputError` on an ill-conditioned Gramian.
        Set False to get the flagged-but-returned value for exploratory work.
    """
    A = np.asarray(A, dtype=float)
    x0 = np.asarray(x0, dtype=float).ravel()
    xf = np.asarray(xf, dtype=float).ravel()
    n = A.shape[0]
    if x0.shape[0] != n or xf.shape[0] != n:
        raise ResilienceInputError("x0 and xf must match the dimension of A.")

    W = controllability_gramian(A, B, T, n_steps)
    delta = xf - expm(A * T) @ x0

    # Eigendecomposition of the (symmetric PSD) Gramian: cleaner rank handling
    # than a bare solve, and it gives the reachable/unreachable split directly.
    eigvals, eigvecs = np.linalg.eigh(W)
    eigvals = np.clip(eigvals, 0.0, None)
    top = float(eigvals.max()) if eigvals.size else 0.0
    if top <= 0.0:
        raise ResilienceInputError(
            "Controllability Gramian is identically zero; the system is uncontrollable.",
            detail={"dimension": n},
        )

    tol = top * max(n, 1) * np.finfo(float).eps
    keep = eigvals > tol
    rank = int(keep.sum())
    smallest_kept = float(eigvals[keep].min())
    cond = float(top / smallest_kept) if smallest_kept > 0 else float("inf")

    # Project the requested displacement onto the retained subspace.
    coords = eigvecs.T @ delta
    total_sq = float(coords @ coords)
    dropped_sq = float(coords[~keep] @ coords[~keep]) if (~keep).any() else 0.0
    unreachable = (dropped_sq / total_sq) if total_sq > 0 else 0.0

    energy = float(np.sum(coords[keep] ** 2 / eigvals[keep]))
    well_conditioned = bool(cond <= cond_limit and rank == n)

    if rank < n:
        interp = (
            f"Gramian is rank-deficient ({rank}/{n}); {unreachable:.1%} of the requested "
            "displacement lies in an uncontrollable direction and was truncated. The "
            "reported energy covers only the reachable component."
        )
    elif not well_conditioned:
        interp = (
            f"Gramian is ill-conditioned (cond = {cond:.3g} > {cond_limit:.3g}). The system "
            "is near-uncontrollable in at least one direction and the energy is not "
            "numerically trustworthy; treat it as a lower bound at best."
        )
    else:
        interp = (
            f"Well-conditioned (cond = {cond:.3g}). Energy is comparative — meaningful "
            "against other targets on the same network, not as an absolute quantity."
        )

    if strict and not well_conditioned:
        raise ResilienceInputError(
            "Controllability Gramian is too ill-conditioned for a trustworthy control "
            "energy. Re-run with strict=False to obtain the flagged estimate.",
            detail={
                "condition_number": cond,
                "cond_limit": cond_limit,
                "rank": rank,
                "dimension": n,
                "unreachable_fraction": unreachable,
            },
        )

    return ControlEnergyResult(
        control_energy=energy,
        horizon=float(T),
        gramian_condition_number=cond,
        gramian_rank=rank,
        dimension=n,
        unreachable_fraction=unreachable,
        truncated_modes=int(n - rank),
        well_conditioned=well_conditioned,
        interpretation=interp,
    )


def control_energy(
    A: np.ndarray,
    B: np.ndarray,
    x0: np.ndarray,
    xf: np.ndarray,
    T: float = 1.0,
    n_steps: int = 400,
) -> float:
    """Minimum energy to steer the linear system from ``x0`` to ``xf`` by time T.

    Raises :class:`ResilienceInputError` when the Gramian is too ill-conditioned
    for the answer to mean anything. Use :func:`control_energy_detailed` when you
    need the conditioning diagnostics or want the flagged estimate anyway.
    """
    return control_energy_detailed(A, B, x0, xf, T, n_steps).control_energy
