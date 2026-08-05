"""Cox proportional-hazards regression and Harrell's C, in numpy and scipy only.

Why hand-rolled rather than lifelines or statsmodels: neither is a dependency of
this package, and adding a heavyweight one to compute a partial likelihood that
fits in eighty lines would make the mortality path harder to install than the
rest of the repo. The cost of hand-rolling is that the estimator has to earn
trust, so it is checked three ways in ``tests/test_survival.py``:

1. against an example small enough to differentiate by hand,
2. by finite-difference agreement between the analytic gradient/Hessian and the
   log partial likelihood they claim to be derivatives of,
3. by recovering a planted log hazard ratio from simulated survival data —
   the same "does the estimator find an effect known to be there" argument the
   synthetic CSD fixture exists to make.

Ties are handled by Breslow's approximation. Efron's is more accurate when many
events share a time; here follow-up is recorded in whole months over ~20 years,
so ties are common and the choice is not cosmetic. Breslow is used because it is
the one that can be verified against a hand computation, and its bias attenuates
hazard ratios toward the null — an estimator that understates an effect is the
safer failure mode for a repo whose headline result is a null.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from ..exceptions import GeroQueryError

# Newton-Raphson stops when the largest coefficient step falls below this.
CONVERGENCE_TOL = 1e-9
MAX_ITERATIONS = 100

# Ridge added to the Hessian diagonal only if it is singular. Collinear
# predictors (a clock and its own age acceleration, say) otherwise produce a
# LinAlgError that reads like a bug rather than a modelling mistake.
_RIDGE = 1e-8


class SurvivalInputError(GeroQueryError):
    """Malformed survival input: shape mismatch, no events, negative times."""

    code = "survival_input_error"
    http_status = 422


@dataclass(frozen=True)
class CoxResult:
    """A fitted Cox model.

    ``hazard_ratio`` is per unit of the predictor as supplied. When
    :func:`cox_regression` is called with ``standardize=True`` (the default),
    that unit is one standard deviation — which is the only way a hazard ratio
    for an epigenetic clock in years and one for a variance slope in arbitrary
    units can be put beside each other.
    """

    covariates: list[str]
    coefficients: list[float]
    standard_errors: list[float]
    hazard_ratios: list[float]
    ci_low: list[float]
    ci_high: list[float]
    z_scores: list[float]
    p_values: list[float]
    n_subjects: int
    n_events: int
    log_likelihood: float
    concordance: float
    converged: bool
    n_iterations: int
    standardized: bool
    scale: list[float] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "covariates": self.covariates,
            "coefficients": [round(c, 6) for c in self.coefficients],
            "standard_errors": [round(s, 6) for s in self.standard_errors],
            "hazard_ratios": [round(h, 4) for h in self.hazard_ratios],
            "ci_low": [round(c, 4) for c in self.ci_low],
            "ci_high": [round(c, 4) for c in self.ci_high],
            "z_scores": [round(z, 3) for z in self.z_scores],
            "p_values": self.p_values,
            "n_subjects": self.n_subjects,
            "n_events": self.n_events,
            "log_likelihood": round(self.log_likelihood, 4),
            "concordance": round(self.concordance, 4),
            "converged": self.converged,
            "n_iterations": self.n_iterations,
            "standardized": self.standardized,
            "scale": [round(s, 6) for s in self.scale],
        }

    def summary_rows(self) -> list[dict]:
        """One row per covariate, for tabular display."""
        return [
            {
                "covariate": name,
                "coefficient": round(self.coefficients[i], 4),
                "hazard_ratio": round(self.hazard_ratios[i], 4),
                "ci_low": round(self.ci_low[i], 4),
                "ci_high": round(self.ci_high[i], 4),
                "p_value": self.p_values[i],
                # The interval, not the p-value, is what says "we cannot tell".
                "excludes_null": bool(self.ci_low[i] > 1.0 or self.ci_high[i] < 1.0),
            }
            for i, name in enumerate(self.covariates)
        ]


def _validate(x: np.ndarray, time: np.ndarray, event: np.ndarray) -> None:
    if x.ndim != 2:
        raise SurvivalInputError(f"Design matrix must be 2-D, got shape {x.shape}.")
    if not (len(time) == len(event) == x.shape[0]):
        raise SurvivalInputError(
            "time, event and the design matrix must have the same number of rows.",
            detail={"n_rows": x.shape[0], "n_time": len(time), "n_event": len(event)},
        )
    if x.shape[0] == 0:
        raise SurvivalInputError("No subjects.")
    if np.any(time < 0):
        raise SurvivalInputError("Follow-up times must be non-negative.")
    if not np.all(np.isin(event, (0, 1))):
        raise SurvivalInputError("Event indicator must be 0 or 1.")
    if event.sum() == 0:
        raise SurvivalInputError("No events observed; a Cox model is not identified.")
    if not np.all(np.isfinite(x)):
        raise SurvivalInputError("Design matrix contains non-finite values.")


def _breslow_terms(
    x: np.ndarray, time: np.ndarray, event: np.ndarray, beta: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    """Log partial likelihood, gradient, and negative Hessian, Breslow ties.

    Walks unique times from the largest down, accumulating the risk-set sums.
    Every subject still under observation at time t is in R(t), so processing
    downward means the risk set only ever grows and each subject is added once.
    """
    n, p = x.shape
    eta = x @ beta
    # Subtract the max before exponentiating: eta can reach several hundred
    # during a bad Newton step and exp() would overflow to inf, turning a
    # recoverable step into a nan that poisons every later iteration.
    eta_shift = eta - eta.max()
    w = np.exp(eta_shift)

    order = np.argsort(time, kind="mergesort")
    xs, ws, ts, es = x[order], w[order], time[order], event[order]
    etas = eta_shift[order]

    loglik = 0.0
    grad = np.zeros(p)
    hess = np.zeros((p, p))

    s0 = 0.0
    s1 = np.zeros(p)
    s2 = np.zeros((p, p))

    i = n - 1
    while i >= 0:
        t = ts[i]
        j = i
        while j >= 0 and ts[j] == t:
            j -= 1
        block = slice(j + 1, i + 1)

        wb = ws[block]
        xb = xs[block]
        s0 += wb.sum()
        s1 += wb @ xb
        s2 += (xb * wb[:, None]).T @ xb

        eb = es[block].astype(bool)
        d = int(eb.sum())
        if d:
            loglik += float(etas[block][eb].sum() - d * np.log(s0))
            mean = s1 / s0
            grad += xb[eb].sum(axis=0) - d * mean
            hess += d * (s2 / s0 - np.outer(mean, mean))
        i = j
    return loglik, grad, hess


def concordance_index(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """Harrell's C: the share of comparable pairs the risk score orders correctly.

    A pair is comparable when the earlier of the two is a death — otherwise
    censoring means we do not know who died first. Ties in risk count as half,
    which is what makes C = 0.5 the value of a coin flip rather than of a
    constant predictor.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)

    concordant = 0.0
    comparable = 0.0
    for i in np.flatnonzero(event == 1):
        # Subject i died at time[i]; anyone still observed after that is
        # comparable regardless of their own eventual status.
        later = time > time[i]
        if not later.any():
            continue
        r = risk[later]
        comparable += r.size
        concordant += float((risk[i] > r).sum()) + 0.5 * float((risk[i] == r).sum())
    if comparable == 0:
        return float("nan")
    return concordant / comparable


def cox_regression(
    x: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    covariates: list[str] | None = None,
    *,
    standardize: bool = True,
    alpha: float = 0.05,
) -> CoxResult:
    """Fit a Cox proportional-hazards model by Newton-Raphson.

    Args:
        x: design matrix, subjects by covariates.
        time: follow-up time. Units set the time scale only; hazard ratios are
            unitless, so years and months give identical coefficients.
        event: 1 if the follow-up ended in death, 0 if censored.
        covariates: names, for the result. Defaults to ``x0, x1, ...``.
        standardize: divide each column by its standard deviation so a hazard
            ratio is per SD. On by default because the whole point here is to
            compare predictors measured in different units, and a per-year HR
            for a clock is not comparable to a per-unit HR for a variance slope.
        alpha: two-sided level for the confidence interval.

    Raises:
        SurvivalInputError: malformed input, or no events.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    _validate(x, time, event)

    names = covariates or [f"x{i}" for i in range(x.shape[1])]
    if len(names) != x.shape[1]:
        raise SurvivalInputError(
            f"Got {len(names)} covariate names for {x.shape[1]} columns.",
            detail={"names": names},
        )

    # Centre always (does not change the partial likelihood, improves
    # conditioning); scale only when asked.
    scale = np.ones(x.shape[1])
    if standardize:
        sd = x.std(axis=0, ddof=1)
        # A constant column carries no information; leaving its scale at 1
        # keeps the fit from dividing by zero and lets the singular-Hessian
        # branch report it as unidentified rather than as nan.
        sd[sd == 0] = 1.0
        scale = sd
    xw = (x - x.mean(axis=0)) / scale

    beta = np.zeros(x.shape[1])
    converged = False
    iterations = 0
    loglik = float("nan")
    hess = np.eye(x.shape[1])

    while iterations < MAX_ITERATIONS:
        iterations += 1
        loglik, grad, hess = _breslow_terms(xw, time, event, beta)
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.solve(hess + _RIDGE * np.eye(hess.shape[0]), grad)
        # Backtrack if the step made the likelihood worse. Without this a
        # near-separated covariate sends beta to infinity in one jump and the
        # fit reports a hazard ratio of 1e17 rather than failing to converge.
        factor = 1.0
        for _ in range(20):
            trial = beta + factor * step
            trial_ll, _, _ = _breslow_terms(xw, time, event, trial)
            if trial_ll >= loglik:
                break
            factor /= 2.0
        beta = beta + factor * step
        if np.max(np.abs(factor * step)) < CONVERGENCE_TOL:
            converged = True
            break

    loglik, _, hess = _breslow_terms(xw, time, event, beta)
    try:
        covariance = np.linalg.inv(hess)
    except np.linalg.LinAlgError:
        covariance = np.linalg.inv(hess + _RIDGE * np.eye(hess.shape[0]))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    z_crit = float(stats.norm.ppf(1.0 - alpha / 2.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(se > 0, beta / se, np.nan)
    p = 2.0 * stats.norm.sf(np.abs(z))

    risk = xw @ beta
    return CoxResult(
        covariates=list(names),
        coefficients=[float(v) for v in beta],
        standard_errors=[float(v) for v in se],
        hazard_ratios=[float(v) for v in np.exp(beta)],
        ci_low=[float(v) for v in np.exp(beta - z_crit * se)],
        ci_high=[float(v) for v in np.exp(beta + z_crit * se)],
        z_scores=[float(v) for v in z],
        p_values=[float(v) for v in p],
        n_subjects=int(x.shape[0]),
        n_events=int(event.sum()),
        log_likelihood=float(loglik),
        concordance=float(concordance_index(time, event, risk)),
        converged=converged,
        n_iterations=iterations,
        standardized=standardize,
        scale=[float(v) for v in scale],
    )


def likelihood_ratio_test(full: CoxResult, reduced: CoxResult) -> dict:
    """Does the larger model explain more than the smaller one?

    This is the test behind the third cross-layer question — whether resilience
    adds anything *over* an epigenetic clock. Comparing two hazard ratios side
    by side cannot answer that; nesting the models and testing the difference
    in log partial likelihood can.

    The caller is responsible for the models actually being nested and fitted on
    the same rows. A likelihood ratio between models fitted on different
    complete-case subsets is meaningless, and nothing here can detect it, so
    :func:`crosslayer_report` builds one design matrix and slices it.
    """
    df = len(full.covariates) - len(reduced.covariates)
    if df <= 0:
        raise SurvivalInputError(
            "The full model must have more covariates than the reduced model.",
            detail={"full": full.covariates, "reduced": reduced.covariates},
        )
    if full.n_subjects != reduced.n_subjects or full.n_events != reduced.n_events:
        raise SurvivalInputError(
            "Nested models must be fitted on identical rows.",
            detail={
                "full": {"n": full.n_subjects, "events": full.n_events},
                "reduced": {"n": reduced.n_subjects, "events": reduced.n_events},
            },
        )
    statistic = 2.0 * (full.log_likelihood - reduced.log_likelihood)
    p_value = float(stats.chi2.sf(statistic, df))
    return {
        "statistic": round(statistic, 4),
        "df": df,
        "p_value": p_value,
        "added": [c for c in full.covariates if c not in reduced.covariates],
        "concordance_full": round(full.concordance, 4),
        "concordance_reduced": round(reduced.concordance, 4),
        "concordance_gain": round(full.concordance - reduced.concordance, 4),
    }
