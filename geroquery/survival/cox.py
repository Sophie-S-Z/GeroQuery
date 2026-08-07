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

**Survey weights.** NHANES is a stratified multistage probability sample, so an
unweighted fit answers a question about the 2,517 people who were measured. Pass
``weights`` and the estimator maximizes Binder's (1992) weighted pseudo-
likelihood instead, which answers the question about the population they were
drawn to represent. Pass ``strata`` and ``psu`` as well and the variance becomes
the design-based (Taylor-linearized) one, which is wider because subjects inside
a primary sampling unit are correlated. Both matter and they matter separately:
the weights move the estimate, the design moves the interval.

One consequence has to be stated rather than absorbed: a likelihood ratio test
is not valid against a pseudo-likelihood. :func:`likelihood_ratio_test` refuses a
weighted fit, and :func:`wald_test` — which uses the design-based covariance —
is the nested test that replaces it.
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

    # Survey design. ``variance`` says which of the three estimators produced
    # ``standard_errors``, because an interval is only interpretable if you know
    # what it was computed against: "model" assumes independent observations,
    # "robust" drops the model but not the independence, "design" drops both.
    weighted: bool = False
    variance: str = "model"
    population_size: float = 0.0
    population_events: float = 0.0
    n_strata: int = 0
    n_psu: int = 0
    lonely_psu_strata: int = 0
    # Per covariate: the design-based variance over the variance the same
    # weighted fit would have if every subject were their own PSU. Above 1 means
    # clustering cost information; below 1 means stratification bought some.
    design_effect: list[float] = field(default_factory=list)
    covariance: list[list[float]] = field(default_factory=list)

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
            "weighted": self.weighted,
            "variance": self.variance,
            "population_size": round(self.population_size, 1),
            "population_events": round(self.population_events, 1),
            "n_strata": self.n_strata,
            "n_psu": self.n_psu,
            "lonely_psu_strata": self.lonely_psu_strata,
            "design_effect": [round(d, 3) for d in self.design_effect],
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
    x: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    beta: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Log partial likelihood, gradient, and negative Hessian, Breslow ties.

    Walks unique times from the largest down, accumulating the risk-set sums.
    Every subject still under observation at time t is in R(t), so processing
    downward means the risk set only ever grows and each subject is added once.

    ``weights`` are survey weights, and they enter in two places: the risk-set
    sums and the event contributions. Putting them in only one of the two gives a
    fit that converges and is wrong, which is why
    ``test_a_weight_of_two_is_the_same_as_the_row_appearing_twice`` exists.
    """
    n, p = x.shape
    sw = np.ones(n) if weights is None else weights
    eta = x @ beta
    # Subtract the max before exponentiating: eta can reach several hundred
    # during a bad Newton step and exp() would overflow to inf, turning a
    # recoverable step into a nan that poisons every later iteration.
    eta_shift = eta - eta.max()
    w = sw * np.exp(eta_shift)

    order = np.argsort(time, kind="mergesort")
    xs, ws, ts, es = x[order], w[order], time[order], event[order]
    etas, sws = eta_shift[order], sw[order]

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
        if eb.any():
            # d is the *weighted* number of events at this time; it reduces to
            # the count when the weights are all one.
            dw = sws[block][eb]
            d = float(dw.sum())
            loglik += float(dw @ etas[block][eb] - d * np.log(s0))
            mean = s1 / s0
            grad += dw @ xb[eb] - d * mean
            hess += d * (s2 / s0 - np.outer(mean, mean))
        i = j
    return loglik, grad, hess


def _score_residuals(
    x: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    beta: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Per-subject score contributions, summing to the gradient. Shape (n, p).

    These are the martingale-residual form of the score,

        U_i = w_i [ d_i (x_i - xbar(T_i))
                    - e^{eta_i} * sum_{t_k <= T_i} (x_i - xbar(t_k)) dLambda0(t_k) ]

    with Breslow's baseline increment. They are the only ingredient of the
    sandwich variance that is not already computed by the fit, and the identity
    ``sum_i U_i == gradient`` is what pins them to the likelihood.
    """
    n, p = x.shape
    sw = np.ones(n) if weights is None else weights
    eta = x @ beta
    eta_shift = eta - eta.max()
    risk = sw * np.exp(eta_shift)

    order = np.argsort(time, kind="mergesort")
    inverse = np.empty(n, dtype=int)
    inverse[order] = np.arange(n)
    xs, ts, es = x[order], time[order], event[order]
    rs, sws = risk[order], sw[order]

    # First pass, descending: xbar(t) and the baseline increment at each event
    # time, since both need the risk set.
    xbar = np.zeros((n, p))
    increment = np.zeros(n)
    s0 = 0.0
    s1 = np.zeros(p)
    i = n - 1
    while i >= 0:
        t = ts[i]
        j = i
        while j >= 0 and ts[j] == t:
            j -= 1
        block = slice(j + 1, i + 1)
        s0 += rs[block].sum()
        s1 += rs[block] @ xs[block]
        eb = es[block].astype(bool)
        if eb.any():
            xbar[block] = s1 / s0
            increment[block] = float(sws[block][eb].sum()) / s0
        i = j

    # Second pass, ascending: the cumulative terms every subject still at risk
    # has accrued by their own exit time.
    cumulative_increment = 0.0
    cumulative_weighted_xbar = np.zeros(p)
    expected = np.zeros((n, p))
    i = 0
    while i < n:
        t = ts[i]
        j = i
        while j < n and ts[j] == t:
            j += 1
        block = slice(i, j)
        if es[block].any():
            cumulative_increment += increment[i]
            cumulative_weighted_xbar = cumulative_weighted_xbar + increment[i] * xbar[i]
        # A subject exiting at t is at risk *through* t, so this time's
        # increment counts for them.
        expected[block] = cumulative_increment * xs[block] - cumulative_weighted_xbar
        i = j

    observed = np.where(es[:, None].astype(bool), xs - xbar, 0.0)
    residuals = sws[:, None] * observed - rs[:, None] * expected
    return residuals[inverse]


def _design_meat(
    scores: np.ndarray, strata: np.ndarray | None, psu: np.ndarray | None
) -> tuple[np.ndarray, int, int, int]:
    """Taylor-linearized variance of the total score, for a stratified design.

    Within each stratum the PSU totals are treated as an independent sample of
    size n_h, so the contribution is the usual with-replacement estimator

        (n_h / (n_h - 1)) * sum_a (u_ha - ubar_h)(u_ha - ubar_h)'

    A stratum with one PSU has no within-stratum degrees of freedom and
    contributes nothing. That is standard (a "certainty" PSU), but it is a
    silently narrowed interval if nobody is told, so the count is returned and
    carried onto the result rather than logged.
    """
    if strata is None and psu is None:
        # No design: every observation is its own cluster, which is the
        # Lin-Wei robust sandwich.
        return scores.T @ scores, 0, len(scores), 0

    n = len(scores)
    strata_arr = np.zeros(n, dtype=object) if strata is None else np.asarray(strata)
    psu_arr = np.arange(n) if psu is None else np.asarray(psu)

    meat = np.zeros((scores.shape[1], scores.shape[1]))
    lonely = 0
    total_psu = 0
    unique_strata = np.unique(strata_arr)
    for h in unique_strata:
        in_stratum = strata_arr == h
        # PSU ids are only unique within a stratum — NHANES reuses 1 and 2 in
        # every SDMVSTRA — so they are always grouped inside the stratum loop.
        units = np.unique(psu_arr[in_stratum])
        total_psu += len(units)
        if len(units) < 2:
            lonely += 1
            continue
        totals = np.array([scores[in_stratum & (psu_arr == a)].sum(axis=0) for a in units])
        centred = totals - totals.mean(axis=0)
        meat += (len(units) / (len(units) - 1)) * (centred.T @ centred)
    return meat, len(unique_strata), total_psu, lonely


def concordance_index(
    time: np.ndarray,
    event: np.ndarray,
    risk: np.ndarray,
    weights: np.ndarray | None = None,
) -> float:
    """Harrell's C: the share of comparable pairs the risk score orders correctly.

    A pair is comparable when the earlier of the two is a death — otherwise
    censoring means we do not know who died first. Ties in risk count as half,
    which is what makes C = 0.5 the value of a coin flip rather than of a
    constant predictor.

    ``weights`` makes each pair count ``w_i * w_j``, so the result describes the
    population the sample represents rather than the sample. Without it a
    C-index computed on an unequal-probability sample is a statement about who
    happened to be selected.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    w = np.ones(len(time)) if weights is None else np.asarray(weights, dtype=float)

    concordant = 0.0
    comparable = 0.0
    for i in np.flatnonzero(event == 1):
        # Subject i died at time[i]; anyone still observed after that is
        # comparable regardless of their own eventual status.
        later = time > time[i]
        if not later.any():
            continue
        r, pair = risk[later], w[i] * w[later]
        comparable += float(pair.sum())
        concordant += float(pair[risk[i] > r].sum()) + 0.5 * float(pair[risk[i] == r].sum())
    if comparable == 0:
        return float("nan")
    return concordant / comparable


def cox_regression(
    x: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    covariates: list[str] | None = None,
    *,
    weights: np.ndarray | None = None,
    strata: np.ndarray | None = None,
    psu: np.ndarray | None = None,
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
        weights: survey weights. Supplying them switches the objective to
            Binder's weighted pseudo-likelihood and the variance to a sandwich,
            because a model-based variance under weighting is not the variance
            of anything.
        strata, psu: design variables (NHANES ``SDMVSTRA``, ``SDMVPSU``).
            Supplying them switches the sandwich's meat to the Taylor-linearized
            design variance, which accounts for clustering and is wider.
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

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if weights.shape != time.shape:
            raise SurvivalInputError("weights must have one entry per subject.")
        if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
            raise SurvivalInputError("Survey weights must be positive and finite.")
    if (strata is not None or psu is not None) and weights is None:
        raise SurvivalInputError(
            "Design variables were given without weights. A clustered variance "
            "around an estimate that ignores selection probability mixes two "
            "different populations; supply weights too."
        )

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
        loglik, grad, hess = _breslow_terms(xw, time, event, beta, weights)
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
            trial_ll, _, _ = _breslow_terms(xw, time, event, trial, weights)
            if trial_ll >= loglik:
                break
            factor /= 2.0
        beta = beta + factor * step
        if np.max(np.abs(factor * step)) < CONVERGENCE_TOL:
            converged = True
            break

    loglik, _, hess = _breslow_terms(xw, time, event, beta, weights)
    try:
        bread = np.linalg.inv(hess)
    except np.linalg.LinAlgError:
        bread = np.linalg.inv(hess + _RIDGE * np.eye(hess.shape[0]))

    covariance = bread
    variance_kind = "model"
    n_strata = n_psu = lonely = 0
    deff = np.ones(x.shape[1])

    if weights is not None:
        # Under a pseudo-likelihood the inverse information is not the variance
        # of anything, so the sandwich is not an option here — it is the only
        # correct estimator.
        scores = _score_residuals(xw, time, event, beta, weights)
        independent, _, n_psu, _ = _design_meat(scores, None, None)
        covariance = bread @ independent @ bread
        variance_kind = "robust"
        if strata is not None or psu is not None:
            clustered, n_strata, n_psu, lonely = _design_meat(scores, strata, psu)
            design_covariance = bread @ clustered @ bread
            # The design effect: how much the interval widens once the design is
            # respected, holding the weights fixed. Comparing against the
            # model-based variance instead would fold in the weighting and give a
            # number that is not about clustering at all.
            with np.errstate(divide="ignore", invalid="ignore"):
                deff = np.where(
                    np.diag(covariance) > 0,
                    np.diag(design_covariance) / np.diag(covariance),
                    np.nan,
                )
            covariance = design_covariance
            variance_kind = "design"

    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    z_crit = float(stats.norm.ppf(1.0 - alpha / 2.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(se > 0, beta / se, np.nan)
    p = 2.0 * stats.norm.sf(np.abs(z))

    risk = xw @ beta
    total = float(len(time)) if weights is None else float(weights.sum())
    events = float(event.sum()) if weights is None else float(weights[event == 1].sum())
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
        concordance=float(concordance_index(time, event, risk, weights)),
        converged=converged,
        n_iterations=iterations,
        standardized=standardize,
        scale=[float(v) for v in scale],
        weighted=weights is not None,
        variance=variance_kind,
        population_size=total,
        population_events=events,
        n_strata=n_strata,
        n_psu=n_psu,
        lonely_psu_strata=lonely,
        design_effect=[float(v) for v in deff],
        covariance=[[float(v) for v in row] for row in covariance],
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
    if full.weighted or reduced.weighted:
        raise SurvivalInputError(
            "A likelihood ratio test is not valid against a weighted "
            "pseudo-likelihood: twice the difference is not chi-squared, so the "
            "p-value would come from the wrong reference distribution. Use "
            "wald_test, which uses the design-based covariance.",
            detail={"weighted": True},
        )
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
        "test": "likelihood_ratio",
        "statistic": round(statistic, 4),
        "df": df,
        "p_value": p_value,
        "added": [c for c in full.covariates if c not in reduced.covariates],
        "concordance_full": round(full.concordance, 4),
        "concordance_reduced": round(reduced.concordance, 4),
        "concordance_gain": round(full.concordance - reduced.concordance, 4),
    }


def wald_test(fit: CoxResult, covariates: list[str]) -> dict:
    """Are these coefficients jointly zero, judged by the fit's own covariance?

    The nested test that survives survey weighting. Because it reads the
    covariance matrix rather than the likelihood, it inherits whichever variance
    the fit used — model, robust, or design — and so a weighted, clustered fit
    gets a test that accounts for both. That is the whole reason it exists here.
    """
    if not fit.covariance:
        raise SurvivalInputError("This fit carries no covariance matrix.")
    missing = [c for c in covariates if c not in fit.covariates]
    if missing or not covariates:
        raise SurvivalInputError(
            f"Cannot test {missing or 'an empty set of'} covariates.",
            detail={"available": fit.covariates},
        )
    index = [fit.covariates.index(c) for c in covariates]
    beta = np.asarray(fit.coefficients)[index]
    block = np.asarray(fit.covariance)[np.ix_(index, index)]
    try:
        statistic = float(beta @ np.linalg.solve(block, beta))
    except np.linalg.LinAlgError as exc:
        raise SurvivalInputError(
            "The covariance block for these covariates is singular; they are "
            "collinear and cannot be tested jointly.",
            detail={"covariates": covariates},
        ) from exc
    return {
        "test": "design_wald" if fit.variance == "design" else "wald",
        "statistic": round(statistic, 4),
        "df": len(index),
        "p_value": float(stats.chi2.sf(statistic, len(index))),
        "tested": list(covariates),
        "variance": fit.variance,
    }
