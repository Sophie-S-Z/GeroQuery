"""Adaptive shrinkage across the corpus, and the local false sign rate.

Every other estimator in this repo looks at one gene at a time. This one looks
at all 41,983 at once, and that is the whole point: the pooled estimates are not
41,983 unrelated questions, they are one distribution, and most of it is at
zero. Knowing that is information about any single gene, and the pooling
throws it away.

**What it produces, and why each is worth having.**

``lfsr`` — the *local false sign rate* (Stephens 2017): the posterior
probability that the **sign** of a gene's effect is wrong. This is the question
the tool is actually asked. Someone typing a gene name wants to know whether it
goes up or down with age and how much to trust that; a p-value answers "is the
effect exactly zero", which nobody believes and nobody asked. The lfsr is
defined for every gene, including the 38,276 whose interval crosses zero and
which today get a verdict and nothing else — "we cannot tell" is honest but it
is not the same as "we cannot tell, and here is how close we are".

``posterior_mean`` — the shrunken effect. An estimate that is extreme *and*
imprecise is mostly noise, and the corpus knows what the noise looks like. This
is the fix for a defect in the shipped UI rather than a refinement: the
landscape table ranks genes by pooled *g* among those whose interval excludes
zero, which is selecting on significance and then ranking on an extreme of a
noisy statistic. That is the winner's curse, and the top of that table is
enriched for exactly the genes that will shrink on replication. Ranking on the
posterior mean asks for "large *and* well measured" instead of "large".

**The model.** Effects are drawn from a mixture of zero-centred normals — a
point mass at zero plus a geometric grid of standard deviations wide enough to
cover the observed spread. Mixture weights are fit by EM over the whole corpus;
each gene's posterior is then analytic. The unimodal-at-zero assumption is not
an assumption of convenience here: that most genes do not change with age is the
corpus's own headline finding.

**What this is not.** The lfsr is not a p-value and must never be labelled one;
it is a posterior probability under a prior estimated from the corpus. It is
also *corpus-relative* — the same gene gets a different lfsr if the panel
changes — so it sits behind the estimator version in the evidence baseline like
every other pooled quantity.

Reference: M. Stephens, "False discovery rates: a new deal", *Biostatistics*
18(2):275-294, 2017.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

# EM stops when the penalised log-likelihood improves by less than this. The
# quantity being fit is a set of mixture weights used only to shrink; a tighter
# tolerance buys digits nobody reads.
_TOLERANCE = 1e-7
_MAX_ITERATIONS = 2000

# Dirichlet concentration on the null component. Stephens' default, and it is
# deliberately asymmetric: it biases the fit *toward* calling things null, so
# the shrinkage errs on the conservative side. In a repo whose headline result
# is a null, a prior that has to be argued into declaring an effect is the right
# way round.
_NULL_PRIOR_WEIGHT = 10.0


@dataclass(frozen=True)
class Shrunk:
    """Corpus-level posterior for a set of estimates.

    Arrays are parallel to the inputs and in the same order.
    """

    posterior_mean: np.ndarray
    posterior_sd: np.ndarray
    lfsr: np.ndarray
    #: P(effect > 0 | data) and P(effect < 0 | data). They do not sum to 1 --
    #: the remainder is the posterior mass on exactly zero.
    prob_positive: np.ndarray
    prob_negative: np.ndarray
    #: Fitted weight on the **exact** point mass at zero.
    #:
    #: This is *not* the fraction of genes with no age effect, and reading it as
    #: one will understate the null badly -- it came out 0.66 on a corpus
    #: planted at 0.90. The grid's narrowest components sit at a standard
    #: deviation an order of magnitude below the smallest standard error, so
    #: "exactly zero" and "very nearly zero" fit the data almost identically and
    #: the mixture has no way to apportion mass between them. The point mass is
    #: identified only up to that near-null neighbourhood. Use
    #: :meth:`prior_mass_below`, which asks the question that *is* identified.
    null_proportion: float
    #: The fitted grid, for reproducibility and for tests that need to see it.
    sd_grid: np.ndarray
    weights: np.ndarray

    def __len__(self) -> int:
        return int(self.posterior_mean.shape[0])

    def prior_mass_below(self, threshold: float) -> float:
        """Fitted prior probability that a gene's true effect is under
        ``threshold`` in magnitude.

        The honest version of "how much of this corpus is null". Where
        :attr:`null_proportion` splits hairs the data cannot split, this
        integrates the whole fitted prior over an interval the caller names, so
        the answer is stable no matter how the near-zero mass happened to be
        distributed across the narrow components.

        On a corpus planted with 90% null effects it recovers 0.906 at a
        threshold of 0.1, against a true 0.915.
        """
        if threshold <= 0:
            raise ValueError("threshold must be positive.")
        nonzero = self.sd_grid > 0
        within = 2.0 * stats.norm.cdf(threshold / np.where(nonzero, self.sd_grid, 1.0)) - 1.0
        return float(self.weights[0] + np.sum(self.weights[1:] * within[1:]))


def _sd_grid(effects: np.ndarray, standard_errors: np.ndarray) -> np.ndarray:
    """A geometric grid of prior SDs, from "smaller than any SE" to "wider than
    any observed effect".

    Stephens' construction. The grid only has to be dense enough that the
    mixture can approximate the true prior; the EM decides which parts of it
    matter, so being generous costs a little arithmetic and nothing else. The
    leading 0.0 is the point mass at exactly no effect.
    """
    se_min = float(np.min(standard_errors))
    # The largest effect that is not explicable as noise, doubled for headroom.
    excess = float(np.max(effects**2 - standard_errors**2))
    sd_max = 2.0 * np.sqrt(excess) if excess > 0 else 8.0 * se_min
    sd_min = se_min / 10.0
    if sd_max <= sd_min:
        sd_max = sd_min * 16.0
    n_steps = int(np.ceil(np.log2(sd_max / sd_min) / np.log2(np.sqrt(2.0)))) + 1
    grid = sd_min * np.sqrt(2.0) ** np.arange(n_steps)
    return np.concatenate([[0.0], grid])


def _fit_weights(likelihood: np.ndarray) -> np.ndarray:
    """EM for the mixture weights, with a Dirichlet penalty on the null.

    ``likelihood`` is (n_observations, n_components). The penalty enters as
    pseudo-counts, which is what keeps the null weight from collapsing to zero
    just because the grid can always fit a little better with a wider component.
    """
    n_components = likelihood.shape[1]
    prior_counts = np.ones(n_components)
    prior_counts[0] = _NULL_PRIOR_WEIGHT
    weights = np.full(n_components, 1.0 / n_components)

    previous = -np.inf
    for _ in range(_MAX_ITERATIONS):
        weighted = likelihood * weights
        total = weighted.sum(axis=1, keepdims=True)
        # A row can underflow to zero if every component is a terrible fit.
        # Falling back to the current weights leaves that row uninformative
        # rather than propagating a nan through the whole fit.
        degenerate = total[:, 0] <= 0
        if degenerate.any():
            weighted[degenerate] = weights
            total[degenerate, 0] = 1.0
        responsibility = weighted / total

        counts = responsibility.sum(axis=0) + prior_counts - 1.0
        weights = np.maximum(counts, 0.0)
        weights /= weights.sum()

        objective = float(np.sum(np.log(total[:, 0]))) + float(
            np.sum((prior_counts - 1.0) * np.log(np.maximum(weights, 1e-300)))
        )
        if objective - previous < _TOLERANCE:
            break
        previous = objective

    return weights


def adaptive_shrinkage(
    effects: np.ndarray | list[float],
    standard_errors: np.ndarray | list[float],
) -> Shrunk:
    """Fit the corpus-wide prior and return each estimate's posterior.

    Parameters
    ----------
    effects, standard_errors:
        Pooled effect sizes and their standard errors, one pair per gene. Same
        length, all standard errors strictly positive.

    Notes
    -----
    Vectorised over genes and components -- the corpus is ~42,000 x ~25, which
    is one matrix, not a loop.
    """
    y = np.asarray(effects, dtype=float)
    se = np.asarray(standard_errors, dtype=float)
    if y.shape != se.shape or y.ndim != 1:
        raise ValueError("effects and standard_errors must be 1-D and equal length.")
    if y.size == 0:
        raise ValueError("Need at least one estimate to shrink.")
    if np.any(se <= 0):
        raise ValueError("All standard errors must be positive.")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(se)):
        raise ValueError("effects and standard_errors must all be finite.")

    grid = _sd_grid(y, se)
    # Marginal likelihood of each estimate under each component: the prior
    # variance and the sampling variance simply add.
    marginal_sd = np.sqrt(se[:, None] ** 2 + grid[None, :] ** 2)
    likelihood = stats.norm.pdf(y[:, None], loc=0.0, scale=marginal_sd)

    weights = _fit_weights(likelihood)

    posterior_weight = likelihood * weights
    row_total = posterior_weight.sum(axis=1, keepdims=True)
    row_total[row_total <= 0] = 1.0
    posterior_weight /= row_total

    # Per-component posterior: the usual normal-normal conjugate update. The
    # first column is the point mass, whose posterior mean and variance are both
    # exactly zero -- the expressions below give that for free at grid == 0.
    shrink_factor = grid[None, :] ** 2 / marginal_sd**2
    component_mean = y[:, None] * shrink_factor
    component_var = shrink_factor * se[:, None] ** 2

    posterior_mean = np.sum(posterior_weight * component_mean, axis=1)
    # Law of total variance: average within-component variance plus the spread
    # of the component means. Dropping the second term is the classic way to
    # ship an over-confident posterior.
    second_moment = np.sum(posterior_weight * (component_var + component_mean**2), axis=1)
    posterior_sd = np.sqrt(np.maximum(second_moment - posterior_mean**2, 0.0))

    # Sign probabilities. The point mass contributes to neither tail: it is mass
    # on exactly zero, which is why prob_positive + prob_negative < 1.
    nonzero = grid > 0
    component_sd = np.sqrt(np.maximum(component_var, 0.0))
    safe_sd = np.where(component_sd > 0, component_sd, 1.0)
    tail_positive = np.where(nonzero[None, :], stats.norm.sf(0.0, component_mean, safe_sd), 0.0)
    prob_positive = np.sum(posterior_weight * tail_positive, axis=1)
    prob_negative = np.sum(posterior_weight * (nonzero[None, :] - tail_positive), axis=1)

    # lfsr: the chance of being wrong about the sign if forced to call one. Mass
    # sitting exactly on zero counts against *both* directions, which is what
    # makes this quantity conservative where an ordinary tail probability is not.
    null_mass = np.clip(1.0 - prob_positive - prob_negative, 0.0, 1.0)
    lfsr = np.minimum(prob_positive, prob_negative) + null_mass

    return Shrunk(
        posterior_mean=posterior_mean,
        posterior_sd=posterior_sd,
        lfsr=np.clip(lfsr, 0.0, 1.0),
        prob_positive=prob_positive,
        prob_negative=prob_negative,
        null_proportion=float(weights[0]),
        sd_grid=grid,
        weights=weights,
    )
