"""Per-gene young-vs-old contrasts from an expression matrix.

This is the step between "a study's expression matrix" and "a row in the
signature table": given a matrix and two disjoint sample groups, produce one
standardized effect size per gene.

The maths is the same Hedges' *g* as :mod:`geroquery.harmonize.effect_size`,
evaluated over all genes at once instead of one at a time — a 55,000-probe array
with 100 samples is 55,000 scalar calls otherwise, per contrast, per study.
``test_matrix_hedges_g_matches_scalar_implementation`` asserts the two agree, so
the vectorized copy cannot drift from the reference.

Three decisions are made here that a reader should be able to find and argue
with, so they are constants rather than inline literals:

``LINEAR_SCALE_PERCENTILE`` / ``LINEAR_SCALE_THRESHOLD``
    GEO's declared ``value_type`` is not reliable enough to decide whether a
    matrix is already log-scaled — "count" covers both raw MAS5 intensities and
    normalized values, and some "transformed count" sets are linear. A
    standardized mean difference computed on a linear intensity scale is
    dominated by whichever genes happen to be brightly expressed, so the scale is
    inferred from the data instead of trusted from the header.

``LOG_FLOOR``
    Background-subtracted intensities can be zero or negative, which has no
    logarithm. Flooring compresses the bottom of the range; dropping those probes
    instead would silently delete the lowest-expressed genes, which is a worse
    bias than compressing them.

``MIN_GROUP_SIZE``
    Below this, the pooled SD is too unstable for the effect size to mean
    anything, and Hedges' correction is being asked to rescue a sample it cannot.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# A matrix whose high percentile is above this is treated as linear intensity.
# Log2-scaled microarray data essentially never exceeds ~16-20; linear MAS5
# intensities routinely run into the thousands.
LINEAR_SCALE_PERCENTILE = 99.0
LINEAR_SCALE_THRESHOLD = 50.0

LOG_FLOOR = 1.0

MIN_GROUP_SIZE = 3


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted p-values, order-preserving.

    NaN inputs stay NaN and are excluded from the ranking rather than being
    treated as p=1, which would inflate the denominator.
    """
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    n = int(finite.sum())
    if n == 0:
        return out

    idx = np.flatnonzero(finite)
    order = idx[np.argsort(p[idx], kind="stable")]
    ranks = np.arange(1, n + 1, dtype=float)
    scaled = p[order] * n / ranks
    # Enforce monotonicity from the largest p downwards, then clip to 1.
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    out[order] = np.minimum(adjusted, 1.0)
    return out


def looks_linear(values: np.ndarray) -> bool:
    """True if the matrix appears to be on a linear intensity scale."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return False
    return bool(np.percentile(finite, LINEAR_SCALE_PERCENTILE) > LINEAR_SCALE_THRESHOLD)


def to_log2(frame: pd.DataFrame) -> pd.DataFrame:
    """Log2-transform, flooring non-positive values at :data:`LOG_FLOOR`."""
    return pd.DataFrame(
        np.log2(np.maximum(frame.to_numpy(dtype=float), LOG_FLOOR)),
        index=frame.index,
        columns=frame.columns,
    )


def collapse_probes(frame: pd.DataFrame, genes: pd.Series) -> pd.DataFrame:
    """Reduce a probe-level matrix to one row per gene.

    Where several probes measure the same gene, the one with the highest mean
    signal wins (the "MaxMean" rule). Averaging instead would pull a gene's
    profile towards its worst-performing probe: on these platforms a large share
    of probes for any given gene sit at background and carry no signal, so their
    contribution is noise with the same weight as the informative probe.
    """
    labels = genes.reindex(frame.index)
    values = frame.to_numpy(dtype=float)
    keep = labels.notna().to_numpy() & np.isfinite(values).all(axis=1)
    if not keep.any():
        return frame.iloc[:0]

    sub = frame.loc[keep]
    means = sub.to_numpy(dtype=float).mean(axis=1)
    ranked = pd.DataFrame({"gene": labels.loc[keep].to_numpy(), "mean": means}, index=sub.index)
    # Stable sort on (gene, -mean); the first row per gene is the MaxMean probe.
    winners = ranked.sort_values(["gene", "mean"], ascending=[True, False], kind="stable")
    winners = winners[~winners["gene"].duplicated()]
    out = sub.loc[winners.index]
    out.index = pd.Index(winners["gene"].to_numpy(), name="gene")
    return out


def hedges_g_matrix(young: np.ndarray, old: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Row-wise Hedges' *g* and its standard error.

    Args:
        young: ``(n_genes, n_young)`` values for the young group.
        old: ``(n_genes, n_old)`` values for the old group.

    Returns:
        ``(g, se)``, each of shape ``(n_genes,)``. Positive *g* means higher in
        the old group. Genes with zero pooled variance get ``g = 0`` rather than
        an infinity — no difference is detectable, which is a finite statement.
    """
    n1, n2 = young.shape[1], old.shape[1]
    df = n1 + n2 - 2
    if df < 1:
        raise ValueError("A contrast needs at least 3 samples in total.")

    v1 = young.var(axis=1, ddof=1)
    v2 = old.var(axis=1, ddof=1)
    pooled_sd = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / df)
    delta = old.mean(axis=1) - young.mean(axis=1)
    d = np.divide(delta, pooled_sd, out=np.zeros_like(delta), where=pooled_sd > 0)

    j = 1.0 - (3.0 / (4.0 * df - 1.0))
    g = j * d
    se = np.sqrt((n1 + n2) / (n1 * n2) + (g * g) / (2.0 * (n1 + n2)))
    return g, se


def contrast_effects(
    frame: pd.DataFrame,
    young_samples: list[str],
    old_samples: list[str],
    *,
    genes: pd.Series | None = None,
    already_log: bool = False,
) -> pd.DataFrame:
    """One standardized age effect per gene for a single young-vs-old contrast.

    Args:
        frame: probe-by-sample (or gene-by-sample) expression matrix.
        young_samples: column labels forming the young group.
        old_samples: column labels forming the old group.
        genes: optional probe -> gene label mapping; when given, probes are
            collapsed to genes before the contrast.
        already_log: skip the linear-scale check and log2 transform.

    Returns:
        A frame indexed by gene with ``effect_size``, ``standard_error``,
        ``p_value``, ``q_value``, ``direction``, ``n_young``, ``n_old``.

    Raises:
        ValueError: either group is smaller than :data:`MIN_GROUP_SIZE`, or the
            groups overlap.
    """
    from scipy import stats

    young_samples = [c for c in young_samples if c in frame.columns]
    old_samples = [c for c in old_samples if c in frame.columns]
    overlap = set(young_samples) & set(old_samples)
    if overlap:
        raise ValueError(f"Young and old groups share samples: {sorted(overlap)}")
    if len(young_samples) < MIN_GROUP_SIZE or len(old_samples) < MIN_GROUP_SIZE:
        raise ValueError(
            f"Each group needs at least {MIN_GROUP_SIZE} samples; "
            f"got young={len(young_samples)}, old={len(old_samples)}."
        )

    used = frame[young_samples + old_samples]
    if not already_log and looks_linear(used.to_numpy(dtype=float)):
        used = to_log2(used)
    if genes is not None:
        used = collapse_probes(used, genes)
    used = used.loc[np.isfinite(used.to_numpy(dtype=float)).all(axis=1)]
    if used.empty:
        return pd.DataFrame(
            columns=[
                "effect_size",
                "standard_error",
                "p_value",
                "q_value",
                "direction",
                "n_young",
                "n_old",
            ]
        )

    young = used[young_samples].to_numpy(dtype=float)
    old = used[old_samples].to_numpy(dtype=float)
    g, se = hedges_g_matrix(young, old)
    # The same normal approximation the random-effects pool downstream assumes.
    # Using a t-test here and a z there would make the per-study p-values
    # inconsistent with the pooled one.
    p = 2.0 * stats.norm.sf(np.abs(np.divide(g, se, out=np.zeros_like(g), where=se > 0)))

    return pd.DataFrame(
        {
            "effect_size": g,
            "standard_error": se,
            "p_value": p,
            "q_value": benjamini_hochberg(p),
            "direction": np.where(g >= 0, "up", "down"),
            "n_young": len(young_samples),
            "n_old": len(old_samples),
        },
        index=used.index,
    )
