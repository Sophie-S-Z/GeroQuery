"""Batch-effect handling.

A deliberately simple, transparent location adjustment: for each feature, remove
the per-batch mean shift while preserving the grand mean and any biological
covariate structure supplied. This is the mean-only special case of ComBat /
limma's ``removeBatchEffect``; the production ETL swaps in full ComBat for bulk
and Harmony/scVI for single-cell. Kept simple here because a correct, legible
adjustment beats a black box for a portfolio core.
"""

from __future__ import annotations

import numpy as np


def remove_batch_effect(
    matrix: np.ndarray, batch: np.ndarray, biological_group: np.ndarray | None = None
) -> np.ndarray:
    """Return ``matrix`` with additive per-batch shifts removed.

    Parameters
    ----------
    matrix: (n_features, n_samples)
    batch: (n_samples,) batch label per sample
    biological_group: optional (n_samples,) label to *preserve* (e.g. young/old).
        When given, batch means are estimated within each biological group and
        then re-centered, so genuine biological differences are not flattened.
    """
    X = np.asarray(matrix, dtype=float)
    batch = np.asarray(batch)
    if X.ndim != 2:
        raise ValueError("matrix must be 2-D (features x samples).")
    if batch.shape[0] != X.shape[1]:
        raise ValueError("batch length must equal number of samples (columns).")

    corrected = X.copy()
    grand_mean = X.mean(axis=1, keepdims=True)

    if biological_group is None:
        for b in np.unique(batch):
            cols = batch == b
            batch_mean = X[:, cols].mean(axis=1, keepdims=True)
            corrected[:, cols] += grand_mean - batch_mean
        return corrected

    # Preserve biological structure: within each biological group, align batches
    # to that group's mean.
    biological_group = np.asarray(biological_group)
    for g in np.unique(biological_group):
        g_cols = biological_group == g
        g_mean = X[:, g_cols].mean(axis=1, keepdims=True)
        for b in np.unique(batch[g_cols]):
            cols = g_cols & (batch == b)
            if not cols.any():
                continue
            batch_mean = X[:, cols].mean(axis=1, keepdims=True)
            corrected[:, cols] += g_mean - batch_mean
    return corrected


def batch_variance(matrix: np.ndarray, batch: np.ndarray) -> float:
    """Mean across features of the variance of per-batch feature means.

    A scalar summary of how much batch structure remains; used in tests to show
    correction shrinks it.
    """
    X = np.asarray(matrix, dtype=float)
    batch = np.asarray(batch)
    per_feature = []
    for f in range(X.shape[0]):
        means = [X[f, batch == b].mean() for b in np.unique(batch)]
        per_feature.append(np.var(means))
    return float(np.mean(per_feature))
