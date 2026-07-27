"""Dynamic Network Biomarker (DNB / l-DNB) composite index.

DNB theory (Chen et al., 2012, *Sci Rep*; single-sample landscape l-DNB, Liu et
al., 2017, *PLoS Comput Biol*) identifies a pre-transition ("tipping") state from
three signals within a co-fluctuating group of markers:

  * rising standard deviation *within* the group,
  * rising correlation *within* the group, and
  * falling correlation *between* the group and the background.

The composite is ``I = SD_in * PCC_in / PCC_out``. It is a more structured
early-warning construct than raw variance because it separates a coordinated
"leading" module from the rest of the system.

This is a **pure-Python, cross-sectional adaptation** for the small clinical
marker panels GeroQuery handles. Two honest caveats:

  * With no external background (a single small panel), ``PCC_out`` is undefined,
    so the index reduces to the numerator ``SD_in * PCC_in`` — clearly the
    "no-background" reduction, not the full three-term index.
  * When a module/background split *can* be formed (>= 2 markers each), we derive
    it heuristically from the correlation structure and compute the full index.

Like the age-stratified dispersion indicator, this is exploratory on
cross-sectional data and is not a validated tipping-point detector.
"""

from __future__ import annotations

import numpy as np


def _abs_corr(X: np.ndarray) -> np.ndarray:
    """Absolute Pearson correlation matrix of columns; NaNs -> 0."""
    if X.shape[0] < 3 or X.shape[1] < 2:
        return np.zeros((X.shape[1], X.shape[1]))
    c = np.abs(np.corrcoef(X, rowvar=False))
    return np.nan_to_num(c, nan=0.0)


def auto_module(X: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    """Split markers into a co-fluctuating module vs background from correlation.

    Markers whose mean absolute correlation to the others is above the median are
    the module; the rest are background. Returns ``None`` if either side would
    have fewer than two markers (so the full three-term index is not meaningful).
    """
    p = X.shape[1]
    if p < 4:
        return None
    c = _abs_corr(X)
    np.fill_diagonal(c, 0.0)
    connectivity = c.sum(axis=1) / (p - 1)
    thresh = float(np.median(connectivity))
    module = np.where(connectivity > thresh)[0]
    background = np.where(connectivity <= thresh)[0]
    if module.size < 2 or background.size < 2:
        return None
    return module, background


def dnb_index(X: np.ndarray, partition: tuple[np.ndarray, np.ndarray] | None = None) -> float:
    """DNB composite for one block of samples ``X`` (n_samples, n_markers).

    If ``partition`` (module_idx, background_idx) is given the full three-term
    index is computed; otherwise the no-background reduction ``SD_in * PCC_in``.
    """
    X = np.asarray(X, dtype=float)
    if X.shape[0] < 3 or X.shape[1] < 2:
        return float("nan")

    if partition is None:
        module = np.arange(X.shape[1])
        background = None
    else:
        module, background = partition

    sd_in = float(np.mean(np.std(X[:, module], axis=0, ddof=1)))
    cm = _abs_corr(X[:, module])
    iu = np.triu_indices_from(cm, k=1)
    pcc_in = float(np.mean(cm[iu])) if iu[0].size else 0.0

    if background is None or len(background) == 0:
        return sd_in * pcc_in

    cfull = _abs_corr(X)
    pcc_out = float(np.mean(cfull[np.ix_(module, background)]))
    pcc_out = max(pcc_out, 1e-6)  # guard division
    return sd_in * pcc_in / pcc_out


def dnb_per_block(
    X: np.ndarray, blocks: list[np.ndarray], partition: tuple[np.ndarray, np.ndarray] | None
) -> list[float]:
    """DNB index for each block (e.g. age stratum), using a fixed partition."""
    return [dnb_index(X[idx], partition) for idx in blocks]
