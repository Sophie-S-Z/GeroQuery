"""Per-gene contrast statistics: the step from an expression matrix to a signature row.

The load-bearing test here is
``test_matrix_hedges_g_matches_scalar_implementation``. ``differential.py`` holds
a vectorized copy of the effect-size maths that already exists in
``effect_size.py``; without a test tying them together, a fix to one would
silently leave the other wrong.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from geroquery.harmonize.differential import (
    LOG_FLOOR,
    MIN_GROUP_SIZE,
    benjamini_hochberg,
    collapse_probes,
    contrast_effects,
    hedges_g_matrix,
    looks_linear,
    to_log2,
)
from geroquery.harmonize.effect_size import hedges_g

# ---- effect size ----------------------------------------------------------


def test_matrix_hedges_g_matches_scalar_implementation():
    rng = np.random.default_rng(11)
    young = rng.normal(0.0, 1.0, (40, 9))
    old = rng.normal(0.7, 1.3, (40, 12))

    g, se = hedges_g_matrix(young, old)
    for i in range(young.shape[0]):
        reference = hedges_g(young[i], old[i])
        assert g[i] == pytest.approx(reference.effect, rel=1e-12)
        assert se[i] == pytest.approx(reference.standard_error, rel=1e-12)


def test_matrix_hedges_g_sign_convention_is_positive_when_higher_in_old():
    young = np.zeros((1, 4))
    young[0] = [1.0, 1.1, 0.9, 1.0]
    old = np.array([[3.0, 3.1, 2.9, 3.0]])
    g, _ = hedges_g_matrix(young, old)
    assert g[0] > 0


def test_zero_variance_gene_gets_a_finite_zero_effect():
    """A gene with the same value everywhere has no detectable difference, which
    is a finite statement — not an infinity that would then dominate a pool."""
    constant = np.full((1, 5), 7.0)
    g, se = hedges_g_matrix(constant, constant.copy())
    assert g[0] == 0.0
    assert np.isfinite(se[0])


# ---- multiple testing -----------------------------------------------------


def test_benjamini_hochberg_matches_worked_example():
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    q = benjamini_hochberg(p)
    assert q == pytest.approx([0.05, 0.05, 0.05, 0.05, 0.05])


def test_benjamini_hochberg_is_monotone_and_bounded():
    rng = np.random.default_rng(3)
    p = rng.random(500)
    q = benjamini_hochberg(p)
    assert (q >= p - 1e-12).all()
    assert (q <= 1.0).all()
    order = np.argsort(p)
    assert (np.diff(q[order]) >= -1e-12).all()


def test_benjamini_hochberg_excludes_nan_from_the_denominator():
    """Treating a missing p-value as 1.0 would inflate n and shrink every q."""
    with_nan = benjamini_hochberg(np.array([0.01, 0.02, np.nan]))
    without = benjamini_hochberg(np.array([0.01, 0.02]))
    assert np.isnan(with_nan[2])
    assert with_nan[:2] == pytest.approx(without)


# ---- scale detection ------------------------------------------------------


def test_linear_intensities_are_detected_and_log_scaled():
    linear = pd.DataFrame(np.array([[10.0, 4000.0], [50.0, 9000.0]]))
    assert looks_linear(linear.to_numpy())
    assert not looks_linear(np.array([[2.0, 14.0], [3.0, 11.0]]))


def test_log2_floors_non_positive_values_rather_than_dropping_them():
    frame = pd.DataFrame([[-5.0, 0.0, 4.0]])
    out = to_log2(frame)
    assert out.shape == frame.shape
    assert out.iloc[0, 0] == pytest.approx(np.log2(LOG_FLOOR))
    assert out.iloc[0, 2] == pytest.approx(2.0)


# ---- probe collapse -------------------------------------------------------


def test_collapse_probes_keeps_the_highest_mean_probe_per_gene():
    frame = pd.DataFrame(
        [[1.0, 1.0], [9.0, 11.0], [5.0, 5.0]],
        index=["p_dim", "p_bright", "p_other"],
    )
    genes = pd.Series({"p_dim": "A", "p_bright": "A", "p_other": "B"})
    out = collapse_probes(frame, genes)
    assert list(out.index) == ["A", "B"]
    assert out.loc["A"].tolist() == [9.0, 11.0]


def test_collapse_probes_drops_unannotated_probes():
    frame = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=["p1", "p2"])
    out = collapse_probes(frame, pd.Series({"p1": "A"}))
    assert list(out.index) == ["A"]


# ---- the whole contrast ---------------------------------------------------


def _matrix(rng, n_genes=30, n_samples=10, shift=0.0):
    data = rng.normal(0.0, 1.0, (n_genes, n_samples))
    return pd.DataFrame(data + shift, columns=[f"GSM{i}" for i in range(n_samples)])


def test_contrast_effects_recovers_a_planted_shift():
    rng = np.random.default_rng(5)
    # Large groups on purpose: at the n=3-15 the real GEO panel supplies, the
    # sampling error on a single gene's g is of the same order as the effect —
    # which is the finding in docs/RESULTS_GEO_SIGNATURES.md, not something a
    # unit test should have to tolerate.
    frame = _matrix(rng, n_genes=1, n_samples=4000)
    young = list(frame.columns[:2000])
    old = list(frame.columns[2000:])
    frame.loc[:, old] += 1.5  # a clean 1.5 SD increase in the old group

    out = contrast_effects(frame, young, old, already_log=True)
    assert out.loc[0, "effect_size"] == pytest.approx(1.5, abs=0.1)
    assert out.loc[0, "direction"] == "up"
    assert out.loc[0, "p_value"] < 1e-6
    assert out.loc[0, "n_young"] == 2000
    assert out.loc[0, "n_old"] == 2000


def test_contrast_effects_refuses_overlapping_groups():
    rng = np.random.default_rng(6)
    frame = _matrix(rng)
    with pytest.raises(ValueError, match="share samples"):
        contrast_effects(frame, list(frame.columns[:5]), list(frame.columns[4:]))


def test_contrast_effects_refuses_groups_below_the_minimum():
    rng = np.random.default_rng(7)
    frame = _matrix(rng)
    small = list(frame.columns[: MIN_GROUP_SIZE - 1])
    rest = list(frame.columns[MIN_GROUP_SIZE:])
    with pytest.raises(ValueError, match="at least"):
        contrast_effects(frame, small, rest)


def test_contrast_effects_ignores_non_sample_columns():
    """The GDS table carries annotation columns beside the sample values; naming
    the samples explicitly must be enough to exclude them."""
    rng = np.random.default_rng(8)
    frame = _matrix(rng, n_genes=5, n_samples=8)
    frame.insert(0, "IDENTIFIER", ["G1", "G2", "G3", "G4", "G5"])
    out = contrast_effects(
        frame, list(frame.columns[1:5]), list(frame.columns[5:]), already_log=True
    )
    assert len(out) == 5
    assert np.isfinite(out["effect_size"]).all()
