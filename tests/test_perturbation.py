"""Perturbation alignment — planted-truth tests for the scoring math.

The layer's whole claim is that a knockdown's downstream response can be
compared with the pooled aging direction. These tests plant responses whose
answer is known by construction, and pin the two corrections the real data
forced: excluding the perturbed gene's own readout, and removing the
response-magnitude confound.
"""

import numpy as np
import pytest

from geroquery.harmonize.perturbation import (
    adjust_for_magnitude,
    align,
    align_many,
)


def test_a_response_that_copies_the_aging_direction_aligns_at_one():
    direction = np.array([1.0, -2.0, 0.5, 3.0, -1.5])
    assert align(direction * 2.0, direction).alignment == pytest.approx(1.0)


def test_the_opposite_response_aligns_at_minus_one():
    direction = np.array([1.0, -2.0, 0.5, 3.0, -1.5])
    assert align(-direction, direction).alignment == pytest.approx(-1.0)


def test_an_unrelated_response_aligns_near_zero():
    rng = np.random.default_rng(0)
    direction = rng.normal(size=4000)
    assert abs(align(rng.normal(size=4000), direction).alignment) < 0.1


def test_alignment_ignores_the_scale_of_the_response():
    """Cosine, not covariance: a weak and a strong copy agree in direction."""
    direction = np.array([1.0, -2.0, 0.5, 3.0])
    weak = align(direction * 0.01, direction)
    strong = align(direction * 100.0, direction)
    assert weak.alignment == pytest.approx(strong.alignment)
    assert strong.magnitude > weak.magnitude  # but magnitude is kept, not lost


def test_a_flat_response_has_no_direction_rather_than_a_fake_one():
    direction = np.array([1.0, -2.0, 0.5, 3.0])
    assert np.isnan(align(np.ones(4), direction).alignment)


def test_a_flat_aging_direction_is_a_programmer_error_not_a_score():
    with pytest.raises(ValueError, match="no direction"):
        align(np.array([1.0, 2.0, 3.0]), np.ones(3))


def test_the_perturbed_gene_is_dropped_from_its_own_score():
    """Knocking down gene G drives G's own readout down.

    Scoring that entry would let a gene's own aging estimate answer a question
    about its downstream effect, so the index is excluded. Here the response is
    flat apart from the perturbed gene, so excluding it leaves no direction.
    """
    direction = np.array([2.0, 1.0, -1.0, -2.0])
    response = np.array([-5.0, 0.0, 0.0, 0.0])  # only the knocked-down gene moves
    assert align(response, direction).alignment < 0  # self entry drives the score
    assert np.isnan(align(response, direction, skip=0).alignment)


def test_align_many_matches_scoring_each_row_on_its_own():
    rng = np.random.default_rng(1)
    direction = rng.normal(size=50)
    responses = rng.normal(size=(7, 50))
    skips = [0, 3, None, 49, None, 12, 7]
    many = align_many(responses, direction, skip=skips)
    for row, skip, got in zip(responses, skips, many, strict=True):
        assert got.alignment == pytest.approx(align(row, direction, skip=skip).alignment)
        assert got.magnitude == pytest.approx(align(row, direction, skip=skip).magnitude)


def test_adjustment_removes_the_magnitude_confound():
    """On the real screens, alignment correlates ~0.4 with response magnitude.

    A knockdown that disturbs the cell more scores as more aging-like whatever
    it hit, so ranking on the raw score ranks damage. The adjusted score must
    carry none of that.
    """
    rng = np.random.default_rng(2)
    magnitude = np.exp(rng.normal(size=800))
    # Alignment driven entirely by magnitude, plus gene-specific noise.
    alignment = 0.3 * np.log(magnitude) + rng.normal(scale=0.05, size=800)
    z = adjust_for_magnitude(alignment, magnitude)
    assert abs(np.corrcoef(z, np.log(magnitude))[0, 1]) < 1e-8
    assert z.std() == pytest.approx(1.0)
    assert z.mean() == pytest.approx(0.0, abs=1e-8)


def test_adjustment_keeps_a_signal_that_is_not_magnitude():
    rng = np.random.default_rng(3)
    magnitude = np.exp(rng.normal(size=600))
    planted = rng.normal(size=600)
    z = adjust_for_magnitude(0.2 * np.log(magnitude) + planted, magnitude)
    assert np.corrcoef(z, planted)[0, 1] > 0.95


def test_adjustment_needs_enough_perturbations_to_fit_a_slope():
    with pytest.raises(ValueError, match="at least"):
        adjust_for_magnitude(np.array([0.1, 0.2]), np.array([1.0, 2.0]))


def test_adjustment_rejects_a_non_positive_magnitude():
    with pytest.raises(ValueError, match="positive"):
        adjust_for_magnitude(np.zeros(40), np.zeros(40))
