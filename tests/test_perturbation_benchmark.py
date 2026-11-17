"""The knockdown-prediction benchmark, on planted data with a known answer.

The benchmark's job is to report that a trained model barely beats "predict the
average knockdown". That conclusion is only worth anything if the benchmark can
see a gain when one is really there, so both cases are planted here.
"""

from __future__ import annotations

import numpy as np
import pytest

from geroquery.etl.perturbation_benchmark import benchmark_matrix, run, trainable
from geroquery.sources.perturbseq import Screen


def _planted(n=160, k=12, readout=90, signal=0.0, seed=0):
    """Targets = one shared average profile + ``signal`` x a gene-specific part."""
    rng = np.random.default_rng(seed)
    features = rng.normal(size=(n, k))
    average = rng.normal(scale=2.0, size=readout)
    coefficients = rng.normal(size=(k, readout))
    targets = (
        average + signal * (features @ coefficients) + rng.normal(scale=0.2, size=(n, readout))
    )
    return features, targets


def test_the_benchmark_sees_a_gain_when_the_response_really_is_gene_specific():
    base, model, _, _ = benchmark_matrix(*_planted(signal=1.0), seed=1)
    assert model.mean() > base.mean() + 0.10


def test_the_benchmark_reports_no_gain_when_every_knockdown_does_the_same_thing():
    base, model, _, _ = benchmark_matrix(*_planted(signal=0.0), seed=1)
    assert abs(model.mean() - base.mean()) < 0.05
    # And the shared profile alone already gets most directions right, which is
    # why a bare accuracy figure says nothing about gene-specific learning.
    assert base.mean() > 0.8


def test_the_baseline_never_sees_the_gene_it_is_scored_on():
    """A leaked average would make the baseline unbeatable and the comparison
    meaningless, so the held-out gene must not be in its own training mean.

    Four genes, two folds: one gene's profile is the loud negation of the other
    three. Trained without it, the baseline predicts the shared profile and gets
    that gene's directions exactly wrong. Were it in its own training mean, it
    would dominate those two rows and get them exactly right.
    """
    rng = np.random.default_rng(0)
    shared = rng.normal(size=60)
    targets = np.vstack([-10.0 * shared, shared, shared, shared])
    base, _, _, _ = benchmark_matrix(rng.normal(size=(4, 3)), targets, folds=2, seed=0)
    assert base[0] == 0.0


def test_a_screen_needs_genes_that_were_both_knocked_down_and_read_out():
    screen = Screen(
        key="replogle_rpe1",
        cell_line="RPE1",
        responses=np.ones((2, 3)),
        readout_genes=np.array(["A", "B", "C"]),
        readout_symbols=np.array(["A", "B", "C"]),
        perturbed_genes=np.array(["X", "Y"]),
        symbols=np.array(["X", "Y"]),
        has_effect=np.ones(2, dtype=bool),
        is_control=np.zeros(2, dtype=bool),
    )
    with pytest.raises(ValueError, match="knocked down and read out"):
        trainable(screen)


def test_too_few_genes_to_cross_validate_is_an_error_not_a_number():
    with pytest.raises(ValueError, match="at least"):
        benchmark_matrix(*_planted(n=6), folds=5)


def test_run_reports_the_gain_and_whether_the_model_won():
    rng = np.random.default_rng(3)
    readout = np.array([f"G{i}" for i in range(80)])
    responses = rng.normal(size=(60, 80))
    screen = Screen(
        key="replogle_rpe1",
        cell_line="RPE1",
        responses=responses,
        readout_genes=readout,
        readout_symbols=readout,
        perturbed_genes=readout[:60],
        symbols=readout[:60],
        has_effect=np.ones(60, dtype=bool),
        is_control=np.zeros(60, dtype=bool),
    )
    result = run(screen, folds=3)
    assert result.n_genes == 60
    assert 0.0 <= result.baseline_accuracy <= 1.0
    assert result.model_wins == (result.gain_pp > 0)
