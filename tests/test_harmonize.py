"""M3 harmonize — scientific correctness on planted-truth synthetic data."""

import numpy as np
import pytest

from geroquery.harmonize import (
    batch_variance,
    hedges_g,
    random_effects,
    remove_batch_effect,
)


def test_hedges_g_recovers_planted_effect():
    rng = np.random.default_rng(0)
    young = rng.normal(0.0, 1.0, 400)
    old = rng.normal(2.0, 1.0, 400)  # true standardized diff ~= 2
    es = hedges_g(young, old)
    assert es.effect == pytest.approx(2.0, abs=0.2)
    assert es.direction == "up"
    assert es.standard_error > 0


def test_hedges_g_sign_flips_for_down_regulation():
    rng = np.random.default_rng(1)
    young = rng.normal(2.0, 1.0, 200)
    old = rng.normal(0.0, 1.0, 200)
    es = hedges_g(young, old)
    assert es.effect < 0
    assert es.direction == "down"


def test_random_effects_homogeneous_studies():
    pooled = random_effects([0.5, 0.5, 0.5, 0.5], [0.1, 0.1, 0.1, 0.1])
    assert pooled.pooled_effect == pytest.approx(0.5, abs=1e-6)
    assert pooled.i2 == pytest.approx(0.0, abs=1e-6)
    assert pooled.tau2 == pytest.approx(0.0, abs=1e-9)
    assert pooled.ci_low < 0.5 < pooled.ci_high
    assert pooled.n_studies == 4


def test_random_effects_detects_heterogeneity():
    pooled = random_effects([0.1, 0.9, 0.5, 0.0], [0.1, 0.1, 0.1, 0.1])
    assert pooled.tau2 > 0
    assert pooled.i2 > 50  # substantial heterogeneity


def test_random_effects_single_study_passthrough():
    pooled = random_effects([0.7], [0.2])
    assert pooled.pooled_effect == pytest.approx(0.7)
    assert pooled.standard_error == pytest.approx(0.2)
    assert pooled.i2 == 0.0


def test_random_effects_rejects_nonpositive_se():
    with pytest.raises(ValueError):
        random_effects([0.5, 0.5], [0.1, 0.0])


def test_batch_correction_reduces_batch_signal_and_preserves_biology():
    rng = np.random.default_rng(3)
    n_features, per = 30, 10
    # 40 samples: 2 biological groups x 2 batches.
    group = np.array(["young"] * (2 * per) + ["old"] * (2 * per))
    batch = np.array((["A"] * per + ["B"] * per) * 2)
    base = rng.normal(0, 1, (n_features, 40))
    group_effect = np.where(group == "old", 3.0, 0.0)  # real biology
    batch_effect = np.where(batch == "B", 5.0, 0.0)  # nuisance
    X = base + group_effect + batch_effect

    before = batch_variance(X, batch)
    corrected = remove_batch_effect(X, batch, biological_group=group)
    after = batch_variance(corrected, batch)
    assert after < before * 0.2  # batch signal largely removed

    # Biological difference (old - young) is preserved.
    diff_before = X[:, group == "old"].mean() - X[:, group == "young"].mean()
    diff_after = corrected[:, group == "old"].mean() - corrected[:, group == "young"].mean()
    assert diff_after == pytest.approx(diff_before, abs=0.3)
