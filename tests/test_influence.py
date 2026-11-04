"""Influence diagnostics for a random-effects pool (§3.3 of the visual plan).

These are the estimator behind two charts: the leave-one-out strip beside the
forest plot (§4.4) and the corpus-influence bar chart (§4.5). The tests here are
the §3.3 contract — a planted outlier that leave-one-out must find, Cook's
distance must rank first, and whose leave-one-out estimate must equal an
independent from-scratch pool of the same subset. "It runs" is not the bar.
"""

from __future__ import annotations

import numpy as np

from geroquery.harmonize import leave_one_out, random_effects
from geroquery.harmonize.influence import corpus_influence


def test_leave_one_out_matches_an_independent_pool_of_the_same_subset():
    """The LOO estimate for study i equals ``random_effects`` on the other k-1.

    This is the clause that catches an influence measure that is merely
    plausible: the strip in §4.4 draws ``estimates[i]`` as "the pooled effect
    with contrast i removed", so it must literally be that pool, not an
    approximation of it.
    """
    effects = [0.40, 0.44, 0.38, 0.42, 0.95]  # the last one is the outlier
    ses = [0.10, 0.11, 0.09, 0.10, 0.12]

    loo = leave_one_out(effects, ses)

    for i in range(len(effects)):
        subset_e = [e for j, e in enumerate(effects) if j != i]
        subset_s = [s for j, s in enumerate(ses) if j != i]
        expected = random_effects(subset_e, subset_s).pooled_effect
        assert loo.estimates[i] == np.float64(expected) or np.isclose(
            loo.estimates[i], expected, atol=1e-12
        ), f"LOO[{i}]={loo.estimates[i]} != independent pool {expected}"


def test_the_planted_outlier_is_the_most_influential_contrast():
    """Leave-one-out finds the discordant contrast; Cook's distance ranks it #1.

    k-1 concordant contrasts near g=0.4 plus one deliberately discordant one at
    g=0.95. Removing the outlier must move the pooled estimate more than removing
    any concordant contrast, and Cook's distance (the squared standardised shift)
    must put it first.
    """
    effects = [0.40, 0.44, 0.38, 0.42, 0.95]
    ses = [0.10, 0.11, 0.09, 0.10, 0.12]
    outlier = 4

    loo = leave_one_out(effects, ses)

    assert int(np.argmax(np.abs(loo.shifts))) == outlier
    assert int(np.argmax(loo.cooks)) == outlier
    # Removing the outlier pulls the pooled effect down (the outlier was high),
    # so the full estimate exceeds the estimate without it.
    assert loo.full_effect > loo.estimates[outlier]


def test_cooks_distance_is_the_standardised_squared_shift():
    """D_i = (full - without_i)^2 / se_full^2 for an intercept-only pool.

    A defining identity, one derivative out from the leave-one-out estimates, so
    it fails if Cook's distance is computed from the wrong variance.
    """
    effects = [0.2, 0.5, -0.1, 0.3, 0.9]
    ses = [0.12, 0.10, 0.15, 0.11, 0.13]

    full = random_effects(effects, ses)
    loo = leave_one_out(effects, ses)

    expected = loo.shifts**2 / full.standard_error**2
    assert np.allclose(loo.cooks, expected)


def test_a_single_discordant_contrast_dominates_the_corpus_influence_ranking():
    """One contrast present with a discordant effect in many genes ranks first.

    This is the §4.5 claim made testable: aggregate each contrast's absolute
    leave-one-out shift across the whole (gene, study) table, and the contrast
    that is systematically bending the corpus rises to the top. "STUDY_BAD"
    carries a high effect in every gene while the others are concordant and low.
    """
    rows = []
    for g in range(30):
        base = 0.30 + 0.01 * (g % 5)
        for study, eff in (
            ("STUDY_A", base),
            ("STUDY_B", base + 0.02),
            ("STUDY_C", base - 0.02),
            ("STUDY_BAD", base + 0.90),
        ):
            rows.append(
                {
                    "gene_id": f"G{g}",
                    "study_id": study,
                    "effect_size": eff,
                    "standard_error": 0.12,
                }
            )

    ranking = corpus_influence(rows)

    assert ranking[0]["study_id"] == "STUDY_BAD"
    assert ranking[0]["total_abs_shift"] > ranking[1]["total_abs_shift"]
    assert ranking[0]["n_genes"] == 30
    # The ranking covers every contrast and nothing else.
    assert {r["study_id"] for r in ranking} == {"STUDY_A", "STUDY_B", "STUDY_C", "STUDY_BAD"}
