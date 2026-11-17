"""The perturbation exporter: scoring, replication gating and coverage.

Planted screens with a known answer, plus consistency checks against the
committed artifact.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from geroquery.etl.build_perturbation import (
    EXTREME_QUANTILE,
    aging_direction,
    combine,
    score_screen,
)
from geroquery.sources.perturbseq import Screen

RNG = np.random.default_rng(0)


def _pooled(**rows) -> pd.DataFrame:
    """A pooled frame with the columns the direction needs."""
    gene_ids = list(rows)
    return pd.DataFrame(
        {
            "gene_id": gene_ids,
            "species": ["human"] * len(gene_ids),
            "g": [rows[g][0] for g in gene_ids],
            "verdict": [rows[g][1] for g in gene_ids],
        }
    )


def _screen(responses, perturbed, *, readout, key="replogle_rpe1", controls=0) -> Screen:
    responses = np.asarray(responses, dtype=float)
    n = len(responses)
    is_control = np.array([False] * (n - controls) + [True] * controls)
    perturbed = np.array(list(perturbed) + [""] * controls)
    return Screen(
        key=key,
        cell_line="RPE1" if key == "replogle_rpe1" else "K562",
        responses=responses,
        readout_genes=np.array(readout),
        readout_symbols=np.array(readout),
        perturbed_genes=perturbed,
        symbols=np.array([g or "non-targeting" for g in perturbed]),
        has_effect=np.ones(n, dtype=bool),
        is_control=is_control,
    )


def test_the_aging_direction_uses_resolved_human_genes_only():
    pooled = _pooled(
        A=(0.5, "rises_with_age"),
        B=(-0.4, "falls_with_age"),
        C=(0.02, "no_evidence"),
    )
    pooled.loc[len(pooled)] = ["M", "mouse", 0.9, "rises_with_age"]
    direction = aging_direction(pooled)
    assert list(direction.index) == ["A", "B"]
    assert direction["A"] == pytest.approx(0.5)


def _wide_screen(n_genes=400, n_pert=60, key="replogle_rpe1", seed=1):
    rng = np.random.default_rng(seed)
    readout = [f"G{i}" for i in range(n_genes)]
    perturbed = [f"G{i}" for i in range(n_pert)]
    return readout, perturbed, rng


def test_a_knockdown_that_copies_the_aging_direction_scores_above_one_that_opposes_it():
    readout, perturbed, rng = _wide_screen()
    direction = pd.Series(rng.normal(size=len(readout)), index=readout)
    responses = rng.normal(scale=0.1, size=(len(perturbed), len(readout)))
    responses[0] += 3.0 * direction.to_numpy()  # planted: moves with aging
    responses[1] -= 3.0 * direction.to_numpy()  # planted: moves against it
    scored = score_screen(_screen(responses, perturbed, readout=readout), direction)
    assert scored.scores.loc["G0", "alignment"] > 0.8
    assert scored.scores.loc["G1", "alignment"] < -0.8
    assert scored.scores.loc["G0", "z"] > scored.scores.loc["G1", "z"]


def test_the_perturbed_genes_own_readout_cannot_drive_its_score():
    """A knockdown that moves nothing but its own target must not be scored as
    aging-like merely because that gene happens to rise with age."""
    readout, perturbed, rng = _wide_screen()
    direction = pd.Series(rng.normal(size=len(readout)), index=readout)
    responses = rng.normal(scale=0.1, size=(len(perturbed), len(readout)))
    responses[0][0] = -8.0  # G0's own readout drops hard; everything else is noise
    scored = score_screen(_screen(responses, perturbed, readout=readout), direction)
    assert abs(scored.scores.loc["G0", "alignment"]) < 0.2


def test_a_knockdown_with_no_downstream_movement_gets_no_score_at_all():
    """Its own target is excluded, so nothing is left to compare. That is an
    absent answer, not an alignment of zero."""
    readout, perturbed, rng = _wide_screen()
    direction = pd.Series(rng.normal(size=len(readout)), index=readout)
    responses = rng.normal(scale=0.1, size=(len(perturbed), len(readout)))
    responses[0] = 0.0
    responses[0][0] = -8.0
    scored = score_screen(_screen(responses, perturbed, readout=readout), direction)
    assert "G0" not in scored.scores.index


def test_scoring_needs_enough_shared_aging_genes():
    readout, perturbed, rng = _wide_screen(n_genes=50, n_pert=30)
    direction = pd.Series(rng.normal(size=10), index=readout[:10])
    with pytest.raises(ValueError, match="aging direction"):
        score_screen(
            _screen(rng.normal(size=(30, 50)), perturbed, readout=readout),
            direction,
        )


def _two_scored():
    """Two screens over the same genes, with planted agreement and discord."""
    readout, perturbed, rng = _wide_screen(n_genes=400, n_pert=60)
    direction = pd.Series(rng.normal(size=len(readout)), index=readout)
    out = []
    for key, seed in (("replogle_rpe1", 2), ("replogle_k562_gwps", 3)):
        r = np.random.default_rng(seed).normal(scale=0.1, size=(len(perturbed), len(readout)))
        r[0] += 4.0 * direction.to_numpy()  # agrees: toward, in both
        r[1] -= 4.0 * direction.to_numpy()  # agrees: away, in both
        sign = 1.0 if key == "replogle_rpe1" else -1.0
        r[2] += sign * 4.0 * direction.to_numpy()  # discordant
        out.append(score_screen(_screen(r, perturbed, readout=readout, key=key), direction))
    return out


def test_only_genes_both_screens_rank_extreme_and_agree_on_replicate():
    rpe1, k562 = _two_scored()
    combined = combine([rpe1, k562]).set_index("gene_id")
    assert combined.loc["G0", "replicates"]
    assert combined.loc["G0", "direction"] == "toward_aging"
    assert combined.loc["G1", "replicates"]
    assert combined.loc["G1", "direction"] == "away_from_aging"
    assert not combined.loc["G2", "replicates"]
    assert combined.loc["G2", "direction"] == "discordant"


def test_a_gene_only_one_screen_scored_is_reported_as_such_not_as_disagreement():
    rpe1, k562 = _two_scored()
    k562.scores.drop(index="G0", inplace=True)
    combined = combine([rpe1, k562]).set_index("gene_id")
    assert combined.loc["G0", "n_screens"] == 1
    assert not combined.loc["G0", "replicates"]
    assert combined.loc["G0", "direction"] == ""


def test_a_readout_gene_that_was_never_perturbed_carries_no_score():
    rpe1, k562 = _two_scored()
    combined = combine([rpe1, k562]).set_index("gene_id")
    assert combined.loc["G399", "n_screens"] == 0
    assert combined.loc["G399", "rpe1_readout"]
    assert not combined.loc["G399", "rpe1_perturbed"]
    assert pd.isna(combined.loc["G399", "rpe1_z"])


def _committed():
    root = Path(__file__).resolve().parents[1] / "frontend" / "public" / "data"
    meta = root / "meta.json"
    table = root / "perturbation.parquet"
    if not table.exists():
        return None, None
    return json.loads(meta.read_text())["perturbation"], pd.read_parquet(table)


def test_the_committed_coverage_partitions_every_resolved_gene():
    meta, _ = _committed()
    if meta is None:  # pragma: no cover - only before the layer is built
        pytest.skip("no perturbation layer committed")
    c = meta["coverage"]
    assert (
        c["resolved_scored"]
        + c["resolved_no_effect"]
        + c["resolved_readout_only"]
        + c["resolved_absent"]
        == c["resolved_genes"]
    )


def test_the_committed_table_agrees_with_its_own_summary():
    meta, table = _committed()
    if meta is None:  # pragma: no cover - only before the layer is built
        pytest.skip("no perturbation layer committed")
    assert int(table["replicates"].sum()) == meta["n_replicating"]
    assert int((table["n_screens"] == 2).sum()) == meta["n_scored_in_both"]
    # Every replicating gene is extreme in both screens and names its direction.
    replicating = table[table["replicates"]]
    assert set(replicating["direction"]) <= {"toward_aging", "away_from_aging"}
    for line in ("rpe1", "k562"):
        scored = table[table[f"{line}_scored"]][f"{line}_z"]
        cut = scored.abs().quantile(1 - EXTREME_QUANTILE)
        assert (replicating[f"{line}_z"].abs() >= cut).all()


def test_the_committed_screens_spread_wider_than_their_own_controls():
    """Guides that target nothing still produce an alignment. Knockdowns must
    spread wider, or the score is reading noise."""
    meta, _ = _committed()
    if meta is None:  # pragma: no cover - only before the layer is built
        pytest.skip("no perturbation layer committed")
    for screen in meta["screens"]:
        assert screen["knockdown_alignment_sd"] > 1.4 * screen["control_alignment_sd"]
