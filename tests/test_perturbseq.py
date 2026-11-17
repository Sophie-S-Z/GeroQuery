"""Perturb-seq screen loading and its selection rule, offline on a planted file.

Also pins the two facts the layer's credibility rests on, against the committed
artifact: knockdowns lower their own target, and the screens disagree about the
ranking they produce. Both skip when the artifact is not present.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from geroquery.exceptions import SourceError
from geroquery.sources import perturbseq

h5py = pytest.importorskip("h5py", reason="perturbation extra not installed")


def _plant(path: Path, *, responses, readout, index, effect_p) -> Path:
    with h5py.File(path, "w") as handle:
        handle["X"] = np.asarray(responses, dtype=np.float32)
        handle["var/gene_id"] = np.array([g.encode() for g in readout])
        handle["var/gene_name"] = np.arange(len(readout), dtype=np.int16)
        handle["var/__categories/gene_name"] = np.array(
            [f"SYM{i}".encode() for i in range(len(readout))]
        )
        handle["obs/gene_transcript"] = np.array([n.encode() for n in index])
        handle["obs/energy_test_p_value"] = np.asarray(effect_p, dtype=float)
    return path


def _simple(tmp_path: Path, **overrides) -> Path:
    kwargs = {
        "responses": [
            [-2.0, 2.0, 3.0],  # AAA: its own readout goes down, as it should
            [0.5, -1.0, -1.0],  # BBB: the energy test finds no effect
            [2.0, 2.0, 1.0],  # CCC: perturbed but not read out
            [0.1, -0.1, 0.2],  # non-targeting control
            [1.0, 1.0, -1.0],  # DDD: no Ensembl id upstream
        ],
        "readout": ["ENSG001", "ENSG002", "ENSG003"],
        "index": [
            "1_AAA_P1P2_ENSG001",
            "2_BBB_P1_ENSG002",
            "3_CCC_P1P2_ENSG009",
            "4_non-targeting_non-targeting_non-targeting",
            "5_DDD_P1P2_nan",
        ],
        "effect_p": [0.001, 0.5, 0.01, 0.001, 0.001],
    }
    kwargs.update(overrides)
    return _plant(tmp_path / "screen.h5ad", **kwargs)


def test_a_screen_carries_both_identities_of_the_perturbed_gene(tmp_path):
    screen = perturbseq.read_screen(_simple(tmp_path), "replogle_rpe1")
    assert screen.cell_line == "RPE1"
    assert list(screen.perturbed_genes) == ["ENSG001", "ENSG002", "ENSG009", "", ""]
    assert list(screen.symbols) == ["AAA", "BBB", "CCC", "non-targeting", "DDD"]
    assert list(screen.readout_genes) == ["ENSG001", "ENSG002", "ENSG003"]


def test_the_energy_test_decides_which_knockdowns_are_scored(tmp_path):
    screen = perturbseq.read_screen(_simple(tmp_path), "replogle_rpe1")
    assert list(screen.has_effect) == [True, False, True, True, True]
    # Of the five rows only AAA and CCC are real, joinable knockdowns that did
    # something: BBB had no effect, one row is a control, one has no gene id.
    assert list(screen.scorable) == [True, False, True, False, False]
    assert screen.n_scored == 2


def test_a_readout_gene_that_is_infinite_anywhere_is_dropped_everywhere(tmp_path):
    """Dropping per-perturbation would give each knockdown its own gene set,
    and scores computed over different genes are not comparable."""
    path = _simple(
        tmp_path,
        responses=[
            [1.0, np.inf, 3.0],
            [0.5, 0.0, -1.0],
            [2.0, 1.0, 2.0],
            [0.1, -0.1, 0.2],
            [1.0, 1.0, -1.0],
        ],
    )
    screen = perturbseq.read_screen(path, "replogle_rpe1")
    assert list(screen.readout_genes) == ["ENSG001", "ENSG003"]
    assert screen.responses.shape == (5, 2)
    assert np.isfinite(screen.responses).all()


def test_on_target_uses_only_scorable_knockdowns_that_were_also_read_out(tmp_path):
    screen = perturbseq.read_screen(_simple(tmp_path), "replogle_rpe1")
    # AAA alone qualifies: CCC is not read out, and the rest are not scorable.
    assert perturbseq.on_target(screen) == pytest.approx([-2.0])


def test_an_unknown_screen_key_names_the_known_ones(tmp_path):
    with pytest.raises(SourceError, match="replogle_rpe1"):
        perturbseq.read_screen(_simple(tmp_path), "not_a_screen")


def test_an_unparseable_perturbation_name_fails_loudly(tmp_path):
    path = _simple(
        tmp_path,
        index=["1_AAA_P1P2_ENSG001", "oops", "3_C_P1_ENSG009", "4_D_P1_ENSG010", "5_E_P1_nan"],
    )
    with pytest.raises(SourceError, match="guides"):
        perturbseq.read_screen(path, "replogle_rpe1")


def test_a_screen_with_no_joinable_target_is_an_error_not_an_empty_answer(tmp_path):
    path = _simple(
        tmp_path,
        index=[f"{i}_non-targeting_non-targeting_non-targeting" for i in range(5)],
    )
    with pytest.raises(SourceError, match="Ensembl"):
        perturbseq.read_screen(path, "replogle_rpe1")


def _committed() -> dict | None:
    path = Path(__file__).resolve().parents[1] / "frontend" / "public" / "data" / "meta.json"
    if not path.exists():
        return None
    return json.loads(path.read_text()).get("perturbation")


def test_the_committed_screens_pass_their_own_positive_control():
    """Knocking a gene down must lower that gene. If this fails, the layer's
    negative results are a loading bug rather than a finding."""
    meta = _committed()
    if not meta:  # pragma: no cover - only when the layer has not been built
        pytest.skip("no perturbation layer committed")
    for screen in meta["screens"]:
        assert screen["on_target_median"] < -0.3
        assert screen["on_target_fraction_down"] > 0.95


def test_the_committed_layer_records_that_the_screens_disagree():
    """The honest headline, pinned: the full ranking does not replicate."""
    meta = _committed()
    if not meta:  # pragma: no cover - only when the layer has not been built
        pytest.skip("no perturbation layer committed")
    assert 0.4 < meta["cross_screen_sign_agreement"] < 0.6
    assert meta["n_replicating"] < 0.05 * meta["n_scored_in_both"]
