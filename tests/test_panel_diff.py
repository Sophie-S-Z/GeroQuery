"""Tests for the living-evidence loop.

What matters here is not that a diff is produced but that it says the right
thing about a corpus that changed. The classification is the product: a rebuild
that reports "data refreshed" is worthless, and one that reports "FOXO3 now
excludes zero, and GDS1803 left the panel" is the thing nothing else in aging
publishes.
"""

from __future__ import annotations

import json

import pytest
import yaml

from geroquery.etl.panel_diff import (
    CHANGE_ORDER,
    EFFECT_MOVE_THRESHOLD,
    ESTIMATOR,
    MIN_CONTRASTS,
    Snapshot,
    diff_snapshots,
    render_markdown,
)

WORKFLOW = ".github/workflows/living-evidence.yml"


def _estimate(gene: str, g: float, low: float, high: float, k: int = 8) -> dict:
    return {
        "gene_id": gene,
        "species": "human",
        "omic_layer": "transcriptome",
        "g": g,
        "ci_low": low,
        "ci_high": high,
        "k": k,
        "i2": 20.0,
    }


def _snapshot(estimates: dict, accessions: list[str], taken_on: str = "2026-01-01") -> Snapshot:
    return Snapshot(
        taken_on=taken_on,
        manifest_version="2026.3",
        checksums_verified_on="2026-08-05",
        data_version="test+0000",
        accessions=accessions,
        series=accessions,
        n_contrasts=len(accessions),
        estimates=estimates,
    )


# --- classification ----------------------------------------------------------


@pytest.mark.parametrize(
    ("before", "after", "expected"),
    [
        # The headline event: a claim that could not be made now can be.
        ((0.20, -0.10, 0.50), (0.40, 0.12, 0.68), "newly_excludes_zero"),
        # And the one that matters just as much and is far more tempting to bury.
        ((0.40, 0.12, 0.68), (0.20, -0.10, 0.50), "no_longer_excludes_zero"),
        # A sign flip only counts when at least one side was actually a claim.
        ((0.40, 0.12, 0.68), (-0.40, -0.68, -0.12), "direction_reversed"),
        ((-0.02, -0.30, 0.26), (0.03, -0.25, 0.31), None),
        # Magnitude movement above the noise floor, verdict unchanged.
        ((0.40, 0.12, 0.68), (0.60, 0.30, 0.90), "effect_moved"),
        # Below the floor, and k unchanged: not news.
        ((0.40, 0.12, 0.68), (0.44, 0.16, 0.72), None),
    ],
)
def test_a_gene_is_labelled_by_what_actually_changed(before, after, expected):
    old = _snapshot({"G|human|transcriptome": _estimate("G", *before)}, ["GDS1"])
    new = _snapshot({"G|human|transcriptome": _estimate("G", *after)}, ["GDS1"])
    diff = diff_snapshots(old, new)
    if expected is None:
        assert diff.changes == []
    else:
        assert [c["kind"] for c in diff.changes] == [expected]


def test_more_contrasts_with_the_same_verdict_is_still_reported():
    old = _snapshot({"G|human|transcriptome": _estimate("G", 0.4, 0.12, 0.68, k=8)}, ["GDS1"])
    new = _snapshot({"G|human|transcriptome": _estimate("G", 0.4, 0.12, 0.68, k=11)}, ["GDS1"])
    (change,) = diff_snapshots(old, new).changes
    assert change["kind"] == "evidence_added"


def test_a_gene_losing_contrasts_is_reported_not_silently_dropped():
    old = _snapshot({"G|human|transcriptome": _estimate("G", 0.4, 0.12, 0.68, k=11)}, ["GDS1"])
    new = _snapshot({"G|human|transcriptome": _estimate("G", 0.4, 0.12, 0.68, k=8)}, ["GDS1"])
    (change,) = diff_snapshots(old, new).changes
    assert change["kind"] == "evidence_removed"


def test_the_effect_move_threshold_is_the_documented_one():
    """Pinned so the noise floor cannot drift without someone deciding to."""
    assert EFFECT_MOVE_THRESHOLD == pytest.approx(0.10)
    assert MIN_CONTRASTS == 3


# --- panel membership --------------------------------------------------------


def test_datasets_leaving_the_panel_are_reported_as_prominently_as_arrivals():
    """A withdrawal shrinks the evidence behind every gene in that dataset. A
    diff that only listed additions would hide the corpus getting smaller."""
    old = _snapshot({}, ["GDS1", "GDS2", "GDS3"])
    new = _snapshot({}, ["GDS1", "GDS3", "GDS4"])
    diff = diff_snapshots(old, new)

    assert diff.accessions_added == ["GDS4"]
    assert diff.accessions_removed == ["GDS2"]

    report = render_markdown(diff)
    assert "GDS4" in report and "GDS2" in report
    assert "shrinks the evidence" in report


def test_new_and_dropped_genes_are_distinguished_from_moved_ones():
    old = _snapshot({"A|human|transcriptome": _estimate("A", 0.3, -0.1, 0.7)}, ["GDS1"])
    new = _snapshot({"B|human|transcriptome": _estimate("B", 0.3, -0.1, 0.7)}, ["GDS1"])
    kinds = {c["kind"] for c in diff_snapshots(old, new).changes}
    assert kinds == {"new_gene", "dropped_gene"}


# --- the report --------------------------------------------------------------


def test_an_unchanged_rebuild_says_so_rather_than_producing_an_empty_report():
    same = {"G|human|transcriptome": _estimate("G", 0.4, 0.12, 0.68)}
    diff = diff_snapshots(_snapshot(same, ["GDS1"]), _snapshot(same, ["GDS1"]))
    assert diff.is_empty
    assert "Nothing moved" in render_markdown(diff)


def test_report_refuses_to_truncate_silently():
    """A section that says 'showing 2 of 40' tells a reader the corpus moved a
    lot. A quietly cut list tells them the opposite."""
    old = _snapshot(
        {f"G{i}|human|transcriptome": _estimate(f"G{i}", 0.2, -0.1, 0.5) for i in range(40)},
        ["GDS1"],
    )
    new = _snapshot(
        {f"G{i}|human|transcriptome": _estimate(f"G{i}", 0.4, 0.12, 0.68) for i in range(40)},
        ["GDS1"],
    )
    report = render_markdown(diff_snapshots(old, new), max_rows_per_kind=2)
    assert "Showing 2 of 40" in report


def test_report_forbids_auto_merge_and_hand_picking_in_words():
    """The guardrails live in the PR body because that is what a reviewer reads."""
    old = _snapshot({}, ["GDS1"])
    new = _snapshot({}, ["GDS1", "GDS2"])
    report = render_markdown(diff_snapshots(old, new))
    assert "must not be auto-merged" in report
    assert "rule-selected" in report
    assert "Do not hand-edit the panel" in report


def test_every_change_kind_has_a_label_and_an_ordering():
    from geroquery.etl.panel_diff import CHANGE_LABELS

    assert set(CHANGE_LABELS) == set(CHANGE_ORDER)
    # Verdict changes must outrank bookkeeping, or the report buries the news.
    assert CHANGE_ORDER.index("newly_excludes_zero") < CHANGE_ORDER.index("evidence_added")
    assert CHANGE_ORDER.index("no_longer_excludes_zero") < CHANGE_ORDER.index("effect_moved")


# --- round trip --------------------------------------------------------------


@pytest.mark.parametrize("suffix", [".json", ".json.gz"])
def test_snapshot_round_trips_through_disk(tmp_path, suffix):
    original = _snapshot({"G|human|transcriptome": _estimate("G", 0.4, 0.12, 0.68)}, ["GDS1"])
    path = original.write(tmp_path / f"snap{suffix}")
    restored = Snapshot.read(path)
    assert restored.to_dict() == original.to_dict()


def test_gzip_is_much_smaller(tmp_path):
    estimates = {
        f"G{i}|human|transcriptome": _estimate(f"G{i}", 0.1 * i, -0.2, 0.5) for i in range(500)
    }
    snapshot = _snapshot(estimates, ["GDS1"])
    plain = snapshot.write(tmp_path / "s.json").stat().st_size
    packed = snapshot.write(tmp_path / "s.json.gz").stat().st_size
    assert packed < plain / 4


# --- the workflow that drives it --------------------------------------------


def test_workflow_never_merges_and_only_ever_opens_a_pull_request():
    """The single most important property of the loop.

    An auto-updating corpus that can merge itself is a corpus nobody reviewed.
    """
    raw = open(WORKFLOW, encoding="utf-8").read()
    workflow = yaml.safe_load(raw)

    steps = workflow["jobs"]["rebuild-and-diff"]["steps"]
    uses = [step.get("uses", "") for step in steps]
    assert any("create-pull-request" in u for u in uses)
    assert not any("merge" in u.lower() for u in uses)
    assert "--auto" not in raw and "gh pr merge" not in raw


def test_workflow_runs_on_a_schedule_and_can_be_triggered_by_hand():
    workflow = yaml.safe_load(open(WORKFLOW, encoding="utf-8"))
    # PyYAML parses a bare `on:` key as the boolean True.
    triggers = workflow.get("on", workflow.get(True))
    assert "schedule" in triggers
    assert "workflow_dispatch" in triggers


def test_workflow_rebuilds_through_the_rule_rather_than_editing_the_panel():
    raw = open(WORKFLOW, encoding="utf-8").read()
    assert "build_signatures" in raw
    assert "panel_diff" in raw
    assert "rule-selected" in raw


def test_diff_json_is_serializable():
    """The workflow writes it to a file and attaches it to the PR."""
    old = _snapshot({"A|human|transcriptome": _estimate("A", 0.2, -0.1, 0.5)}, ["GDS1"])
    new = _snapshot({"A|human|transcriptome": _estimate("A", 0.5, 0.2, 0.8)}, ["GDS1", "GDS2"])
    payload = json.dumps(diff_snapshots(old, new).to_dict())
    assert "newly_excludes_zero" in payload


def test_the_diff_refuses_to_compare_across_estimators():
    """A method change must never be published as evidence moving.

    Swapping DerSimonian-Laird for the Hartung-Knapp interval retracted 1,164
    claims across this corpus without a byte of new data. If the loop had run
    across that change it would have opened a pull request reporting all of them
    as findings, and every one would have been wrong about why.
    """
    before = _snapshot({"A|human|transcriptome": _estimate("A", 0.5, 0.2, 0.8)}, ["GDS1"])
    after = _snapshot({"A|human|transcriptome": _estimate("A", 0.5, -0.1, 1.1)}, ["GDS1"])
    object.__setattr__(before, "estimator", "dersimonian_laird")
    object.__setattr__(after, "estimator", "hartung_knapp_modified")

    with pytest.raises(ValueError, match="not the evidence"):
        diff_snapshots(before, after)


def test_a_snapshot_records_which_estimator_produced_it():
    snapshot = _snapshot({"A|human|transcriptome": _estimate("A", 0.5, 0.2, 0.8)}, ["GDS1"])
    assert snapshot.to_dict()["estimator"] == ESTIMATOR
    # And a snapshot written before the field existed reads back as the
    # estimator that was in use then, rather than as unknown.
    legacy = dict(snapshot.to_dict())
    legacy.pop("estimator")
    assert Snapshot.from_dict(legacy).estimator == "dersimonian_laird"
