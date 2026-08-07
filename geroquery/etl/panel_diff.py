"""The living-evidence loop: what changed in the evidence, not just in the data.

Run:
    python -m geroquery.etl.panel_diff --snapshot evidence/latest.json
    python -m geroquery.etl.panel_diff --before evidence/latest.json \\
                                       --after  evidence/candidate.json \\
                                       --markdown evidence/CHANGELOG.md

A scheduled rebuild that reports "the data was refreshed" is a cron job. This
reports **which pooled estimates moved and by how much**, which is the thing a
reader of an aging database has never been able to see:

    CDKN2A   g +0.07 [-0.20, +0.35] k=14  ->  g +0.11 [-0.14, +0.36] k=16   still null
    FOXO3    g -0.31 [-0.62, +0.01] k=7   ->  g -0.38 [-0.65, -0.11] k=9    newly excludes zero

Curated databases publish *entries*. This publishes *how the evidence changed*,
with the interval, which is the only honest way to say that a claim strengthened.

Two rules the loop must not be allowed to break, both enforced here rather than
left to a reviewer's memory:

1. **The panel stays rule-selected.** Datasets enter because they satisfy the
   GEO DataSets query and the contrast rules, never because someone liked the
   result. :func:`diff_snapshots` reports panel membership changes but has no
   mechanism to include or exclude a dataset, and the workflow that calls it
   opens a pull request rather than committing.
2. **Datasets leaving are as newsworthy as datasets arriving.** A GDS being
   withdrawn or reclassified upstream silently shrinks the evidence behind every
   gene in it. A diff that only reported additions would hide that.
"""

from __future__ import annotations

import argparse
import glob as glob_module
import json
import traceback
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

from ..harmonize import random_effects
from ..sources.manifest import MANIFEST_VERSION, VERIFIED_ON
from ..store import GeroStore

# Process exit codes. "The evidence moved" and "something broke" must never
# share one, because the workflow branches on it: exit 1 promotes the candidate
# over the baseline, and a crash reaching that branch overwrites the baseline
# with a snapshot that was never diffed.
EXIT_UNCHANGED = 0
EXIT_CHANGED = 1
EXIT_ERROR = 2

# Genes pooled from fewer contrasts than this are omitted from the snapshot.
# Not a significance filter: with one or two contrasts the random-effects
# between-study variance is not estimable, so the interval is not a summary of
# anything and comparing it across rebuilds would report noise as news.
MIN_CONTRASTS = 3

# Identifies the interval a snapshot's estimates were computed with. Bump it
# whenever the estimator changes, and reset the baseline in the same commit —
# diffing across estimators reports a change of method as a change of evidence.
ESTIMATOR = "hartung_knapp_modified"

# A pooled effect must move by at least this much (in Hedges' g) before it is
# reported as a magnitude change. Below it, the movement is rebuild jitter and
# reporting it would bury the real changes.
EFFECT_MOVE_THRESHOLD = 0.10

# Change kinds, most newsworthy first. The ordering is the report ordering.
CHANGE_ORDER: tuple[str, ...] = (
    "newly_excludes_zero",
    "no_longer_excludes_zero",
    "direction_reversed",
    "effect_moved",
    "evidence_added",
    "evidence_removed",
    "new_gene",
    "dropped_gene",
)

CHANGE_LABELS: dict[str, str] = {
    "newly_excludes_zero": "Interval now excludes zero (a claim strengthened)",
    "no_longer_excludes_zero": "Interval now includes zero (a claim weakened)",
    "direction_reversed": "Pooled direction reversed",
    "effect_moved": f"Pooled effect moved by more than {EFFECT_MOVE_THRESHOLD}",
    "evidence_added": "More contrasts, same verdict",
    "evidence_removed": "Fewer contrasts, same verdict",
    "new_gene": "Newly measurable (reached the contrast minimum)",
    "dropped_gene": "No longer measurable (fell below the contrast minimum)",
}


def _excludes_zero(low: float, high: float) -> bool:
    return low > 0 or high < 0


@dataclass(frozen=True)
class Snapshot:
    """Pooled estimates plus the panel that produced them, at one point in time."""

    taken_on: str
    manifest_version: str
    checksums_verified_on: str
    data_version: str
    accessions: list[str]
    series: list[str]
    n_contrasts: int
    estimates: dict[str, dict[str, Any]]
    # Which interval produced these estimates. Without it a change of estimator
    # is indistinguishable from evidence moving: swapping DerSimonian-Laird for
    # Hartung-Knapp retracted 1,164 claims across this corpus, and a diff that
    # did not know would have published every one of them as a finding.
    estimator: str = ESTIMATOR

    def to_dict(self) -> dict[str, Any]:
        return {
            "taken_on": self.taken_on,
            "manifest_version": self.manifest_version,
            "checksums_verified_on": self.checksums_verified_on,
            "data_version": self.data_version,
            "accessions": self.accessions,
            "series": self.series,
            "n_contrasts": self.n_contrasts,
            "min_contrasts": MIN_CONTRASTS,
            "n_estimates": len(self.estimates),
            "estimator": self.estimator,
            "estimates": self.estimates,
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> Snapshot:
        return cls(
            taken_on=raw.get("taken_on", "unknown"),
            manifest_version=raw.get("manifest_version", "unknown"),
            checksums_verified_on=raw.get("checksums_verified_on", "unknown"),
            data_version=raw.get("data_version", "unknown"),
            accessions=list(raw.get("accessions", [])),
            series=list(raw.get("series", [])),
            n_contrasts=int(raw.get("n_contrasts", 0)),
            estimates=dict(raw.get("estimates", {})),
            # Snapshots written before the estimator was recorded were all
            # DerSimonian-Laird; naming that is better than calling it unknown.
            estimator=raw.get("estimator", "dersimonian_laird"),
        )

    def write(self, path: Path) -> Path:
        """Write the snapshot, gzipping when the path says ``.gz``.

        The baseline is ~42,000 pooled estimates — 9.9 MB as plain JSON, about
        1 MB gzipped. The loop is designed to **overwrite** one baseline rather
        than accumulate a snapshot per run, so repository growth is bounded by
        the changelog (a few KB a month) and not by the corpus.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(self.to_dict(), indent=1, sort_keys=True)
        if path.suffix == ".gz":
            import gzip

            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(payload)
        else:
            path.write_text(payload, encoding="utf-8")
        return path

    @classmethod
    def read(cls, path: Path) -> Snapshot:
        path = Path(path)
        if path.suffix == ".gz":
            import gzip

            with gzip.open(path, "rt", encoding="utf-8") as handle:
                return cls.from_dict(json.load(handle))
        return cls.from_dict(json.loads(path.read_text(encoding="utf-8")))


def take_snapshot(store: GeroStore | None = None, *, today: str | None = None) -> Snapshot:
    """Pool every gene with at least :data:`MIN_CONTRASTS` contrasts.

    Reads the partitioned Parquet directly rather than looping the service one
    gene at a time: the panel carries ~46,000 genes, and a per-gene round trip
    through the resolver would make the loop too slow to run on a schedule,
    which is the one thing it has to do.
    """
    store = (store or GeroStore()).ensure_built()
    # recursive=True is required for ``**`` to cross the species/omic_layer
    # partition directories; without it this silently matches nothing.
    files = sorted(
        glob_module.glob((store.sig_dir / "**" / "*.parquet").as_posix(), recursive=True)
    )
    if not files:
        raise FileNotFoundError(f"No signature parquet under {store.sig_dir}.")

    # Read the dataset root, not the individual files: `species` and
    # `omic_layer` are hive partition keys living in the directory names, so a
    # per-file read loses exactly the two columns the grouping needs.
    frame = pd.read_parquet(store.sig_dir)
    for column in ("species", "omic_layer"):
        if column not in frame.columns:
            raise KeyError(
                f"Partition column {column!r} was not reconstructed from "
                f"{store.sig_dir}. Read the dataset root, not the parquet files."
            )
    frame = frame.dropna(subset=["effect_size", "standard_error"])
    frame["species"] = frame["species"].astype(str)
    frame["omic_layer"] = frame["omic_layer"].astype(str)

    estimates: dict[str, dict[str, Any]] = {}
    grouped = frame.groupby(["gene_id", "species", "omic_layer"], observed=True)
    for (gene_id, species, omic_layer), rows in grouped:
        if len(rows) < MIN_CONTRASTS:
            continue
        pooled = random_effects(rows["effect_size"].tolist(), rows["standard_error"].tolist())
        estimates[f"{gene_id}|{species}|{omic_layer}"] = {
            "gene_id": gene_id,
            "species": species,
            "omic_layer": omic_layer,
            "g": round(pooled.pooled_effect, 4),
            "ci_low": round(pooled.ci_low, 4),
            "ci_high": round(pooled.ci_high, 4),
            "k": int(pooled.n_studies),
            "i2": round(pooled.i2, 1),
        }

    studies = store.list_studies()
    return Snapshot(
        taken_on=today or date.today().isoformat(),
        manifest_version=MANIFEST_VERSION,
        checksums_verified_on=VERIFIED_ON,
        data_version=store.version(),
        accessions=sorted({s.study_id for s in studies if s.study_id}),
        series=sorted({s.series_id for s in studies if s.series_id}),
        n_contrasts=len(studies),
        estimates=estimates,
    )


@dataclass
class PanelDiff:
    """What changed between two snapshots."""

    before: Snapshot
    after: Snapshot
    accessions_added: list[str] = field(default_factory=list)
    accessions_removed: list[str] = field(default_factory=list)
    changes: list[dict[str, Any]] = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        return not (self.changes or self.accessions_added or self.accessions_removed)

    def by_kind(self) -> dict[str, list[dict[str, Any]]]:
        out: dict[str, list[dict[str, Any]]] = {kind: [] for kind in CHANGE_ORDER}
        for change in self.changes:
            out[change["kind"]].append(change)
        return out

    def to_dict(self) -> dict[str, Any]:
        counts = {kind: len(rows) for kind, rows in self.by_kind().items() if rows}
        return {
            "before": {"taken_on": self.before.taken_on, "data_version": self.before.data_version},
            "after": {"taken_on": self.after.taken_on, "data_version": self.after.data_version},
            "accessions_added": self.accessions_added,
            "accessions_removed": self.accessions_removed,
            "n_changes": len(self.changes),
            "counts": counts,
            "changes": self.changes,
        }


def _classify(old: dict[str, Any], new: dict[str, Any]) -> str | None:
    """Which single label best describes this gene's movement.

    Ordered by newsworthiness: a gene whose interval crossed zero has changed
    what it claims, which matters more than that its k rose. Only one label is
    assigned, so the report reads as a list of distinct events rather than the
    same gene appearing five times.
    """
    was, now = _excludes_zero(old["ci_low"], old["ci_high"]), _excludes_zero(
        new["ci_low"], new["ci_high"]
    )
    if now and not was:
        return "newly_excludes_zero"
    if was and not now:
        return "no_longer_excludes_zero"
    # A direction reversal only means something if at least one side was a claim.
    if (was or now) and (old["g"] > 0) != (new["g"] > 0):
        return "direction_reversed"
    if abs(new["g"] - old["g"]) >= EFFECT_MOVE_THRESHOLD:
        return "effect_moved"
    if new["k"] > old["k"]:
        return "evidence_added"
    if new["k"] < old["k"]:
        return "evidence_removed"
    return None


def diff_snapshots(before: Snapshot, after: Snapshot) -> PanelDiff:
    """Compare two snapshots and classify every gene that moved.

    Refuses to compare across estimators. The whole output of this loop is a
    claim that the *evidence* moved, and a diff between a DerSimonian-Laird
    baseline and a Hartung-Knapp rebuild would publish 1,164 retracted claims
    that no new data caused. Changing the estimator means resetting the baseline
    in the same commit, deliberately, with a changelog entry that says so.
    """
    if before.estimator != after.estimator:
        raise ValueError(
            f"Refusing to diff a {before.estimator!r} baseline against a "
            f"{after.estimator!r} rebuild: the difference would be the method, "
            f"not the evidence. Reset the baseline instead."
        )
    diff = PanelDiff(before=before, after=after)

    old_acc, new_acc = set(before.accessions), set(after.accessions)
    diff.accessions_added = sorted(new_acc - old_acc)
    # Reported with equal prominence: an upstream withdrawal silently shrinks
    # the evidence behind every gene in that dataset.
    diff.accessions_removed = sorted(old_acc - new_acc)

    for key, new in after.estimates.items():
        old = before.estimates.get(key)
        if old is None:
            diff.changes.append({"kind": "new_gene", "key": key, "after": new, "before": None})
            continue
        kind = _classify(old, new)
        if kind:
            diff.changes.append({"kind": kind, "key": key, "before": old, "after": new})

    for key, old in before.estimates.items():
        if key not in after.estimates:
            diff.changes.append({"kind": "dropped_gene", "key": key, "before": old, "after": None})

    order = {kind: i for i, kind in enumerate(CHANGE_ORDER)}
    diff.changes.sort(key=lambda c: (order[c["kind"]], c["key"]))
    return diff


def _format_estimate(row: dict[str, Any] | None) -> str:
    if row is None:
        return "—"
    return f"g {row['g']:+.2f} [{row['ci_low']:+.2f}, {row['ci_high']:+.2f}] k={row['k']}"


def render_markdown(diff: PanelDiff, *, max_rows_per_kind: int = 25) -> str:
    """The changelog, as a pull-request body.

    Truncation is reported rather than silent: a section that says "showing 25
    of 812" tells a reader the corpus moved a lot, which a quietly cut list
    would hide.
    """
    after, before = diff.after, diff.before
    lines = [
        "# Evidence changelog",
        "",
        f"**{before.taken_on} → {after.taken_on}** · manifest "
        f"{before.manifest_version} → {after.manifest_version} · data version "
        f"`{before.data_version[:12]}` → `{after.data_version[:12]}`",
        "",
        f"Panel: {before.n_contrasts} → **{after.n_contrasts} contrasts** over "
        f"{len(after.accessions)} GEO DataSets / {len(after.series)} Series. "
        f"Pooled estimates with k ≥ {MIN_CONTRASTS}: "
        f"{len(before.estimates)} → **{len(after.estimates)}**.",
        "",
    ]

    if diff.is_empty:
        lines += ["Nothing moved. The corpus and every pooled estimate are unchanged.", ""]
        return "\n".join(lines)

    if diff.accessions_added or diff.accessions_removed:
        lines += ["## Panel membership", ""]
        if diff.accessions_added:
            lines.append(
                f"**Added ({len(diff.accessions_added)}):** "
                + ", ".join(f"`{a}`" for a in diff.accessions_added)
            )
        if diff.accessions_removed:
            lines += [
                "",
                f"**Removed ({len(diff.accessions_removed)}):** "
                + ", ".join(f"`{a}`" for a in diff.accessions_removed),
                "",
                "> A dataset leaving the panel shrinks the evidence behind every gene "
                "in it. Confirm this is an upstream withdrawal or reclassification and "
                "not a fetch failure before merging.",
            ]
        lines.append("")

    grouped = diff.by_kind()
    for kind in CHANGE_ORDER:
        rows = grouped[kind]
        if not rows:
            continue
        lines += [f"## {CHANGE_LABELS[kind]} — {len(rows)}", ""]
        lines += ["| gene | species | before | after |", "|---|---|---|---|"]
        for change in rows[:max_rows_per_kind]:
            row = change["after"] or change["before"]
            lines.append(
                f"| `{row['gene_id']}` | {row['species']} | "
                f"{_format_estimate(change['before'])} | "
                f"{_format_estimate(change['after'])} |"
            )
        if len(rows) > max_rows_per_kind:
            lines.append(
                f"\n*Showing {max_rows_per_kind} of {len(rows)}. "
                f"The full diff is in the attached JSON.*"
            )
        lines.append("")

    lines += [
        "---",
        "",
        "**This PR must not be auto-merged.** The panel is rule-selected — datasets "
        "enter because they satisfy the GEO DataSets query and the contrast rules, "
        "never because of their result. Review that the membership change above is "
        "explained by upstream, then merge. Do not hand-edit the panel.",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--snapshot", type=Path, help="take a snapshot and write it here")
    parser.add_argument("--before", type=Path, help="baseline snapshot JSON")
    parser.add_argument("--after", type=Path, help="candidate snapshot JSON")
    parser.add_argument("--markdown", type=Path, help="write the changelog here")
    parser.add_argument("--json", type=Path, help="write the full diff as JSON here")
    args = parser.parse_args(argv)

    if args.snapshot:
        snapshot = take_snapshot()
        snapshot.write(args.snapshot)
        print(f"snapshot: {args.snapshot}")
        print(f"estimates: {len(snapshot.estimates)}")
        print(f"contrasts: {snapshot.n_contrasts}")
        print(f"data_version: {snapshot.data_version}")
        return 0

    if not (args.before and args.after):
        parser.error("provide --snapshot, or both --before and --after")

    diff = diff_snapshots(Snapshot.read(args.before), Snapshot.read(args.after))
    report = render_markdown(diff)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(report, encoding="utf-8")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(diff.to_dict(), indent=1), encoding="utf-8")
    print(report)
    # Exit 1 when something moved, so a workflow can gate PR creation on it
    # without parsing the report.
    return EXIT_CHANGED if not diff.is_empty else EXIT_UNCHANGED


def cli() -> int:
    """`main` with failures mapped to their own exit code.

    "The evidence moved" is exit 1, and so is any uncaught Python exception —
    which meant the workflow's `code -eq 1` test read a crash as a change. It
    then promoted a candidate that had never been diffed, overwriting the
    baseline, and died on the missing changelog. The estimator guard, whose
    whole purpose is refusing to publish a method change as an evidence change,
    was the most likely thing to trip it.

    Errors get their own code so the two can never be confused again. The
    traceback still goes to stderr; nothing is swallowed.
    """
    try:
        return main()
    except Exception:
        traceback.print_exc()
        return EXIT_ERROR


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(cli())
