# evidence/ — the living evidence loop

This directory is what makes GeroQuery a resource rather than a snapshot.

| File | What it is |
|---|---|
| `baseline.json.gz` | Every pooled estimate with k ≥ 3, plus the panel that produced them, as of the last accepted rebuild. **Overwritten**, not accumulated |
| `CHANGELOG.md` | Append-only record of what moved, newest first. This is the product |

## What the loop does

`.github/workflows/living-evidence.yml` runs monthly. It re-runs the GEO
DataSets selection query, rebuilds the panel from whatever now satisfies it,
pools every gene again, and diffs the result against `baseline.json.gz`.

The output is not "the data was refreshed". It is a changelog of **effect sizes
that moved, with their intervals**:

```
FOXO3   g -0.31 [-0.62, +0.01] k=7   ->  g -0.38 [-0.65, -0.11] k=9   newly excludes zero
CDKN2A  g +0.07 [-0.20, +0.35] k=14  ->  g +0.11 [-0.14, +0.36] k=16  still null
```

Curated aging databases publish *entries*. This publishes *how the evidence
changed*, which is the only honest way to say that a claim strengthened.

## Three rules the loop must not break

**1. Nothing auto-merges.** The workflow opens a pull request and stops. A human
approves every panel change. `test_workflow_never_merges_and_only_ever_opens_a_pull_request`
fails the build if a merge step is ever added.

**2. Datasets leaving are reported as prominently as datasets arriving.** A GDS
withdrawn or reclassified upstream silently shrinks the evidence behind every
gene in it. A diff that only listed additions would hide the corpus shrinking.

**3. The panel stays rule-selected.** Datasets enter because they satisfy

```
age[Subset Variable Type] AND gds[Filter]
  AND (Homo sapiens[Organism] OR Mus musculus[Organism])
```

and the contrast rules in `sources/geo.py` — never because someone liked the
result. `panel_diff` observes; it has no mechanism to include or exclude a
dataset. **Do not hand-edit the panel.** That rule is the only thing keeping
selection bias out.

## Running it by hand

```bash
# Bootstrap or refresh the baseline
python -m geroquery.etl.panel_diff --snapshot evidence/baseline.json.gz

# Compare a candidate rebuild against it
python -m geroquery.etl.panel_diff \
    --before evidence/baseline.json.gz \
    --after  evidence/candidate.json.gz \
    --markdown evidence/changelog-latest.md
```

Exit code is 1 when something moved, 0 when nothing did, so a workflow can gate
on it without parsing the report.

## Reading a changelog entry

`k` is the number of contrasts behind the pooled estimate. A change in `k`
without a change in the interval means new data arrived and agreed. A change in
the interval without a change in `k` means the same data was reprocessed —
which should be rare, and is worth asking about.

The section headings are ordered by newsworthiness: intervals that crossed zero
first, bookkeeping last. A gene appears under exactly one heading.
