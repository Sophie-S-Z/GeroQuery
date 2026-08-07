# Evidence changelog

Newest first. Each entry is one accepted rebuild of the GEO aging panel, and
records which pooled estimates moved — not that data was refreshed.

See [`README.md`](README.md) for how the loop works and the three rules it
must not break.

---

## 2026-08-06 — baseline reset: the estimator changed, not the evidence

**No new data.** The panel is byte-identical to 2026-08-05 — same 32 contrasts,
same manifest, same data version. What changed is the interval.

The reported confidence interval is now **Hartung-Knapp-Sidik-Jonkman**, widened
to the DerSimonian-Laird interval wherever HK comes out narrower. DL substitutes
an *estimated* between-study variance into the weights and then proceeds as if it
were known, which makes the interval too narrow — and a quarter of this corpus
pools fewer than ten contrasts, where that bias is worst.

| | DerSimonian-Laird | Reported now |
|---|---|---|
| Genes with an interval excluding zero | 4,875 | **3,711** |
| Median interval width | 0.843 | 0.928 |
| Median 95% prediction interval width | — | 1.831 |

**1,164 claims retracted. None created.** A wider interval can only ever take a
verdict away, and `test_widening_the_interval_can_only_retract_a_verdict_never_create_one`
asserts that over the corpus rather than trusting it.

This is the one kind of change the loop must *not* publish as a diff, because
"1,164 genes stopped excluding zero" would read as evidence moving when no
evidence moved. Snapshots now record the estimator that produced them and
`diff_snapshots` refuses to compare across two. The baseline is reset here
deliberately; the next scheduled run compares against it.

Also new in every snapshot's source estimate: a **95% prediction interval**. The
confidence interval says where the *mean* effect is; the prediction interval says
where the *next study's* effect would fall. Under this corpus's heterogeneity the
second is about twice as wide, and conflating them is the most common misreading
of a forest plot.

---

## 2026-08-05 — baseline established

First snapshot. Nothing to compare against yet.

| | |
|---|---|
| Panel | 32 contrasts over 32 GEO DataSets / 27 Series |
| Pooled estimates with k ≥ 3 | **41,983** |
| Manifest | 2026.3, checksums verified 2026-08-05 |
| Data version | `2024.1+6650ee535bcb` |

The next scheduled run compares against this. From then on, this file records
what changed.
