# Real gene aging signatures — GEO DataSets panel

**Date:** 2026-08-03
**Data:** 31 checksum-pinned GEO DataSets → 32 young-vs-old contrasts → 485,905 gene-level effect sizes
**Code path:** `sources/geo.py` → `harmonize/differential.py` → `harmonize/meta.py`
**Reproduce:** `make signatures` (needs network, ~313 MB download, ~25 min)

---

## 0. What this replaced

Until this build, `signatures.csv` was generated. Twelve genes, hand-assigned
"true" effects (`CDKN2A: +1.20`, `LMNB1: -1.00`), Gaussian scatter around them,
fake `GSE9xxxx` accessions. The directions encoded real biology, which made it
convincing and made it worse: a reader who checked whether p16 went up with age
in GeroQuery got a yes, and the yes meant nothing.

This document is what the same question returns now that the answer comes from
data.

## 1. The panel, and how it was chosen

The panel is **not hand-picked**. It is the output of a rule applied to a query.

1. Ask GEO DataSets for every curated dataset that declares an **age subset
   variable** in human or mouse: 189 records.
2. Keep those that yield an adult young-vs-old contrast under the bands in
   `sources/geo.py` — human young 18–40 y / old ≥60 y, mouse young 2–8 mo / old
   ≥18 mo, with the middle deliberately left empty.
3. Restrict to the control arm of every other subset variable (disease state,
   agent, genotype, protocol, stress, time). A dataset whose confounder has **no
   recognisable control level is dropped**, not guessed at.
4. Split, rather than pool, across tissue and cell-type variables.
5. Require ≥3 samples per group after all restrictions.

**31 datasets survive, producing 32 contrasts from 27 independent GEO Series.**
Change a rule and the panel changes; nothing in it was selected because of what
it showed.

| Tissue | Human contrasts | Mouse contrasts |
|---|---|---|
| Skeletal muscle | 10 | 4 |
| Brain | 4 | 1 |
| Bone marrow / HSC | 1 | 2 |
| Blood | 1 | — |
| Kidney | — | 2 |
| Thymus | — | 2 |
| Adipose | — | 2 |
| Liver | — | 2 |
| Heart | — | 1 |
| Lung | — | 1 |
| **Total** | **16** | **16** |

### Why GEO DataSets and not GEO Series

A GSE series matrix carries age as free text — `age: 67`, `Age (yrs): 67`, `67y`,
or buried in a sample title, differently in every series. A GDS is the
curator-built view: age is a declared subset variable with named levels and
explicit sample lists, and the matrix ships with platform annotation joined, so
every probe carries a gene symbol and an Entrez id. **Both fragile steps of GEO
ingestion — deciding who is old, and deciding what a probe measures — are
answered by GEO's curators rather than by a regex of ours.**

The cost is stated plainly in §5: GEO stopped producing GDS records, so this is a
microarray-era panel with small per-study groups.

## 2. Method

Per contrast:

1. Detect the value scale from the data (99th percentile > 50 ⇒ linear) and
   log2-transform, flooring at 1. GEO's declared `value_type` is not reliable
   enough to decide this — "count" covers both raw MAS5 intensities and
   normalized values.
2. Drop probes that map to more than one gene (`1029///1030`): a probe measuring
   two genes cannot be attributed to either.
3. Collapse probes to genes by **MaxMean** — the highest-mean probe per gene wins.
   Averaging would pull each gene towards its background-level probes.
4. Hedges' *g* (old vs young), with the small-sample correction, and its SE.
5. Benjamini–Hochberg within contrast.
6. Resolve Entrez → Ensembl through mygene.info; unresolvable numeric ids keep an
   `ENTREZ:` id, unresolvable bare symbols are dropped.

Across contrasts: DerSimonian–Laird random effects (`harmonize/meta.py`,
unchanged — this is the first time it has seen a real effect size).

## 3. Result

### 3.1 The panel is honest about its power

| | Human | Mouse |
|---|---|---|
| Genes with ≥3 contrasts | 20,850 | 21,133 |
| Median contrasts per gene | 12 | 10 |
| Median \|pooled *g*\| | 0.143 | 0.206 |
| Median I² | 21.0% | 38.2% |
| **Genes at BH q < 0.05** | **689 (3.3%)** | **44 (0.2%)** |
| — of which up with age | 310 | 13 |
| — of which down with age | 379 | 31 |

Median 95% CI width across all pooled estimates: **0.84**. That is the number to
keep in mind. With groups of 3–15 per study, an effect of *g* = 0.3 — a real,
biologically meaningful shift — sits comfortably inside the confidence interval
of "nothing". The mouse side is weaker still: more tissues, fewer samples each.

**This panel can detect large, consistent, cross-tissue effects. It cannot rule
out moderate ones.** A null here means "not detected by this instrument", never
"absent".

### 3.2 What does come through

The strongest human signals are recognisable aging biology, which is the check
that the pipeline works at all:

| Symbol | pooled *g* | 95% CI | q | k | I² |
|---|---|---|---|---|---|
| DMRT2 | −1.54 | [−2.37, −0.72] | 0.022 | 11 | 84% |
| ARHGAP36 | +1.34 | [+0.61, +2.08] | 0.028 | 10 | 81% |
| CYCS | −1.13 | [−1.56, −0.69] | 3.0e−4 | 16 | 63% |
| **CDKN1A** (p21) | **+1.07** | [+0.55, +1.59] | 0.008 | 12 | 67% |
| MT1F | +1.04 | [+0.59, +1.48] | 0.0017 | 14 | 61% |
| MT2A | +1.02 | [+0.58, +1.46] | 0.0018 | 14 | 61% |
| UQCRFS1 | −1.02 | [−1.41, −0.63] | 2.2e−4 | 14 | 50% |
| SELENOM | +1.01 | [+0.65, +1.36] | 4.1e−5 | 11 | 29% |
| DLAT | −1.01 | [−1.46, −0.55] | 0.0036 | 14 | 63% |

Two coherent programmes fall out without being looked for: **metallothioneins up**
(MT1F, MT2A, MT1HL1) and **mitochondrial / OXPHOS machinery down** (CYCS,
UQCRFS1, DLAT, MTLN). Both are among the most replicated transcriptional
signatures of aging in the literature. Nobody told the pipeline to find them.

### 3.3 The textbook genes — four replicate, six do not

| Symbol | pooled *g* (human) | 95% CI | q | k | Verdict |
|---|---|---|---|---|---|
| CDKN1A (p21) | +1.07 | [+0.55, +1.59] | 0.008 | 12 | **replicates** |
| FOXO3 | +0.82 | [+0.43, +1.21] | 0.007 | 16 | **replicates** |
| SIRT1 | +0.58 | [+0.30, +0.86] | 0.008 | 12 | **replicates** |
| IGF1 | −0.53 | [−0.79, −0.27] | 0.008 | 16 | **replicates** |
| MTOR | −0.27 | [−0.52, −0.02] | 0.358 | 14 | suggestive only |
| GDF15 | +0.27 | [−0.01, +0.55] | 0.457 | 15 | null |
| **CDKN2A (p16)** | **+0.07** | **[−0.20, +0.35]** | **0.925** | **14** | **null** |
| LMNB1 | +0.11 | [−0.23, +0.45] | 0.892 | 14 | null |
| TP53 | +0.04 | [−0.30, +0.38] | 0.970 | 14 | null |
| TERT | −0.04 | [−0.35, +0.27] | 0.969 | 13 | null |
| B2M | +0.09 | [−0.15, +0.33] | 0.872 | 16 | null |

**CDKN2A/p16 — the single most-cited transcriptional marker of aging — does not
replicate.** Fourteen contrasts, six up and eight down, pooled effect
indistinguishable from zero, and low heterogeneity (I² = 14%), so the studies
agree with each other about there being nothing to see.

The synthetic slice this replaced had p16 planted at *g* = +1.20. The old API
test asserted every CDKN2A direction was "up", and passed.

Two things are worth separating here.

- **This is partly an instrument problem, and it is a known one.** p16INK4a is
  expressed at very low levels in bulk tissue, in a small subpopulation of cells.
  The canonical p16-rises-with-age results come from qPCR and single-cell assays,
  not bulk arrays. Several `CDKN2A` probes on these platforms also cannot
  distinguish the p16INK4a and p14ARF transcripts, which are different proteins
  from overlapping reading frames. A bulk microarray panel is close to the worst
  available instrument for this specific transcript.
- **It is not only an instrument problem**, because the same instrument found
  p21, IGF1, the metallothioneins, and the OXPHOS decline with tight intervals.
  The estimator works. p16 specifically is not measurable this way.

LMNB1's human null is similar in character; in mouse it trends in the expected
direction (*g* = −0.42, p = 0.075, k = 14) without reaching significance.

## 4. What this changes about GeroQuery's claims

**Before:** "query a gene, get its multi-omic aging signature" — over generated
numbers.

**Now:** "query a gene, get a random-effects pooled estimate over up to 16 real
published contrasts per species, with a confidence interval and a heterogeneity
statistic wide enough to tell you when the answer is *we cannot say*."

The second is a smaller claim and a true one. It also makes the API's
`meta_signatures` block worth reading: `ci_low`, `ci_high`, `heterogeneity_i2`,
and `n_studies` were decorative when the inputs were manufactured to agree.

## 5. Limitations — read before citing anything above

1. **Microarray-era only.** GEO stopped curating GDS records around 2016. There
   is no RNA-seq in this panel. Probe-level detection limits, cross-hybridization,
   and compressed dynamic range apply to every number in it.
2. **Small groups.** 3–15 samples per group. Median CI width 0.84. Moderate
   effects are invisible.
3. **Tissue coverage is lopsided.** Half the human contrasts are skeletal muscle;
   there is one blood contrast (n = 9) and no proteome or methylome at all. The
   `omic_layer` column reads `transcriptome` for every row.
4. **Cross-study processing heterogeneity is uncorrected.** Twelve platforms,
   three decades of normalization practice. `harmonize/batch.py` exists but is not
   applied across studies — this is exactly the confound recount3 would remove.
5. **Some studies share subjects.** 32 contrasts come from 27 GEO Series: GEO
   splits some experiments across two array halves (GDS287/GDS288,
   GDS355/GDS356, GDS472/GDS473, GDS2961/GDS2962). A random-effects pool assumes
   independent studies. `series_id` is carried on every `Study` row so a caller
   can check this rather than trust the study count — but the pooling **does not
   currently account for it**. In practice the paired halves carry mostly disjoint
   probe sets, so few genes are double-counted; "few" is not "none".
6. **Sex is mostly unrecorded.** Only 8 of 32 contrasts come from a dataset that
   declares a sex variable. The rest are `unspecified` — which means GEO never
   said, not that the contrast mixes sexes.
7. **Pooling across strains and sexes inflates within-group SD**, which biases *g*
   towards zero. Conservative, but it is a bias.
8. **Age bands are thresholds, not biology.** Human 40/60 and mouse 8/18 mo are
   defensible conventions, not measurements. Datasets with only middle-aged
   groups are excluded entirely rather than stretched to fit.
9. **~22% of probe rows were dropped** at identifier resolution (141,156 of
   627,155), almost all retired probes on the oldest platforms with symbols
   mygene.info no longer recognises. Coverage is worse on GPL75/GPL81/GPL91 than
   on GPL570.
10. **No multiple-testing correction across the panel-level pool.** The BH
    q-values in §3 are computed across genes within a species — appropriate for
    "which genes move" but not for "is this specific gene, chosen after looking,
    real".

## 6. What would move this forward

- **recount3 or ARCHS4** for uniformly-processed RNA-seq. This removes limitation
  4 (the biggest one) and 1, and brings blood and modern tissue coverage.
- **A methylation panel** (e.g. GSE40279) would make `omic_layer` mean something
  and would let the 236 wired clocks run on real data — currently the largest
  unexercised capability in the repo.
- **A p16-specific check** by qPCR-derived data, to confirm §3.3's reading that
  the null is instrumental.
- **Nested random effects by series** to handle limitation 5 properly.

---

## Provenance

Every dataset in the panel is declared in `geroquery/sources/manifest.py` with its
exact URL, byte count, and SHA-256, verified on 2026-08-03. Nothing enters the
store unverified. `pytest -m live` re-checks that the URLs still serve those bytes
and that GEO's age-subset query still returns a comparable population.
