# Visualisation and statistics — what to build, and why

**Written:** 2026-08-06 · **Status:** specification. Items marked **[built]** landed
in this session; everything else is scoped for later and is ordered.

This is the design document behind the next layer of GeroQuery's charts and
estimators. It exists because the front end has four views and roughly one chart
per view, and the gap is not that the charts are plain — it is that **each view
answers one question and the data supports several.**

Every recommendation below is checked against the corpus rather than taken from a
methods paper, because several standard recommendations turn out not to apply
here, and one of them is currently top of the roadmap. §2 is the measurement; §7
is the list of things *not* to build and is the most important section.

---

## 1. The bar a new chart has to clear

Three rules, inherited from the repo and not negotiable:

1. **JavaScript renders numbers, it never derives them.** Every statistic is
   computed in Python by the estimator the API serves and shipped in the Parquet.
   A second implementation in the browser is a second estimator, and two
   estimators drift.
2. **A chart must answer a question the table cannot.** Prettier is not a reason.
   If the reader could get the same answer by reading three numbers, ship the
   three numbers.
3. **A new statistic ships with a test that fails when the statistic is wrong** —
   a defining property, an identity one derivative out, or a planted answer whose
   value is known by construction. Not "it runs".

And one added here, from the visualisation review:

4. **Colour is validated, not chosen.** Categorical and diverging palettes go
   through `validate_palette.js` in *both* themes before they ship. §6 records
   what that found in the palette already shipping.

---

## 2. What the corpus actually is, measured

This determines which methods are informative, and it is unusual enough that the
textbook ordering is wrong here. All figures computed from the shipped
`pooled.parquet` / `contrasts.parquet` on 2026-08-06.

| Property | Value | Consequence |
|---|---|---|
| Pooled estimates | 41,983 (20,850 human · 21,133 mouse) | Corpus-scale views are viable |
| Contrasts per gene, *k* | median **11**, IQR 10–14, max 16, min 3 | **Not** a literature meta-analysis |
| Genes with k ≥ 10 | **31,976** (18,320 human · 13,656 mouse) | k-hungry diagnostics apply to 76% of the corpus, not a fringe |
| I² | median **29.2%**, 25% of genes at exactly 0 | Heterogeneity is real but not dominant |
| Median SE across the 32 contrasts | 0.34 → 0.88, a **2.6× spread** | A funnel has genuine precision spread |
| Symbols present in both species | **14,733** | The cross-species view is a 14.7k-point scatter, not an anecdote |
| Human contrasts by tissue | skeletal muscle **10 of 16**, brain 4, blood 1, bone marrow 1 | Tissue is a confounder, not a facet |
| Mouse contrasts by tissue | 9 tissues, 1–3 contrasts each | Mouse can be facetted; human cannot |

### 2.1 The property that changes everything: this panel is near-balanced

In a conventional meta-analysis, *k* varies per outcome because someone searched
the literature and found what was published. Here **every gene is measured in
every contrast that carries a probe for it.** *k* varies from 3 to 16 because
microarray platforms differ in coverage — not because of what any study found.

Two consequences, and they point in opposite directions from the current roadmap:

- **Selective outcome reporting is impossible by construction.** You cannot
  cherry-pick an outcome when the outcome is the entire transcriptome and the
  panel is rule-selected. Classic publication bias is not the threat here.
- **Every gene is pooled from very nearly the same 32 contrasts.** So a
  per-gene influence diagnostic is not just answering "which study carries this
  gene" — aggregated across 41,983 genes it answers **"which of our 32 contrasts
  is bending the whole corpus"**, which is a question about the panel itself and
  which nothing in the repo can currently answer.

That second point is the single highest-value item in this document, and it is
not in the current roadmap at all.

---

## 3. Statistical measures, ranked

### 3.1 Local false sign rate (lfsr) — **the best single addition** **[built]**

**What it answers.** "What is the probability that the *direction* is wrong?"

This is the question the tool is actually asked. A user typing a gene wants to
know whether it goes up or down with age and how much to trust that. A p-value
answers "is the effect exactly zero", which nobody believes and nobody asked.
The lfsr (Stephens 2017) is the posterior probability of getting the *sign*
wrong, and it is defined for every gene — including the 38,276 currently labelled
`no_evidence`, which today get a verdict and nothing else.

**Why it fits this corpus specifically.** Adaptive shrinkage needs a large set of
estimates with standard errors, drawn from a common unimodal-at-zero prior. That
is exactly 41,983 pooled effects per species: most genes genuinely do not change
with age, so the "mostly null" assumption is not an assumption of convenience, it
is the finding. And an empirical-Bayes posterior directly fixes §3.2.

**Cost.** One pass over the pooled table per species. Pure numpy: fit a mixture of
zero-centred normals over a fixed grid of variances by EM, take posterior
means and sign probabilities. No new dependency.

**Test.** Plant a corpus with a known mixture (say 90% null, 10% at g = 0.8) and
assert the recovered null proportion and that lfsr calibrates — among genes with
lfsr ≤ 0.05, at most 5% have the wrong sign. That is a property that fails if the
posterior is wrong and passes only if it is right.

**Pitfall to record.** lfsr is *not* a p-value and must not be labelled one. It is
also a corpus-relative quantity: the same gene gets a different lfsr if the panel
changes, so it belongs behind the estimator version in the evidence baseline like
everything else.

### 3.2 Shrunken effects, and the winner's curse in the table we already ship **[built]**

> **Measured, after building it.** The defect was worse than argued above. Of the
> 800 mouse genes the table can draw from, **719 have lfsr > 0.2**, and the row
> the table led with — `Kif21a`, *g* = +1.154, 95% CI [0.180, 2.127], k = 10 —
> carries a **46% posterior probability that its direction is wrong**. `Rrad`
> (0.445) and `Matr3` (0.435) sat directly beneath it. Human is better and not
> clean: the former top row `ARHGAP36` has lfsr 0.150, while `SELENOM` (0.000)
> and `MT2A` (0.005) ranked below it.
>
> Corpus-wide, **3,707 genes have an interval excluding zero and only 916 have
> lfsr ≤ 0.05.** The two disagree because an interval is a per-gene statement
> and the table is a corpus-wide selection over 41,983 of them; the lfsr carries
> that multiplicity and the interval cannot. Ranking a corpus-wide selection by
> a per-gene interval is the error, and it was shipping.

**The problem, in the current UI.** The Landscape view's "Strongest replicated
effects" table ranks by pooled *g* among genes whose interval excludes zero. That
is ranking on an extreme of a noisy statistic after selecting on significance,
which is the textbook winner's curse: the top of that table is enriched for genes
whose effect was overestimated, and the largest estimated effect is the *most*
likely to shrink on replication. The table's own caption warns that rank order is
not a claim about importance — which is honest, and is not a fix.

**The fix.** Ship the posterior mean from §3.1 alongside the raw pooled *g*, and
rank on it. Shrinkage pulls imprecise extremes toward zero and leaves
well-measured ones nearly untouched, so the ranking becomes "large *and*
well-measured" rather than "large".

**Test.** Simulate from a known prior, rank by raw and by shrunken estimate, and
assert the shrunken ranking has lower mean squared error against the truth. This
is a property with a known answer, not a smoke test.

### 3.3 Per-gene influence diagnostics, and the corpus-level view they unlock

**Per gene:** leave-one-out estimates, Cook's distance, and the two Baujat
coordinates (contribution to Q on x, influence on the pooled estimate on y), from
Viechtbauer & Cheung (2010). For a gene pooled from 14 contrasts, "which one
contrast is carrying this result" is the first question a sceptical reader asks
and the forest plot cannot answer it — a wide row and an influential row look the
same.

**Aggregated over the corpus — the part that is novel.** Because the panel is
near-balanced (§2.1), summing each contrast's influence across all 41,983 genes
gives a per-*contrast* influence score. That answers: *is any one of our 32
contrasts systematically bending the published corpus?* If GDS5226's two
adipocyte contrasts turn out to drive a disproportionate share of mouse verdicts,
that is a finding about the panel, and the panel is the thing this repo asks to be
trusted on.

**Cost.** k extra pooled fits per gene: ~460k fits, which is a vectorised numpy
pass, not a loop. Precompute for the ~3,707 genes carrying a verdict at minimum;
the corpus aggregate needs all of them, so budget the full pass once in the ETL.

**Test.** A planted outlier: pool k−1 concordant contrasts plus one deliberately
discordant one, and assert leave-one-out identifies it, that Cook's distance
ranks it first, and that removing it moves the pooled estimate by the amount the
leave-one-out fit predicted. The last clause is the one that catches an
influence measure that is merely plausible.

### 3.4 Restricted mean survival time on the mortality view

**What it answers.** "How much earlier do people in the top dysregulation
quartile die?" — in **years**, which a hazard ratio cannot say.

The mortality result currently reports HR/SD and Harrell's C. HR is a ratio on a
scale nobody has intuition for, and it is only interpretable at all if
proportional hazards holds over 20 years of follow-up, which §3.6 has never
checked. RMST is model-free: the area under the survival curve to a horizon *t*,
and the difference between two groups is a difference in event-free years. It
stays valid under non-proportional hazards, which is precisely where a single HR
becomes an average of two different truths.

**Fit here.** The design-based machinery already exists — `survival/cox.py` has
weighted estimation and a Taylor-linearised design variance, so a weighted RMST
with a design-based interval is an extension of an estimator that is already
tested, not a new one. Report at t = 10 and t = 15 years.

**Test.** Under exponential survival RMST has a closed form,
`(1 − exp(−λt)) / λ`; assert the estimator recovers it. And assert that with
weights all equal to 1 the weighted estimator equals the unweighted one exactly —
the same defining property the weighted Cox path is already tested on.

### 3.5 Three-level random effects: contrast within series

32 contrasts come from 27 Series, so some contrasts share subjects and the
two-level pool slightly overstates independence. A three-level model
(contrast nested in series) is the correct fix and is the last correctness debt
in the pooling.

**Kept at this rank deliberately.** It is real, and its effect will be small: only
5 of 27 series contribute more than one contrast, so the induced correlation is
confined to a handful of rows. It changes published numbers, which means an
estimator-version bump and a baseline reset, so it should land in a commit of its
own and not be bundled with a visualisation change.

### 3.6 Proportional-hazards diagnostics

Scaled Schoenfeld residuals against time, per covariate, with a smoother and a
±2 SE band; plus the global test for a non-zero slope. Currently listed as a
limitation with nothing behind it, and it is the assumption most likely to fail
over a 20-year follow-up. This is also the gate on §3.4: if PH holds, the HR is
fine and RMST is a convenience; if it does not, RMST becomes the headline and the
HR needs a caveat.

**Caveat to record.** The standard Schoenfeld test is defined for an unweighted
partial likelihood. Under the survey-weighted pseudo-likelihood the reference
distribution is not the usual one — the same reason `likelihood_ratio_test`
already raises on a weighted fit. Report the residual *plot* under weighting and
the *test* on the unweighted fit, and say which is which.

### 3.7 Small-study effects — reinterpreted, and demoted

The roadmap lists Egger's regression and a funnel plot for genes with k ≥ 10 as
priority 3. Two findings move it down and change what it means.

**First, Egger's test is miscalibrated on Hedges' *g*.** The standard error of an
SMD is a function of the effect size itself, so effect and SE are correlated
*by construction* and the regression finds asymmetry where there is none
(Pustejovsky & Rodgers 2019). The entire corpus is Hedges' *g*. A naive Egger
here would produce thousands of false positives and they would look like a
finding. If it is implemented at all it must use a sample-size-based predictor
(`sqrtninv`) or the modified SE, and the docs must say why.

**Second, and more fundamental: publication bias is not the mechanism here**
(§2.1). Nothing was selectively reported. So funnel asymmetry in this corpus
cannot mean what it means in a literature synthesis.

**What it can legitimately be used for.** Reframed as a *design-confounding*
diagnostic: if a gene's effect correlates with contrast precision, and precision
correlates with tissue and platform (it does — the 2.6× SE spread tracks study
size and array generation), then asymmetry is evidence that the pooled estimate
is confounded by study design. That is a real and interesting claim. It is also a
different claim from the one a funnel plot is normally read as making, so the
label on the chart has to carry it.

---

## 4. Visualisations, ranked

Each entry names the form, and the form is chosen from the data's job — magnitude,
identity, polarity, distribution, or agreement — not from what looks impressive.

### 4.1 Volcano over the whole corpus — with the double-filter trap avoided **[built]**

**Job:** the joint distribution of effect and certainty, which the current
effect-size histogram hides. Today the Landscape bins *g* alone, so a reader
cannot see that the largest effects are also the least certain — which is the
corpus's most important structural fact.

**Form:** scatter, effect on x, −log₁₀ of the certainty measure on y, coloured by
verdict using the existing diverging polarity palette. 41,983 points needs
density handling: bin to a hex or square grid and draw counts, with individual
points only in the sparse tails where they are individually meaningful.

**The trap, and it is a real one.** Selecting features by effect size *and*
significance — the two thresholds a volcano plot invites — does not control FDR.
Benjamini-Hochberg bounds the error rate over all features, not over the subset
you carved out with a second filter, and the feature with the largest estimated
effect is a very likely false positive. So: y is the **lfsr** from §3.1, not a
p-value, and the quadrant guides are labelled as *reading* aids with no
thresholding behaviour attached. The chart must not offer a "significant and
large" selector, because that selector is the bug.

### 4.2 Cross-species scatter — the view most likely to produce a finding

**Job:** agreement between two measurements of the same thing.

Human pooled *g* against mouse pooled *g*, one point per ortholog group. Both are
computed today and nothing shows the comparison as a comparison. **14,733
symbols** are present in both species, so this is a real density plot.

**Form:** scatter with a y = x identity line — *not* a regression line as the
primary mark. The question is agreement with the identity, not slope. Where a
fitted line is shown it must be **Deming/orthogonal**, because both axes carry
measurement error and ordinary least squares is biased toward zero slope when the
x variable is noisy; OLS here would systematically understate conservation and
the understatement would look like a result. Quadrant counts (concordant up,
concordant down, discordant) as a small summary beside it.

**Prerequisite:** the exporter must carry `ortholog_group`. Symbol-case matching
is a proxy and will silently mismatch paralogues.

**Companion, not replacement:** a Bland-Altman panel (mean vs difference) answers
"does disagreement grow with effect size", which the scatter cannot show cleanly.

### 4.3 Weighted survival curves by dysregulation quartile

**Job:** the chart every reader of a mortality result expects and this one does
not have. Hazard ratios with no survival curve is a table asking to be trusted.

**Form:** four survey-weighted Kaplan-Meier curves (`svykm`-equivalent: weights in
the risk sets), with a **number-at-risk table beneath aligned to the x ticks** —
non-optional, it is what makes the tail readable — and RMST differences (§3.4)
called out as the summary. Curves ordered by quartile so the legend order matches
the visual order.

**Colour:** this is ordered magnitude, not identity — a sequential single-hue ramp
light→dark, never four categorical hues.

### 4.4 Leave-one-out strip beside the forest plot

**Job:** turn the forest plot from a display into an audit.

**Form:** a narrow companion panel sharing the forest plot's y-axis (one row per
contrast, same order, same labels). Each row shows the pooled estimate **with
that contrast removed**, against a vertical reference line at the full pooled
estimate. A row that jumps is a contrast carrying the result.

Sharing the y-axis is what makes it read in one glance; a separate Baujat scatter
is the same information and requires the reader to join two charts by eye.

### 4.5 Corpus influence — which of our 32 contrasts bends the answer

**Job:** provenance, which is this repo's whole argument.

**Form:** a horizontal bar chart, 32 rows, ranked by aggregate influence
(§3.3), annotated with tissue, platform and *n*. It sits in the Panel view, which
currently lists the contrasts without ever saying what they do.

This is the chart that has no equivalent in any comparable tool, because no
comparable tool re-derives its own numbers.

### 4.6 Tissue small-multiples for one gene

**Job:** decomposition — the pooled number hides that half the human evidence is
one tissue.

**Form:** small multiples, one mini forest per tissue, shared x scale, sorted by
contrast count. **Small multiples, not colour**: mouse has 9 tissues, above the
8-hue ceiling for a categorical scale, and a 9th hue is never generated. Facetting
is the correct answer and it is also the clearer one.

**Honest limit to render, not hide:** for human this will show four panels of
which one holds 10 contrasts and two hold one each. That asymmetry *is* the
limitation section made visible, and the empty-ish panels should be drawn, not
dropped.

### 4.7 Prediction interval on the forest plot — **[shipped, previous session]**

Kept in the list because it is the reference standard for the others: it made
visible, in one dashed rule, that most genes' next study could land anywhere.

**Worth surfacing further:** only **1,280** of 41,983 genes have a *prediction*
interval excluding zero, against 3,707 for the confidence interval. The stronger
and more honest headline number is currently computed and never displayed.

---

## 5. Interaction and accessibility requirements

Applies to every chart above; these are not per-chart decisions.

- **Hover by default.** Every mark gets a tooltip; line/area get a crosshair.
  A rendered SVG chart is interactive and shipping one that is not is a choice
  nobody makes deliberately.
- **Identity is never colour alone.** ≥ 2 series means a legend is present; ≤ 4
  series are also direct-labelled. The existing forest-plot legend already does
  this and is the pattern to copy.
- **A table view exists for every chart.** The Landscape already pairs its
  histogram with a table; that pairing is the standard, not an accident of that
  view.
- **Text wears text tokens.** Values and labels stay in ink/ink-soft/ink-faint;
  the coloured mark beside them carries the identity. No orange numerals.
- **Dark mode is selected, not flipped** — see §6, where it currently is not.
- **Wide charts scroll inside their own container.** The 393px Playwright project
  exists to catch the page scrolling sideways instead; new charts join it.

---

## 6. The palette audit, and a defect it found

The shipped palette is a **diverging** scale — `--up` warm, `--down` cool, `--null`
neutral gray — which is the correct encoding for polarity and is already right in
structure. Running `validate_palette.js` against it in both themes:

| Theme | Pair | Normal-vision ΔE | CVD worst | Verdict |
|---|---|---|---|---|
| Light `#e9e7e2` | `--null` ↔ `--down` | **20.7** | 19.0 protan | pass |
| Light | `--up` ↔ `--down` | 26.1 | 19.5 protan | pass |
| **Dark `#101116`** | **`--null` ↔ `--down`** | **12.8** | 12.7 deutan | **FAIL — below the 15 floor** |

**In dark mode, "interval crosses zero" and "interval below zero" are hard to tell
apart even with full colour vision.** Those are two of the three verdicts, they
are the fill colours on the forest plot and the landscape histogram, and the
distinction between them is the product's entire thesis. Light mode is fine; dark
mode lightened the neutral to `#9a968f` and that is what closed the gap.

**Fix, found by search rather than by eye:** `--null: #7d7568` in dark mode —
normal-vision ΔE **19.1**, CVD 11.8 protan / 18.9 tritan, contrast ≥ 3:1 against
the dark surface. It is essentially the light theme's warm gray, which is the
tell: the neutral should not have been lightened with the rest of the ramp.

The two remaining validator FAILs on this palette (lightness band, chroma floor)
are the categorical checks applied to a diverging scale — a low-chroma neutral
midpoint is *required* for diverging and correctly flagged as gray. The validator
says so itself: *scope: categorical palettes only.* Recorded here so the next
person does not "fix" the midpoint by saturating it.

**New palettes still needed:** none for the charts above. §4.3 needs a sequential
ramp (single hue, light→dark) and §4.6 uses facetting rather than hues
specifically to avoid needing a 9-colour categorical scale.

---

## 7. What not to build, and why

The most useful part of this document.

| Not this | Why |
|---|---|
| **Egger's test as shipped in every meta-analysis package** | Miscalibrated on Hedges' *g* by construction (§3.7). It would produce thousands of confident false positives across 31,976 eligible genes, and they would be indistinguishable from findings. |
| **A funnel plot labelled "publication bias"** | The mechanism does not exist in this corpus (§2.1). The plot may be worth drawing; that label would be false. |
| **Trim-and-fill** | Imputes studies that were never suppressed, for the same reason. It would invent data, in a repo whose first guarantee is that nothing is asserted that was not measured. |
| **Fail-safe N** | Deprecated in the methods literature and answers a question nobody asked. |
| **A "significant and large" selector on the volcano** | Double filtering breaks the FDR guarantee (§4.1). The convenience feature *is* the bug. |
| **A dual-axis chart anywhere** | Two y-scales invite any correlation the author wants. Two measures → two charts or a shared index. |
| **A 9-hue categorical tissue scale** | Above the 8-hue ceiling; small multiples are both the rule and the clearer answer (§4.6). |
| **Any statistic computed in JavaScript** | Standing repo rule. Applies to the new ones: lfsr, shrinkage, leave-one-out and RMST are all Python, shipped in the Parquet. |
| **Ranking by raw pooled *g*** | Winner's curse (§3.2). Keep the column, rank on the shrunken one. |

---

## 8. Two corrections to published numbers, found while measuring §2

Both are the repo's characteristic failure mode — nothing errors, and the wrong
number looks exactly like the right one.

1. **`HANDOFF_2026-08-06.md` reports 3,711 genes excluding zero; the shipped
   artifact has 3,707.** The published figure was computed from the p-value at
   full precision, and the previous "4,875 under DerSimonian-Laird" is likewise
   the count of p < 0.05 — the DL *interval* excludes zero for 4,871. The
   difference is four genes whose bound rounds onto zero at the fourth decimal.
   This is bug #20's family again: **a tool whose thesis is the interval reported
   its own headline from a p-value.** The counts in §2 are computed from the
   intervals as they ship.

2. **Bug #20 was fixed on four of six interval columns.** The `+ 0.0` that
   collapses negative zero was applied to `ci_low`, `ci_high`, `ci_low_dl` and
   `ci_high_dl`, and not to `pi_low`/`pi_high`, which were added in the same
   commit. Two genes — human `APPL1` and mouse `Ckap2l` — shipped
   `pi_high = -0.0`, printing as "−0.000" and comparing as *not less than zero*
   in JavaScript. Fixed this session by routing every bound through one `_bound`
   helper and asserting over a `BOUND_COLUMNS` tuple rather than over two
   hand-named columns — *the original test named its columns literally, which is
   why it could only ever catch the bug already known.*

---

## 9. Build order

**Phase 1 is done.** `harmonize/shrinkage.py` with 13 tests, the `Certainty`
view, the re-ranked landscape table carrying an `lfsr` column, the dark-mode
neutral, and the negative-zero fix. Gate: **436 Python tests** (was 423) ·
**34 Playwright** (was 26) · ruff, black, mypy clean.

| Phase | Items | Why here |
|---|---|---|
| **1 — done** | §8.2 negative-zero fix · §6 dark-mode neutral · §4.1 volcano on lfsr · §3.1 lfsr + §3.2 shrunken effects | The two fixes are defects. lfsr is the enabling statistic for the volcano, and both change what the Landscape view claims. |
| **2** | §3.3 influence + §4.4 leave-one-out strip + §4.5 corpus influence | One estimator, two charts, and the only genuinely novel view in the document. |
| **3** | §4.2 cross-species scatter (needs `ortholog_group` in the exporter) | Highest chance of an actual finding; blocked on one export column. |
| **4** | §3.6 PH diagnostics → §3.4 RMST → §4.3 survival curves | Ordered: the diagnostic decides whether RMST is a convenience or the headline. |
| **5** | §3.5 three-level model | Correctness debt. Own commit, estimator bump, baseline reset — never bundled. |
| **6** | §4.6 tissue small-multiples | Cheap, and honest about the muscle-heavy panel. |

---

## 10. Sources

- Kossmeier, Tran & Voracek, *Charting the landscape of graphical displays for
  meta-analysis and systematic reviews: a comprehensive review, taxonomy, and
  feature analysis*, BMC Med Res Methodol 2020 — the 11-class taxonomy, and the
  finding that diagnostic plots are the most underused class —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7006175/
- Viechtbauer & Cheung, *Outlier and influence diagnostics for meta-analysis*,
  Res Synth Methods 2010 — leave-one-out, Cook's distance, Baujat coordinates —
  https://onlinelibrary.wiley.com/doi/10.1002/jrsm.11
- Stephens, *False discovery rates: a new deal*, Biostatistics 2017 — adaptive
  shrinkage and the local false sign rate —
  https://academic.oup.com/biostatistics/article/21/1/15/5050477
- Pustejovsky & Rodgers, *Testing for funnel plot asymmetry of standardized mean
  differences*, Res Synth Methods 2019 — why Egger's test is miscalibrated on
  SMDs — https://onlinelibrary.wiley.com/doi/abs/10.1002/jrsm.1332
- metafor `regtest` — sample-size-based predictors (`sqrtninv`) as the remedy —
  https://wviechtb.github.io/metafor/reference/regtest.html
- *Inflated false discovery rate due to volcano plots: problem and solutions*,
  Brief Bioinform 2021 — why double filtering breaks the BH guarantee —
  https://academic.oup.com/bib/article/22/5/bbab053/6184412
- Royston & Parmar, *Restricted mean survival time: an alternative to the hazard
  ratio*, BMC Med Res Methodol 2013 —
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3922847/
- Cole & Hernán, *Adjusted survival curves with inverse probability weights*,
  Comput Methods Programs Biomed 2004 — weighted KM as direct standardisation —
  https://pubmed.ncbi.nlm.nih.gov/15158046/
- Van den Noortgate et al., *Three-level meta-analysis of dependent effect
  sizes* — contrast-within-series nesting —
  https://www.um.es/metaanalysis/pdf/5044.pdf
- Lex, Gehlenborg et al., *UpSet: Visualization of Intersecting Sets*, IEEE
  TVCG 2014 — the 4–30 set band, for any later panel-overlap view —
  https://vdl.sci.utah.edu/publications/2014_infovis_upset/
- Deming regression / errors-in-both-variables, for §4.2 —
  https://en.wikipedia.org/wiki/Deming_regression
- The winner's curse and shrinkage —
  https://associationofanaesthetists-publications.onlinelibrary.wiley.com/doi/10.1111/anae.16161
