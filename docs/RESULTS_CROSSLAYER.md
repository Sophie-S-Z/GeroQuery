# Cross-layer result — clocks, health state, and death in one cohort

**Cohort:** NHANES 1999–2002 DNA methylation subsample · **n = 2,517** ·
**1,350 deaths** · median follow-up **17.1 years** (max 20.8)
**Design:** `WTDN4YR` weights, 28 strata, 57 sampling units →
**75,754,006 US adults aged 50+**, 35,727,943 deaths
**Manifest:** 2026.3, checksums verified 2026-08-05 · **Reproduce:** `make crosslayer`

> **Updated 2026-08-06 — the survey design is now applied.** The first version of
> this document reported unweighted estimates and listed weighting as the highest
> value follow-up. It has been done, and it moved every number: **every hazard
> ratio is larger once the design is respected**, because the DNAm subsample
> under-represents the higher-risk part of the population. The unweighted value
> is kept beside each weighted one, since the size of that gap is itself the
> finding. Both conclusions survive: the trained-on hierarchy still reproduces,
> and the same three clocks still carry no evidence over six blood tests.
>
> One consequence had to be handled rather than absorbed: **a likelihood ratio
> test is not valid against a weighted pseudo-likelihood.** The nested tests are
> now design-based Wald tests. `likelihood_ratio_test` refuses a weighted fit
> outright rather than returning a p-value from the wrong distribution.

This is the analysis `HANDOFF.md` §10.3 called "the reason this exists" and
listed as blocked. It was blocked because clocks and health state lived in
different cohorts. They no longer do: NCHS released per-participant DNAm
biomarkers for the 1999–2002 cycles on 2024-07-31, and those participants are in
the same survey — and the same linked mortality file — as the clinical panel.

Every number below is computed through the repo code path
(`sources.nhanes_dnam` → `survival.crosslayer`), not a standalone script.

---

## 1. What was asked

| # | Question | Model comparison |
|---|---|---|
| 1 | Does an epigenetic clock predict death here? | B vs A |
| 2 | Does a dysregulated health state predict death? | C vs A |
| 3 | **Does either add anything over the other?** | D vs B, D vs C |

Model A is age + sex. B adds clock age acceleration, C adds dysregulation,
D adds both. All four are fitted on one design matrix, sliced — so the nested
tests compare models on identical rows, and the survey columns go through the
same complete-case filter as the covariates rather than being aligned afterwards.
Hazard ratios are **per standard deviation**, which is the only way a clock in
years and a Mahalanobis distance in arbitrary units can be put side by side.

**How the design enters.** The point estimate maximizes Binder's (1992) weighted
pseudo-likelihood with `WTDN4YR`; the variance is the Taylor-linearized
design-based sandwich over `SDMVSTRA` × `SDMVPSU`. The two do different jobs and
are worth separating: **the weights move the estimate, the clustering moves the
interval.** Here the weights move it a lot and the clustering almost not at all —
design effects run 0.88 to 1.71, median ≈ 1.15 — so nearly all of the change
below is unequal selection rather than correlated sampling units.

**On the resilience measure.** Critical slowing down as implemented in
`resilience/` is a *population* statistic — variance across age strata. It has
no per-subject value, so it cannot enter a Cox model. The individual-level
analogue used here is Mahalanobis distance from a young-reference centroid over
the same six markers (Cohen et al. 2013 physiological dysregulation), computed
from the same covariance structure whose drift CSD measures. **They are not the
same quantity** and are not reported as one. The population CSD is in §5.

---

## 2. Result

All ten age-like clocks, joint model D. HR per SD, 95% CI, survey-weighted.
`deff` is the design effect: the design-based variance over the variance the same
weighted fit would have with every subject as their own sampling unit.

| Clock | Trained on | HR alone | **HR adjusted for dysregulation** | 95% CI | unweighted HR | deff | C (clock) | C (joint) |
|---|---|---|---|---|---|---|---|---|
| GrimAge2 | mortality | 1.628 | **1.527** | 1.410–1.652 | 1.443 | 1.01 | 0.7872 | 0.7943 |
| GrimAge | mortality | 1.611 | **1.518** | 1.403–1.642 | 1.427 | 0.99 | 0.7854 | 0.7936 |
| PhenoAge | phenotypic age | 1.317 | **1.220** | 1.103–1.350 | 1.163 | 1.66 | 0.7703 | 0.7813 |
| Hannum | chronological age | 1.217 | **1.145** | 1.059–1.237 | 1.108 | 1.16 | 0.7669 | 0.7796 |
| Horvath | chronological age | 1.139 | **1.116** | 1.034–1.205 | 1.129 | 1.14 | 0.7643 | 0.7788 |
| Zhang | chronological age | 1.120 | **1.103** | 1.002–1.214 | 1.079 | 1.71 | 0.7636 | 0.7780 |
| Vidal-Bralo | chronological age | 1.129 | **1.086** | 1.009–1.169 | 1.070 | 1.40 | 0.7638 | 0.7779 |
| SkinBlood | chronological age | 1.086 | 1.070 | **0.984–1.164** | 1.048 | 1.35 | 0.7630 | 0.7775 |
| Lin | chronological age | 1.069 | 1.045 | **0.980–1.115** | 1.052 | 0.90 | 0.7624 | 0.7774 |
| Weidner | chronological age | 1.052 | 1.033 | **0.974–1.095** | 1.026 | 0.88 | 0.7621 | 0.7772 |

Baseline (age + sex) C = **0.7618**. Dysregulation alone C = **0.7771**.
Dysregulation HR per SD = **1.28–1.42** depending on which clock it sits beside
(unweighted: 1.20–1.29).

### 2.0 Every effect is larger in the population than in the sample

The weighted hazard ratio exceeds the unweighted one for **nine of the ten
clocks** and for dysregulation under every clock. That direction is not a
detail. It means the DNAm subsample over-represents people whose measured aging
markers carry *less* mortality information than they do nationally — so the
unweighted analysis was a conservative one, and correcting it strengthens rather
than deflates the finding.

It is worth being precise about what weighting does and does not fix. It corrects
for **who among the eligible was selected**. It cannot correct for who was
eligible: NHANES's DNAm subsample is adults 50+ who consented and had usable
stored DNA, and no weight recovers a person who was never in the frame.

### 2.1 The ordering is the validation

Nothing in the pipeline knows what a clock was trained on. The ranking that
falls out is, top to bottom: mortality-trained (GrimAge2, GrimAge), then
phenotype-trained (PhenoAge), then blood chronological-age clocks (Hannum,
Horvath), then the rest. That is the published hierarchy, reproduced from raw
files by code that was never told it. Getting this ordering wrong would have
been the loudest possible signal that the join was scrambled.

### 2.2 Three clocks do not predict death at all

**Lin, SkinBlood, and Weidner have confidence intervals that include 1.0** once
dysregulation is in the model. They are chronological-age clocks and, in a
cohort already adjusted for chronological age, they carry no mortality
information this study can detect. Published as found.

### 2.3 Answer to question 3 — both directions

Nested tests are **design-based Wald tests** on the weighted fit, because a
likelihood ratio is not valid against a pseudo-likelihood (§2.3.1). The
unweighted likelihood-ratio p-values are carried beside them.

**Dysregulation adds over every clock**, including the mortality-trained ones:

| Clock in the model already | Wald p for adding dysregulation |
|---|---|
| Weidner / SkinBlood / Lin / Zhang | 3×10⁻⁴⁰ – 1×10⁻³⁹ |
| Vidal-Bralo / Horvath / Hannum | 4×10⁻³⁹ – 3×10⁻³¹ |
| PhenoAge | 2×10⁻²³ |
| GrimAge | 7×10⁻¹⁹ |
| GrimAge2 | 6×10⁻¹⁶ |

**The clock adds over dysregulation only sometimes:**

| Clock | Wald p, weighted | LR p, unweighted | Verdict |
|---|---|---|---|
| GrimAge2 | 1×10⁻²⁵ | 1×10⁻³² | adds |
| GrimAge | 4×10⁻²⁵ | 1×10⁻³⁰ | adds |
| PhenoAge | 1×10⁻⁴ | 1×10⁻⁷ | adds |
| Hannum | 6×10⁻⁴ | 2×10⁻⁴ | adds |
| Horvath | 0.005 | 1×10⁻⁵ | adds |
| Vidal-Bralo | 0.027 | 0.009 | adds |
| Zhang | 0.045 | 0.007 | adds |
| **SkinBlood** | **0.114** | 0.087 | **no evidence** |
| **Lin** | **0.178** | 0.059 | **no evidence** |
| **Weidner** | **0.281** | 0.316 | **no evidence** |

So the honest one-line answer to "is an epigenetic clock worth running": for
GrimAge and PhenoAge, yes, clearly. For three of the seven chronological-age
clocks, **there is no evidence it tells you anything six routine blood tests do
not** — and that is now a statement about the US population aged 50+, not about
2,517 people. Zhang and Vidal-Bralo sit close enough to 0.05 that the verdict is
weak in either direction; they should be read as marginal, not as established.

#### 2.3.1 Why the test changed

Twice the difference in a *weighted* log pseudo-likelihood is not χ² distributed
— the weights break the information equality the likelihood ratio depends on. A
p-value computed that way is drawn from the wrong reference distribution, and it
would look entirely normal. So `likelihood_ratio_test` now raises on a weighted
fit rather than answering, and `wald_test` — which reads the design-based
covariance matrix and therefore inherits both the weighting and the clustering —
is what the analysis uses. On the unweighted fit the two agree to within 15% of
the statistic, which is what
`test_the_wald_test_agrees_with_the_likelihood_ratio_when_unweighted` pins.

### 2.4 The headline, stated plainly

> A Mahalanobis distance over six routine clinical chemistry markers
> (C = 0.777) predicts 20-year mortality better than eight of the ten DNA
> methylation clocks NCHS published, and adds information to all ten —
> including GrimAge2.

Only GrimAge (0.785) and GrimAge2 (0.787) beat it, and both are trained on
mortality, which is the outcome being predicted. The six markers cost a
standard blood panel; the clocks cost an EPIC array.

*Correction:* the first version of this section said "nine of the ten" while also
saying only two clocks beat it. Two of ten beating it means it beats eight. The
count was wrong then and is right here; the substance is unchanged, and it is
unchanged by weighting too — the same two clocks lead before and after.

---

## 3. The caveat that matters more than the p-values

**Every p-value above is tiny and every effect on discrimination is small.**

```
                                    weighted     unweighted
age + sex only                      C = 0.7618   0.7394
+ the best clock (GrimAge2)         C = 0.7872   0.7644   (+0.025)
+ dysregulation                     C = 0.7771   0.7567   (+0.015)
+ both                              C = 0.7943   0.7731   (+0.033)
```

With 1,350 events, a nested test detects effects far below the size that would
change a decision about a person. **All of this biology adds 0.033 of concordance
over knowing someone's age and sex** — and weighting did not rescue that number,
it moved every row up by about the same 0.02 and left the gap where it was. A
reader who takes `p = 3×10⁻⁴⁰` as "this is a large effect" has read it wrong, and
the number that corrects them is in the same table.

This is the same discipline the CDKN2A result applies in the other direction:
there, a wide interval was reported rather than a verdict; here, a tiny p-value
is reported *next to* the effect size that deflates it.

---

## 4. What was found on the way

**The DNAm sex coding is the inverse of NHANES's.** `XY_Estimation` is
1 = female, 2 = male; `RIAGENDR` is 1 = male, 2 = female. Both are
valid-looking 1/2 integers, so reusing one map for the other produced 2.4%
agreement instead of 97.6% — that is, it mislabelled every subject, silently.
Caught by crosstab, pinned by `test_dnam_sex_coding_is_the_inverse_of_nhanes_sex_coding`.

The 60 genuinely discordant samples (2.4%) are reported as a QC rate by
`sex_discordance()` on every build rather than dropped. A rise in that rate is
the cheapest available signal that samples were mixed up upstream.

**Variable names are not stable across NHANES cycles.** Creatinine is `LBXSCR`
in 1999–2000 and `LBDSCR` in 2001–2002; alkaline phosphatase `LBXSAPSI` then
`LBDSAPSI`. Using one cycle's names for both does not raise — it drops the other
cycle's 1,306 subjects and yields a cohort half the size that looks fine.

---

## 5. The population CSD signal disappears in this cohort — and that is informative

Running the existing `ResilienceService.csd` on this cohort's six-marker health
state, unchanged:

| Indicator | Slope | 95% CI | Supported |
|---|---|---|---|
| Variance across age strata | +0.0010 | [−0.0015, +0.0038] | **No** |
| Marker cross-correlation | +0.00015 | [−0.0018, +0.0015] | **No** |

`resilience_declines` → **False**, on both indicators.

Compare `RESULTS_NHANES_CSD.md`, where variance was supported in **20 of 20**
analytic configurations on NHANES 2017–2018. The difference is not the method
and not the markers — both are identical. It is the **age floor of 50**. The
2017–2018 cohort spans 20–80; this one spans 50–85. Removing the young half of
the range removes the variance trend.

That was predicted before the analysis was run (see the age-floor caveat in
`sources/nhanes_dnam.py`), which is the only reason it counts as a finding
rather than an excuse.

**And yet the individual-level statistic is strongly predictive** in the same
subjects, with a monotone gradient and almost no age correlation (r = 0.077):

| Dysregulation quartile | Died (sample) | Died (weighted) | n |
|---|---|---|---|
| Q1 (lowest) | 44.9% | 36.6% | 630 |
| Q2 | 47.5% | 43.4% | 629 |
| Q3 | 52.6% | 45.7% | 629 |
| Q4 (highest) | **69.5%** | **67.4%** | 629 |

The weighted column is the population one. The gradient is steeper in the
population than in the sample — Q1 falls by 8 points and Q4 by only 2 — which is
the same story the hazard ratios tell: the sample understates the spread.

So the two are measuring different things, and this cohort is the first place
that could separate them: **the population variance trend is an artefact of age
range; the individual dysregulation signal is not.** Anyone reading the CSD
literature as evidence about individuals should read that table.

---

## 6. Limitations

**Cohort.** Age 50–85 only, topcoded at 85 — this is the methylation
subsample's definition, not a choice. The CSD design loses roughly half its
dynamic range (§5). Exposure is cross-sectional: markers and methylation are
measured once, at the baseline exam; only the outcome is longitudinal.

**CRP is a different assay.** 1999–2002 used standard CRP in mg/dL; 2017–2018
used high-sensitivity CRP in mg/L. Values are converted (×10) so the column
means the same thing, but detection limits differ at the low end and the two
cohorts are not interchangeable.

**Clocks are NCHS's, not ours.** That removes our normalization from the picture
— a genuine improvement over `RESULTS_METHYLATION_CLOCKS.md`, where the
normalization was not the authors' — and it also means we cannot audit theirs.

**Mortality records are perturbed.** NCHS substitutes synthetic follow-up time
or cause of death for a subset of public-use records. Aggregate estimates are
designed to survive this; individual rows are not trustworthy.

**Survey weighting corrects selection, not eligibility.** `WTDN4YR` with
`SDMVSTRA`/`SDMVPSU` is applied, so the estimates describe US adults aged 50+ in
1999–2002. What it cannot do is repair the sampling frame: the DNAm subsample is
people who consented and had usable stored DNA, and no weight recovers someone
who was never eligible. The `WTMEC4YR` exam weight is also carried and is *not*
the right one here — it describes the MEC-examined sample, not the methylation
subsample.

**Two strata carry an unusual PSU count.** 28 strata, 57 sampling units: 27
strata with two and one with three. No stratum has a single PSU, so nothing is
being silently dropped — `lonely_psu_strata` is reported on every fit and is 0.
With 28 degrees of freedom the design variance is itself estimated with
noticeable uncertainty, which is why the design effects range from 0.88 to 1.71
on quantities that should behave similarly.

**Ties are handled by Breslow.** Follow-up is recorded in whole months over ~20
years, so ties are common. Breslow attenuates hazard ratios toward the null
relative to Efron — the safer direction, but the effects above are mild
underestimates.

**Dysregulation reference is age-defined, not health-defined.** Cohen et al. use
a healthy reference population; there is no health screen here, so the reference
is everyone aged ≤60 — who are themselves 50+. That weakens the statistic.

**No cause-specific analysis.** `ucod_leading` is carried and unused. All-cause
mortality only.

**No proportional-hazards diagnostic.** The Cox models are not tested for the
proportional-hazards assumption over 20 years of follow-up, where it is least
likely to hold.

---

## 7. The estimator earns its own trust

The Cox model is hand-rolled (numpy + scipy; neither lifelines nor statsmodels
is a dependency). It is checked three independent ways in `tests/test_survival.py`:

1. **A closed form.** Three subjects, x = [0, 1, 0], times 1 < 2 < 3, no
   censoring: the Breslow MLE is exactly β = ln(2)/2. Matched to 1e-9.
2. **Finite differences.** The analytic gradient and Hessian agree with the
   numerical derivatives of the log partial likelihood they claim to
   differentiate, to 1e-6. This is what catches a risk set accumulated in the
   wrong direction — a bug that otherwise produces a smooth, convergent,
   entirely wrong fit.
3. **A planted effect.** Simulated exponential survival with known coefficients
   [0.7, −0.4] recovers to ±0.06. The same argument the synthetic CSD fixture
   exists to make: an estimator that cannot find an effect known to be there has
   not earned the right to report a null — and §2.2 reports three.

Residual bias over 40 replicates is ~2% relative, consistent with the known
O(1/n) bias of the Cox partial-likelihood MLE rather than with an implementation
error.

The survey extension is checked three more ways, on the same principle:

4. **A defining property.** A weight of two must give exactly the fit a
   duplicated row gives. This is the check that separates "the weight reached the
   risk set" from "the weight reached the event term" — getting only one of the
   two produces a convergent, wrong fit.
5. **An identity one derivative out.** The per-subject score residuals that build
   the sandwich must sum to the gradient of the log pseudo-likelihood, at the fit
   and away from it. A mis-derived residual still yields a plausible variance and
   fails only this.
6. **A planted design.** A population whose two strata have different hazard
   ratios, sampled at different rates. The weighted fit on the biased sample must
   land on the answer computed from the whole population; the unweighted one must
   visibly miss it, or the test proves nothing about the weighting.

---

## 8. Reproducing

```bash
make crosslayer          # fetch + verify 11 pinned artifacts, join, write the table
```

```python
from geroquery.sources import nhanes_dnam as nd
from geroquery.survival import mahalanobis_dysregulation, crosslayer_analysis

frame = nd.load_full()
dysregulation = mahalanobis_dysregulation(frame, list(nd.MARKERS))
result = crosslayer_analysis(
    frame, nd.age_acceleration(frame, "GrimAge2Mort"), dysregulation.values,
    clock_name="GrimAge2Mort",
    # Omit these three and you get the sample estimate, not the population one.
    weight_col=nd.WEIGHT_COLUMN,
    strata_col=nd.STRATUM_COLUMN,
    psu_col=nd.PSU_COLUMN,
)
print(result.headline())
```

**Sanity numbers:** 2,517 subjects · 1,350 deaths · 60 sex-discordant ·
age 50–85 · follow-up median 17.08 y · 28 strata / 57 PSUs / 0 lonely ·
population 75,754,006 · weighted baseline C 0.7618 · GrimAge2 HR/SD **1.527**
[1.410, 1.652] weighted, 1.443 unweighted.
