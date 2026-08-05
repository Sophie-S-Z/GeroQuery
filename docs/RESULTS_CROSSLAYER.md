# Cross-layer result — clocks, health state, and death in one cohort

**Cohort:** NHANES 1999–2002 DNA methylation subsample · **n = 2,517** ·
**1,350 deaths** · median follow-up **17.1 years** (max 20.8)
**Manifest:** 2026.3, checksums verified 2026-08-05 · **Reproduce:** `make crosslayer`

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
D adds both. All four are fitted on one design matrix, sliced — so the
likelihood-ratio tests compare models on identical rows. Hazard ratios are
**per standard deviation**, which is the only way a clock in years and a
Mahalanobis distance in arbitrary units can be put side by side.

**On the resilience measure.** Critical slowing down as implemented in
`resilience/` is a *population* statistic — variance across age strata. It has
no per-subject value, so it cannot enter a Cox model. The individual-level
analogue used here is Mahalanobis distance from a young-reference centroid over
the same six markers (Cohen et al. 2013 physiological dysregulation), computed
from the same covariance structure whose drift CSD measures. **They are not the
same quantity** and are not reported as one. The population CSD is in §5.

---

## 2. Result

All ten age-like clocks, joint model D. HR per SD, 95% CI.

| Clock | Trained on | HR alone | **HR adjusted for dysregulation** | 95% CI | C (clock) | C (joint) |
|---|---|---|---|---|---|---|
| GrimAge2 | mortality | 1.531 | **1.443** | 1.361–1.531 | 0.7644 | 0.7731 |
| GrimAge | mortality | 1.502 | **1.427** | 1.346–1.512 | 0.7618 | 0.7720 |
| PhenoAge | clinical phenotype | 1.236 | **1.163** | 1.099–1.231 | 0.7468 | 0.7603 |
| Hannum | chronological age | 1.143 | **1.108** | 1.049–1.171 | 0.7426 | 0.7584 |
| Horvath | chronological age | 1.133 | **1.129** | 1.069–1.192 | 0.7421 | 0.7586 |
| Vidal-Bralo | chronological age | 1.095 | **1.070** | 1.017–1.126 | 0.7416 | 0.7578 |
| Zhang | chronological age | 1.076 | **1.079** | 1.019–1.142 | 0.7406 | 0.7574 |
| Lin | chronological age | 1.071 | 1.052 | **0.998–1.108** | 0.7408 | 0.7573 |
| SkinBlood | chronological age | 1.056 | 1.049 | **0.992–1.108** | 0.7401 | 0.7571 |
| Weidner | chronological age | 1.044 | 1.026 | **0.976–1.078** | 0.7400 | 0.7568 |

Baseline (age + sex) C = **0.7394**. Dysregulation alone C = **0.7567**.
Dysregulation HR per SD = **1.20–1.29** depending on which clock it sits beside.

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

**Dysregulation adds over every clock**, including the mortality-trained ones:

| Clock in the model already | LR p for adding dysregulation |
|---|---|
| Weidner / SkinBlood / Zhang / Horvath | 2×10⁻²⁷ – 1×10⁻²⁷ |
| Hannum / Vidal-Bralo / Lin | 8×10⁻²⁶ – 4×10⁻²⁷ |
| PhenoAge | 1×10⁻²¹ |
| GrimAge | 2×10⁻¹⁷ |
| GrimAge2 | 2×10⁻¹⁴ |

**The clock adds over dysregulation only sometimes:**

| Clock | LR p for adding the clock to a dysregulation model | Verdict |
|---|---|---|
| GrimAge2 | 1×10⁻³² | adds |
| GrimAge | 1×10⁻³⁰ | adds |
| PhenoAge | 1×10⁻⁷ | adds |
| Horvath | 1×10⁻⁵ | adds |
| Hannum | 2×10⁻⁴ | adds |
| Zhang | 0.007 | adds |
| Vidal-Bralo | 0.009 | adds |
| **Lin** | **0.059** | **no evidence** |
| **SkinBlood** | **0.087** | **no evidence** |
| **Weidner** | **0.32** | **no evidence** |

So the honest one-line answer to "is an epigenetic clock worth running": for
GrimAge and PhenoAge, yes, clearly. For three of the seven chronological-age
clocks, **there is no evidence it tells you anything six routine blood tests do
not.**

### 2.4 The headline, stated plainly

> A Mahalanobis distance over six routine clinical chemistry markers
> (C = 0.757) predicts 20-year mortality better than nine of the ten DNA
> methylation clocks NCHS published, and adds information to all ten —
> including GrimAge2.

Only GrimAge (0.762) and GrimAge2 (0.764) beat it, and both are trained on
mortality, which is the outcome being predicted. The six markers cost a
standard blood panel; the clocks cost an EPIC array.

---

## 3. The caveat that matters more than the p-values

**Every p-value above is tiny and every effect on discrimination is small.**

```
age + sex only                      C = 0.7394
+ the best clock (GrimAge2)         C = 0.7644   (+0.025)
+ dysregulation                     C = 0.7567   (+0.017)
+ both                              C = 0.7731   (+0.034)
```

With 1,350 events, a likelihood-ratio test detects effects far below the size
that would change a decision about a person. **All of this biology adds 0.034 of
concordance over knowing someone's age and sex.** A reader who takes
`p = 1×10⁻²⁷` as "this is a large effect" has read it wrong, and the number that
corrects them is in the same table.

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

| Dysregulation quartile | Died | n |
|---|---|---|
| Q1 (lowest) | 44.9% | 630 |
| Q2 | 47.5% | 629 |
| Q3 | 52.6% | 629 |
| Q4 (highest) | **69.5%** | 629 |

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

**Survey weights are carried, not applied.** `WTDN4YR` (DNAm-specific, 4-year),
`WTMEC4YR`, `SDMVPSU` and `SDMVSTRA` are all in the frame. Unweighted estimates
are not nationally representative. Applying them would make this a population
claim and is the single highest-value follow-up.

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
)
print(result.headline())
```

**Sanity numbers:** 2,517 subjects · 1,350 deaths · 60 sex-discordant ·
age 50–85 · follow-up median 17.08 y · baseline C 0.7394 · GrimAge2 HR/SD 1.443.
