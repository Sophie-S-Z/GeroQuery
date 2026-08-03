# Critical slowing down in NHANES 2017–2018: one signal replicates, one does not

**Headline.** In 4,895 real NHANES adults aged 20–80, the **variance** of the
multi-marker health state rises with age robustly across every analytic choice we
tried (20/20 configurations, bootstrap 95% CI excludes zero in all 20). The
**cross-correlation** among those same markers shows **no evidence of a trend**
(0/20 configurations). Critical slowing down predicts both. We find one.

We report this because it is what the data says. A synthetic fixture in this repo
shows both signals cleanly — and that fixture was built with the effect planted in
it, so it was never evidence about human biology. Keeping the null visible next to
it is the point.

---

## 1. Data

| | |
|---|---|
| Source | NHANES 2017–2018 (release "J"), CDC/NCHS |
| Files | `DEMO_J`, `BIOPRO_J`, `HSCRP_J`, `CBC_J` — pinned by SHA-256 in `geroquery/sources/manifest.py` |
| Join | `SEQN` (inner) |
| Filter | age ≥ 20, complete cases on all six markers (no imputation) |
| n | **4,895** subjects |
| Age range | 20–80 (**topcoded at 80**; 354 subjects sit at the topcode) |

Markers: albumin, creatinine, glucose (BIOPRO_J), high-sensitivity CRP (HSCRP_J),
lymphocyte percent, RDW (CBC_J).

Reproduce with:

```bash
make data                       # fetch + SHA-256 verify + harmonize
python -m geroquery.etl.build_data
```

Nothing here was computed from a generator. Every row is a real measurement on a
real participant.

## 2. Method

Cross-sectional proxy for critical slowing down, as used in the aging literature:
subjects are sorted by age into *K* equal-count strata, and within each stratum we
compute

* **variance** of the health state (mean of globally z-scored markers), and
* **cross-correlation**: mean absolute off-diagonal Pearson correlation among markers.

Two choices do real work here:

**Within-stratum age detrending.** Markers drift with age, so age heterogeneity
*inside* a bin inflates that bin's variance and can manufacture a rising-variance
trend with no critical slowing down present. We regress age out within each stratum
before computing either indicator, and report the undetrended series alongside.

**Subject-level bootstrap CIs, not Kendall τ.** A trend over ~6 points has almost
no resolution. Empirically, with a *perfect* monotone increase:

| K | best achievable Kendall p | can ever reach p < 0.05? |
|---|---|---|
| 4 | 0.083 | **no** |
| 6 | 0.0028 | yes, but that is the floor |
| 8+ | < 0.001 | yes |

At K = 4 the test cannot return a significant result no matter how clean the trend
is. Conversely, a bare `slope > 0` check fires on pure noise roughly half the time.
So significance is established by resampling **subjects** — the actual sampling unit —
800 times and taking a percentile interval on the slope across strata. A signal
counts as supported only when its 95% CI lies entirely above zero.

## 3. Result

Default configuration (K = 6, detrended, raw scale):

| Stratum midpoint (yrs) | n | Variance | Cross-correlation |
|---|---|---|---|
| 25.3 | 816 | 0.0984 | 0.0967 |
| 36.2 | 816 | 0.0851 | 0.1132 |
| 47.4 | 816 | 0.1854 | 0.1150 |
| 57.5 | 816 | 0.1166 | 0.1230 |
| 65.6 | 816 | 0.1452 | 0.1516 |
| 76.9 | 815 | 0.1651 | 0.1187 |

* **Variance**: slope **+0.00128 / year**, bootstrap 95% CI **[+0.00045, +0.00219]** — excludes zero. **Supported.**
* **Cross-correlation**: slope +0.00064 / year, 95% CI **[−0.00039, +0.00123]** — includes zero. **Not supported.**
* `resilience_declines` → **False** (it requires *both*).

Note the variance column is not monotone (the 47.4 stratum overshoots, 57.5 dips).
The trend is real but noisy; this is why the claim rests on an interval and not on
the shape of six numbers.

## 4. Sensitivity: 20 configurations

`log-transform of skewed markers (CRP, creatinine, glucose) ∈ {raw, log}` ×
`within-stratum age detrending ∈ {no, yes}` × `K ∈ {4, 6, 8, 10, 12}`,
800 subject-level bootstrap resamples each.

| scale | detrend | K | var slope | var 95% CI | var supported | xcorr slope | xcorr 95% CI | xcorr supported |
|---|---|---|---|---|---|---|---|---|
| raw | no | 4 | +0.00111 | [+0.00038, +0.00215] | **yes** | +0.00089 | [-0.00021, +0.00155] | no |
| raw | no | 6 | +0.00128 | [+0.00044, +0.00221] | **yes** | +0.00065 | [-0.00036, +0.00123] | no |
| raw | no | 8 | +0.00125 | [+0.00060, +0.00206] | **yes** | +0.00032 | [-0.00043, +0.00091] | no |
| raw | no | 10 | +0.00124 | [+0.00057, +0.00206] | **yes** | +0.00014 | [-0.00057, +0.00074] | no |
| raw | no | 12 | +0.00137 | [+0.00062, +0.00214] | **yes** | +0.00025 | [-0.00053, +0.00082] | no |
| raw | yes | 4 | +0.00112 | [+0.00039, +0.00216] | **yes** | +0.00084 | [-0.00027, +0.00145] | no |
| raw | yes | 6 | +0.00128 | [+0.00045, +0.00219] | **yes** | +0.00064 | [-0.00039, +0.00123] | no |
| raw | yes | 8 | +0.00126 | [+0.00061, +0.00206] | **yes** | +0.00032 | [-0.00046, +0.00089] | no |
| raw | yes | 10 | +0.00124 | [+0.00058, +0.00205] | **yes** | +0.00014 | [-0.00056, +0.00075] | no |
| raw | yes | 12 | +0.00136 | [+0.00063, +0.00214] | **yes** | +0.00026 | [-0.00052, +0.00084] | no |
| log | no | 4 | +0.00088 | [+0.00053, +0.00126] | **yes** | +0.00030 | [-0.00041, +0.00091] | no |
| log | no | 6 | +0.00096 | [+0.00057, +0.00134] | **yes** | +0.00010 | [-0.00059, +0.00066] | no |
| log | no | 8 | +0.00090 | [+0.00057, +0.00127] | **yes** | -0.00009 | [-0.00067, +0.00054] | no |
| log | no | 10 | +0.00091 | [+0.00057, +0.00128] | **yes** | -0.00012 | [-0.00073, +0.00049] | no |
| log | no | 12 | +0.00095 | [+0.00057, +0.00130] | **yes** | -0.00005 | [-0.00071, +0.00054] | no |
| log | yes | 4 | +0.00089 | [+0.00054, +0.00128] | **yes** | +0.00023 | [-0.00049, +0.00083] | no |
| log | yes | 6 | +0.00097 | [+0.00057, +0.00132] | **yes** | +0.00008 | [-0.00062, +0.00066] | no |
| log | yes | 8 | +0.00091 | [+0.00057, +0.00128] | **yes** | -0.00009 | [-0.00070, +0.00053] | no |
| log | yes | 10 | +0.00090 | [+0.00056, +0.00128] | **yes** | -0.00012 | [-0.00072, +0.00049] | no |
| log | yes | 12 | +0.00095 | [+0.00058, +0.00130] | **yes** | -0.00004 | [-0.00072, +0.00054] | no |

**Tally: variance supported 20/20. Cross-correlation supported 0/20** (and its
point slope is not even consistently positive — 14/20 — flipping sign on the log
scale at K ≥ 8).

Two things worth noting:

* **Log-transforming shrinks the variance slope by ~25% and tightens the CI sharply**
  (+0.00128 → +0.00097 at K=6; CI width 0.0017 → 0.0008). Some of the raw-scale
  variance growth is driven by a few extreme CRP/glucose values. The effect survives
  the transform, but it is smaller than the raw scale suggests.
* **Detrending changed almost nothing** (+0.00128 → +0.00128 at K=6). The
  within-stratum age confound was not driving this result. The detrending fix is
  still correct — it just is not what mattered here.

## 5. Comparison with the synthetic fixture

Same estimator, same settings (K = 6, detrended, 800 bootstraps), run on
`clinical_synthetic_csd` — the generated fixture in this repo:

| Stratum midpoint | Variance (synthetic) | Cross-correlation (synthetic) | Cross-correlation (real) |
|---|---|---|---|
| ~25 | 0.0481 | 0.4254 | 0.0967 |
| ~35 | 0.0669 | 0.4708 | 0.1132 |
| ~45 | 0.0980 | 0.7083 | 0.1150 |
| ~56 | 0.1530 | 0.7660 | 0.1230 |
| ~68 | 0.1736 | 0.7854 | 0.1516 |
| ~79 | 0.1900 | 0.8153 | 0.1187 |

Synthetic: variance slope +0.00284 CI [+0.00218, +0.00357] **supported**;
cross-correlation slope +0.00769 CI [+0.00613, +0.00908] **supported**;
`resilience_declines` → **True**.

The synthetic cross-correlation nearly doubles across the age range (0.43 → 0.82).
The real one moves from 0.097 to 0.119 and is not monotone.

**Why the fixture shows both.** It is generated from a single shared latent factor
whose variance is hard-coded to grow with age (`0.4 + 1.8 * aging` in
`etl/build_fixtures.py`). One growing common factor inflates the state variance
*and* the marker cross-correlation simultaneously, by construction. It could not
have come out any other way.

So the synthetic result is **method validation** — evidence that the estimator
recovers an effect known to be present — and nothing more. It is registered in the
store as a separate dataset (`clinical_synthetic_csd`, description prefixed
`SYNTHETIC`) precisely so it cannot be mistaken for the real one
(`clinical_nhanes_slice`, prefixed `REAL`).

## 6. What this does and does not show

**Does show.** Between-individual dispersion in a six-marker health state increases
with age in a large real human cohort, robustly to stratification resolution,
skew correction, and within-stratum age adjustment.

**Does not show.**

1. **Not a tipping point.** Rising variance is *consistent with* critical slowing
   down, but it is also consistent with plain accumulating heterogeneity — people
   diverging in health with age for reasons having nothing to do with a bifurcation.
   The cross-correlation signal is what would have distinguished these, and it is
   absent.
2. **Not within-individual dynamics.** NHANES is cross-sectional. Age strata are a
   proxy for trajectories; no relaxation time is observed. Discriminating the two
   explanations above needs longitudinal data (see caveat 2).
3. **Not nationally representative.** Survey weights (`WTMEC2YR`) are carried in the
   data as `survey_weight` but **not applied**. All estimates here are unweighted.

### Caveats carried in code

These are in `geroquery/sources/nhanes.py:CAVEATS`, not only in this document:

1. **Age topcoded at 80.** This compresses exactly the oldest stratum — the one where
   any critical slowing down should be strongest. The effect is measured on a
   truncated age range and could plausibly be underestimated.
2. **Cross-sectional design.** No follow-up, so no within-individual relaxation time.
3. **Survey weights not applied.**
4. **No adjustment** for fasting status, medication, or comorbidity. Some of the
   variance growth is likely disease prevalence rather than a dynamical property of
   aging.

## 7. Honest summary

A method that reproduces both textbook early-warning signals on synthetic data
reproduces one of two on real data. The variance signal is solid. The
cross-correlation signal is not there at this sample size and design. We are not
going to hunt across specifications until it appears — the 20-configuration sweep
above is the whole search, reported in full.

`resilience_declines` returns `False` on this dataset. That is the correct answer.

---

*Generated through the repository code path: `geroquery.sources.nhanes` →
`geroquery.resilience.ResilienceService.csd`. Manifest version 2018.0, checksums
verified 2026-08-03.*
