# The aging clocks, run on real DNA methylation

**Date:** 2026-08-04
**Data:** two checksum-pinned GEO 450K blood series — GSE64495 (n=106 controls, ages 2.3–73.7) and GSE30870 (n=40, newborns vs nonagenarians)
**Code path:** `sources/methylation.py` → `clocks/library.py` + `clocks/pyaging_clocks.py`
**Reproduce:** `GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.fetch_artifacts GSE64495 GSE30870`, then run the clocks in the Python 3.12 environment

---

## 0. What this closes

Every previous handoff carried the same line: *236 clocks are wired and have
never run on real data.* They were validated against each other and against the
shape of their own metadata, which proves they do not crash. It does not prove
they are right.

That was the largest unexercised capability in the repository, and — unlike most
gaps — it is *checkable*, because published MAEs exist for these clocks and
because one of the two datasets ships the authors' own Horvath output per sample.

**150 of 209 attempted clocks ran. Horvath's clock reproduces the authors'
published per-sample values to r = 0.998, MAE 1.39 years.**

## 1. The strongest available check: agreement with published output

GSE64495's metadata carries `dna methylation age` — the Horvath clock as the
paper's authors ran it. Comparing against that is qualitatively better than
comparing against chronological age: it isolates *our wrapper* from the clock's
own accuracy. A clock that is 4 years off chronology should be 4 years off in the
same direction for us too; disagreement there is our bug.

| Clock | r vs published Horvath age | MAE vs published |
|---|---|---|
| **biolearn:Horvathv1** | **0.9984** | **1.39 y** |
| **pyaging:horvath2013** | **0.9984** | **1.40 y** |
| pyaging:pipekfilteredh | 0.9907 | 2.04 y |
| pyaging:pipekelasticnet | 0.9826 | 2.33 y |
| biolearn:Horvathv2 / pyaging:skinandblood | 0.9826 | 2.52 y |

Two independent implementations of Horvath 2013 — different libraries, different
matrix orientations, different normalization internals — land within 0.01 years
of each other and within 1.4 years of the published values. The residual gap is
what you would expect from normalization differences (the authors used Horvath's
own BMIQ-based pipeline; biolearn and pyaging each do their own).

**This is the result that says the wrappers are correct rather than merely
non-crashing.**

## 2. Accuracy against chronological age

48 clocks declare `chronological_age` in years and so have a meaningful MAE.

| Clock | MAE (y) | r | Published MAE |
|---|---|---|---|
| biolearn:Horvathv2 / pyaging:skinandblood | **2.46** | 0.987 | ~2.5 (skin & blood) |
| biolearn:AltumAge / pyaging:altumage | 2.50 | 0.986 | ~2.8 |
| pyaging:intrinclock | 2.50 | 0.984 | — |
| pyaging:pipekfilteredh | 3.02 | 0.978 | — |
| biolearn:Horvathv1 / pyaging:horvath2013 | 3.61 | 0.973 | **~3.6** |

Horvath 2013's published mean absolute error across tissues is about 3.6 years.
We measure **3.61**. Median MAE across all 48 age clocks is 7.04 years, which is
the expected spread once clocks trained on other tissues, other age ranges, and
other outcomes are included.

## 3. The extreme-contrast check

Correlation is easy to get right for the wrong reason when age spans 0–100. So
the second dataset is deliberately bimodal: cord blood at birth versus
nonagenarians, nothing in between. Every one of the 18 headline clocks separates
them in the correct direction.

| Clock | newborn | nonagenarian | MAE (y) | r |
|---|---|---|---|---|
| biolearn:Horvathv1 | **0.23** | 79.7 | 7.39 | 0.989 |
| biolearn:Horvathv2 | 0.22 | 82.0 | 5.90 | 0.997 |
| biolearn:Hannum | −4.84 | 92.4 | 5.03 | 0.995 |
| biolearn:AltumAge | 1.77 | 72.0 | 11.46 | 0.994 |
| biolearn:YingCausAge | −24.7 | 83.7 | 17.46 | 0.995 |
| biolearn:PhenoAge | −38.9 | 81.8 | 25.79 | 0.992 |
| biolearn:DunedinPACE | 1.04 | 1.17 | — (pace) | 0.464 |
| biolearn:DunedinPoAm38 | 0.75 | 0.98 | — (pace) | 0.837 |
| biolearn:DNAmTL | −7.49 | −9.37 | — (kb) | −0.967 |
| biolearn:PEDBE | 3.86 | 16.4 | 40.28 | 0.975 |

**Horvath's clock puts cord blood at 0.23 years.** That is the single most
convincing number here: it is not a correlation, it is an absolute prediction on
an absolute scale, and it is right.

Three readings that matter, and that a bare MAE column would hide:

- **PhenoAge and Hannum go negative on cord blood.** They are linear models fit
  on adults; extrapolating below their training range produces negative ages.
  That is a property of the clocks, not a defect in the wrapper, and it is why
  MAE alone is a poor summary for them.
- **PEDBE has the worst MAE (40 y) and is behaving correctly.** It is a
  *paediatric* clock. Asked about 90-year-olds, it saturates. A pipeline that
  ranked clocks by MAE and quietly dropped the worst would be discarding the one
  clock that is right about children.
- **DunedinPACE's weak correlation (0.46) is not a failure.** It predicts pace of
  aging, a rate, not an age — and it still moves in the right direction
  (1.04 → 1.17). Its `predicted_outcome` is `pace_of_aging`, which is exactly why
  that field is never defaulted.

`DNAmTL` (telomere length) correlates **−0.967** with age. Telomeres shorten. The
sign is the check.

## 4. What did not run, and why

59 of 209 attempts failed. Every failure is catalogued rather than swallowed, and
they fall into four kinds — none of which is "the wrapper is broken":

| Cause | n | Example |
|---|---|---|
| Clock is not a methylation clock at all | ~25 | `pyaging:bitage` (transcriptomic), `pyaging:camilloh3k27ac` and the other histone-ChIP clocks |
| Needs EPIC-array probes absent from 450K | ~8 | `biolearn:DeconvoluteBloodEPIC`, `TwelveCellDeconvoluteBloodEPIC` |
| Needs specific CpGs this series does not carry | ~20 | `biolearn:DownSyndrome`, `biolearn:Pasta`, `biolearn:REG` |
| Optional dependency unavailable | 6 | the `GPAge*` family needs `GPy` |

The `GPy` case is worth recording because it is a genuine, unresolvable
conflict rather than a missing install: **GPy requires numpy < 2, and biolearn's
`cvxpy` dependency requires numpy ≥ 2.** Installing GPy to gain 6 clocks breaks
the 63 biolearn clocks. It was attempted, verified to break the stack, and
reverted. Those six stay unavailable until GPy supports numpy 2.

## 5. The bug this exposed

The first full run failed **236 times out of 236**, with every clock reporting
`requires features not present in the input`. That message is a lie the code told
itself; there were two real causes.

1. **Orientation.** The adapter stores betas CpG-by-sample, which is how GEO ships
   them. Every clock wrapper takes samples-as-rows. `MethylationDataset.matrix()`
   now does the conversion in one place, and a test pins it.
2. **The preparation was deleting the clock's own coefficients.** `read_series`
   originally dropped every CpG with any missing value — 7,510 of them in
   GSE30870, and three were Horvath coefficients. The clock then correctly
   reported that its features were absent. Scattered gaps are now mean-imputed
   per CpG and the imputed fraction reported: **0.04% of cells in GSE30870**,
   which is small enough to be uninteresting and is printed anyway so it cannot
   quietly stop being small.

Both are the same class of mistake as the CDKN1A identifier bug in the expression
panel: a *data-preparation* fault wearing the costume of an incompatible input.

## 6. Limitations

1. **Two series, both blood, both 450K.** Nothing here says anything about clock
   performance in other tissues, on EPIC arrays, or on sequencing-based
   methylation.
2. **n = 106 and n = 40.** Small. The MAEs have real sampling error; treat the
   third decimal as noise.
3. **GSE64495's age range stops at 73.7.** The oldest stratum, where clocks
   typically degrade, is represented only by GSE30870's nonagenarians.
4. **Normalization is not the authors'.** Each library applies its own; we did not
   reimplement Horvath's BMIQ pipeline. The 1.39-year gap in §1 is mostly this.
5. **The `sex` covariate is approximate.** GSE30870 does not record sex, so it is
   passed as a constant; covariate-dependent clocks (GrimAge) are therefore
   reported with correlation but their absolute values should not be read.
6. **150 of 209, not 236 of 236.** The full pyaging sweep on GSE64495 was still
   running when this was written; the remaining clocks are the slow neural models.
   Results are written incrementally to `clock_results.csv` and the run resumes.
7. **No mortality or morbidity outcome.** Neither series has follow-up, so the
   mortality-predicting clocks (GrimAge, PhenoAge) can be checked for internal
   consistency and direction but not for what they were built to predict.

## 7. What this makes possible next

The blocking dependency for the project's central claim is now gone. Age
acceleration — clock age minus chronological age — can be computed on real
samples, and the next step is to correlate it against the resilience metrics
computed on NHANES. That cross-layer link is the reason GeroQuery exists and is
still the one headline claim with no evidence behind it in either direction.

---

## Provenance

Both series are declared in `geroquery/sources/manifest.py` with exact URL, byte
count, and SHA-256, verified 2026-08-04. `sources/methylation.py` validates each
series' declared characteristic keys and array platform at parse time, so an
upstream format change fails loudly rather than producing a frame of NaN ages.
