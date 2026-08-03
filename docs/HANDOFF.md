# GeroQuery — Handoff

**Date:** 2026-08-03 (updated same day — Tier 1.1 complete)
**Supersedes:** [`HANDOFF_TIER0_REVIEW.md`](HANDOFF_TIER0_REVIEW.md) — the original "tested skeleton, no real data" review, kept verbatim for context.
**Status:** Tier 0 and Tier 1.1 complete. **No synthetic data layer remains.**

---

## 0. TL;DR

The first handoff's verdict was *"the pipes are real, the water is fake."*
The second said the water was real everywhere except the gene signatures.

**There is now no layer of GeroQuery that is not backed by a checksum-pinned
upstream.** The only generated table left is `clinical_synthetic_csd`, which is
deliberate: it has critical slowing down planted in it so the estimator can be
tested against a known answer.

| | Original | After Tier 0 | Now |
|---|---|---|---|
| Clinical data | 720 synthetic subjects | 4,895 real NHANES adults | unchanged |
| Aging clocks | 3 hand-written demos | 236 real published clocks | unchanged |
| **Gene signatures** | synthetic | **still synthetic** | **485,905 real Hedges' *g* over 46,091 genes**, 31 pinned GEO DataSets |
| **Curated knowledge** | 12 hand-written genes | 12 hand-written genes | **2,720 real HAGR assertions** across 5 databases |
| **Interventions** | 7 hand-written | 7 hand-written | **1,340 real DrugAge / GenDR records** |
| Pinned artifacts | 0 | 4 | **41** (~320 MB, 5 upstreams) |
| Tests | 77 | 184 offline + 9 live, 82% | **295 offline + 12 live, 84%** |

Two headline scientific results, both published as found:

1. **On real NHANES**, health-state variance rises with age robustly (20/20
   configurations) and marker cross-correlation does not (0/20) — a partial null.
   [`RESULTS_NHANES_CSD.md`](RESULTS_NHANES_CSD.md)
2. **On the real GEO panel**, CDKN2A/p16 — the most-cited transcriptional marker
   of aging — does not replicate (*g* = +0.07, CI [-0.20, +0.35], k=14, I²=14%),
   while p21/CDKN1A, IGF1, the metallothioneins, and a mitochondrial/OXPHOS
   decline do. [`RESULTS_GEO_SIGNATURES.md`](RESULTS_GEO_SIGNATURES.md)

Everything in one place — sources, architecture, use cases, limitations:
[`OVERVIEW.md`](OVERVIEW.md).

---

## 1. Environment setup (read this first)

**Python 3.14 will not work for the clock layer.** Use 3.10–3.12.

The repo itself runs fine on 3.14; only the optional `clocks` extra is constrained (`ecos`, a transitive biolearn dependency, has no 3.14 wheel).

### The exact sequence that works — all free, no accounts, no tokens

```bash
# 1. Get a 3.12 toolchain (uv provisions it; ~20 MB)
uv python install 3.12
uv venv --python 3.12 .venv312

# 2. Install the project + both clock libraries
VIRTUAL_ENV=.venv312 uv pip install -e ".[dev,clocks]"

# 3. torch, CPU-only build (~200 MB instead of ~2.5 GB for the CUDA build)
VIRTUAL_ENV=.venv312 uv pip install torch --index-url https://download.pytorch.org/whl/cpu

# 4. Verify — expect 63 biolearn clocks, then 173 pyaging clocks
./.venv312/Scripts/python.exe -c "from geroquery.clocks.library import load_library_clocks; print(load_library_clocks()[1].to_dict())"
GEROQUERY_ALLOW_NETWORK=1 ./.venv312/Scripts/python.exe -c "from geroquery.clocks.pyaging_clocks import load_pyaging_clocks; print(load_pyaging_clocks()[1].to_dict())"
```

Windows: `.venv312/Scripts/python.exe`. macOS/Linux: `.venv312/bin/python`.

### Four upstream landmines — all hit, all now pinned

These cost real debugging time. They are encoded in `pyproject.toml`'s `clocks` extra with comments, but know them:

1. **biolearn under-declares its dependencies.** `biolearn.model` imports `torch` and `seaborn` at module scope but lists neither in `install_requires`. A clean `pip install biolearn` gives you a package that raises `ModuleNotFoundError` the moment you touch `ModelGallery`.
2. **biolearn 0.9.1 is not pandas-3 compatible.** DunedinPACE quantile-normalizes in place, writing into the DataFrame's numpy buffer. pandas 3 made copy-on-write mandatory, so that buffer is always read-only → `assignment destination is read-only`. **Nothing outside biolearn can work around this**; the extra pins `pandas>=2.0,<3`.
3. **pyaging's `@progress` decorator does `logger = args[-1]`.** So `load_clock_metadata(dir=..., logger=..., indent_level=...)` raises `IndexError`, and `load_clock_metadata("d", logger, 1)` makes it treat the int `1` as the logger. Exactly two positional args, logger last.
4. **pyaging raises a bare, message-less `NameError`** when 100% of a clock's features are missing. Unwrapped it reads as an internal bug rather than "your columns aren't CpG ids".

### Free-tier notes

- Everything above is free. No API keys, no accounts.
- pyaging downloads its catalog and per-clock artifacts from **Hugging Face unauthenticated**. It warns about rate limits; `HF_TOKEN` is optional and only speeds downloads.
- On Windows, HF caching warns about symlinks. Harmless (uses more disk). Silence with `HF_HUB_DISABLE_SYMLINKS_WARNING=1`, or enable Developer Mode.

---

## 2. What is real and what is not

**This table is the single most important thing in this document. Keep it accurate.**

| Layer | Status | Where |
|---|---|---|
| **Gene aging signatures** | **REAL** — 31 pinned GEO DataSets, 32 contrasts, 27 Series, 485,905 effect sizes | `sources/geo.py`, `harmonize/differential.py`, `etl/build_signatures.py` |
| **Curated knowledge** | **REAL** — GenAge, CellAge, LongevityMap, GenDR (2,720 mammalian assertions) | `sources/hagr.py` |
| **Interventions** | **REAL** — DrugAge + GenDR (1,340 records, 54 NIA ITP) | `sources/hagr.py` |
| Clinical / resilience | **REAL** — NHANES 2017-2018, n=4,895, SHA-256 pinned | `sources/nhanes.py`, `sources/manifest.py` |
| Aging clocks | **REAL** — 63 biolearn + 173 pyaging published clocks | `clocks/library.py`, `clocks/pyaging_clocks.py` |
| Gene ID resolution | **REAL** — batched mygene.info + disk cache (request path and bulk ETL path) | `idmap/mygene.py`, `idmap/bulk.py` |
| Tissue expression | **REAL** — GTEx Portal v2, median TPM + UBERON | `sources/gtex.py` |
| Meta-analysis maths | **REAL** — DerSimonian-Laird, and now fed real effect sizes for the first time | `harmonize/meta.py` |
| `clinical_synthetic_csd` | **SYNTHETIC ON PURPOSE** — CSD planted; method validation only | `etl/build_fixtures.py`, separate dataset id |
| recount3, ARCHS4, CELLxGENE, UK Biobank, GTEx-protected | **STUBS** — declare licence/capabilities, no fetch body | `sources/federated.py` |

Enforcement, not just documentation:
- Real and synthetic clinical data are **separate dataset ids** (`clinical_nhanes_slice` vs `clinical_synthetic_csd`), with `REAL`/`SYNTHETIC` prefixes in their descriptions, and a test asserts they never merge.
- `docs/DATA_SOURCES.md` has a **status column**: live / bundled / stub / planned.
- The React panel tags the resilience card `real data` and the signature card `demonstration slice`.

---

## 3. What was built this session

### 3.1 Data layer — pinned, verified, download-on-demand
- `sources/manifest.py` — every artifact declared with exact URL, byte count, SHA-256, licence, attribution.
- `sources/fetch.py` — **no path returns unverified bytes.** Cache hits are re-verified rather than trusted; downloads land on a temp file and are only renamed after the digest matches; failed verification deletes rather than leaving a valid-looking file. Offline raises `NetworkDisabledError` instead of silently reaching out.
- `sources/nhanes.py` — joins 4 XPT files on `SEQN` → the six-marker frame. 600-row **real** sample committed for offline CI; full cohort via `make data`. Returns `(frame, mode)` so a sample-derived number is never reported as the cohort result.

> **The NHANES URL pattern is a trap.** `wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT` returns HTML, not XPORT. The working pattern is `wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt`. A test pins this.

### 3.2 Resilience — defects fixed, result published
- `recovery.py` — negative AR(1) coefficients were being clipped into the decaying regime and reported as "fast recovery — high resilience". Now classified into four regimes with CIs.
- `csd.py` — within-stratum age detrending + **subject-level bootstrap CIs**. `resilience_declines` now requires both CIs to exclude zero. Previously a positive slope alone sufficed, which fires on noise roughly half the time per indicator.
- `control.py` — rank-truncated eigendecomposition; reports `cond(W)`, rank, unreachable fraction; refuses above `cond > 1e10` unless `strict=False`.
- `service.py` — **had a live `NameError` on every `control_energy` call** because nothing tested it. Fixed; that file is now at 100% coverage.

### 3.3 Clocks — two real libraries, correctly labelled
- **63 biolearn clocks**: Horvath, Hannum, PhenoAge, GrimAge, DunedinPACE.
- **173 pyaging clocks**, registered metadata-only (artifacts download lazily — eagerly fetching 173 would be hundreds of MB just to answer `GET /v1/clocks`).
- **`predicted_outcome` is preserved, never defaulted.** 19 of the 63 biolearn models are *not* aging clocks — BMI, smoking, cholesterol, prostate cancer, depression, sex, cell proportions. Defaulting unmapped outputs to `chronological_age` would present a prostate-cancer classifier as an aging clock. Same guard on the pyaging side (grip strength, VO2max, smoking exposure).
- **Cross-library validation:** fed byte-identical input, biolearn's `Horvathv1` returns 53.22 and pyaging's `horvath2013` returns 53.21. Two independent implementations, two different matrix orientations, same answer. A live test asserts agreement within 0.05.

Bugs found only by running the real libraries — a stub would have passed all of these:
- The transpose handed biolearn a **read-only** buffer → DunedinPACE failed.
- biolearn's `LinearModel.predict` does `matrix_data.loc["intercept"] = 1`, **mutating the caller's DataFrame**. Now copied.
- GrimAge needs **age and sex metadata**; the wrapper passed an empty frame, making every covariate-dependent clock unusable. `apply_clock(..., sample_metadata=...)` now plumbs it through.
- `library_clocks()` **swallowed every exception**, so "biolearn is installed but broken" was indistinguishable from "not installed". Now `LibraryStatus` reports the reason, surfaced at `GET /v1/clocks` under `tiers`. This is what diagnosed landmines 1 and 3 above.
- pyaging **imputes missing features from reference values** rather than failing, so a prediction can be produced from almost no real data. The wrapper now refuses above 50% imputation (`MAX_IMPUTED_FRACTION`) rather than returning a confident-looking artifact of pyaging's reference cohort.

### 3.4 Store / API defects
DuckDB connection reuse; glob-level partition pruning (`species=…/omic_layer=…` resolved in the path, because a bound `?` parameter defeats pruning); cached `ensure_built`; `except Exception` narrowed to real lookup failures; `InterventionNotFoundError` moved to `exceptions.py`; CI `workflow_dispatch` + weekly schedule so `smoke-live` can actually fire.

---

## 4. What to do next — prioritized

### TIER 1.1 — DONE

Real gene signatures shipped. See [`RESULTS_GEO_SIGNATURES.md`](RESULTS_GEO_SIGNATURES.md)
for the method, the panel-selection rule, the results, and the ten limitations.
Three things worth carrying forward from doing it:

- **The panel is rule-selected, not hand-picked.** Every GEO DataSet declaring an
  `age` subset variable (189 records) is run through the contrast rules in
  `sources/geo.py`; 31 survive. Change a rule and the panel changes. Resist the
  urge to add or drop an accession by hand — that is where selection bias enters.
- **The handoff predicted "real effect sizes will be much noisier and some
  textbook genes will not replicate."** Both happened. Median 95% CI width is
  0.84 and p16 is a clean null. It was handled the way the CSD null was.
- **The next source matters more than more of this one.** recount3 or ARCHS4
  would remove the biggest confound (cross-study processing heterogeneity) and the
  biggest coverage gap (no RNA-seq, no modern blood) in one move.

### TIER 1.1b — Real methylation, to close the clock loop *(now the highest leverage)*

236 clocks are wired and have never seen real data. That is the largest
unexercised capability in the repo, and unlike most gaps it is *checkable*:
published MAEs exist, so running Horvath / Hannum / PhenoAge / DunedinPACE on a
GEO methylation series with age labels (GSE40279 is the obvious one) tells you
whether the wrappers are correct rather than merely non-crashing.

The GEO ingestion layer built for Tier 1.1 mostly transfers: `sources/geo.py`
parses SOFT and builds contrasts, `sources/fetch.py` verifies, `manifest.py`
pins. What is new is that methylation series are GSE-only (no GDS), so the
free-text age parsing this build deliberately avoided has to be written for one
or two specific series — which is tractable when it is two series rather than
thirty.

Then compute age acceleration and correlate it against the resilience metrics.
That cross-layer claim is the reason this project exists.

### TIER 1.2 — Longitudinal resilience *(closes the biggest scientific caveat)*

The CSD result's central limitation is that NHANES is cross-sectional: age strata proxy for within-individual trajectories, and no relaxation time is observed. That is precisely what would distinguish critical slowing down from ordinary accumulating heterogeneity — the thing the current null cannot resolve.

Candidates, all free or application-based:
- **NHANES III linked mortality files** — free, adds a survival outcome to the existing pipeline.
- **ELSA / SHARE / HRS** — longitudinal aging cohorts, free registration.
- **UK Biobank** — the ideal dataset; requires an application and a fee. The federate-only adapter and licence gate are already in place.

`resilience/recovery.py` (AR(1) relaxation) is implemented and tested but has **never run on real longitudinal data.** It is the module with the most theoretical value and the least empirical exercise.

### TIER 1.3 — (folded into 1.1b above)

236 clocks are wired, but the repo has no real methylation matrix to run them on. Closing that loop is high-visibility and mostly plumbing:
- Pull a GEO methylation series with age labels (e.g. GSE40279, Hannum's cohort).
- Run Horvath / Hannum / PhenoAge / DunedinPACE; report MAE against chronological age.
- **This is a real, checkable validation** — published MAEs exist, so you find out whether the wrappers are correct rather than merely non-crashing.
- Then compute age acceleration and correlate it against the resilience metrics. That cross-layer claim is the reason this project exists.

### TIER 2 — Worth doing, not blocking

- **Survey weights.** `WTMEC2YR` is carried but not applied, so nothing is nationally representative. Applying it makes the CSD result a population claim.
- **Age topcoding.** NHANES topcodes at 80, compressing the stratum where CSD should be strongest. Pooling multiple NHANES cycles widens the effective range.
- **Confounder adjustment.** No fasting/medication/comorbidity adjustment. Some of the variance growth is likely disease prevalence rather than a dynamical property of aging.
- **GEO/recount3 adapters** to replace the remaining `FederatedStub`s.
- **Streamlit dashboard** has not been repointed at the real data (the React app has).

---

## 5. Open questions that need a human decision

1. **Is the CSD null the headline, or a stepping stone?** The README currently leads with it. That is honest and distinctive, but it means the marquee result is a partial negative. Leading with the *infrastructure* (verified data layer + 236 clocks) is also defensible. A positioning call, not a technical one.

2. ~~How much does the synthetic signature layer cost you?~~ **Resolved by Tier 1.1.** The replacement question is smaller and better: the marquee gene-level result is now a *non-replication* of p16. That is a real, defensible, and slightly contrarian finding. Lead with it, or lead with the infrastructure? Same positioning call as question 1.

3. **UK Biobank application?** The single dataset that would resolve the cross-sectional caveat properly. Weeks of latency and a fee, but the adapter seam already exists.

4. **Do the 173 pyaging clocks earn their keep?** They are metadata-only until used, so they cost nothing at rest — but they inflate `GET /v1/clocks` from 3 entries to 239. Consider a `?library=` filter or a curated default subset.

---

## 6. Gotchas for whoever picks this up

- **`GEROQUERY_ALLOW_NETWORK` is off by default.** Tests and CI are hermetic by construction. Set it to `1` for anything touching an upstream.
- **`get_registry()` memoizes a module-level registry.** Tests that fake a clock library must reset `clocks.registry._REGISTRY`, or they pass alone and fail in a full run. This actually happened.
- **Two interpreters, on purpose.** 3.14 runs the core suite without the clock libraries; 3.12 runs everything including them. Both must stay green — that is what proves the optional tier is genuinely optional.
- **Live tests are excluded by default** (`-m 'not live'`). Run with `GEROQUERY_ALLOW_NETWORK=1 pytest -m live`. They are the only thing that catches upstream drift: GTEx changing its API, CDC moving a file, a checksum going stale.
- **The GTEx age-bracket trap.** The open v2 API *accepts* `attributeSubset=ageBracket` and returns HTTP 200 — with one undivided group. Donor age is dbGaP-controlled. A live test pins this so the claim is rechecked automatically rather than from memory.
- **`.venv312/` should be gitignored** before committing.
- **`make signatures` takes ~25 minutes on a cold cache**, dominated by parsing 313 MB of gzipped SOFT and ~71k mygene.info lookups. Both cache; a rebuild is minutes.
- **`fetch.cache_path` is prefixed with the manifest key.** URL basenames collide: HAGR serves both DrugAge and GenDR from `/dataset.zip`. The checksum layer caught it, which is what it is for, but do not "simplify" the prefix away.
- **mygene's `ensembl` field is sometimes a list containing other species' ids.** A mouse gene's list can lead with `ENSFALG...` (flycatcher). `_ensembl_gene` filters by the taxid's expected prefix and falls back to `ENTREZ:`. This was a live bug that silently put flycatcher and ferret identifiers into the canonical id space.
- **`signatures_full.csv` is gitignored; `signatures_curated.csv` is committed.** `LocalSignatureSource` prefers the full one when present, so a developer who has run `make data` sees different counts from CI. Tests pin `prefer_full=False` for exactly this reason — do the same in any new test.

---

## 7. Reproducing the current state

```bash
make data          # fetch + verify every upstream, build every table (network, ~320 MB, ~25 min)
make signatures    # just the GEO + HAGR half of it
make data-offline  # same, from the committed real samples             (no network)
make test          # offline suite
GEROQUERY_ALLOW_NETWORK=1 pytest -m live     # live suite
make ci            # lint + typecheck + tests
```

Gate at time of writing: **ruff, black, mypy clean; 295 passed / 1 skipped on Python 3.14; 7 of 12 live tests pass on 3.14 (5 need the clock extra); 84% coverage.**
