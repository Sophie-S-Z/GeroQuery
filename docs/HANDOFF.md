# GeroQuery — complete handoff

**Date:** 2026-08-04 · **Branch:** `main` · **Manifest:** 2026.2 (checksums verified 2026-08-04)
**Status:** every data layer is backed by a checksum-pinned upstream. One synthetic table remains, deliberately.

This is the authoritative document. Read §1 and §3 if you read nothing else.

---

## 1. Where the project stands

The first review of this repository called it *"the pipes are real, the water is
fake."* It was a well-architected, fully-tested system over generated data:
twelve genes with hand-assigned effect sizes, fabricated GEO accessions
(`GSE90001`–`GSE90228`), synthetic p-values, and a clinical cohort labelled
"NHANES" that was a random draw.

Everything below is now measured or ingested.

| Layer | Then | Now |
|---|---|---|
| Gene aging signatures | 12 genes, invented effect sizes and accessions | **485,905 Hedges' *g* estimates over 46,091 genes** from 31 pinned GEO DataSets → 32 contrasts → 27 independent GEO Series |
| Clinical / resilience | 720 synthetic subjects | **4,895 real NHANES 2017–2018 adults**, SHA-256 pinned |
| Curated knowledge | 12 hand-written genes | **2,720 assertions over 1,880 genes** from five HAGR databases |
| Interventions | 7 hand-written | **1,340 records** from DrugAge + GenDR, 54 flagged NIA ITP |
| Aging clocks | 3 hand-written demos | **240 clocks**, and **436 clock-dataset runs on real methylation** |
| DNA methylation | none | **2 pinned GEO 450K blood series** (n=106 + n=40) |
| Gene identifier table | 24 hand-picked genes | **2,560 records / 1,835 ortholog groups**, generated from HAGR + mygene |
| Species lifespans | 9 hand-transcribed | parsed from the pinned AnAge release |
| Pinned artifacts | 0 | **44** (~679 MB, 5 upstreams) |
| Tests | 77 | **321 offline + 12 live**, 78% coverage |

**The only synthetic table left is `clinical_synthetic_csd`**, and it exists on
purpose — see §5.2.

### The three scientific results, all published as found

1. **On real NHANES**, health-state variance rises with age in 20/20 analytic
   configurations and marker cross-correlation in 0/20. A partial null.
   → [`RESULTS_NHANES_CSD.md`](RESULTS_NHANES_CSD.md)
2. **On the real GEO panel**, CDKN2A/p16 — the most-cited transcriptional marker
   of aging — does not replicate (*g* = +0.07, CI [−0.20, +0.35], k=14, I²=14%),
   while p21/CDKN1A, IGF1, the metallothioneins and a mitochondrial decline do.
   → [`RESULTS_GEO_SIGNATURES.md`](RESULTS_GEO_SIGNATURES.md)
3. **On real methylation**, Horvath's clock reproduces the study authors' own
   published per-sample ages to **r = 0.998, MAE 1.39 y**, matches its published
   3.6-year MAE against chronological age to three significant figures, and puts
   cord blood at **0.23 years**.
   → [`RESULTS_METHYLATION_CLOCKS.md`](RESULTS_METHYLATION_CLOCKS.md)

---

## 2. Environment setup (read this before touching the clock layer)

**Python 3.14 will not work for the clock layer.** Use 3.10–3.12. The repo itself
runs fine on 3.14; only the optional `clocks` extra is constrained (`ecos`, a
transitive biolearn dependency, has no 3.14 wheel).

### The exact sequence that works — all free, no accounts, no tokens

```bash
uv python install 3.12
uv venv --python 3.12 .venv312
VIRTUAL_ENV=.venv312 uv pip install -e ".[dev,clocks]"
VIRTUAL_ENV=.venv312 uv pip install torch --index-url https://download.pytorch.org/whl/cpu
```

Windows: `.venv312/Scripts/python.exe`. macOS/Linux: `.venv312/bin/python`.
Expect 63 biolearn clocks and 173 pyaging clocks.

### Five upstream landmines — all hit, all pinned

1. **biolearn under-declares its dependencies.** `biolearn.model` imports `torch`
   and `seaborn` at module scope but lists neither in `install_requires`. A clean
   `pip install biolearn` raises `ModuleNotFoundError` the moment you touch
   `ModelGallery`.
2. **biolearn 0.9.1 is not pandas-3 compatible.** DunedinPACE quantile-normalizes
   in place, writing into the DataFrame's numpy buffer; pandas 3 made
   copy-on-write mandatory, so that buffer is always read-only. **Nothing outside
   biolearn can work around this**; the extra pins `pandas>=2.0,<3`.
3. **pyaging's `@progress` decorator does `logger = args[-1]`.** Exactly two
   positional args, logger last.
4. **pyaging raises a bare, message-less `NameError`** when 100% of a clock's
   features are missing. Unwrapped it reads as an internal bug rather than "your
   columns aren't CpG ids".
5. **GPy and biolearn cannot coexist.** GPy requires numpy < 2; biolearn's `cvxpy`
   requires numpy ≥ 2. Installing GPy to gain the 6 `GPAge*` clocks breaks all 63
   biolearn clocks. Attempted, verified, reverted. **Do not retry** until GPy
   supports numpy 2.

### Free-tier notes

Everything is free — no API keys, no accounts. pyaging downloads its catalog and
per-clock artifacts from Hugging Face unauthenticated; `HF_TOKEN` is optional and
only speeds downloads. On Windows, HF warns about symlinks — harmless; silence
with `HF_HUB_DISABLE_SYMLINKS_WARNING=1`.

---

## 3. Complete data inventory

Every row is fetched from a named upstream and verified against a SHA-256 pinned
in `geroquery/sources/manifest.py`. Nothing enters the store unverified.

### 3.1 Pinned artifacts — 44 files, ~679 MB

| Group | Files | Size | What it is |
|---|---|---|---|
| **GEO expression** | 31 | 313.0 MB | GDS full SOFT files — the aging-signature panel |
| **GEO methylation** | 2 | 358.2 MB | 450K series matrices — the clock validation panel |
| **NHANES 2017–2018** | 4 | 7.2 MB | DEMO_J, BIOPRO_J, HSCRP_J, CBC_J (XPORT) |
| **HAGR** | 6 | 0.2 MB | GenAge human + models, CellAge, LongevityMap, DrugAge, GenDR |
| **AnAge** | 1 | 0.2 MB | Maximum-lifespan estimates across species |

### 3.2 Live-queried, not pinned

| Source | Role |
|---|---|
| **mygene.info** | Identifier harmonization → canonical Ensembl. ~71k identifiers resolved at build; batched POST, disk-cached |
| **GTEx Portal v2** | Tissue expression context, median TPM + UBERON, ~54 tissues/gene |
| **biolearn** (63 clocks) | Published epigenetic + clinical clocks |
| **pyaging** (173 clocks) | Published clocks, metadata-registered, artifacts lazy |

### 3.3 Derived, committed artifacts

Regenerable, not hand-written. Committed because the offline path must work.

| File | Size | Generated by |
|---|---|---|
| `idmap/data/genes.json` | 0.88 MB | `etl/build_idmap.py` |
| `idmap/data/anage.json` | 2 KB | `etl/build_idmap.py` |
| `sources/data/signatures_curated.csv` | 5.05 MB | `etl/build_signatures.py` |
| `sources/data/curated_knowledge.csv` | 2.09 MB | `etl/build_signatures.py` |
| `sources/data/interventions.csv` | 0.19 MB | `etl/build_signatures.py` |
| `sources/data/studies.csv` | 10 KB | `etl/build_signatures.py` |
| `sources/data/nhanes_2017_2018_sample.csv` | 50 KB | `etl/build_data.py --write-sample` |
| `sources/data/signatures_full.csv` | 60.7 MB | **git-ignored** — reproducible from the manifest |
| `sources/data/clinical_nhanes_full.csv` | 0.43 MB | **git-ignored** — same reason |

`idmap/data/tissues.json` (12 UBERON terms) is a controlled vocabulary, not data.

### 3.4 Contracted but not implemented

`sources/federated.py` declares these with capabilities and licence so the store,
API, and licence tests can reason about them. Querying one raises.

recount3 · ARCHS4 · CELLxGENE Census · GTEx protected tier · UK Biobank · dbGaP ·
ADNI · MIMIC · HRS-restricted

---

## 4. How everything works right now

### 4.1 Architecture

Nine modules, strictly layered. Lower layers never import higher ones.

```
ui  ->  api  ->  ( store, clocks, resilience, harmonize )  ->  ( sources, idmap, knowledge )
                    M4      M5        M6          M3              M2       M1
```

| Module | Files | Responsibility |
|---|---|---|
| **idmap** | `resolver.py`, `mygene.py`, `bulk.py` | Any identifier → canonical Ensembl; UBERON tissues; fractional age via AnAge |
| **sources** | `manifest.py`, `fetch.py`, `geo.py`, `methylation.py`, `nhanes.py`, `hagr.py`, `gtex.py`, `local_fixture.py`, `federated.py`, `base.py` | Uniform adapters + the checksum-verifying fetch layer + the cache-vs-federate gate |
| **harmonize** | `effect_size.py`, `differential.py`, `meta.py`, `batch.py` | Hedges' *g*, per-gene contrasts + BH FDR, DerSimonian–Laird pooling |
| **store** | `store.py` | Partitioned Parquet via DuckDB + SQLite metadata + content-hash versioning |
| **clocks** | `registry.py`, `library.py`, `pyaging_clocks.py`, `phenoage.py`, `service.py` | 240 clocks; predicted-outcome metadata; age acceleration; mortality risk |
| **resilience** | `csd.py`, `recovery.py`, `control.py`, `service.py` | Critical slowing down, AR(1) relaxation, network control energy |
| **knowledge** | `references.py`, `hallmarks.py` | Citation layer + López-Otín vocabulary over measured data |
| **api** | `app.py`, `service.py`, `schemas.py`, `cache.py` | FastAPI `/v1`, error envelope, pagination, format switch, versioned LRU |
| **etl** | `build_data.py`, `build_signatures.py`, `build_idmap.py`, `build_fixtures.py`, `fetch_artifacts.py` | Offline batch. Never on the request path |
| **ui** | `streamlit_app.py`, `frontend/` | Streamlit dashboard + React showcase |

### 4.2 The signature pipeline, end to end

```
GEO DataSets query (age subset variable)  ->  189 candidate records
        |  contrast rules in sources/geo.py
        v
31 datasets -> 32 contrasts (27 GEO Series)
        |  harmonize/differential.py
        v
per-gene Hedges' g + BH q, per contrast
        |  idmap/bulk.py  (Entrez -> Ensembl)
        v
485,905 rows -> partitioned Parquet (species / omic_layer)
        |  harmonize/meta.py  at query time
        v
random-effects pooled estimate with CI, I², k
```

**The panel is rule-selected, not hand-picked.** Ask GEO for every curated
dataset declaring an `age` subset variable in human or mouse (189 records); keep
those yielding an adult young-vs-old contrast (human 18–40 vs ≥60, mouse 2–8 mo
vs ≥18 mo, middle deliberately empty); restrict to the control arm of every other
subset variable, dropping the dataset if no control level is recognisable; split
rather than pool across tissue; require ≥3 per group. Change a rule and the panel
changes. **Do not add or drop an accession by hand — that is where selection bias
enters.**

### 4.3 The clock layer

| Tier | Count | Availability |
|---|---|---|
| biolearn | 63 | Needs the `clocks` extra (Python 3.10–3.12) |
| pyaging | 173 | Same; artifacts download lazily from Hugging Face |
| **PhenoAge** (real Levine coefficients) | 1 | **Always** — no optional dependency |
| Transparent `*_demo` reference clocks | 3 | Always; test instruments, not clinical |

`predicted_outcome` is **preserved, never defaulted**. 19 of the 63 biolearn
models predict BMI, smoking, cholesterol, prostate cancer, depression, or sex.
Defaulting them to `chronological_age` would present a cancer classifier as an
aging clock.

`clocks/service.py` also surfaces `mortality_risk_10yr` for any clock exposing it
— currently PhenoAge, which yields a biological age and a calibrated 10-year
mortality risk from one fit.

### 4.4 Data-layer guarantees, enforced in code

1. **Nothing enters the store unverified.** `fetch_artifact` returns bytes whose
   SHA-256 matches the manifest, or raises. **Cache hits are re-verified, not
   trusted.** Downloads land on a temp file and are renamed only after the digest
   matches; a failed verification deletes rather than leaving a valid-looking file.
2. **Offline is hermetic by construction.** With `GEROQUERY_ALLOW_NETWORK` unset
   (the default), a cache miss raises `NetworkDisabledError`. Tests and CI cannot
   silently depend on an upstream.
3. **Licence gate before persistence.** The store calls `assert_cacheable()`
   before writing; controlled sources refuse and a test proves it.
4. **Mode is returned, not logged.** The NHANES and signature adapters report
   whether a number came from the full table or the committed slice, so a
   slice-derived estimate can never be reported as the full-cohort result.
5. **Real and synthetic never merge.** Separate dataset ids, `SYNTHETIC`
   description prefix, and a test.
6. **The data version is a content digest**, so identical data always reports the
   same version.

### 4.5 API surface (`/v1`)

| Endpoint | Returns |
|---|---|
| `GET /gene/{id}/signature` | Per-study signatures + pooled estimates (CI, I², k) |
| `GET /gene/{id}/card` | Signatures + curated flags + linked interventions |
| `POST /geneset/signature` | Batch, with unresolved genes reported separately |
| `GET /clocks` | All registered clocks with `predicted_outcome` + library status |
| `POST /clock/apply` · `/clock/compare` | Run clocks; age acceleration; mortality risk |
| `GET /intervention/{name}` | Every organism tested, mammals first |
| `POST /resilience/csd` · `/resilience/recovery` | Resilience metrics with bootstrap CIs |
| `GET /studies` · `/sources` · `/datasets` · `/version` | Provenance |
| `GET /healthz` | Liveness |

Service-layer methods not yet exposed over HTTP: `gene_report`,
`list_curated_genes`, `references` — the Streamlit dashboard calls these directly.

All list endpoints support `limit`/`offset`; data endpoints support
`format=json|csv|parquet`. Every error uses one envelope:
`{"error": {"code", "message", "detail"}}`.

### 4.6 Datasets in the store

| Dataset id | Rows | What |
|---|---|---|
| `clinical_nhanes_slice` | 4,895 | Real NHANES, six resilience markers |
| `clinical_nhanes_phenoage` | 4,894 | Real NHANES, the nine markers PhenoAge needs |
| `clinical_synthetic_csd` | 720 | **Synthetic on purpose** — see §5.2 |

---

## 5. What was built, in detail

### 5.1 New modules

**`sources/geo.py`** — GEO DataSets adapter. Uses GEO *DataSets* rather than
*Series* because a GDS is the curator-built view: age is a declared subset
variable with explicit sample lists, and the matrix ships with platform
annotation joined. Both fragile steps of GEO ingestion — deciding who is old,
deciding what a probe measures — are answered by GEO's curators rather than by a
regex of ours.

**`sources/methylation.py`** — GEO Series adapter for 450K methylation. Uses
*Series* because no GDS exists for methylation, so free-text age parsing is
unavoidable — acceptable for two hand-verified series in a way it would not be
for thirty. Validates each series' declared characteristic keys and array
platform at parse time, so an upstream format change fails loudly instead of
producing a frame of NaN ages.

**`sources/hagr.py`** — the five HAGR databases. Two boundaries enforced in code:
non-mammalian rows are parsed, counted, and dropped rather than arriving labelled
as if mammalian; and **DrugAge records no gene targets, so none are asserted**
(the old table linked rapamycin → MTOR by hand — that edge was never in the
data). LongevityMap's non-significant entries are kept: a gene list built only
from hits cannot tell you a gene was tested and did not replicate.

**`harmonize/differential.py`** — vectorized per-gene Hedges' *g* + BH FDR. Three
judgement calls are named constants rather than inline literals: the
scale-detection threshold (GEO's declared `value_type` is not reliable enough to
decide whether a matrix is log-scaled), the log floor, and the minimum group
size. `test_matrix_hedges_g_matches_scalar_implementation` pins it
element-by-element against `effect_size.hedges_g` so the fast copy cannot drift.

**`idmap/bulk.py`** — bulk resolution for the ETL. Keeps one consolidated map per
species instead of ~50,000 sub-kilobyte cache files. Unresolvable identifiers are
stored as `None` so a rebuild does not re-ask a free public service about the same
thousands of failures. Unlike the request path, an offline cache miss is **fatal**
— a build that silently dropped unresolved genes would write a plausible-looking
but incomplete table.

**`clocks/phenoage.py`** — the real published Levine coefficients, plus a
column-alias layer so GeroQuery's marker names and the coefficients' unit-suffixed
names can talk without either side renaming.

**`etl/build_signatures.py`**, **`etl/build_idmap.py`** — the two ETLs producing
everything derived.

**`knowledge/`** — `references.py` (28 real papers, never an invented PMID) and
`hallmarks.py` (the López-Otín taxonomy). A citation layer *over* measured data,
not a source of it.

### 5.2 The one synthetic table

`clinical_synthetic_csd` — 720 subjects with critical slowing down planted by
construction (a latent factor whose variance grows with age). **It is not a
stand-in for missing data.** It is the only way to ask whether the estimator
recovers an effect that is *known* to be there, which no real dataset can answer,
because in real data you do not know the answer. It has its own dataset id, a
`SYNTHETIC` description prefix, and a test asserting it never merges with NHANES.

### 5.3 A design decision worth preserving

The three PhenoAge-only markers (MCV, WBC, ALP) went into a **separate dataset
id** rather than three more columns on `clinical_nhanes_slice`. The resilience
service infers its biomarker list by *excluding* known non-biomarker columns, so
widening the clinical table in place would have pulled three markers into the
health state and silently moved the published CSD numbers.

**Never widen `nhanes.MARKER_MAP`.** Extra markers go in
`PHENOAGE_ONLY_MARKER_MAP` and get their own dataset id.

---

## 6. Every bug found

Each was exposed by running real data. Each would have passed a stub.

| # | Bug | Consequence |
|---|---|---|
| 1 | **`idmap/mygene.py` wrote other species' identifiers into the canonical id space.** mygene's `ensembl` field is sometimes a list of homology cross-references; taking `[0]` gave mouse genes `ENSFALG…` (collared flycatcher) and `ENSMPUG…` (ferret) ids | Those became the join key for a mouse gene's signatures. Affected the **live request path**, not just the ETL |
| 2 | **Symbol-built gene tables silently fail to join.** mygene resolves the symbol `CDKN1A` to `ENSG00000205643` — a copy on an alternative MHC haplotype scaffold — while Entrez `1026` gives the primary `ENSG00000124762` | **p21 had twelve real contrasts and the API reported none.** `build_idmap` now resolves by Entrez wherever HAGR supplies one |
| 3 | **Truncating a sorted alias list drops the tail.** Capping aliases at 12 lost `p16` from human CDKN2A | The offline resolver answered "p16" with the mouse ortholog |
| 4 | **The download cache collided on URL basename.** HAGR serves both DrugAge and GenDR from a path ending `/dataset.zip` | They shared one cache entry. The checksum check turned it into a loud failure rather than a silent wrong answer — which is what it is for — but `cache_path` is now prefixed with the manifest key |
| 5 | **The data version was not reproducible.** It hashed Parquet filenames and byte sizes, and DuckDB writes a partitioned COPY in parallel | Two builds of identical data could report different versions. Now a content digest |
| 6 | **236 clock failures that were one data-prep bug.** Betas handed over in the stored CpG-by-sample orientation, and `read_series` deleting every CpG with any gap — including three Horvath coefficients | Every clock reported "requires features not present in the input" — a data-preparation fault wearing the costume of an incompatible input |
| 7 | **`/v1/intervention/{name}` returned `matches[0]`** | Asking about rapamycin returned the *C. elegans* result (+19.8%) instead of mouse (+13.0%) — a number about nematodes presented as the answer about a drug |
| 8 | **A test that passed only when you had not run `make data`** | It relied on the git-ignored full NHANES table being absent |

Pre-existing defects also fixed earlier in this work: negative AR(1) coefficients
were being clipped into the decaying regime and reported as *high resilience*;
CSD accepted a positive slope alone (fires on noise ~half the time per indicator)
and now requires bootstrap CIs excluding zero; `resilience/service.py` had a live
`NameError` on every `control_energy` call because nothing tested it.

---

## 7. The divergence with `origin/main`, and how it was resolved

While this work was in progress, a separate session merged **PR #1** into
`origin/main` — an independent revamp toward the same goal, "no fabricated data",
by the opposite route. Its environment had biomedical APIs blocked (HTTP 403), so
it *deleted* the fabricated signatures and replaced them with a hand-encoded,
hand-cited knowledge module. This branch ingested and pinned the real sources.

Integration was chosen over discarding either side.

**Carried forward from `main`:** `clocks/phenoage.py`; the ~1,000-line Streamlit
redesign and `.streamlit/config.toml`; `knowledge/references.py` and the hallmark
vocabulary; `clocks/service.py`'s `mortality_risk_10yr` (its bare
`except Exception` narrowed to `ClockInputError`); `CanonicalGene.ortholog_group`.

**Deliberately not carried:** `knowledge/aging_knowledge.py`, whose hand-written
per-gene evidence asserts what the ingested panel now measures — and disagrees
with for CDKN2A. Shipping both would mean two answers to one question, one of
them unfalsifiable. Also `main`'s `api/app.py`, which drops
`/gene/{id}/signature`, `/geneset/signature` and `/studies`; and
`example_cohort_simulated.csv`, superseded by real NHANES.

> ⚠️ **A trap worth knowing.** Git merged `api/app.py`, `clocks/service.py`,
> `idmap/resolver.py`, `conftest.py` and `test_clocks.py` **without conflict**,
> because this branch had not touched the same lines — so it silently took
> `main`'s versions wholesale, including the removal of three endpoints.
> **"No conflict" is not "no decision."** Every file the merge touched was
> reviewed individually.

---

## 8. Quality state

```
ruff        clean
black       clean
mypy        clean — 53 source files
pytest      321 passed, 1 skipped (12 live deselected)
coverage    78%
```

Live suite: `GEROQUERY_ALLOW_NETWORK=1 pytest -m live` → 7 pass, 5 skipped on
Python 3.14 (they need the `clocks` extra).

### Coverage by area

| Well covered | Not covered, and why |
|---|---|
| `sources/geo.py` 97%, `store.py` 95%, `resilience/*` 89–100%, `methylation.py` 94%, `differential.py` 94%, `models.py` 99%, `manifest.py` 100% | `ui/streamlit_app.py` **0%** — presentation, manually exercised |
| `idmap/resolver.py` 95%, `mygene.py` 88%, `bulk.py` 83% | `etl/build_idmap.py`, `build_fixtures.py`, `fetch_artifacts.py` **0%** — offline scripts run by `make`, not on any request path |
| `hagr.py` 86%, `build_signatures.py` 71%, `local_fixture.py` 97% | `nhanes.py` 70% — the uncovered part is the network fetch path, exercised by the live suite |

### 17 test files

`test_api` · `test_clocks` · `test_clocks_library` · `test_crosscutting` ·
`test_data_layer` · `test_differential` · `test_etl_signatures` · `test_geo` ·
`test_hagr` · `test_harmonize` · `test_idmap` · `test_methylation` ·
`test_network_paths` · `test_resilience` · `test_sources` · `test_store` ·
`conftest`

Notable tests, because they encode judgement rather than mechanics:
- `test_p16_does_not_replicate_across_the_real_panel` **and**
  `test_p21_does_replicate_across_the_real_panel` — the pair is what stops a null
  being dismissed as a broken estimator.
- `test_matrix_hedges_g_matches_scalar_implementation` — pins the vectorized copy
  to the reference implementation.
- `test_ensembl_cross_reference_list_does_not_leak_another_species_id` — bug #1.
- `test_all_missing_cpgs_are_dropped_but_partial_ones_are_imputed` — bug #6.
- `test_live_geo_soft_urls_still_serve_the_pinned_artifacts` — catches upstream
  drift a checksum cannot (a moved URL).
- `test_curated_genes_are_browsable_without_a_network` — guards the regression
  where the offline resolver silently returned 22 of 1,880 genes.

---

## 9. Running it

```bash
pip install -e ".[dev,ui]"          # add ",clocks" on Python 3.10-3.12

make data          # every upstream, verified          (~679 MB, ~30 min cold)
make signatures    # just the GEO + HAGR half
make idmap         # regenerate the derived identifier tables
make data-offline  # from the committed real samples   (no network)

make api           # http://localhost:8000/docs
make dashboard     # Streamlit
make ci            # ruff + mypy + 321 tests

GEROQUERY_ALLOW_NETWORK=1 pytest -m live     # upstream drift check
```

Docker: `docker compose up` (API :8000, dashboard :8501).
React: `cd frontend && npm install && npm run dev`.

### Sanity numbers for a rebuild

31 GEO DataSets / 32 contrasts / 27 Series · 485,905 signature rows over 46,091
genes · 40,585-row committed slice · 2,720 curated assertions over 1,880 genes ·
1,340 interventions · NHANES 4,895 (resilience) and 4,894 (PhenoAge) · 44 pinned
artifacts · rapamycin in *Mus musculus* = +13.0% · Horvath vs published = MAE 1.39 y.

---

## 10. Limitations

### Gene signatures
1. **Microarray-era only.** GEO stopped curating DataSets ~2016; no RNA-seq.
2. **Small groups — 3–15 per study; median 95% CI width 0.84.** This panel detects
   large, consistent, cross-tissue effects and cannot rule out moderate ones.
   **A null means "not detected by this instrument", never "absent".**
3. **Lopsided tissue coverage** — half the human contrasts are skeletal muscle;
   one blood contrast (n=9).
4. **`omic_layer` is `transcriptome` everywhere.** The multi-omic schema is real;
   the multi-omic *data* is one layer deep.
5. **Cross-study processing heterogeneity uncorrected** — 12 platforms, three
   decades of normalization practice.
6. **32 contrasts from 27 Series.** Paired array halves share subjects; the pool
   assumes independence. `series_id` makes this checkable but the pooling does not
   account for it.
7. **Sex is `unspecified` for 24 of 32 contrasts** — GEO never recorded it.
   Deliberately distinct from `both`.
8. **~22% of probe rows dropped** at identifier resolution, concentrated on the
   oldest platforms.

### Methylation and clocks
9. **Two blood 450K series only.** Nothing established about other tissues, EPIC
   arrays, or sequencing-based methylation.
10. **n = 106 and n = 40.** Treat the third decimal as noise.
11. **GSE30870 is swept to 200 of 236**; GSE64495 is complete. The outstanding 36
    are the slowest neural models. Nothing in the results depends on them.
12. **`GPAge*` cannot run** — see §2, landmine 5.
13. **Normalization is not the authors'.** The 1.39-year gap in the Horvath
    comparison is mostly this.
14. **No mortality outcome.** Neither series has follow-up, so mortality clocks
    can be checked for direction and internal consistency but not for what they
    were built to predict.

### Clinical / resilience
15. **NHANES is cross-sectional.** Age strata proxy for within-individual
    trajectories; no relaxation time is observed — which is precisely what would
    distinguish critical slowing down from accumulating heterogeneity.
16. **Age is topcoded at 80**, compressing the stratum where CSD should be
    strongest. 354 subjects sit at the cap.
17. **Survey weights carried but not applied.** Estimates are not nationally
    representative.
18. **No fasting, medication, or comorbidity adjustment.**

### Curated knowledge
19. **Only mammalian rows loaded** — 2,276 of 5,029 HAGR assertions are yeast,
    worm, or fly and are dropped. Ortholog mapping would recover them.
20. **DrugAge gene links are absent because DrugAge has none.**
21. **Gene→hallmark annotation was dropped** with `aging_knowledge.py`.
    `gene_report` keeps an empty `hallmarks` key so the shape is stable.

### Operational
22. **HAGR pins go stale by design** — an upstream update fails the build loudly.
    Intended, but it needs periodic maintenance.
23. **GTEx open cannot be age-stratified.** The API accepts
    `attributeSubset=ageBracket` and returns HTTP 200 with one undivided group;
    donor age is dbGaP-controlled. A live test pins this.
24. **The React frontend has not been repointed** at the newest fields.

---

## 11. Next steps, prioritized

### 11.1 Correlate age acceleration against resilience — *the reason this project exists*

Now unblocked and the highest-value remaining work. Clock age minus chronological
age is computable on real samples; resilience metrics are computable on real
NHANES. **This cross-layer claim is still the one headline claim with no evidence
behind it in either direction.**

The obstacle is that they live in different cohorts — GEO methylation samples and
NHANES participants are different people. In increasing order of effort: find a
GEO series carrying both methylation and clinical chemistry; pursue NHANES DNA
methylation subsamples; or apply for a cohort with both.

### 11.2 recount3 or ARCHS4

The biggest single improvement to the signature panel: removes limitation 1
(microarray-only) and 5 (uncorrected processing heterogeneity) in one move, and
brings modern blood coverage. **ARCHS4 is the lighter path** — one HDF5, no
R/Bioconductor bridge.

### 11.3 Longitudinal resilience

`resilience/recovery.py` (AR(1) relaxation) is implemented, tested, and has
**never run on real longitudinal data**. It has the most theoretical value and the
least empirical exercise. NHANES III linked mortality files are free;
ELSA/SHARE/HRS need registration; UK Biobank is ideal and the federate-only
adapter exists. This resolves limitation 15.

### 11.4 Finish the GSE30870 sweep

36 clocks short. The runner writes incrementally and resumes from
`clock_results.csv`. Letting it run, not work.

### 11.5 Smaller

- **Nested random effects by series** — handles limitation 6 properly.
- **Expose `gene_report`, `list_curated_genes`, `references` over HTTP** — the
  service methods exist; only Streamlit calls them.
- **Re-add a sourced gene→hallmark mapping.**
- **Apply survey weights** — makes the CSD result a population claim.
- **Ortholog mapping for HAGR** — recovers 2,276 assertions.
- **`?library=` filter on `GET /v1/clocks`** — 240 entries is a lot unfiltered.
- **Pool multiple NHANES cycles** to get past the age-80 topcode.

---

## 12. Gotchas

- **`GEROQUERY_ALLOW_NETWORK` is off by default.** Set `=1` for anything touching
  an upstream.
- **Two interpreters, on purpose.** Python 3.14 runs the core suite; `.venv312`
  runs the clock libraries. Both staying green is what proves the optional tier is
  genuinely optional.
- **Do not install GPy.** See §2, landmine 5.
- **`signatures_full.csv` is gitignored; `signatures_curated.csv` is committed.**
  `LocalSignatureSource` prefers the full one when present, so a developer who has
  run `make data` sees different counts from CI. **Tests pin `prefer_full=False`
  for exactly this reason — do the same in any new test.**
- **Never widen `nhanes.MARKER_MAP`.** See §5.3.
- **Do not hand-edit the GEO panel.** See §4.2.
- **`get_registry()` memoizes** a module-level registry. Tests faking a clock
  library must reset `clocks.registry._REGISTRY`, or they pass alone and fail in a
  full run.
- **Live tests are excluded by default** (`-m 'not live'`). They are the only
  thing that catches upstream drift: a moved URL, a renamed file inside a zip, a
  changed GEO query population.
- **biolearn's gallery init takes ~2 minutes.** Not a hang.
- **`fetch.cache_path` is prefixed with the manifest key.** Do not "simplify" that
  away — see bug #4.
- **pyaging imputes missing features from reference values.** The wrapper refuses
  above 50% imputation rather than return a confident-looking artifact of
  pyaging's reference cohort.
- **The NHANES URL pattern is a trap.** `wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT`
  returns HTML, not XPORT. The working pattern is
  `wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt`. A test pins it.

---

## 13. Commit history for this work

```
ef1dfd9  docs: correct the methylation sweep coverage
bb9a6da  docs: session handoff for the integration and methylation work
cd2ef49  merge: reconcile with the PR #1 revamp on main
10189b2  docs: methylation clock validation results
171b840  test: methylation parsing and the PhenoAge clock
6c91cf0  feat: integrate the main-branch revamp; ingest real methylation;
         source the last hand-written tables
4757514  docs: session handoff for the real-data ingestion work
e0c3301  feat: replace every synthetic data layer with checksum-pinned real sources
```

All authored by Sophie Zhang. Branch `real-data-ingestion` is retained at
`e0c3301` as the pre-merge state.

---

## 14. Document map

| Document | What it is |
|---|---|
| **`HANDOFF.md`** | This document. Authoritative |
| `OVERVIEW.md` | Sources, architecture, use cases, limitations — the reference table |
| `DATA_SOURCES.md` | Cache-vs-federate contract and licence table |
| `RESULTS_NHANES_CSD.md` | The clinical resilience result (partial null) |
| `RESULTS_GEO_SIGNATURES.md` | The gene signature result (p16 does not replicate) |
| `RESULTS_METHYLATION_CLOCKS.md` | The clock validation result |
| `HANDOFF_2026-08-03.md`, `HANDOFF_2026-08-04.md` | Session records |
| `HANDOFF_TIER0_REVIEW.md`, `HANDOFF_TIER1_SIGNATURES.md`, `HANDOFF_PR1_REVAMP.md` | Historical |
| `WORKING_NOTES.md` | Scratch: verified URLs, variable maps, rebuild numbers |
