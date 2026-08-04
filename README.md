# 🧬 GeroQuery

**An open-source, multi-omic + clinical aging-data aggregator with a dynamical-systems resilience layer no other aging portal exposes.**

Type a gene, get its harmonized, multi-omic, cross-species, clock-aware, resilience-aware relationship to aging — in one place, behind one API, reproducibly.

[![CI](https://github.com/praneel/geroquery/actions/workflows/ci.yml/badge.svg)](https://github.com/praneel/geroquery/actions) · License: MIT · Python ≥3.10

---

## The problem

Researchers and longevity companies spend months cleaning and harmonizing public aging datasets before they can use them. Transcriptomics lives in GEO/recount3/ARCHS4, methylation in scattered GEO series, proteomics in PRIDE and supplementary tables, single-cell in CELLxGENE Census, clinical aging in NHANES, curated knowledge in HAGR/OpenGenes, interventions in DrugAge and the NIA ITP — each with different gene identifiers, tissue labels, age representations, and access methods.

Existing databases are each partial: HAGR/OpenGenes are curated catalogs, not live expression engines; the clock ecosystem (pyaging, biolearn) computes biological age but doesn't connect clocks back to datasets or interventions. **Nothing unifies transcriptome + methylome + proteome + clinical biological-age + intervention data behind one clean, harmonized, cross-species API — and nothing exposes the dynamical-systems "resilience / criticality" view of aging.**

## The solution

GeroQuery is three layers over one harmonized core:

1. **Harmonization + federation backend** — a FastAPI service over a DuckDB/Parquet analytical store. Gene IDs are harmonized via a canonical Ensembl space (mygene.info in production), tissues via UBERON, and species ages normalized to **fractional lifespan** using AnAge maximum-lifespan values, so an "old mouse" and an "old human" are comparable. Open, redistributable data is cached; large/controlled sources are federated live and **never re-hosted**.
2. **Clock + biological-age service** — wraps real `biolearn` clocks (with transparent reference clocks bundled for offline use), always surfacing *what each clock predicts* (chronological age vs mortality vs pace), because those decouple — and refusing to relabel a non-age predictor as an age clock.
3. **Resilience / criticality module (the differentiator)** — critical-slowing-down indicators (rising variance and cross-correlation across age strata), a DOSI-style recovery-rate metric, and an optional network-control-energy view.

A Streamlit dashboard and a React showcase page both call the same service.

## The process

The build is organized around **eight deep modules** (Ousterhout's sense) — each hides substantial complexity behind a small, stable interface and is tested in isolation:

```
ui  ->  api  ->  ( store, clocks, resilience, harmonize )  ->  ( sources, idmap )
                       M4     M5        M6         M3            M2       M1
```

Lower modules never import higher ones, so every layer is testable with the layer below mocked.

| Module | Responsibility | Correctness is… |
|---|---|---|
| **M1 `idmap`** | resolve any gene ID → canonical; UBERON tissues; AnAge fractional age | scientific |
| **M2 `sources`** | uniform source adapters; `capabilities()` + `license()`; cache-vs-federate gate | plumbing |
| **M3 `harmonize`** | Hedges' g effect sizes; per-gene contrasts + BH FDR (`differential.py`); DerSimonian-Laird random-effects meta-analysis; batch correction | scientific |
| **M4 `store`** | partitioned Parquet (DuckDB) + SQLite metadata + versioning | plumbing |
| **M5 `clocks`** | clock registry with predicted-outcome metadata; validation; age acceleration | scientific |
| **M6 `resilience`** | CSD indicators, recovery rate, control energy, documented fallback | scientific |
| **M7 `api`** | FastAPI `/v1`, error envelope, pagination, `format=json\|csv\|parquet` | plumbing |
| **M8 `ui`** | Streamlit dashboard + React showcase | manual QA |

## What is real data and what is not

This distinction is the first thing you should read, not a footnote.

| Layer | Status |
|---|---|
| **Clinical / resilience** (`clinical_nhanes_slice`) | **Real.** NHANES 2017–2018, 4,895 complete cases aged 20–80, fetched from CDC and verified against pinned SHA-256 checksums. A 600-row real sample is committed so tests run offline. |
| **Aging clocks** | **Real, and now validated on real data.** 236 published clocks wrapped from `biolearn` (63) and `pyaging` (173), plus the real Levine **PhenoAge** coefficients implemented directly. Run against two checksum-pinned GEO 450K blood series: **436 clock-dataset runs; Horvath reproduces the authors' own published per-sample ages to r=0.998, MAE 1.39 y**, and predicts **0.23 years for cord blood**. See [`docs/RESULTS_METHYLATION_CLOCKS.md`](docs/RESULTS_METHYLATION_CLOCKS.md). |
| **DNA methylation** | **Real.** GSE64495 (450K, n=106 controls, ages 2.3–73.7) and GSE30870 (newborns vs nonagenarians), checksum-pinned. |
| **Gene identifier resolution** | **Real.** Bundled canonical table, with live batched mygene.info resolution and an on-disk cache behind it. |
| **Tissue expression** (`gtex-open`) | **Real.** GTEx Portal API v2, queried live: median TPM per tissue with UBERON ids. |
| **Gene aging signatures** | **Real.** 31 checksum-pinned [GEO DataSets](docs/RESULTS_GEO_SIGNATURES.md) -> 32 young-vs-old contrasts from 27 independent GEO Series -> **485,905 Hedges' *g* estimates over 46,091 genes**, 10 tissues, human + mouse. A 40,585-row slice (curated aging genes plus their orthologs, both species) is committed so tests run offline. |
| **Curated knowledge & interventions** | **Real.** The five HAGR databases, checksum-pinned: GenAge (307 human + 2,205 model-organism genes), CellAge (949), LongevityMap (1,325 assertions including nulls), DrugAge (3,423 experiments -> 1,334 compound x organism records), GenDR (214). |
| **`clinical_synthetic_csd`** | **Synthetic, deliberately.** Critical slowing down is planted in it. Exists to check that estimators recover a known effect. Registered as a separate dataset with a `SYNTHETIC` description so it cannot be confused with the NHANES one. |

The honest one-line summary: **every layer is now backed by a checksum-pinned upstream.** The one synthetic table left is `clinical_synthetic_csd`, which exists on purpose: it has critical slowing down planted in it so the estimator can be checked against a known answer.

Full source table, architecture, use cases, and limitations: [`docs/OVERVIEW.md`](docs/OVERVIEW.md).

## Headline result

Running the resilience module on real NHANES data, critical slowing down predicts
two early-warning signals. **One replicates and one does not:**

| Indicator | Result across 20 analytic configurations |
|---|---|
| Health-state **variance** rises with age | **Supported in 20/20** (bootstrap 95% CI excludes zero in all) |
| Marker **cross-correlation** rises with age | **Supported in 0/20** |

`resilience_declines` returns **`False`** on this dataset, because it requires
both. The synthetic fixture returns `True` — which is exactly what it was built
to do. Full method, sensitivity sweep, and caveats:
[`docs/RESULTS_NHANES_CSD.md`](docs/RESULTS_NHANES_CSD.md).

### And a second one, from the real gene signatures

Pooling 32 real GEO contrasts, **CDKN2A/p16 — the most-cited transcriptional
marker of aging — does not replicate**: *g* = +0.07, 95% CI [-0.20, +0.35], 14
contrasts, I² = 14%. The studies agree with each other that there is nothing to
see. The synthetic slice this replaced had p16 planted at *g* = +1.20, and the
old API test asserted every direction was "up".

The same pipeline does find p21/CDKN1A (*g* = +1.07), IGF1 (-0.53), the
metallothioneins, and a coherent mitochondrial/OXPHOS decline — so the estimator
works and the p16 null is a result, not a bug. Method, top hits, and the ten
limitations that qualify all of it:
[`docs/RESULTS_GEO_SIGNATURES.md`](docs/RESULTS_GEO_SIGNATURES.md).

## Quickstart

```bash
pip install -e ".[dev,ui]"        # add ",clocks" for 236 real published clocks (needs Python 3.10-3.12)
make data                         # fetch + verify every upstream, build every table (~320 MB, ~25 min)
make api                          # http://localhost:8000/docs
make dashboard                    # Streamlit UI
make test                         # 295 tests (+12 live)
```

No network? `make data-offline` builds everything from the committed real samples
— the 600-row NHANES subset and the 40,585-row curated signature slice. Both are
real measurements, not generated stand-ins; they are *samples*, so a number
computed on them is not the full-panel number, and the adapters return which mode
produced it rather than only logging it. Network access is off by default (`GEROQUERY_ALLOW_NETWORK`), so
tests and CI are hermetic by construction rather than by convention.

Docker:

```bash
docker compose up            # API on :8000, dashboard on :8501
```

React showcase:

```bash
cd frontend && npm install && npm run dev   # proxies /v1 to the API on :8000
```

## API (`/v1`)

| Endpoint | What it returns |
|---|---|
| `GET /v1/gene/{id}/signature` | multi-omic aging signature (filters: species, tissue, sex, omic_layer; `format=` switch) |
| `GET /v1/gene/{id}/card` | assembled multi-omic + curated-knowledge + intervention card |
| `POST /v1/geneset/signature` | aggregate signature for a gene set / pathway |
| `GET /v1/clocks` | clocks with predicted-outcome + training-population metadata |
| `POST /v1/clock/apply` · `POST /v1/clock/compare` | apply/compare clocks on a dataset ID or uploaded matrix |
| `GET /v1/intervention/{name}` | intervention record + linked signatures + lifespan data |
| `POST /v1/resilience/csd` · `POST /v1/resilience/recovery` | resilience metrics |
| `GET /v1/studies` · `GET /v1/sources` · `GET /v1/datasets` · `GET /v1/version` | provenance & versioning |

Every error uses one envelope: `{"error": {"code", "message", "detail"}}`.

## Validation against known biology

The test suite asserts *behavior against ground truth*, not internals:

- **The real signature panel** — the pipeline recovers p21/CDKN1A, IGF1, and the mitochondrial decline with intervals excluding zero, and returns an honest null for p16/CDKN2A. Both are pinned as tests: one asserts the effect is found, the other asserts the CI still spans zero. The pair matters because it is what stops a null being dismissed as a broken estimator.
- **The vectorized effect-size path is pinned to the scalar one** — `differential.hedges_g_matrix` is asserted equal to `effect_size.hedges_g` element by element, so the fast copy cannot drift from the reference.
- **Contrast construction** — the GEO rules are tested against the real level vocabulary (`20 - 39 y`, `E16.5`, `P7`, `aged (24 months)`, `ad libitum`), including that a dataset whose confounder has no recognisable control arm is skipped with a stated reason rather than guessed at.
- **Meta-analysis** — recovers a planted pooled effect and detects planted heterogeneity (I², τ²) on synthetic studies.
- **Clocks** — a reference clock recovers known ages to <1.5 yr MAE; a missing CpG/gene yields a precise, actionable error; a mortality clock correctly reports *no* age-acceleration. Wrapped models keep their declared output, so the 19 biolearn models that predict BMI, smoking, cholesterol or cancer status are never labelled age clocks.
- **Cross-library agreement** — fed byte-identical input, biolearn's `Horvathv1` and pyaging's `horvath2013` agree to within 0.05 years. Two independent implementations, two different matrix orientations, one answer: the strongest check available that the wrappers are right and not merely non-crashing.
- **Resilience** — the estimator finds critical slowing down in the fixture that has it planted, and **declines to claim it** on real NHANES where only one of the two signals is supported. That discrimination is the test: an earlier version returned `True` for both, because a positive slope alone was enough.
- **Data integrity** — a truncated or edited cache entry fails its checksum and is refused; an interrupted download never leaves a valid-looking file behind; offline runs use the committed real sample rather than silently reaching out.
- **Reproducibility / provenance / licence** — identical results at a fixed data version; every signature carries source + study + version; controlled sources (UK Biobank, dbGaP, GTEx-protected) are provably federate-only and refuse caching.

## Data sources

See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) for the full cache-vs-federate + licence table. In short: open, redistributable data (GEO summaries, GTEx-open, NHANES, HAGR, OpenGenes, ComputAgeBench) is cached as a harmonized slice; large/controlled sources (recount3 long tail, CELLxGENE Census, UK Biobank, dbGaP, GTEx-protected) are **federated or linked only, never re-hosted** — enforced in code by each adapter's `license()` and a test.

## Limitations (read these)

**On the real data:**

- NHANES age is **topcoded at 80**, compressing exactly the oldest stratum — where any critical slowing down should be strongest. The effect is measured on a truncated age range.
- NHANES is **cross-sectional**. Age strata are a *proxy* for within-individual trajectories (`fallback_used=True` in every response); no relaxation time is observed. Distinguishing critical slowing down from ordinary accumulating heterogeneity needs longitudinal data.
- **Survey weights are carried but not applied.** Results are unweighted and therefore not nationally representative.
- No fasting-status, medication, or comorbidity adjustment.

**On the gene signatures:**

- **Microarray-era only.** GEO stopped curating DataSets around 2016, so there is no RNA-seq in the panel. recount3 or ARCHS4 is the highest-value next source.
- **Small groups — 3-15 samples per study, median 95% CI width 0.84.** This panel detects large, consistent, cross-tissue effects. It cannot rule out moderate ones. A null means "not detected by this instrument", never "absent".
- **Tissue coverage is lopsided** — half the human contrasts are skeletal muscle, there is one blood contrast (n=9), and `omic_layer` reads `transcriptome` for every row. The multi-omic schema is real; the multi-omic *data* is one layer deep.
- **Cross-study processing heterogeneity is uncorrected** (12 platforms, three decades of normalization practice), and **32 contrasts come from 27 Series** — GEO splits some experiments across two array halves, so a few studies share subjects. `series_id` is carried on every study row so this is checkable; the pooling does not yet account for it.
- **Only mammalian HAGR rows are loaded.** 2,276 of 5,029 assertions are yeast, worm, or fly; they are parsed, counted, and dropped, because GeroQuery's canonical gene space is human + mouse.
- **DrugAge records no gene targets, so none are asserted.** The previous hand-written table linked rapamycin -> MTOR; that edge was never in the data.

**On what is still synthetic:**

- `clinical_synthetic_csd` — 720 subjects with critical slowing down planted by construction. A method-validation fixture, deliberately kept, registered as its own dataset id with a `SYNTHETIC` prefix. It is the only generated data in the repository.
- **Reference clocks** (`*_demo`) are transparent linear models for testability, **not** validated clinical instruments.

**On coverage:**

- **The clock libraries need Python 3.10–3.12 and pandas 2.** `ecos` (a biolearn transitive dep) has no 3.14 wheel, and biolearn's DunedinPACE is incompatible with pandas 3's mandatory copy-on-write. biolearn also under-declares `torch` and `seaborn`. All pinned in the `clocks` extra; see `docs/HANDOFF.md` §1.
- **pyaging clocks register metadata-only.** Their model artifacts download lazily from Hugging Face, so `required_features` is empty until first use and the first prediction per clock hits the network.
- **GTEx cannot be stratified by age here.** The open v2 API accepts `attributeSubset=ageBracket` but returns a single undivided group — donor age is in the dbGaP-controlled tier. GTEx contributes tissue context, not an aging signature. A live test pins this so we find out if it ever changes.
- No wet-lab or clinical claims. Metrics are research indicators with stated assumptions.

## Repository layout

```
geroquery/
  idmap/         M1  + bundled reference data (genes, UBERON, AnAge)
  sources/       M2  + pinned data manifest, checksum fetch layer, real GEO
                     DataSets adapter, HAGR adapters, NHANES adapter, GTEx v2
  harmonize/     M3
  store/         M4
  clocks/        M5
  resilience/    M6
  api/           M7
  ui/            M8 (Streamlit)
  clocks/library  real biolearn clock wrappers
  clocks/pyaging_clocks  real pyaging clock wrappers
  etl/           real-data ETL (build_data, build_signatures) + the one
                 synthetic fixture builder
frontend/        M8 (React/Vite showcase)
tests/           295 tests across all modules (+12 live), 84% coverage
notebooks/       end-to-end example
```

## License

Code: **MIT**. Bundled data: redistributable slices of real upstream releases; every upstream keeps its own licence, attributed per adapter and recorded in `sources/manifest.py`. HAGR data is free for non-commercial use with attribution. See [`LICENSE`](LICENSE).
