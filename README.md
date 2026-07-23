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
2. **Clock + biological-age service** — wraps `pyaging` / `biolearn` (with transparent reference clocks bundled for offline use), always surfacing *what each clock predicts* (chronological age vs mortality vs pace), because those decouple.
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
| **M3 `harmonize`** | Hedges' g effect sizes; DerSimonian-Laird random-effects meta-analysis; batch correction | scientific |
| **M4 `store`** | partitioned Parquet (DuckDB) + SQLite metadata + versioning | plumbing |
| **M5 `clocks`** | clock registry with predicted-outcome metadata; validation; age acceleration | scientific |
| **M6 `resilience`** | CSD indicators, recovery rate, control energy, documented fallback | scientific |
| **M7 `api`** | FastAPI `/v1`, error envelope, pagination, `format=json\|csv\|parquet` | plumbing |
| **M8 `ui`** | Streamlit dashboard + React showcase | manual QA |

## Quickstart

```bash
pip install -e ".[dev,ui]"        # or: pip install -r requirements.txt
python -m geroquery.etl.build_fixtures   # (re)generate the bundled data slice
make api                          # http://localhost:8000/docs
make dashboard                    # Streamlit UI
make test                         # 77 tests, ~1s
```

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

- **Cross-species conservation** — CDKN2A/p16 pools **up with age** in both human (g≈1.22) and mouse (g≈1.11) transcriptome; LMNB1 pools **down**.
- **Meta-analysis** — recovers a planted pooled effect and detects planted heterogeneity (I², τ²) on synthetic studies.
- **Clocks** — a reference clock recovers known ages to <1.5 yr MAE; a missing CpG/gene yields a precise, actionable error; a mortality clock correctly reports *no* age-acceleration.
- **Resilience** — on data engineered toward a tipping point, variance **and** cross-correlation rise with age (critical slowing down); recovery rate is higher for a resilient (low-autocorrelation) system; control energy grows with target distance.
- **Reproducibility / provenance / licence** — identical results at a fixed data version; every signature carries source + study + version; controlled sources (UK Biobank, dbGaP, GTEx-protected) are provably federate-only and refuse caching.

## Data sources

See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) for the full cache-vs-federate + licence table. In short: open, redistributable data (GEO summaries, GTEx-open, NHANES, HAGR, OpenGenes, ComputAgeBench) is cached as a harmonized slice; large/controlled sources (recount3 long tail, CELLxGENE Census, UK Biobank, dbGaP, GTEx-protected) are **federated or linked only, never re-hosted** — enforced in code by each adapter's `license()` and a test.

## Limitations (read these)

- The **bundled data is a curated demonstration slice.** Effect *directions* encode established aging biology; magnitudes and per-study scatter are synthetic and labelled as such. The production ETL swaps the generators for live GEO/GTEx/recount3 harvests behind the same schema — the code, interfaces, and tests are the artifact, not the demo numbers.
- **Reference clocks** (`*_demo`) are transparent linear models for testability, **not** validated clinical instruments. Real clocks arrive via `pyaging`/`biolearn` where installed.
- The resilience module's cross-sectional CSD is an **age-stratified proxy** (clearly labelled `fallback_used=True`) unless dense longitudinal data is supplied; it does not observe within-individual recovery.
- Gene/ontology IDs in the bundled slice are internally self-consistent; re-verify against the current mygene.info / UBERON / AnAge releases before publication.
- No wet-lab or clinical claims. Metrics are research indicators with stated assumptions.

## Repository layout

```
geroquery/
  idmap/         M1  + bundled reference data (genes, UBERON, AnAge)
  sources/       M2  + bundled harmonized CSV slice
  harmonize/     M3
  store/         M4
  clocks/        M5
  resilience/    M6
  api/           M7
  ui/            M8 (Streamlit)
  etl/           reproducible fixture builder
frontend/        M8 (React/Vite showcase)
tests/           77 tests across all modules
notebooks/       end-to-end example
```

## License

Code: **MIT**. Bundled data: redistributable demonstration slice; upstream sources retain their own licences (attributed per adapter). See [`LICENSE`](LICENSE).
