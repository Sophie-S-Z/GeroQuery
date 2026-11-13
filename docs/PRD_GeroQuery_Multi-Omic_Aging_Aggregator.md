# PRD: GeroQuery — Open-Source Multi-Omic & Clinical Aging Data Aggregator

**Version:** 1.0
**Status:** Draft for build
**Owner:** Praneel
**Type:** Personal portfolio + open-source research tool
**Working name:** GeroQuery (placeholder; alternatives in §15)
**Primary goal:** Portfolio/learning piece that demonstrates production data-engineering + bioinformatics + domain understanding to longevity companies.
**Secondary goal:** A genuinely useful open tool that fills a real field gap (multi-omic + clinical aging querying with a dynamical-systems resilience layer no competitor exposes).

---

## 0. How to read this document

This PRD follows the standard template (Problem → Solution → User Stories → Implementation Decisions → Testing Decisions → Out of Scope → Further Notes) and extends it with the sections a greenfield build needs: a module decomposition, a phased step-by-step build plan, a data-source contract table, deployment, and a portfolio/launch plan.

A note on terminology used throughout:
- **Federate / federation:** query an external source live at request time and return/cache the result, rather than re-hosting its data. Used for large or licence-restricted sources.
- **Harmonize:** map heterogeneous inputs (gene IDs, tissue labels, age formats, units) onto one internal vocabulary so studies can be compared.
- **Deep module:** a unit that hides substantial complexity behind a small, stable interface (Ousterhout's term). The architecture is organized around these so each can be tested in isolation.
- **Aging signature:** the set of per-gene (or per-feature) age-associated changes (effect size, direction, significance) for a given tissue/species/omic layer.
- **Clock:** a model that maps molecular features to a predicted age or aging rate (e.g., Horvath, PhenoAge, DunedinPACE).
- **Resilience / criticality metrics:** indicators derived from dynamical-systems theory (autocorrelation, variance, recovery rate) that quantify how close a biological system is to a tipping point.

---

## 1. Problem Statement

I want to break into longevity biotech (BioAge, Retro, NewLimit, Gero, Insilico, Altos, Calico) as someone who can do the data engineering and computational biology these companies actually need. To do that I need a portfolio artifact that proves three things at once: I can build production-grade data infrastructure, I understand aging biology, and I can do something the field doesn't already have.

The underlying field problem the artifact should address is real and documented. Researchers and longevity companies spend months cleaning, formatting, and harmonizing public aging datasets before they can use them. The data lives in dozens of incompatible places: transcriptomics in GEO/recount3/ARCHS4, methylation in scattered GEO series, proteomics in PRIDE and supplementary tables, single-cell in CELLxGENE Census and Tabula atlases, clinical/phenotypic aging in NHANES and controlled cohorts, curated gene knowledge in HAGR and Open Genes, and intervention data in DrugAge and the NIA Interventions Testing Program. Each uses different gene identifiers, tissue labels, age representations, and access methods. A surveyed sample of ~400 longevity-biotech participants ranked "more public aging biology datasets" as the single most demanded solution in the field.

The existing aging databases are good but each is partial. HAGR and Open Genes are curated gene/species/drug catalogs, not live expression-query engines. Aging Atlas is multi-omic but hosted in a way that's slow for non-China users and is not a clean redistributable codebase. The clock ecosystem (pyaging, biolearn) computes biological age well but doesn't connect those clocks back to the underlying datasets or to interventions. Nothing unifies transcriptome + methylome + proteome + clinical biological-age + intervention data behind one clean API, harmonized across human and mouse, and nothing exposes the dynamical-systems "resilience / criticality" view of aging that my background is built for.

So the problem is twofold: (a) my career problem — I need a high-signal, end-to-end artifact; and (b) the field problem the artifact should genuinely solve — fragmented, unharmonized aging data with no single queryable multi-omic + clinical interface and no resilience layer.

---

## 2. Solution

Build **GeroQuery**: an open-source, well-documented system with three layers.

1. **A harmonization + federation backend (the real, reusable core).** A FastAPI service over a DuckDB/Parquet analytical store. It ingests a curated, well-harmonized slice of open aging data (transcriptomic aging signatures for human + mouse, a methylation slice, a plasma-proteome slice, curated gene knowledge, and intervention data), and federates live to the large/restricted sources (recount3, CELLxGENE Census, GTEx API) rather than re-hosting them. All gene IDs are harmonized via mygene.info; tissues via UBERON; species ages normalized by fractional lifespan using AnAge maximum-lifespan values.

2. **A clock + biological-age service.** Thin wrappers around `pyaging`, `biolearn`, and the `BioAge` method set, exposed behind a `/clock` endpoint so a user can apply an aging clock to a public dataset (or their own upload) and get predicted age / aging rate, with explicit metadata about *which outcome* each clock predicts (chronological age vs mortality), because those decouple.

3. **A resilience / criticality module (the differentiator).** Given longitudinal or age-stratified data, compute dynamical-systems biomarkers of aging: critical-slowing-down indicators (rising autocorrelation and variance across age strata), a DOSI-style recovery-rate / resilience metric, and optionally a network-control-energy view. This is the part no existing aging portal offers and the part that maps directly to my network-control-theory and dynamical-systems background.

On top of these sits a **Streamlit (or Gradio) dashboard** for fast, scannable exploration — search a gene, see its multi-omic aging signature across tissues and species, apply a clock, view linked interventions, and open the resilience view — plus an optional polished **React showcase page** that calls the same API to demonstrate frontend engineering.

The whole thing ships as a public GitHub repo with a problem→process→solution README, a Dockerized reproducible environment, CI, a live demo on Hugging Face Spaces / Streamlit Community Cloud, and a short bioRxiv preprint describing the harmonization method and the resilience module.

From the user's perspective (a longevity researcher or a hiring manager evaluating me), the solution is: *type a gene or upload a dataset, and immediately get a harmonized, multi-omic, cross-species, clock-aware, resilience-aware view of how it relates to aging — in one place, with one API, reproducibly.*

---

## 3. User Personas / Actors

- **Aging researcher / computational biologist (primary external user):** wants to query a gene or pathway and see its aging signature across omic layers, tissues, and species without writing ingestion code.
- **Longevity-company data scientist (target evaluator):** wants to assess whether I can build the infrastructure they rely on; will look at API design, harmonization rigor, reproducibility, and the resilience novelty.
- **Hiring manager / recruiter (target evaluator):** skims the README, the live demo, and the preprint to judge engineering + domain + communication.
- **Bioinformatics student / open-source contributor (secondary user):** wants to learn from or contribute to a clean, well-documented codebase.
- **Me, the developer/maintainer:** wants a modular codebase I can extend, test in isolation, and talk about fluently in interviews.

---

## 4. User Stories

### Gene / feature query (core)
1. As an aging researcher, I want to search a gene by symbol, Ensembl, or Entrez ID, so that I can find it regardless of which identifier I have.
2. As an aging researcher, I want the system to resolve any gene identifier to a canonical internal ID, so that results from different studies line up on the same gene.
3. As an aging researcher, I want to see a gene's age-associated expression change (effect size, direction, significance) across tissues, so that I understand where in the body it changes with age.
4. As an aging researcher, I want to see a gene's aging signature in both human and mouse side by side, so that I can judge cross-species conservation.
5. As an aging researcher, I want mouse and human ages normalized to fractional lifespan, so that "old mouse" and "old human" are comparable.
6. As an aging researcher, I want to see a gene's methylation-aging relationship where available, so that I can connect transcriptomic and epigenetic aging.
7. As an aging researcher, I want to see a gene's protein-level aging change (plasma proteome) where available, so that I can connect transcript to protein.
8. As an aging researcher, I want a single multi-omic summary card per gene, so that I don't have to assemble layers manually.
9. As an aging researcher, I want to query a pathway or gene set (e.g., a hallmark-of-aging set), so that I can see aggregate aging behavior, not just single genes.
10. As an aging researcher, I want to see which curated databases flag the gene as aging- or longevity-associated (GenAge, Open Genes, CellAge, LongevityMap), so that I get provenance and confidence.
11. As an aging researcher, I want each result annotated with its source study, sample size, and processing method, so that I can judge reliability.
12. As an aging researcher, I want to filter signatures by tissue, species, sex, and age range, so that I can scope to my question.

### Clock / biological-age service
13. As a data scientist, I want to apply a named aging clock to a public dataset by ID, so that I can get predicted ages without rebuilding the clock.
14. As a data scientist, I want to upload my own methylation or expression matrix and apply a clock, so that I can score my own data.
15. As a data scientist, I want each clock annotated with what it predicts (chronological age vs mortality vs pace of aging) and its training population, so that I don't misuse it.
16. As a data scientist, I want to compare several clocks on the same dataset, so that I can see where they agree or diverge.
17. As a data scientist, I want age-acceleration residuals (predicted minus chronological age), so that I can use the standard aging-research readout.
18. As a data scientist, I want the clock service to validate my input format and tell me exactly what's missing (e.g., absent CpG sites), so that I can fix my data quickly.
19. As a data scientist, I want a programmatic endpoint for the clock service, so that I can call it from my own pipelines.

### Resilience / criticality module (differentiator)
20. As a systems-minded researcher, I want to compute critical-slowing-down indicators (autocorrelation, variance) across age strata for a dataset, so that I can see loss of resilience with age.
21. As a systems-minded researcher, I want a DOSI-style recovery-rate / resilience metric where longitudinal data exists, so that I can quantify how fast a system returns to baseline after perturbation.
22. As a systems-minded researcher, I want the module to fall back to age-stratified variance/autocorrelation when longitudinal data is unavailable, so that I still get a usable resilience proxy.
23. As a systems-minded researcher, I want the resilience output tied to a specific dataset and method, so that the result is reproducible and citable.
24. As a systems-minded researcher, I want an optional network-control view (control energy to move an aging network toward a youthful state) on supported datasets, so that I can connect aging to control theory.
25. As a systems-minded researcher, I want clear documentation of the assumptions and limits of each resilience metric, so that I don't over-claim.

### Intervention linkage
26. As a researcher, I want to query an intervention (e.g., rapamycin, caloric restriction) and see which aging signatures it is associated with reversing, so that I can connect treatments to molecular effects.
27. As a researcher, I want to see lifespan-effect data for interventions from DrugAge and the Interventions Testing Program where available, so that I have hard outcome data.
28. As a researcher, I want to go from a gene to the interventions that modulate it, so that I can move from target to candidate.

### Visualization / dashboard
29. As a user, I want a clean search box and a scannable result page, so that I can explore without reading docs.
30. As a user, I want volcano plots, heatmaps, and lifespan/survival curves, so that I can see patterns visually.
31. As a user, I want UMAP embeddings for single-cell views where relevant, so that I can see cell-type-specific aging.
32. As a user, I want to download any result as CSV/Parquet, so that I can use it downstream.
33. As a user, I want shareable permalinks to a query, so that I can send a colleague a specific view.
34. As a user, I want the dashboard to show loading/empty/error states clearly, so that I'm never confused about what happened.

### API / developer experience
35. As a developer, I want a documented REST API with OpenAPI/Swagger, so that I can integrate without reading source.
36. As a developer, I want stable, versioned endpoints, so that my integration doesn't break on updates.
37. As a developer, I want consistent error responses with actionable messages, so that I can debug quickly.
38. As a developer, I want response caching on expensive federated queries, so that repeated calls are fast.
39. As a developer, I want rate-limit-friendly behavior toward upstream sources, so that the tool is a good citizen and doesn't get blocked.
40. As a developer, I want example notebooks showing end-to-end use, so that I can get started in minutes.

### Reproducibility / trust
41. As a user, I want every dataset and version recorded, so that I can reproduce a result months later.
42. As a user, I want the harmonization steps documented and inspectable, so that I can trust the comparisons.
43. As a user, I want explicit limitations stated per feature, so that I can calibrate how much to rely on it.
44. As a maintainer, I want a one-command reproducible environment (Docker/Conda), so that anyone can run it identically.
45. As a maintainer, I want CI that runs tests on every change, so that I don't ship regressions.

### Portfolio / evaluator
46. As a hiring manager, I want a README that states the problem, the data, the method, how to run it, and the limitations, so that I can judge engineering and communication in five minutes.
47. As a hiring manager, I want a live demo I can click without installing anything, so that I can evaluate quickly.
48. As a hiring manager, I want to see tests, CI, and clean module boundaries, so that I can judge production-readiness.
49. As a hiring manager, I want a short preprint or write-up of the novel resilience method, so that I can judge scientific depth.

---

## 5. Module Decomposition (deep modules)

The architecture is organized around deep modules, each with a small stable interface and substantial hidden complexity, each independently testable. This decomposition is the spine of the build and the thing to verify against expectations before writing code.

### M1. `idmap` — Identifier & ontology harmonization
- **Responsibility:** resolve any gene identifier (symbol, alias, Ensembl, Entrez, UniProt) to a canonical internal gene ID; map tissue/cell labels to UBERON / Cell Ontology terms; normalize species ages to fractional lifespan via AnAge.
- **Interface (conceptual):** `resolve_gene(query, species) -> CanonicalGene`; `map_tissue(label) -> UberonTerm`; `fractional_age(age, species) -> float`.
- **Hidden complexity:** alias collisions, deprecated IDs, one-to-many mappings, batch resolution against mygene.info, caching, species-specific ID spaces.
- **Why deep:** the rest of the system depends on stable canonical IDs; this interface should almost never change while the internals (which external service, caching strategy) can.

### M2. `sources` — Source adapters (one per data source, common interface)
- **Responsibility:** fetch and parse data from each external source behind one uniform adapter interface; declare per-source capabilities, licence, and whether data is cached or federated live.
- **Interface:** each adapter implements `fetch_signature(gene, filters)`, `fetch_matrix(query)`, or `fetch_records(query)` plus `capabilities()` and `license()`.
- **Adapters (MVP subset in bold):** **GEO (GEOparse)**, **recount3 (federated, via R bridge or cached extracts)**, ARCHS4 (HDF5), **GTEx (Portal API v2)**, CELLxGENE Census (`cellxgene_census`), Human Protein Atlas, **Open Genes (REST/SQLite)**, **HAGR/GenAge/DrugAge (flat files)**, ComputAgeBench (HF Parquet), **NHANES (clinical, for clocks/resilience)**, NIA ITP (Mouse Phenome Database).
- **Hidden complexity:** wildly different formats (SOFT, HDF5, RSE, JSON, XPT, Parquet), rate limits, retries, pagination, partial failures, licence-driven cache-vs-federate decisions.
- **Why deep:** adding a new source should mean writing one adapter, not touching the rest of the system.

### M3. `harmonize` — Cross-study normalization & meta-analysis
- **Responsibility:** turn raw per-study data into comparable aging signatures: batch-effect handling, effect-size computation, random-effects meta-analysis across studies, age-stratification.
- **Interface:** `to_signature(study_data) -> AgingSignature`; `meta_analyze([signatures]) -> MetaSignature`.
- **Hidden complexity:** ComBat/limma-style batch correction, effect-size standardization, handling missing covariates, avoiding naive matrix concatenation, weighting by sample size/quality.
- **Why deep:** this is the scientific core; correctness here is what separates a real tool from a CSV viewer.

### M4. `store` — Analytical storage & query
- **Responsibility:** persist the harmonized layer as Parquet, expose fast analytical queries via DuckDB, keep relational metadata (genes, studies, ontology maps) queryable; manage dataset versioning.
- **Interface:** `query_signatures(filters)`, `get_gene_card(gene)`, `list_studies()`, `version()`.
- **Hidden complexity:** columnar layout for large matrices, partitioning, predicate pushdown, reading Parquet from local/S3/HF, separating big matrices (DuckDB/Parquet) from small relational metadata (SQLite/Postgres).
- **Why deep:** storage engine choices can change without affecting callers.

### M5. `clocks` — Aging clock & biological-age service
- **Responsibility:** wrap `pyaging`, `biolearn`, and `BioAge` methods; validate inputs; apply clocks; compute age-acceleration; annotate each clock with predicted-outcome and training population.
- **Interface:** `list_clocks()`, `apply_clock(clock_id, matrix) -> ClockResult`, `compare_clocks([ids], matrix)`.
- **Hidden complexity:** input format validation (missing CpGs/genes), per-clock preprocessing, GPU optional path, surfacing the chronological-age-vs-mortality distinction.
- **Why deep:** new clocks are added by the wrapped libraries; our interface stays fixed.

### M6. `resilience` — Dynamical-systems / criticality metrics (differentiator)
- **Responsibility:** compute critical-slowing-down indicators (autocorrelation, variance over age strata or time), DOSI-style recovery-rate/resilience, and an optional network-control-energy view; provide a documented fallback (age-stratified variance/autocorrelation) when longitudinal data is absent.
- **Interface:** `csd_indicators(data, strata)`, `recovery_rate(longitudinal_data)`, `control_energy(network, source_state, target_state)` (stretch).
- **Hidden complexity:** windowing, detrending, stratification, handling sparse longitudinal data, network construction, control-energy computation, clear assumption/limit reporting.
- **Why deep:** this is the novel module; isolating it makes it independently testable and independently publishable.

### M7. `api` — FastAPI service layer
- **Responsibility:** expose REST endpoints, OpenAPI docs, request validation (Pydantic), caching, error formatting, versioning, upstream-rate-limit etiquette.
- **Interface:** HTTP endpoints (see §8 API contract).
- **Hidden complexity:** async orchestration across modules, cache keys, consistent error envelopes, request-size limits for uploads.

### M8. `ui` — Streamlit/Gradio dashboard (+ optional React showcase)
- **Responsibility:** scannable exploration UI calling the API; search, gene cards, plots (volcano/heatmap/survival/UMAP), downloads, permalinks, loading/empty/error states.
- **Interface:** user-facing screens; internally only calls `api`.
- **Hidden complexity:** caching, plot rendering, state, responsive layout; React version decouples entirely via the API.

**Module dependency direction:** `ui` → `api` → (`store`, `clocks`, `resilience`, `harmonize`) → (`sources`, `idmap`). Lower modules never import higher ones. This keeps each layer testable with the layer below mocked.

**Decision to confirm with you before coding (per skill workflow):** Does this 8-module decomposition match your mental model? And which modules do you want unit tests written for first? My recommendation is to prioritize tests for **M1 idmap**, **M3 harmonize**, **M5 clocks**, and **M6 resilience** (the modules where correctness is scientific, not cosmetic), and to rely on lighter integration tests for `sources`, `store`, `api`, and `ui`. Confirm or adjust before Phase 1.

---

## 6. Implementation Decisions

### Architecture
- **Backend-first.** Build the FastAPI + DuckDB/Parquet core before any UI. The API is the reusable, high-signal artifact; the UI is a client of it.
- **Federate vs cache by licence and size.** Redistribute/cache only open, redistributable data (GEO, GTEx-open summaries, NHANES, HAGR, Open Genes [MPL-2.0], Tabula, ComputAgeBench, HPA [CC-BY-SA], CELLxGENE [CC-BY]). Federate live (never re-host) for large or controlled sources (recount3 long tail, CELLxGENE Census, GTEx protected, anything in UK Biobank/dbGaP/ADNI/MIMIC class — which are out of scope for redistribution and only linked).
- **Storage split.** Large harmonized matrices/signatures → Parquet queried by DuckDB. Small relational metadata (genes, studies, ontology maps, clock registry) → SQLite for local/dev, with a clean path to Postgres if hosted. This mirrors the Open Genes pattern (relational metadata) plus a columnar matrix store.
- **Publish the harmonized layer as versioned Parquet on Hugging Face Hub** (free, versioned), the same distribution pattern used by opengenes-mcp and ComputAgeBench. The repo ships the ETL that regenerates it.
- **Python-first, R for ingestion only.** Core in Python (FastAPI, GEOparse, mygene, pandas/polars, DuckDB, pyaging, biolearn). Use R/Bioconductor (recount3, GEOquery, gtexr) only inside offline ETL jobs, bridged via a thin subprocess/`rpy2` boundary or a small R microservice — never on the request path.

### Backend / API
- **FastAPI** for the service layer (async, auto OpenAPI docs, Pydantic validation) — the same choice Open Genes made.
- **REST first.** GraphQL only if query flexibility later demands it (deferred to stretch).
- **Caching:** disk cache for federated responses in MVP; Redis if/when hosted with enough traffic. Cache keys include source + query + data version.
- **Versioned endpoints** under `/v1`. Consistent error envelope `{error: {code, message, detail}}`.

### Harmonization
- **Gene IDs via mygene.info** (batch POST) inside `idmap`, with a local cache of resolved mappings.
- **Tissues via UBERON**, cells via Cell Ontology, resolved through OLS where needed.
- **Cross-species age via fractional lifespan** using AnAge maximum-lifespan values.
- **Batch correction:** ComBat / limma `removeBatchEffect` for bulk; Harmony/scVI for single-cell (single-cell is stretch).
- **Meta-analysis:** random-effects models on standardized effect sizes; never naive concatenation of expression matrices.

### Clocks
- **Wrap, don't rebuild.** Expose `pyaging` + `biolearn` + `BioAge` behind `/clock`. Always surface the predicted-outcome distinction (chronological age vs mortality vs pace) because a clock's age accuracy does not imply mortality-prediction ability.

### Resilience (differentiator)
- **Primary methods:** critical-slowing-down indicators (autocorrelation, variance) computed across age strata and, where available, over longitudinal windows; DOSI-style recovery-rate/resilience following the Gero/Pyrkov methodology as the template.
- **Documented fallback:** when dense longitudinal data is unavailable, compute age-stratified variance/autocorrelation as a resilience proxy and label it as such.
- **Anchor datasets:** NHANES (large, age-stratified, clinical markers), the Lehallier plasma-proteome "crests" dataset (non-linear change around ages 34/60/78), and GTEx age brackets. These are the test beds for the module.
- **Stretch:** network-control-energy view (control energy to move an aging network toward a youthful attractor) on supported datasets — the explicit bridge to network control theory.

### Frontend
- **Streamlit (or Gradio for HF Spaces) for the MVP dashboard** — fastest path to a clickable demo, Python-native, one-click deploy.
- **Optional React showcase page later** that calls the same API — strongest *engineering* signal and reuses existing React experience. Build only after the API is stable.
- **Visualization:** Plotly/Altair for volcano plots, heatmaps, survival/lifespan curves; scanpy for UMAPs in the single-cell stretch.

### Deployment / reproducibility
- **Docker** image; **Conda/Bioconda `environment.yml`** + pinned `requirements.txt`.
- **GitHub Actions** CI (lint, type-check, tests) on every PR.
- **Hosting:** Hugging Face Spaces (best for HF-data integration and an ML-style demo) and/or Streamlit Community Cloud for the dashboard; Render/Fly.io for the API if a separate always-on backend is needed.
- **bioRxiv preprint** describing the harmonization method and the resilience module, paired with a problem→process→solution README.

### Schema (conceptual, not file-level)
- **Gene** (canonical_id, symbol, aliases, species, xrefs).
- **Study** (id, source, omic_layer, species, tissue, sample_size, processing_method, license, url, version).
- **AgingSignature** (gene_id, study_id, tissue, sex, age_range, effect_size, direction, p/q-value, omic_layer).
- **MetaSignature** (gene_id, omic_layer, species, pooled_effect, heterogeneity, n_studies).
- **ClockRegistry** (clock_id, library, predicted_outcome, training_population, input_type, required_features).
- **Intervention** (id, name, type, source [DrugAge/ITP/GenDR], lifespan_effect, organism, linked_gene_ids).
- **ResilienceResult** (dataset_id, method, strata, metric_values, assumptions, version).
- **OntologyMap** (raw_label, canonical_term, ontology, confidence).

### API contract (initial, `/v1`)
- `GET /v1/gene/{id}/signature` — multi-omic aging signature for a gene; filters: species, tissue, sex, age_range, omic_layer.
- `GET /v1/gene/{id}/card` — assembled multi-omic + curated-knowledge + intervention summary card.
- `POST /v1/geneset/signature` — aggregate signature for a gene set/pathway.
- `GET /v1/clocks` — list clocks with predicted-outcome + training-population metadata.
- `POST /v1/clock/apply` — apply a clock to a dataset ID or uploaded matrix; returns predictions + age acceleration.
- `POST /v1/clock/compare` — compare multiple clocks on one input.
- `GET /v1/intervention/{name}` — intervention record + linked signatures + lifespan data.
- `POST /v1/resilience/csd` — critical-slowing-down indicators for a dataset + strata.
- `POST /v1/resilience/recovery` — recovery-rate/resilience for longitudinal input.
- `GET /v1/studies` / `GET /v1/sources` — provenance and capabilities.
- `GET /v1/version` — data + code version.
- All list/query endpoints support pagination and a `format=json|csv|parquet` switch.

---

## 7. Step-by-Step Build Plan (phased)

Each phase ends in a demoable, committable state. Estimated effort assumes solo, part-time; adjust to your schedule.

### Phase 0 — Foundations (scaffold, 1 small slice end-to-end)
1. Initialize repo, licence (MIT or Apache-2.0 for your code; document third-party data licences separately), README skeleton (problem→process→solution), `environment.yml` + pinned `requirements.txt`, pre-commit (ruff/black/mypy), GitHub Actions CI running an empty test suite.
2. Stand up FastAPI with `/v1/version` and `/healthz`. Stand up DuckDB + a tiny Parquet fixture. Confirm one end-to-end request returns real data from one source.
3. Build **M1 `idmap`** first (gene resolution via mygene.info + cache) with unit tests — everything depends on it.
4. **Exit criteria:** `GET /v1/gene/{id}/signature` returns a harmonized result for one tissue from one source (GEO or GTEx), gene IDs resolved through `idmap`, tests green in CI.

### Phase 1 — Transcriptomic core (human + mouse)
5. Implement **M2 `sources`** adapters for GEO (GEOparse) and GTEx (Portal API v2); add recount3 as a federated/cached adapter via the R-bridge ETL.
6. Implement **M3 `harmonize`** for transcriptomic signatures: effect-size computation, batch handling, random-effects meta-analysis; cross-species fractional-age normalization via AnAge.
7. Implement **M4 `store`**: write harmonized signatures to partitioned Parquet; relational metadata to SQLite; DuckDB query layer.
8. Wire `GET /v1/gene/{id}/signature` and `/card` to return real human + mouse transcriptomic signatures with provenance.
9. Add curated-knowledge overlay from **Open Genes + GenAge/CellAge/LongevityMap** so each gene card shows which databases flag it.
10. **Exit criteria:** query any gene, get a harmonized human+mouse transcriptomic aging signature with provenance and curated-knowledge flags; meta-analysis tested on a known aging gene (e.g., CDKN2A/p16) against literature expectation.

### Phase 2 — Clock service
11. Implement **M5 `clocks`**: wrap pyaging + biolearn + BioAge; input validation with precise missing-feature errors; age-acceleration computation; clock registry annotated with predicted-outcome + training population.
12. Endpoints `/v1/clocks`, `/v1/clock/apply`, `/v1/clock/compare`; support dataset-ID input and user upload.
13. Validate against a public methylation dataset with known ages; confirm predictions are sane and the chronological-vs-mortality annotation is surfaced.
14. **Exit criteria:** apply ≥3 clocks to a public dataset and to an uploaded matrix; compare them; correct, well-labeled outputs; tests cover the validation/error paths.

### Phase 3 — Dashboard MVP
15. Implement **M8 `ui`** in Streamlit: search box, gene card, volcano/heatmap/survival plots, clock panel, downloads (CSV/Parquet), permalinks, loading/empty/error states.
16. Deploy to Hugging Face Spaces and/or Streamlit Community Cloud.
17. **Exit criteria:** a stranger can open the live demo, search a gene, apply a clock, and download results without instructions.

### Phase 4 — Resilience / criticality module (the differentiator)
18. Implement **M6 `resilience`**: CSD indicators (autocorrelation/variance across age strata), DOSI-style recovery-rate using the Gero/Pyrkov template, documented age-stratified fallback.
19. Endpoints `/v1/resilience/csd` and `/recovery`; validate on NHANES and the Lehallier proteome dataset; reproduce the qualitative "resilience declines with age" result.
20. Add a resilience view to the dashboard with explicit assumptions/limits text.
21. **Exit criteria:** compute and visualize loss-of-resilience-with-age on at least one real dataset, reproducibly, with limitations stated; this is the headline novelty.

### Phase 5 — Intervention linkage + proteomics/methylation breadth
22. Add **DrugAge + NIA ITP + GenDR** intervention adapters; `/v1/intervention/{name}`; gene→intervention and intervention→signature linkage.
23. Add a plasma-proteome slice (Lehallier) and a methylation slice (ComputAgeBench / GEO) to gene cards for true multi-omic coverage.
24. **Exit criteria:** a gene card shows transcriptome + methylation + proteome + curated knowledge + linked interventions; an intervention query returns lifespan data + associated signatures.

### Phase 6 — Polish, docs, launch
25. Example Jupyter notebooks (end-to-end gene query, clock application, resilience analysis).
26. Finalize README, API docs (OpenAPI/Swagger), architecture diagram, limitations section, contribution guide.
27. Optional **React showcase page** calling the API (frontend-engineering signal).
28. Write and post the **bioRxiv preprint** (harmonization method + resilience module). Enter the **Biomarkers of Aging Challenge** and tag the Biomarkers of Aging Consortium.
29. Launch post (GitHub, relevant communities, LinkedIn), and a short demo video.
30. **Exit criteria:** public repo + live demo + preprint + clean docs; ready to put on a CV and send to longevity companies.

**Stretch (post-launch, optional):** single-cell layer (CELLxGENE Census + Tabula Muris Senis, scanpy UMAPs), network-control-energy resilience view, GraphQL, user accounts/saved queries, additional species.

---

## 8. Testing Decisions

**What makes a good test here:** tests assert *external behavior* of a module through its public interface, not its internals. For scientific modules, "behavior" means *correctness against a known ground truth* (a gene with a well-established aging direction, a clock with published ages, a dataset with a documented resilience trend), plus *robustness* (bad input handled with a precise error). Avoid asserting on internal data structures, exact intermediate values that depend on library versions, or implementation details that will churn.

**Modules to unit-test first (correctness is scientific):**
- **M1 `idmap`:** known ID round-trips (symbol↔Ensembl↔Entrez), alias resolution, deprecated-ID handling, batch resolution, fractional-age math against AnAge values. Mock the external service; test the contract.
- **M3 `harmonize`:** effect-size sign/magnitude on a synthetic study with a planted age effect; meta-analysis pooling on synthetic studies with known pooled effect and heterogeneity; batch-correction reduces a planted batch signal. Plus one real-data sanity check on a canonical aging gene.
- **M5 `clocks`:** clock applied to a fixture with known ages returns predictions within tolerance; missing-feature input yields a precise, actionable error; predicted-outcome metadata is present for every registered clock.
- **M6 `resilience`:** on synthetic data engineered to approach a tipping point, CSD indicators rise as designed; recovery-rate decreases as designed; the fallback path triggers and is labeled when longitudinal data is absent; assumptions are reported.

**Modules to integration-test (behavior is plumbing):**
- **M2 `sources`:** contract tests with recorded fixtures (e.g., VCR-style) so adapters parse real payload shapes without hammering live APIs in CI; a small set of opt-in live "smoke" tests run outside CI.
- **M4 `store`:** round-trip write/read of signatures; DuckDB query returns expected rows; versioning increments correctly.
- **M7 `api`:** endpoint-level tests with the lower layers mocked — status codes, error envelope shape, pagination, `format` switch, OpenAPI schema validity.
- **M8 `ui`:** light smoke tests that the app starts and key callbacks return without error; manual QA for visual states.

**Cross-cutting tests:**
- **Reproducibility test:** the same query at the same data version returns identical results.
- **Provenance test:** every returned signature carries source, study, and version.
- **Licence test:** redistributed data is restricted to the allow-list; a controlled source is never cached/re-hosted (assert the adapter is federate-only).

**Prior art / patterns to mirror:** model adapter/contract testing on how mature bioinformatics portals separate ingestion from serving (cBioPortal's validator + service split is the reference); model the clock-service tests on how pyaging/biolearn validate inputs; model the HF-Parquet distribution + reproducibility checks on ComputAgeBench.

**CI:** GitHub Actions runs lint (ruff), format (black), type-check (mypy), and the unit + recorded-fixture integration tests on every PR. Live smoke tests and full ETL run on a manual/scheduled workflow, not on every PR (to avoid rate-limiting upstream sources).

---

## 9. Data Source Contract (build-time reference)

| Source | Adapter role | Access | Format | Cache vs Federate | Licence note |
|---|---|---|---|---|---|
| NCBI GEO | transcriptome + methylation studies | GEOparse / E-utilities | SOFT / series matrix | Cache harmonized extracts | Public, attribute |
| recount3 | uniform RNA-seq (human+mouse) | R/Bioconductor (ETL) | RSE | Federate/cache extracts | Public |
| GTEx | tissue expression, age brackets | Portal API v2 | JSON/GCT | Cache summaries; federate detail | Open summary; raw protected |
| ARCHS4 | uniform RNA-seq counts | HDF5 | HDF5 | Cache slices | Public |
| CELLxGENE Census | single-cell (stretch) | `cellxgene_census` | SOMA/AnnData | Federate | CC-BY |
| Tabula Muris Senis | mouse single-cell (stretch) | S3/figshare/Census | h5ad | Federate | Public |
| Human Protein Atlas | protein aging views | programmatic/TSV | TSV/JSON | Cache | CC-BY-SA |
| Lehallier plasma proteome | proteome aging slice | supp. tables | XLSX | Cache | Per-paper terms |
| ComputAgeBench | methylation benchmark | Hugging Face | Parquet | Cache | Public |
| NHANES | clinical markers (clocks/resilience) | CDC download | XPT/CSV | Cache | Public |
| Open Genes | curated gene knowledge | REST / SQLite | JSON/SQLite | Cache mirror | MPL-2.0 |
| HAGR (GenAge/CellAge/LongevityMap/DrugAge/GenDR) | curated gene/drug/species | flat-file | CSV/TSV | Cache | Free, attribute |
| NIA ITP | mouse lifespan interventions | Mouse Phenome DB | CSV/portal | Cache | Public; embargo/ack rules |
| mygene.info | gene ID harmonization | REST (batch) | JSON | Cache mappings | Public |
| UBERON / Cell Ontology | tissue/cell harmonization | OBO/OLS | OBO/OWL | Cache | Open |
| **Out of scope (federate/link only, never re-host):** UK Biobank, dbGaP, GTEx-protected raw, ADNI, MIMIC, HRS-restricted | — | controlled-access | — | Link/federate only | Controlled; no redistribution |

---

## 10. Out of Scope

- **Re-hosting controlled-access data** (UK Biobank, dbGaP, GTEx-protected raw, ADNI, MIMIC, HRS-restricted). These are linked or federated only, never cached or redistributed.
- **Rebuilding aging clocks from scratch.** The clock service wraps existing libraries; building new clocks is explicitly out of scope (and would duplicate pyaging/biolearn).
- **Being a universal atlas.** No attempt to ingest all of GEO/recount3. The product is a focused, well-harmonized slice plus live federation for the long tail. Scope creep toward "all aging data everywhere" is the primary failure mode and is excluded by design.
- **Wet-lab or clinical claims.** The tool surfaces public data and computed metrics; it makes no diagnostic or treatment claims, and the resilience metrics are framed as research indicators with stated limits.
- **Real-time streaming / always-on ingestion.** ETL is batch and versioned, not streaming.
- **User authentication / multi-tenant accounts** in the MVP (deferred to stretch).
- **GraphQL** in the MVP (REST only; GraphQL deferred).
- **Single-cell integration** in the MVP (deferred to stretch; the schema and adapters are designed to accommodate it).
- **Mobile-native UI.** Web only.

---

## 11. Risks & Mitigations

- **Harmonization eats the timeline.** *Mitigation:* if ingestion/harmonization exceeds ~40% of total effort, cut omic layers and deepen one; ship transcriptome-only first.
- **Resilience module is data-hungry** (needs dense longitudinal data). *Mitigation:* documented fallback to age-stratified variance/autocorrelation; anchor on NHANES + Lehallier + GTEx brackets where the trend is reproducible.
- **An existing tool already does this.** *Mitigation:* if a live multi-omic + clinical + resilience open API surfaces, narrow to resilience-only or intervention-signature-only, where there is clear white space.
- **Upstream rate limits / API changes.** *Mitigation:* recorded fixtures in CI, caching, polite backoff, federate-only smoke tests off the PR path, pinned source versions.
- **Licence missteps.** *Mitigation:* an explicit allow-list for redistribution and a `license()`-driven cache-vs-federate gate enforced by a test.
- **Scope creep.** *Mitigation:* phase exit criteria are hard gates; stretch items stay in the stretch list until MVP ships.
- **Clock misuse / overclaiming.** *Mitigation:* surface predicted-outcome metadata everywhere; document that chronological-age accuracy ≠ mortality prediction.

---

## 12. Success Metrics

**Portfolio (primary):**
- Public repo with clean modules, tests, CI badge, and a problem→process→solution README.
- A live demo anyone can click (no install).
- A bioRxiv preprint on the harmonization method + resilience module.
- At least one concrete, reproducible "wow" result (loss of resilience with age on a real dataset; cross-species conservation of a known aging gene).
- Evidence of reach: GitHub stars/forks, a launch post, an entry in the Biomarkers of Aging Challenge, and ideally a recruiter/engineer reference to it.

**Product (secondary):**
- Query any gene → harmonized multi-omic + cross-species signature in one call.
- Apply ≥3 clocks to public + uploaded data with correct, labeled outputs.
- Reproduce a published resilience trend.
- Gene card unifies ≥3 omic layers + curated knowledge + interventions.

**Engineering quality:**
- Deterministic results at a fixed data version.
- Every result carries provenance.
- One-command reproducible environment.
- Green CI on every change.

---

## 13. Further Notes

- **Why this is a strong career bet, concretely:** Altos Labs lists open-source track record (a high-impact GitHub repo) as a preferred qualification; the field's own surveyed #1 demand is "more public aging biology datasets"; and there is a documented precedent of an open-source aging package (pyaging) preceding a senior industry ML role. The resilience module maps directly onto the physics-of-aging framing used at Gero and the biomarker focus at BioAge — make that connection explicit in the README and preprint.
- **Lead with one excellent biological question, not the platform.** The most legible demo is: "type a gene, see its harmonized multi-omic aging signature across human and mouse, apply a clock, and view its loss-of-resilience signal." Build the narrative around that.
- **Surface the chronological-vs-mortality clock distinction prominently** — it signals real domain understanding that distinguishes you from a generic data-engineer.
- **Use the README as a portfolio instrument:** state the biological question, the data sources and their licences, the harmonization method, how to run it, validation against known literature, and explicit limitations. Stating limitations is a maturity signal evaluators look for.
- **Verify at build time:** database entry counts and sample sizes grow over time; Aging Atlas and Digital Ageing Atlas maintenance cadence should be checked before depending on them as live sources; confirm current ITP data embargo/acknowledgment terms before redistribution.
- **Sequencing tip:** resist building the React showcase or single-cell layer until the API and resilience module are solid — those are the two things that actually differentiate you, and polish on a thin core reads as a thin core.

---

## 14. Open Questions to Confirm Before Phase 1

1. **Module decomposition:** does the 8-module split in §5 match your mental model, or do you want to merge/split any (e.g., fold `idmap` into `harmonize`, or split `resilience` into CSD vs control-energy submodules)?
2. **Test priorities:** confirm unit-test focus on M1/M3/M5/M6 with integration tests elsewhere, or reprioritize.
3. **Frontend path:** Streamlit-only MVP first with React deferred (recommended), or do you want the React showcase in the MVP for maximum engineering signal?
4. **Scope of MVP omics:** transcriptome-only for Phase 1 (recommended) and add methylation/proteome in Phase 5, or attempt two omic layers from the start?
5. **Name:** keep "GeroQuery" or pick another (see §15)?
6. **Licence for your code:** MIT (maximally permissive, common for portfolio) vs Apache-2.0 (patent grant) — recommend MIT unless you prefer the patent clause.

---

## 15. Appendix: Name candidates

- **GeroQuery** (clear, descriptive)
- **Senescope** / **Senescan** (senescence + scope/scan)
- **AgeAtlas API** (descriptive; note potential confusion with existing Aging Atlas — avoid if so)
- **OmniAge** (multi-omic + age)
- **Resilio** (leans into the resilience differentiator)
- **GeroLake** (data-lake framing)

Pick one before Phase 0 so the repo, package, and demo URLs are consistent.
