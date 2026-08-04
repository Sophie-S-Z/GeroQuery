# GeroQuery — complete overview

**Date:** 2026-08-03 · **Status:** every data layer is real-data-backed

One document covering: where the data comes from, what the system is made of,
what it is for, and what it cannot do.

---

## 1. Data sources

### 1.1 Live and ingested

Everything in this table is fetched from a named upstream, verified against a
SHA-256 pinned in `geroquery/sources/manifest.py`, and reproducible with
`make data`.

| # | Source | What it provides | Scale in GeroQuery | Access | Licence | Adapter |
|---|---|---|---|---|---|---|
| 1 | **NCBI GEO DataSets** | Gene aging signatures — young-vs-old Hedges' *g* per gene | 31 datasets → 32 contrasts (27 GEO Series) → **485,905 effect sizes over 46,091 genes**, 10 tissues, human + mouse | FTP SOFT, checksum-pinned | US public domain, attribute | `sources/geo.py` |
| 2 | **GEO methylation series** | Real 450K beta matrices — the clock validation panel | GSE64495 (n=106 controls, ages 2.3–73.7, ships the authors' Horvath output) + GSE30870 (newborns vs nonagenarians, n=40) | Series matrix, checksum-pinned | US public domain | `sources/methylation.py` |
| 3 | **CDC NHANES 2017–2018** | Clinical chemistry + CBC; resilience *and* the 9 markers PhenoAge needs | **4,895 adults** (6-marker resilience cohort) / **4,894** (9-marker PhenoAge cohort) | 4 XPORT files, checksum-pinned | US public domain (NCHS) | `sources/nhanes.py` |
| 4 | **HAGR — GenAge (human)** | Genes with evidence of a role in human ageing | 307 genes, evidence code preserved | ZIP, checksum-pinned | HAGR, non-commercial + attribution | `sources/hagr.py` |
| 5 | **HAGR — GenAge (models)** | Genes whose manipulation changes lifespan | 2,205 parsed; 132 mammalian loaded | ZIP, pinned | HAGR | `sources/hagr.py` |
| 6 | **HAGR — CellAge** | Genes that induce / inhibit / regulate senescence | 949 genes, with PubMed links | ZIP, pinned | HAGR | `sources/hagr.py` |
| 7 | **HAGR — LongevityMap** | Human variants tested for longevity association, **including nulls** | 1,325 assertions | ZIP, pinned | HAGR | `sources/hagr.py` |
| 8 | **HAGR — DrugAge** | Compounds tested against model-organism lifespan | 3,423 experiments → **1,334 compound × organism records**, 54 flagged NIA ITP | ZIP, pinned | HAGR | `sources/hagr.py` |
| 9 | **HAGR — GenDR** | Genes required for dietary restriction's lifespan effect | 214 genes → 6 dietary-restriction records with real gene links | ZIP, pinned | HAGR | `sources/hagr.py` |
| 10 | **mygene.info** | Gene identifier harmonization → canonical Ensembl | ~71k identifiers resolved at build; batched + disk-cached | REST (batch POST) | Public | `idmap/mygene.py`, `idmap/bulk.py` |
| 11 | **GTEx Portal API v2** | Tissue expression context — median TPM + UBERON | ~54 tissues per gene, queried live | REST | Open summary tier | `sources/gtex.py` |
| 12 | **biolearn** | Published epigenetic / clinical clocks | **63 clocks** (Horvath, Hannum, PhenoAge, GrimAge, DunedinPACE…) | Python package | Library terms | `clocks/library.py` |
| 13 | **pyaging** | Published clocks, metadata-registered, artifacts lazy | **173 clocks** | Hugging Face, unauthenticated | Library terms | `clocks/pyaging_clocks.py` |
| 14 | **UBERON / Cell Ontology** | Tissue label harmonization | 12 coarse terms bundled | OBO extract | Open | `idmap/data/tissues.json` |
| 15 | **AnAge** | Maximum lifespan → fractional-age normalization | bundled | flat file | HAGR | `idmap/resolver.py` |

**Totals across the ingested layers: 44 checksum-pinned artifacts, ~680 MB,
from 5 independent upstreams (NCBI, CDC, HAGR, mygene.info, GTEx).**

Two bundled identifier tables — `idmap/data/genes.json` and `anage.json` —
are *derived* artifacts rather than upstreams: regenerate them with
`python -m geroquery.etl.build_idmap`. They used to be hand-written.

### 1.2 The one synthetic table, and why it exists

| Dataset id | Status | Purpose |
|---|---|---|
| `clinical_synthetic_csd` | **Synthetic on purpose** | 720 subjects with critical slowing down *planted* by construction (a latent factor whose variance grows with age). It answers "does the estimator recover an effect known to be there?" — which no real dataset can answer, because in real data you do not know. Registered as a separate dataset id with a `SYNTHETIC` description prefix; a test asserts it never merges with the NHANES table. |

This is the only generated data in the repository. It is **not** a stand-in for
anything missing.

### 1.3 Contracted but not implemented

Declared in `sources/federated.py` with capabilities and licence so the store,
API, and licence tests can reason about them — but with no fetch body. Querying
one raises.

| Source | Why it is a stub |
|---|---|
| recount3, ARCHS4 | The highest-value next step: uniformly-processed RNA-seq removes cross-study processing heterogeneity, the biggest confound in the current pooling |
| CELLxGENE Census | Single-cell; stretch scope |
| GTEx protected tier | Donor age is dbGaP-controlled — federate/link only, never re-host |
| UK Biobank, dbGaP, ADNI, MIMIC, HRS-restricted | Controlled access; the licence gate refuses to cache them |

---

## 2. Core elements

### 2.1 Architecture

Eight modules, strictly layered. Lower layers never import higher ones, so each
is testable with the layer below mocked.

```
ui  ->  api  ->  ( store, clocks, resilience, harmonize )  ->  ( sources, idmap )
                    M4      M5        M6          M3              M2       M1
```

| Module | Responsibility | Key contents |
|---|---|---|
| **M1 `idmap`** | Any identifier → canonical Ensembl; UBERON tissues; fractional age via AnAge | `resolver.py`, `mygene.py` (request path), `bulk.py` (ETL path) |
| **M2 `sources`** | Uniform adapters, each declaring `capabilities()` and `license()`; the cache-vs-federate gate | `manifest.py`, `fetch.py`, `geo.py`, `nhanes.py`, `hagr.py`, `gtex.py`, `federated.py` |
| **M3 `harmonize`** | The statistics | `effect_size.py` (Hedges' *g*), `differential.py` (vectorized per-gene contrasts, BH), `meta.py` (DerSimonian–Laird), `batch.py` |
| **M4 `store`** | Partitioned Parquet via DuckDB + SQLite metadata + content-hash versioning | `store.py` |
| **M5 `clocks`** | 236 published clocks + 3 transparent reference clocks; predicted-outcome metadata; age acceleration | `registry.py`, `library.py`, `pyaging_clocks.py`, `service.py` |
| **M6 `resilience`** | The differentiator — critical slowing down, recovery rate, control energy | `csd.py`, `recovery.py`, `control.py` |
| **M7 `api`** | FastAPI `/v1`, error envelope, pagination, `format=json\|csv\|parquet`, versioned LRU | `app.py`, `service.py`, `schemas.py` |
| **M8 `ui`** | Streamlit dashboard + React showcase | `ui/streamlit_app.py`, `frontend/` |

### 2.2 The data layer's guarantees

These are enforced in code, not documented and hoped for.

1. **Nothing enters the store unverified.** `fetch_artifact` returns bytes whose
   SHA-256 matches the manifest, or raises. Cache hits are **re-verified**, not
   trusted. Downloads land on a temp file and are renamed only after the digest
   matches. A failed verification deletes rather than leaving a valid-looking file.
2. **Offline is hermetic by construction.** With `GEROQUERY_ALLOW_NETWORK` unset
   (the default) a cache miss raises `NetworkDisabledError` instead of reaching out.
   Tests and CI cannot silently depend on an upstream.
3. **Licence gate before persistence.** The store calls `assert_cacheable()`
   before writing anything; controlled sources refuse.
4. **Full vs sample mode is returned, not logged.** Both the NHANES and signature
   adapters report which table produced a number, so a slice-derived estimate can
   never be reported as the full-cohort result.
5. **Real and synthetic never merge.** Separate dataset ids, `REAL`/`SYNTHETIC`
   description prefixes, and a test.

### 2.3 API surface

| Endpoint | Returns |
|---|---|
| `GET /v1/gene/{id}/signature` | Per-study signatures + random-effects pooled estimates (CI, I², k) |
| `GET /v1/gene/{id}/card` | Signatures + curated flags + linked interventions in one view |
| `POST /v1/geneset/signature` | Batch, with resolved/unresolved reported separately |
| `GET /v1/clocks` | All 239 registered clocks with `predicted_outcome` and library status |
| `POST /v1/clock/apply` · `/clock/compare` | Run clocks on a matrix; age acceleration |
| `GET /v1/intervention/{name}` | Every organism the intervention was tested in, mammals first |
| `POST /v1/resilience/csd` | Variance + cross-correlation trends across age strata, with bootstrap CIs |
| `POST /v1/resilience/recovery` | AR(1) relaxation rate, four-regime classification |
| `GET /v1/studies` · `/sources` · `/datasets` · `/version` | Provenance |

All support `format=json|csv|parquet`; list endpoints support `limit`/`offset`.

### 2.4 Scale

| | |
|---|---|
| Signature rows | 485,905 (60 MB CSV → partitioned Parquet) |
| Genes | 46,091 (19,438 human Ensembl, 21,443 mouse Ensembl, 5,210 `ENTREZ:` fallback) |
| Studies / contrasts | 32, from 31 GEO DataSets, 27 GEO Series |
| Curated assertions | 2,720 mammalian (of 5,029 parsed) across 5 databases, 1,880 genes |
| Interventions | 1,340 |
| Clinical subjects | 4,895 real NHANES adults |
| Clocks | 240 (63 biolearn + 173 pyaging + real Levine PhenoAge + 3 reference) |
| Clocks validated on real methylation | **150 ran**; Horvath reproduces the authors' published per-sample ages to r=0.998, MAE 1.39 y |
| Tests | 295 offline + 12 live, 84% coverage |

---

## 3. Use cases

### 3.1 Works today

| # | Use case | How |
|---|---|---|
| 1 | **"Does gene X change with age?"** with an honest answer | `GET /v1/gene/X/signature` — pooled *g*, 95% CI, I², k. The CI is wide enough to return *we cannot say*, which is the point |
| 2 | **Cross-species conservation** | Human and mouse pooled separately for the same ortholog group; compare direction and magnitude |
| 3 | **Tissue specificity** | Filter by tissue across 10 tissues; a muscle effect and a brain effect are separate estimates |
| 4 | **Gene → curated evidence → intervention, in one call** | `/v1/gene/X/card` joins signatures, GenAge/CellAge/LongevityMap/GenDR assertions, and DrugAge/GenDR interventions on one canonical id |
| 5 | **"How much does compound Y extend lifespan, and in what?"** | `/v1/intervention/Y` returns every organism with a median over DrugAge's *significant* experiments and an experiment count — including `null` for "tested, nothing significant" |
| 6 | **Gene-set enrichment-style summary** | `POST /v1/geneset/signature` pools across a list, reporting unresolved genes separately |
| 7 | **Screen for replicated aging genes** | 689 human genes at BH q<0.05 across ≥3 independent contrasts — a defensible starting list |
| 8 | **Falsify a claimed aging marker** | The CDKN2A result (§ `RESULTS_GEO_SIGNATURES.md`) is the worked example |
| 9 | **Clinical resilience / critical slowing down** | `POST /v1/resilience/csd` on real NHANES with subject-level bootstrap CIs |
| 10 | **Method validation** | Run the same estimator on `clinical_synthetic_csd`, where the answer is known |
| 11 | **Clock catalogue and comparison** | 239 clocks with `predicted_outcome` never defaulted — a prostate-cancer classifier is not presented as an aging clock |
| 12 | **Reproducible data provenance** | `/v1/studies` gives GEO accession, GEO Series, platform, PMID, sample counts, and the exact restrictions applied to build each contrast |
| 13 | **Offline / air-gapped use** | The committed curated slice (40,585 real rows, both species) runs the whole API with no network |

### 3.2 Enabled but not yet exercised

| Use case | What is missing |
|---|---|
| ~~Run the clocks on real data~~ | **Done** — see [`RESULTS_METHYLATION_CLOCKS.md`](RESULTS_METHYLATION_CLOCKS.md) |
| Correlate epigenetic age acceleration against resilience metrics | Same. This cross-layer claim is the reason the project exists |
| AR(1) recovery-rate analysis | Longitudinal data. `resilience/recovery.py` is implemented and tested but has never run on a real trajectory |
| Nationally-representative clinical claims | `WTMEC2YR` is carried but not applied |

---

## 4. Limitations

### 4.1 About the gene signatures

1. **Microarray-era only.** GEO stopped curating DataSets ~2016. No RNA-seq.
2. **Small groups** — 3–15 per study. **Median 95% CI width 0.84.** Moderate
   effects are undetectable; a null means "not detected by this instrument".
3. **Lopsided tissue coverage.** Half the human contrasts are skeletal muscle;
   one blood contrast (n=9); no brain-region resolution.
4. **`omic_layer` is `transcriptome` everywhere.** No proteome, no methylome. The
   multi-omic schema is real; the multi-omic *data* is one layer deep.
5. **Cross-study processing heterogeneity is uncorrected** — 12 platforms, three
   decades of normalization practice.
6. **32 contrasts come from 27 Series.** Paired array halves share subjects; the
   random-effects pool assumes independence. `series_id` exposes this; the pooling
   does not yet account for it.
7. **Sex is `unspecified` for 24 of 32 contrasts** — GEO never recorded it. That
   is deliberately distinct from `both`.
8. **~22% of probe rows dropped** at identifier resolution, concentrated on the
   oldest platforms.
9. **Age bands are conventions** (human 40/60 y, mouse 8/18 mo), not measurements.

### 4.2 About the clinical / resilience layer

10. **NHANES is cross-sectional.** Age strata proxy for within-individual
    trajectories. No relaxation time is observed — which is precisely what would
    distinguish critical slowing down from ordinary accumulating heterogeneity.
11. **Age is topcoded at 80**, compressing the stratum where CSD should be
    strongest. 354 subjects sit at the cap.
12. **Survey weights carried but not applied.** Estimates are unweighted and
    therefore not nationally representative.
13. **No confounder adjustment** — no fasting status, medication, or comorbidity.
    Some variance growth is likely disease prevalence, not a dynamical property.
14. **The headline CSD result is a partial null**: variance rises with age in
    20/20 configurations, cross-correlation in 0/20. See `RESULTS_NHANES_CSD.md`.

### 4.3 About the curated knowledge

15. **Only mammalian rows are loaded.** 2,276 of 5,029 HAGR assertions are yeast,
    worm, or fly and are parsed, counted, and then dropped — GeroQuery's canonical
    gene space is human + mouse. Ortholog mapping would recover them.
16. **DrugAge records no gene targets, so none are asserted.** The previous
    hand-written table linked rapamycin → MTOR; that edge was never in the data.
    Only GenDR, which does record gene-level dependence, populates
    `linked_gene_ids`.
17. **One GenDR row is malformed upstream** and is dropped by a validity check.
18. **HAGR releases are pinned by checksum**, so a HAGR update breaks the build
    loudly. That is intended, but it means the pins need periodic maintenance.

### 4.4 About the clocks

19. **The clocks are now validated on real methylation, with limits.** Two blood
    450K series only; nothing is established about other tissues, EPIC arrays, or
    sequencing-based methylation. 59 of 209 attempts did not run — mostly because
    the clock is not a methylation clock, needs EPIC probes, or needs CpGs the
    series lacks. The `GPAge*` family cannot run at all: GPy requires numpy < 2
    and biolearn's cvxpy requires numpy >= 2.
20. **pyaging imputes missing features from reference values.** The wrapper
    refuses above 50% imputation rather than return a confident-looking artifact
    of pyaging's reference cohort.
21. **The clock tier needs Python 3.10–3.12.** `ecos`, a transitive biolearn
    dependency, has no 3.14 wheel. The core repo runs fine on 3.14.
22. **The 3 `*_demo` reference clocks are transparent test instruments**, not
    clinical ones.

### 4.5 Operational

23. **`make signatures` needs ~313 MB and ~25 minutes**, dominated by SOFT parsing
    and ~71k mygene.info lookups. Both are cached after the first run.
24. **GTEx open cannot be age-stratified.** The API accepts
    `attributeSubset=ageBracket` and returns HTTP 200 with one undivided group —
    donor age is dbGaP-controlled. A live test pins this so it is rechecked
    automatically rather than from memory.
25. **Live tests are excluded by default** (`-m 'not live'`). They are the only
    thing that catches upstream drift.
26. **The Streamlit dashboard has not been repointed** at the real signature data
    (the React app has).

---

## 5. Things to note

- **The most important table in the repo is §1**, and the most important property
  is that every row of it is checksum-pinned. If a pin fails, an upstream changed;
  confirm and bump the pin rather than deleting the check.
- **A null result is a result.** CDKN2A/p16 not replicating, and the CSD
  cross-correlation null, are both reported as findings. The synthetic layers they
  replaced would have said yes to both.
- **Two interpreters on purpose.** Python 3.14 runs the core suite without the
  clock libraries; 3.12 runs everything. Both staying green is what proves the
  optional tier is genuinely optional.
- **`get_registry()` memoizes.** Tests that fake a clock library must reset
  `clocks.registry._REGISTRY` or they pass alone and fail in a full run.
- **Reproduce everything:** `make data` (network) or `make data-offline`
  (committed real samples), then `make ci`.
