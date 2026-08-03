# GeroQuery — Engineering & Strategy Handoff

**Date:** 2026-08-03
**Repo state at review:** branch `main`, clean tree, commit `ed99669` ("Initial commit")
**Verification performed:** full codebase read; `pytest` run locally → **77 passed**, 0 failed (~1s)
**Audience:** whoever picks up GeroQuery next (including future-you)

---

## 0. TL;DR

GeroQuery is a **well-architected, fully-tested skeleton with correct scientific math and no real data**.

The pipes are real. The water is fake.

- ✅ 8 clean deep modules, strict dependency direction, 77 passing tests, CI, Docker, reproducible ETL
- ✅ Genuinely correct implementations: DerSimonian–Laird meta-analysis, critical-slowing-down indicators, AR(1) recovery rate, controllability Gramian
- ❌ **228 synthetic signature rows across 12 genes.** Zero real network calls. Zero real datasets. Three hand-written "demo" clocks.

**The single highest-leverage action is Tier 0 below: load one real dataset end-to-end.** Until that happens, every other improvement is decoration on an unloaded system. After that, three moves (§6) take it from "portfolio project" to "tool a frontier lab actually runs."

---

## 1. What GeroQuery is

A FastAPI + DuckDB/Parquet service answering one question: **"type a gene, get its harmonized aging biology."**

It unifies four things that normally live in four incompatible places:

1. Multi-omic aging signatures (transcriptome / methylome / proteome)
2. Biological-age clocks
3. Intervention & lifespan data
4. A dynamical-systems **resilience / criticality** view — *the differentiator; no other aging portal exposes this*

...behind one versioned `/v1` REST API, with a Streamlit dashboard and a React showcase as clients.

**Stated purpose** (PRD §1): a portfolio artifact aimed at longevity biotech hiring (BioAge, Retro, NewLimit, Gero, Insilico, Altos, Calico), with a real field gap as the secondary target. Both goals are legitimate and they mostly point the same direction — but see §7 for where they diverge.

---

## 2. Architecture as built

```
ui  ->  api  ->  ( store, clocks, resilience, harmonize )  ->  ( sources, idmap )
                    M4     M5        M6         M3              M2       M1
```

Dependency direction is **actually enforced** in the code, not just documented. Lower modules never import higher ones. This is the strongest engineering signal in the repo and should be preserved through every change below.

| Module | Path | Status | Notes |
|---|---|---|---|
| M1 `idmap` | `geroquery/idmap/resolver.py` | Real logic, tiny bundled data | 24 genes (12 ortholog groups, human+mouse), 12 UBERON tissue terms, AnAge lifespans. No live mygene.info. |
| M2 `sources` | `geroquery/sources/` | **Interface only** | `base.py` licence gate is real and load-bearing. `federated.py` is 5 metadata-only stubs with **no fetch bodies**. |
| M3 `harmonize` | `geroquery/harmonize/meta.py` | **Real & correct** | DerSimonian–Laird with τ², I², Cochran's Q. Correctly refuses naive matrix concatenation. |
| M4 `store` | `geroquery/store/store.py` | Real, small-scale | Partitioned Parquet via DuckDB + SQLite metadata + content-hash data versioning. Works; see perf defects §5. |
| M5 `clocks` | `geroquery/clocks/registry.py` | **Scaffolding** | 3 hand-written linear demo clocks. `_library_clocks()` returns `{}` — pyaging/biolearn never actually wrapped. |
| M6 `resilience` | `geroquery/resilience/` | **Real math, synthetic input** | CSD (variance + cross-correlation), AR(1) DOSI recovery, control energy via Gramian. The moat. |
| M7 `api` | `geroquery/api/app.py` | Real, 13 endpoints | One error envelope, pagination, `format=json\|csv\|parquet`, versioned LRU cache. |
| M8 `ui` | `geroquery/ui/`, `frontend/` | Present | Streamlit dashboard + React/Vite showcase. |

**Total:** ~4,169 LOC across Python + JS + tests.

---

## 3. What is real vs. what is scaffolding

This is the central finding. The README is admirably upfront about it (`README.md:100-107` "Limitations (read these)"), but the consequence needs stating plainly:

### Real

- The **module decomposition and dependency discipline**
- **`harmonize/meta.py`** — publication-grade random-effects meta-analysis
- **`resilience/csd.py`, `recovery.py`, `control.py`** — the math is correct and the assumptions are explicitly reported in the result payload (`CSDResult.assumptions`)
- The **licence gate** (`sources/base.py:50-67`) — `assert_cacheable()` raises `LicenseViolationError` on federate-only sources, and a test asserts it. Rule enforced in code, not in a README. This is a genuinely mature pattern.
- **Fractional-lifespan age normalization** via AnAge — the best single idea in the repo
- **Clock `predicted_outcome` metadata** — every clock declares chronological-age vs mortality vs pace. Real domain understanding.
- Determinism: seeded ETL → byte-identical CSVs → content-hashed data version

### Scaffolding

| Claim | Reality |
|---|---|
| "Multi-omic aging signatures" | **228 rows, 12 genes.** Generated by `etl/build_fixtures.py`. Effect *directions* encode real biology (`TRUE_EFFECTS`, line 26); magnitudes are `rng.normal(0, 0.15)`. |
| Study provenance (`GEO:GSE90001`...) | **Fabricated accessions.** URLs point at GEO but the GSE numbers are a counter starting at 90000 (`build_fixtures.py:122`). |
| "NHANES clinical slice" | **720 synthetic subjects** (`build_clinical()`, line 195) generated by a latent factor whose variance is *explicitly* engineered to grow with age (`0.4 + 1.8 * aging`, line 225). |
| Resilience headline result | Currently **proves the generator works**, not that aging biology does. Circular. |
| Federated sources (recount3, GTEx, CELLxGENE, UK Biobank) | `FederatedStub` objects. Declare capabilities + licence correctly. **No `fetch` implementation anywhere.** |
| pyaging / biolearn wrapping | `_library_clocks()` detects installation then returns `{}` (`registry.py:139-154`). |
| mygene.info resolution | Not called. 24 genes hardcoded in `idmap/data/genes.json`. |

> **Bottom line for an evaluator:** the current artifact makes a strong *engineering* claim and a **zero scientific claim**. A data scientist at Gero or NewLimit will find `_library_clocks() -> {}` in about 90 seconds. That's the seam to close first.

---

## 4. Verified inventory

```
signatures.csv            228 rows   (12 ortholog groups × omic × species × tissue × 2-3 studies)
studies.csv               228 rows   (fabricated GSE accessions)
curated_knowledge.csv      50 rows
interventions.csv           7 rows   (rapamycin, CR, metformin, resveratrol, D+Q, NMN, 17a-estradiol)
clinical_nhanes_slice.csv 720 rows   (synthetic)
genes.json                 24 genes  (12 groups, human + mouse)
tissues.json               12 UBERON terms
anage.json                  2 species
```

CI runs ruff + black + mypy + pytest on Python 3.10/3.11/3.12. Green.

---

## 5. Defects found

Ordered by severity. File references are exact.

### HIGH — `recovery.py` silently misreports anti-persistent series
`geroquery/resilience/recovery.py:51`

```python
a_clipped = min(max(a, 1e-6), 1 - 1e-9)  # keep in (0,1) for a decaying process
```

A **negative** AR(1) coefficient (oscillatory / anti-persistent dynamics) is clipped to `1e-6`, producing a huge `recovery_rate` and the interpretation string *"fast recovery — high resilience"* (line 58, `a <= 0.5`). Anti-persistence is not resilience. Real biomarker series with measurement noise or alternating artifacts will hit this and get a confidently wrong, favorable answer.

**Fix:** detect `a < 0` and return an explicit `interpretation` of "oscillatory / AR(1) model inappropriate" rather than coercing into the decaying-process regime.

### HIGH — CSD variance is confounded by within-stratum mean drift
`geroquery/resilience/csd.py:120`

Markers are z-scored **globally** (`_zscore_columns`, line 50), then within-stratum variance of the mean-z state is computed. With 6 strata spanning ages 20–85, each bin covers ~11 years. If markers drift monotonically with age, that within-bin age heterogeneity inflates within-bin variance — **mimicking critical slowing down with no CSD present**. This is the most likely referee objection to the entire differentiator module.

**Fix:** regress out age *within* each stratum before computing variance; report both raw and detrended.

### MEDIUM — CSD trend test is statistically thin
`geroquery/resilience/csd.py:67-77, 40-42`

Kendall τ over `n_strata=6` points. Floor p-value ≈ 0.0028. `resilience_declines` (line 42) is a bare `slope > 0 and slope > 0` with **no significance requirement at all** — it will return `True` on pure noise ~25% of the time.

**Fix:** bootstrap CIs over subject resampling; gate `resilience_declines` on the CI excluding zero.

### MEDIUM — Gramian inversion will silently return garbage
`geroquery/resilience/control.py:61`

`np.linalg.solve(W, delta)` catches only `LinAlgError`. Real biological networks are *near*-uncontrollable: `W` is typically ill-conditioned rather than exactly singular, so `solve` succeeds and returns numerically meaningless energies spanning many orders of magnitude.

**Fix:** compute and report `cond(W)`; use `pinv` with an explicit rank-truncation threshold; refuse (or flag) above a conditioning cutoff.

### MEDIUM — Dead CI job
`.github/workflows/ci.yml:35-36`

```yaml
smoke-live:
  if: github.event_name == 'workflow_dispatch'
```

The `on:` block (lines 3-6) only lists `push` and `pull_request`. **`workflow_dispatch` is never emitted, so this job can never run.** The live-smoke-test path documented in the PRD does not exist in practice.

**Fix:** add `workflow_dispatch:` (and likely `schedule:`) to the `on:` block.

### MEDIUM — Query performance won't survive real data
`geroquery/store/store.py:213`, `geroquery/api/app.py:72-77`

- A fresh `duckdb.connect()` per `query_signatures()` call; no connection reuse.
- `ensure_built()` hits the filesystem on every query (line 194).
- `_paginate` slices **after** full materialization — no predicate/limit pushdown into DuckDB.

Fine at 228 rows. Falls over at the millions of rows real ingestion produces.

### LOW — Swallowed exception
`geroquery/api/service.py:141-143` — bare `except Exception:` in `geneset_signature` silently converts *any* failure (including bugs) into "unresolved gene."

### LOW — Exception class defined inside a function body
`geroquery/api/service.py:194-197` — `InterventionNotFound` is redefined on every call. Move to `exceptions.py`.

### LOW — Missing production hardening
No auth, no rate limiting, no upload size limits — all named in PRD §6 but absent. Acceptable for a demo; blocking for a hosted service.

> Note: `geroquery/_build/` is **correctly gitignored** (`.gitignore:12-15`) and untracked. No action needed there.

---

## 6. THE ROADMAP — how to make this frontier-lab-level

### The trap to avoid

**"Add more data sources" makes it a better aggregator, and every aggregator loses to the next aggregator.** Aging Atlas, HAGR, and Open Genes already own that ground and have more people on them. Competing there is a losing position.

**The reframe:** GeroQuery's differentiator was never *"we have the data."* It is:

> **"We measure resilience, and we can tell you what to perturb."**

The harmonization layer is table stakes that earns the right to make those two claims. The repo currently *leads* with the table stakes and treats the differentiator as a module. **Invert that.**

---

### TIER 0 — Credibility floor
*Non-negotiable. Until one real dataset flows end-to-end, nothing downstream counts.*

#### 0.1 Load one real dataset completely — **START HERE**
Pick **NHANES**. Public, no access agreement, ~50k subjects, and it carries *exactly* the biomarkers `clocks/registry.py:28` already expects (`albumin, creatinine, glucose, crp, lymphocyte_pct, rdw`).

- Replace `build_clinical()` (`etl/build_fixtures.py:195`) with a real XPT ingest (CDC download → pandas `read_sas`).
- Keep the existing schema so `store.py` and the resilience service need **zero** changes.
- Keep the synthetic generator alongside, renamed `build_clinical_synthetic()`, used only by tests.

**Why this specific dataset first:** the moment CSD rises with age in *real* NHANES rather than in your own generator, you stop having a demo and start having a scientific result. Everything else in this roadmap gets easier to justify.

#### 0.2 Wire two real network paths
- **mygene.info** batch resolution in `idmap` (the PRD's own design; cache into the existing `genes.json` shape — the file's `_provenance` field already describes this exact plan).
- **GTEx Portal API v2** as the first live source adapter. Lowest friction: open, REST, age brackets built in. Implement `fetch_signature` on a real `SourceAdapter` subclass rather than `FederatedStub`.

This permanently kills the "no real network calls" objection.

#### 0.3 Actually wrap pyaging / biolearn
`registry.py:139-154`. This is the most visible unfinished seam in the repo. Enumerate the installed library's clocks, wrap each in the existing `ClockInfo` shape, preserve `predicted_outcome` metadata. The interface already exists — only the body is missing.

#### 0.4 Fix the HIGH/MEDIUM defects in §5
Specifically the two resilience defects (5.1, 5.2) **before** publicizing the resilience module. Publishing a CSD result that a referee can attribute to within-stratum mean drift would be worse than not publishing.

**Tier 0 exit criteria:** `GET /v1/gene/{id}/signature` returns real GTEx data for a real gene resolved through real mygene.info; `POST /v1/resilience/csd` on real NHANES shows rising variance + cross-correlation with bootstrap CIs; `GET /v1/clocks` lists real Horvath/PhenoAge/DunedinPACE.

---

### TIER 1 — The three moves that change the category

#### 1.1 Ship a *benchmark*, not just a database

**Highest-leverage single move.** Databases get cited; **benchmarks get adopted**. ComputAgeBench did exactly this for clocks and is now a standard reference in the field.

Build **GeroBench**: held-out aging signatures with known ground truth plus a public leaderboard scoring harmonization / meta-analysis / clock methods.

The scaffolding already exists — **your test suite *is* an eval harness**:
- cross-species conservation (CDKN2A up in human *and* mouse; LMNB1 down)
- planted-effect recovery with known pooled effect and heterogeneity
- engineered tipping points where CSD must rise
- clock MAE against known ages

Promote these from `tests/` into a versioned, public product with a submission format.

> **Why transformational:** a benchmark makes every other lab's work legible *through your interface*. That is how a tool becomes infrastructure instead of a destination. It also converts your biggest current weakness (synthetic data) into a strength — synthetic ground truth is exactly what a benchmark *wants*.

#### 1.2 Go causal — the aggregator → discovery jump

GeroQuery currently answers *"does this gene correlate with age?"* Every existing portal answers that. The frontier question is:

> **"If I perturb this, does the system get younger?"**

Three concrete builds, in order of impact:

**a) Signature-reversal engine (Connectivity-Map style)** — *the single most valuable feature you could add.*
Given an aging signature, rank perturbagens whose transcriptional signature is **anti-correlated** with it. LINCS L1000 is public, ~1.3M profiles.

This is literally the daily workflow at NewLimit, Retro, and Altos. It turns "here is what changes with age" into "here is what to try" — and the second one is what companies pay for. New endpoint: `POST /v1/reversal` taking a signature or gene set, returning ranked perturbagens with connectivity scores.

**b) Mendelian randomization / causal-genetics overlay**
Open Targets Genetics + GWAS Catalog to flag which age-associated genes have **causal human genetic support**. Turns a correlation list into a target list. Fits cleanly as a new field on the existing gene card.

**c) Perturb-seq / reprogramming trajectory ingest**
The dataset class the field is generating *right now*, and which **no existing portal indexes**. First-mover ground is genuinely open here.

#### 1.3 Make resilience longitudinal — claim the whole differentiator

The CSD module is the moat, but it is currently a labelled cross-sectional proxy. The DOSI/Gero result that made the field pay attention required **longitudinal** measurement.

- **Ingest a real longitudinal cohort** so `recovery_rate()` runs on true repeated measures. NHANES has no follow-up — but **wearable / step-count data** does, and Pyrkov's original DOSI work used physical-activity streams exactly this way. MIMIC time series is an alternative.
- **Implement the DOSI decomposition itself:** first PC of the biomarker set; log-variance and autocorrelation as functions of age; extrapolate the age at which relaxation time diverges. That number — *"resilience → 0 around age ~120–150"* — is a headline result that travels on its own.
- **Extend `recovery.py` to panel data:** a multi-subject estimator, not just a single series. Currently it takes one 1-D array.

Then the pitch becomes one sentence a frontier lab understands instantly:

> *"GeroQuery is the only open tool that measures how fast your biology recovers — not just how old it looks."*

---

### TIER 2 — The frontier-lab wrapper

What makes a repo *read* as a lab initiative rather than a side project.

| Move | Why it matters | Effort |
|---|---|---|
| **MCP server** (`geroquery-mcp`) | The current distribution channel. Every LLM agent doing aging research queries you directly. `opengenes-mcp` proved the pattern — out-covering them is a land-grab. Your API is already clean enough to wrap thinly. | Low |
| **Registered prospective predictions** | Publish signed, timestamped predictions **before** the next NIA ITP cohort reports. When one hits, you have what no dashboard can buy: a track record. **Highest-status move on this list.** | Low effort, long payoff |
| **Versioned data releases on HF Hub + DOI** | Makes results citable and reproducible-by-construction. Your content-hash versioning (`store.py:173`) already supports it. | Low |
| **Uncertainty everywhere** | Every effect size, clock prediction, and resilience metric ships with a CI. Frontier labs distrust point estimates; visible calibration signals you've deployed models before. | Medium |
| **Negative-results surface** | Genes/interventions with null or contradictory evidence. Nobody publishes these, everybody needs them — and it's **nearly free** given τ²/I² are already computed in `meta.py`. | Low |
| **`/v1/explain` provenance trace** | Return the full derivation chain for any number: study → effect size → pooling weight → result. This is the auditability bar real biotech data platforms are held to. | Medium |
| **bioRxiv preprint + Biomarkers of Aging Challenge entry** | Already PRD §7 Phase 6. Converts a repo into a citable method. | Medium |

---

### What to cut

- **The React showcase** (`frontend/`). It adds engineering signal you already demonstrate elsewhere and consumes time Tier 0 needs. Your own PRD §13 says it: *"polish on a thin core reads as a thin core."*
- **Breadth of omic layers.** Two *real* layers beat five synthetic ones. Depth is the credibility signal.
- **Chasing the aggregator crown.** Explicitly deprioritize "more sources" against "more causal / more longitudinal."

---

### Recommended sequencing

```
NOW      Tier 0 — real NHANES + real GTEx + real mygene + real pyaging
         + fix resilience defects 5.1 / 5.2
         └─> unlocks the right to make any scientific claim at all

NEXT     Longitudinal resilience + DOSI decomposition
         └─> owns the differentiator outright

THEN     Signature-reversal engine + causal-genetics overlay
         └─> aggregator becomes discovery tool

PARALLEL MCP server + versioned HF releases + negative-results surface
         └─> cheap, compounding distribution

FINALLY  GeroBench leaderboard + preprint + registered predictions
         └─> category ownership
```

---

## 7. Strategic note: where the two goals diverge

The PRD names a primary goal (portfolio) and a secondary goal (real tool). They align through Tier 0 and most of Tier 1 — but they **split at Tier 2**.

- **Pure portfolio optimum:** polish, demo video, launch post, React page. Fast, visible, shallow.
- **Real-tool optimum:** benchmark, MCP server, registered predictions. Slower, compounding, and *far* more impressive to the exact people the portfolio targets.

Recommendation: **optimize for the real tool.** A hiring manager at Gero cannot distinguish a good portfolio project from a great one, but they can instantly recognize a tool their own team would install. The real-tool path dominates the portfolio path on portfolio outcomes.

---

## 8. Concrete first-week task list

1. Add `workflow_dispatch:` to `.github/workflows/ci.yml:3` — 1 line, unblocks live smoke tests.
2. Fix `recovery.py:51` negative-AR(1) handling; add a regression test with an anti-persistent series.
3. Add within-stratum age detrending to `csd.py`; add bootstrap CIs; gate `resilience_declines` on the CI.
4. Write `sources/nhanes.py` — a real `SourceAdapter` that downloads and parses NHANES XPT.
5. Rename `build_clinical` → `build_clinical_synthetic`; point the ETL at the real adapter; keep the synthetic path test-only.
6. Re-run the resilience test suite against **real** NHANES. **Record whatever happens — including if the effect is weaker than the synthetic version.** That honesty is itself the portfolio asset.
7. Implement `_library_clocks()` against pyaging.
8. Update `README.md` §"Limitations" to reflect the new real/synthetic boundary.

---

## 9. Open questions for the owner

1. **Longitudinal data source** for real DOSI recovery — wearables, MIMIC, or a public cohort with follow-up? This gates Tier 1.3 and is the only genuinely hard sourcing problem in the roadmap.
2. **Benchmark scope** — clocks only (competes with ComputAgeBench) or harmonization + meta-analysis (uncontested)? Recommend the latter.
3. **Hosting** — is an always-on API needed, or is HF Spaces + the Python package sufficient? Determines whether auth/rate-limiting move from LOW to blocking.
4. **LINCS licensing** — confirm redistribution terms before wiring the reversal engine into the cacheable path; the `assert_cacheable` gate will need a correct `License` entry.

---

## 10. File reference map

| Concern | File |
|---|---|
| Synthetic data generation | `geroquery/etl/build_fixtures.py` |
| Meta-analysis (the good part) | `geroquery/harmonize/meta.py` |
| CSD indicators | `geroquery/resilience/csd.py` |
| DOSI recovery rate | `geroquery/resilience/recovery.py` |
| Control energy | `geroquery/resilience/control.py` |
| Clock registry + the `{}` stub | `geroquery/clocks/registry.py:139` |
| Licence gate | `geroquery/sources/base.py:50` |
| Federated stubs (no bodies) | `geroquery/sources/federated.py` |
| Orchestration | `geroquery/api/service.py` |
| HTTP layer | `geroquery/api/app.py` |
| Storage + versioning | `geroquery/store/store.py` |
| Reference gene/tissue/lifespan data | `geroquery/idmap/data/*.json` |
| Product spec | `docs/PRD_GeroQuery_Multi-Omic_Aging_Aggregator (1).md` |
