# Code review — maintainability & correctness pass (2026-11-03)

A senior-engineer review of the whole codebase for dead code, duplication, unused UI, over-
complexity, legacy code, redundant calls, abandoned files, and technical debt — followed by a formal
lint, targeted fixes, and a safe cleanup. The goal was to simplify and harden without breaking
anything.

## Baseline

The codebase was already in good health: `ruff`, `black`, and `mypy` were clean, and **511 tests
passed**. The only red was environmental, not code: the container ships packages newer than the
project's pins (`mcp 2.3.0` vs the pinned `mcp<2`; `pandas 3.0.5` vs the `clocks` extra's `<3`), and
the default `python` (3.11) has an empty site-packages while the dependencies live under
`python3.13`. So this was a cleanup-and-harden pass, not a rescue.

## What changed in this pass

| Change | Type | File(s) |
|---|---|---|
| Harden the MCP transport test to skip cleanly on an incompatible SDK version | test robustness | `tests/test_mcp_transport.py` |
| Fix the "About & data" adapter table (every Source `—`, every Licence a raw dict) + regression test | **bug fix** | `ui/streamlit_app.py`, `tests/test_dashboard.py` |
| `clock_tiers()` now computes the built-in count from the registry (was hardcoded `3`, ignoring PhenoAge) + reconciliation test | inconsistency | `api/service.py`, `tests/test_api.py` |
| Cox `n_psu` no longer reports the subject count for a weighted *unclustered* fit (now `0`, matching the unweighted path) + test | correctness | `survival/cox.py`, `tests/test_survival.py` |
| Remove the network-control-energy module (zero reach from any surface) + its service method, exports, tests, and doc mentions | dead-code removal | `resilience/control.py` (deleted), `resilience/service.py`, `resilience/__init__.py`, `tests/test_resilience.py`, `README.md`, `docs/HANDOFF.md`, `docs/OVERVIEW.md` |
| Remove the unused `rbridge`/`rpy2` optional extra (no importer, no `.R` files, not in Makefile/CI) | unused dep | `pyproject.toml` |
| Add a mypy override so the optional `mcp` SDK's version drift can't break type-checking | tech debt | `pyproject.toml` |

All `ruff`/`black`/`mypy`/`pytest` stayed green after each step.

## 1. Dead code

- **`resilience/control.py` (removed).** Network control energy (controllability Gramian, minimum
  control energy) — ~190 lines, correct and tested, but reachable from **nothing**: no API route, no
  Streamlit tab, no MCP tool, no ETL step. *Impact:* smaller resilience surface, one fewer "stretch"
  that implied a capability the product never exposed. *Risk:* low — removal was contained (service
  method, exports, tests, and three doc lines updated together); the live metrics (CSD, AR(1)
  recovery) are untouched.
- **`rbridge`/`rpy2` extra (removed).** Declared but never imported; no R scripts; not wired to
  Makefile or CI. *Impact:* one fewer phantom install path. *Risk:* none.
- **`gene_report`'s empty `knowledge`/`hallmarks` slots (retained, documented).** These read as dead
  but are an *intentional* forward-compat placeholder — the docstring (`api/service.py`) explains the
  shape is kept stable for a future sourced gene→hallmark mapping. Left as-is deliberately.

## 2. Duplicate logic

- **Scalar vs matrix Hedges' g** (`harmonize/effect_size.py` ↔ `harmonize/differential.py`) — an
  *intentional* pair: the scalar version is the documented oracle that a test uses to pin the
  vectorized one. Not consolidated (removing it would orphan the test oracle).
- **CSV reading** split between `csv.DictReader` (`local_fixture.py`, `hagr.py`) and `pd.read_csv`
  (NHANES/GEO/store adapters). A low-priority consolidation opportunity; left alone this pass (each
  adapter's read is small and correct; a shared helper is a refactor, not a fix).

## 3. Unused UI components

None found. The one UI defect was a *rendering bug*, not an unused component: the About-tab adapter
table was built against the pre-merge `sources()` shape (`source_name` + a flat `license`), so every
Source cell rendered `—` and every Licence cell printed a raw Python dict. The dashboard test only
asserted "no exception", so it shipped green. **Fixed, with a regression test that asserts real
content and the absence of the dict/dash artifacts.**

## 4. Overly complex implementations

Nothing was refactored for its own sake. The heavy statistical machinery (random-effects + Hartung-
Knapp meta-analysis, corpus-wide shrinkage/lfsr, weighted Cox with the design-based sandwich, the
CSD bootstrap) was independently reviewed and found **correct and appropriately scoped** — complexity
that matches the problem, not accidental complexity. The one genuinely over-built piece relative to
its use (control energy) was removed (§1).

## 5. Legacy code

- The three `*_demo` reference clocks (`clinical_phenoage_demo`, etc.) read as "demo/legacy" but are
  **load-bearing**: they are the always-available transparent tier, the default across the test suite
  and the shipped notebook, and the fallback when the optional clock libraries aren't installed.
  **Retained.**
- `clock_tiers()` hardcoded `reference n_clocks = 3`, a stale constant that disagreed with
  `list_clocks()` once PhenoAge was always registered (4 built-ins). **Now computed from the registry**
  so it can't drift.

## 6. Redundant database queries / API calls

- No redundant DB/network calls found. `gene_report` does resolve the gene twice (once inside
  `gene_card`, once again for the ortholog list) — but that is an in-memory hashmap lookup on a
  singleton resolver, not a DB/API round-trip, so it is negligible and was left alone per the
  "smallest change" principle.
- The store's DuckDB access is partition-pruned (`_sig_glob`) and the connection is reused across
  queries (verified by an existing test) — already efficient.

## 7. Abandoned / disconnected files

- **`sources/methylation.py`** — a real but *unwired* adapter (not in `all_adapters()`, imported only
  by its test). It encodes genuine GEO-450K/biolearn validation work. *Recommendation:* either wire
  it into a surface or move it under a clearly-labelled validation path; **retained this pass** because
  deleting real validation code is not a clear win and carries more risk than value.
- **Dated/tier handoff docs** (`docs/HANDOFF_2026-08-0*.md`, `HANDOFF_PR1_REVAMP.md`,
  `HANDOFF_TIER0_REVIEW.md`, `HANDOFF_TIER1_SIGNATURES.md`) look superseded, but they are a
  deliberately-maintained session archive: cross-referenced from a code docstring (`api/service.py`),
  the canonical `HANDOFF.md` index, `ROADMAP.md`, and `VISUAL_AND_STATISTICAL_PLAN.md`. **Not deleted**
  — removing them would break live links for marginal gain. A future tidy could move them to
  `docs/archive/` and fix the references in one change.
- `harmonize/batch.py` and `effect_size.py` are library/test-only but are legitimate primitives (the
  latter is the differential test oracle). **Retained.**
- `sources/gtex.py` query methods are dormant (the adapter is a catalog entry in `/v1/sources` but its
  `median_expression`/`sample_expression` are called only in tests). **Retained** as a federated
  catalog entry.

## 8. Technical debt / opportunities

- **Environment vs pins (highest-impact):** the container carries `mcp 2.3.0` and `pandas 3.0.5`,
  both outside the project's deliberate pins. Addressed defensively so neither breaks CI-equivalent
  runs: the transport test now skips cleanly on an incompatible `mcp`, and a mypy override stops the
  SDK's version drift from failing type-checking. The clocks library self-skips under `pandas 3`. The
  code itself targets the pinned versions and is correct there. *The interpreter split (default
  `python` empty; deps under `python3.13`) is a container issue, not a repo one — not "fixed" here.*
- **`MetaSignature.standard_error` vs the Hartung-Knapp CI:** when HK widens the interval, `ci_low/
  ci_high` are HK-based but `standard_error` stays the DL value, so a consumer reconstructing
  `effect ± 1.96·se` gets a narrower band than the shipped CI. Documented in `models.py`; the UI uses
  the CI directly. Left as-is (changing it risks programmatic consumers); flagged here for awareness.
- **Bad-species error behaviour** has three shapes (400 for a metacharacter, 422 for an unknown
  symbol, silent `[]` for a clean-but-absent species). The partition guard itself is sound; unifying
  the UX is a low-severity follow-up.

## Deliberately retained (summary)

`*_demo` reference clocks, `sources/methylation.py`, the dated handoff archive, `harmonize/batch.py`
& `effect_size.py`, `gtex.py` catalog methods, and `gene_report`'s forward-compat slots — each is
either load-bearing, intentional, or cross-referenced, and removing it would cost more than it saves.

---

# Second pass (2026-11-13) — after the data & UI revamp

The revamp added ~3,000 lines (influence, Deming, KM/RMST, seven views, a
palette validator). This pass reviewed the whole tree again with tools rather
than by eye — vulture for unreferenced Python, an import graph against every
entry point (pyproject, Makefile, CI, Dockerfile), a dependency-import audit, an
unused-export / unused-CSS scan of the frontend, an execution of the tutorial
notebook — and an independent adversarial review of everything added since
`b0a24da`. Every deletion below was proven by a repository-wide grep (code,
tests, docs, notebooks, CI) to have zero references, and the full suite was
re-run after each step.

**Gate after the pass:** ruff, black, mypy clean · **533 Python tests** pass (2
environment skips) · **53 Playwright** pass (1 survival skip) · palette
validator passes in all three theme blocks · the production CSP produces zero
violations on all seven views in both themes.

## Bugs found and fixed (each with a test that failed first)

| Issue | Why it was wrong | Fix |
|---|---|---|
| **Cross-species Deming slope** (`build_frontend_data.cross_species`) | Assumed equal error in human and mouse (δ = 1). Mouse pooled estimates carry 3.6× the error variance, so the view headlined a slope of **3.13**; with the measured ratio it is **0.31 [0.15, 0.47]** at r = 0.15. | δ is measured from the pooled SEs and shipped with the correlation; the prose reports both. `deming.py`'s docstring no longer recommends δ = 1 for "same units". |
| Deming line drawn far outside the plot frame | No clip on the fitted line. | `clip: true`; a Playwright test checks it. |
| `deming_regression` on zero covariance | Returned slope 0 with a zero-width interval — a confident answer to an unanswerable question. | Raises `ValueError`. |
| `rmst` past the last follow-up | Silently carried the last survival value forward (all censored by t = 2, τ = 5 → RMST 5.0). Empty input gave a bare numpy error. | Both raise with a clear message. Shipped horizons (10, 15 y) are inside follow-up. |
| **Corpus-influence ranking** | Ranked by *total* shift, which correlates 0.75 with how many genes an array covers (1,007–20,897): it largely ranked array size. Ranked across species, it also split almost perfectly by species, because mouse genes are pooled from fewer, noisier contrasts. | Ranked by mean shift per gene measured, drawn as one panel per species on a shared scale. |
| False prose about the ranking | Panel said shifts were summed over "~42,000" genes (a contrast reaches ≤ 20,897); the Overview said GDS2973 shifts 20,897 genes "more than any other" (six contrasts tie at 20,897). | Prose states the computed coverage range and a within-species comparison. |
| Hard-coded "none of them disagree" | True today, typed rather than read — the failure mode fixed in the Mortality view one commit earlier. | Read from the shipped quadrant counts. |
| Methods text claimed survival curves | The shipped `crosslayer.json` is the full-cohort summary without curves, so the page described a chart it did not show. | The clause appears only when the curves ship. |
| Survival-curve colour ramp | Light Q1 1.39:1 and Q2 2.27:1, dark Q1 2.07:1 against the page — below the 3:1 floor for chart marks. The validator only audited the three verdict colours, so it reported ALL PASS. | New ramps (3.19–10.37:1 light, 4.33–12.60:1 dark); `validate_palette.js` now audits the ramp's contrast and lightness steps and fails the old values. |
| **Tutorial notebook** | Applied the hand-set `clinical_phenoage_demo` to real NHANES adults and printed a mean age acceleration of **+46 years** as a result. Nothing ever executed the notebook. | Uses the real PhenoAge clock on `clinical_nhanes_phenoage` (mean acceleration −0.70, r = 0.906 with age). `tests/test_notebook.py` runs every cell and pins calibration. |

## 1. Dead code — removed

| Item | Why it was dead | Risk checked |
|---|---|---|
| `harmonize.influence.corpus_influence` (+ its test, export) | No production caller; the exporter's `add_influence` re-implemented it and the two **disagreed** (it grouped by gene only, pooling human and mouse contrasts of shared Entrez ids). | The exporter's own test covers the planted-bad-contrast property. |
| `api.schemas.ErrorEnvelope` / `ErrorBody` | Never used as a response model or anywhere else. | Not referenced by OpenAPI routes. |
| `clocks.registry.PHENOAGE_FEATURES` | Unused alias of `phenoage.REQUIRED_FEATURES`. | — |
| `IdResolver.genes_in_group` | No caller, no test. | — |
| `sources.geo.POOLED_VARIABLES` | Enforced nowhere: geo.py pools *every* variable that is not stratified or controlled. Kept as an accurate comment. | Behaviour unchanged. |
| `sources.manifest.GEO_SERIES_BY_ACCESSION`, `GEO_PUBMED_BY_ACCESSION` | Unused derived dicts. | Series ids ship via `studies.json`. |
| `GtexSource.fetch_signature` | Docstring claimed a "federated query path" that does not exist. | — |
| `methylation.get_series`, `MethylationSeries.browser_url`, `GSE_BROWSER` | Unused; the module's real entry point is `read_series`. | The module itself is retained (see §7). |
| Unreachable `standard_error <= 0` guard in `add_influence` | `_pool` → `random_effects` rejects that input first. | — |
| `.result-species` CSS rule | Defined, never referenced. | — |
| R/Bioconductor "ETL bridge" text in `environment.yml` | Describes a bridge removed in the first pass; no R code exists. | — |

## 2. Duplicate logic — consolidated

- `CrossSpecies.jsx` re-implemented the exporter's quadrant rule → each pair now ships its `quadrant`.
- The 95% interval and direction were copied between `ForestPlot` and `TissueFacets` → one `contrastInterval` helper.
- `corpus_influence` vs `add_influence` → one aggregator (above).
- *Retained, intentional:* scalar vs matrix Hedges' g (the scalar is the test oracle).

## 3. Unused UI components

None: every component and every exported frontend function is imported, and every frontend dependency is used. Two classes (`at-risk`, `certainty-summary`) carry no styles by design — they are test hooks.

## 4. Overly complex implementations

- `LeaveOneOut` and `Panel` rebuilt arrays/objects every render, re-running their chart effects on every parent render → memoised.
- *Noted, not changed:* a few views still do display-level arithmetic on shipped numbers (an RMST difference, `c_clock < c_dysregulation`, a `p > 0.05` style threshold). This matches long-standing practice in the codebase and none of it is an estimator; shipping these as fields is a reasonable follow-up.

## 5. Legacy code

- The Streamlit dashboard and its `docker-compose` service stay **demoted, not deleted**, per the owner's decision.
- `*_demo` reference clocks stay: the test suite and the no-library fallback depend on them. They are no longer used to present results on real data.

## 6. Redundant queries / API calls

**One found and fixed.** `db.js` cached the Parquet tables but not its `json()`
helper, so a file read by two views was fetched twice: on an ordinary visit
`meta.json` (masthead and Overview), `studies.json` and `influence_corpus.json`
(Overview and Panel) each went over the network twice. `json()` now caches the
promise the same way the tables do, and never caches a failure. A Playwright test
counts requests and fails on any repeat. Every component was checked for
in-place mutation first, since the parsed objects are now shared; none mutates
them. The server side is unchanged: DuckDB access stays partition-pruned on a
reused connection.

## 7. Abandoned or disconnected files

| File | Verdict |
|---|---|
| `docs/PRD_GeroQuery_Multi-Omic_Aging_Aggregator (1).md` | The original spec, kept; the browser-download ` (1)` suffix is dropped and its two links updated. |
| `docs/WORKING_NOTES.md` | Marks itself "delete once Tier 0 ships", but its verified NHANES variable map and rebuild sanity numbers exist nowhere else, and two docs link to it. **Recommendation:** move those two sections into `DATA_SOURCES.md` / `RESULTS_NHANES_CSD.md`, repoint the links, then delete. Not deleted here: that would lose information. |
| `sources/methylation.py` | Imported only by its test, but it is the documented code path for a published result (Horvath reproducing the authors' ages, r = 0.998). Retained. |
| `evidence/`, `.streamlit/` | Live (the monthly living-evidence workflow; the demoted dashboard). |
| `requirements.txt`, `environment.yml` | Not redundant with `pyproject.toml`: the Dockerfile and the conda path consume them. |
| `pyarrow`, `uvicorn` (no imports) | Used implicitly: pandas' Parquet engine and the server the run command starts. |

## 8. Technical debt still open

- **Vite 5 → 8** to clear the remaining dev-server-only advisories (queued as a separate task; the static build is unaffected).
- **Survival curves on the shipped site** need a full-cohort ETL build.
- **§3.6 Schoenfeld diagnostics** remain unbuilt.
- **`WORKING_NOTES.md` migration** (above).
