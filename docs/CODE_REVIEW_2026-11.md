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
