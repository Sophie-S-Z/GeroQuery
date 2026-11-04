# Handoff — data & UI revamp (2026-11)

A record of the revamp that built out the unbuilt phases of
`VISUAL_AND_STATISTICAL_PLAN.md`, surfaced the honest headline, added the palette
validator, and gave the static React app an overview landing. Everything below is
shipped on `claude/geroquery-revamp-data-ui-459dby` and is green under
`ruff`/`black`/`mypy`, the Python suite, the palette validator, and the Playwright
suite.

## What this revamp was

The product is the **static React explorer in `frontend/`** — React 18 + Vite 5 +
Observable Plot + hyparquet, reading two Parquet files and a handful of JSON in
the browser, no backend. The Python side is an ETL that turns real upstreams
(GEO transcriptome contrasts, NHANES clinical + DNAm, HAGR curated data,
pyaging/biolearn clocks) into precomputed statistics. The design spec
(`VISUAL_AND_STATISTICAL_PLAN.md`) had shipped only its Phase 1. This revamp
built the rest, under one unbroken rule: **every statistic is computed in Python
by the estimator the API serves and shipped in the Parquet/JSON; the browser
renders it and never derives it.** A new chart adds its numbers to the exporter
with a failing-if-wrong test *first*, then the JSX that draws them.

## New estimators (Python, test-first)

| Module | What it adds | Tests |
|---|---|---|
| `geroquery/harmonize/influence.py` | Leave-one-out for a random-effects pool (Viechtbauer & Cheung 2010): per-contrast pooled estimate with that contrast removed, Cook's distance, Baujat coords; `corpus_influence` aggregates each contrast's shift across every gene | `tests/test_influence.py` — planted-outlier identity, Cook's ranking, LOO = independent pool |
| `geroquery/harmonize/deming.py` | Deming / orthogonal (errors-in-both-variables) regression with a jackknife slope interval | `tests/test_deming.py` — exact recovery on a line, slope(y~x)·slope(x~y)=1, steeper-than-OLS under x-noise |
| `geroquery/survival/km.py` | Survey-weighted Kaplan-Meier and restricted mean survival time | `tests/test_km.py` — exponential closed form, unit-weight identity, KM = 1−ECDF, RMST = area under curve |

## Exporter (`geroquery/etl/build_frontend_data.py`)

Added, all pure functions of the committed inputs (no corpus rebuild needed):

- `ortholog_group` and `uniprot` on `pooled.parquet` (join on the bundled curated
  gene table) — the cross-species pairing and the AlphaFold out-link.
- `loo_g` / `loo_shift` on `contrasts.parquet`, and `influence_corpus.json` (32
  contrasts ranked by total leave-one-out shift).
- `crossspecies.json` — ortholog pairs, the Deming fit, concordance quadrants.
- `prediction_interval` in `meta.json` — **1,280** genes clear the prediction
  interval against **3,707** for the confidence interval (the honest headline).
- Per-quartile weighted survival curves + RMST in `crosslayer.json` (needs the
  full per-subject cohort; see "Rebuilding" below).

**`--enrich`**: `python -m geroquery.etl.build_frontend_data --enrich` layers the
derived fields onto the already-committed full-corpus artifacts in place, since a
plain `build` here would only produce the 40,585-row curated slice. A full build
produces the same fields natively, so the two paths cannot drift.
`tests/test_etl_signatures.py` checks both the unit behaviour and the committed
artifacts as data.

## New / changed frontend views

- **Overview** (new default tab) — leads with the prediction-interval headline, a
  certainty funnel (pooled → mean-excludes-zero → next-study-too), scope tiles,
  three clickable findings, a Methods panel, and "Running it". `Overview.jsx`.
- **Gene** — a **leave-one-out audit strip** beside the forest plot
  (`LeaveOneOut.jsx`), **tissue small-multiples** (`TissueFacets.jsx`), and an
  **AlphaFold structure out-link** via `uniprot`.
- **Cross-species** (new tab) — human vs mouse scatter with a y=x identity line
  and the precomputed Deming fit; the honest finding is that only 6 of 677
  orthologs resolve a direction in both species, and none disagree.
  `CrossSpecies.jsx`.
- **Panel** — a **corpus-influence** bar chart above the provenance table, with a
  "corpus shift" column joining the two. `CorpusInfluence.jsx`.
- **Mortality** — **weighted survival curves + RMST** by dysregulation quartile
  (`SurvivalCurves.jsx`), rendered when the cohort carries them.

## Palette validator

`frontend/validate_palette.js` is the `§1 rule 4` / `§6` validator the spec named
but never had: it reads every theme block in `styles.css`, computes OKLab ΔE under
normal and simulated protan/deutan/tritan vision (Machado 2009) plus WCAG
contrast, hard-gates the confusable verdict pairs at the §6 floors, and runs in CI
(`npm run validate:palette`). Writing it surfaced a real bug — the dark `--null`
neutral was fixed only on the toggle path, not the `@media(prefers-color-scheme)`
path, so system-dark readers still got the collapsed verdict colours; now fixed.

## Legacy Streamlit dashboard — demoted, not deleted

`geroquery/ui/streamlit_app.py` + `theme.py` are a second, divergent design
system on stale data. They are retained for local exploration and their tests but
marked deprecated in the module docstring and the README; **the static React site
is the product and new UI work goes to `frontend/`.**

## Rebuilding and running

```bash
# Python side (the estimators + API + MCP)
python3.13 -m pytest                    # 529 passed, survival/clocks skips
python3.13 -m ruff check geroquery tests && python3.13 -m black --check geroquery tests
python3.13 -m mypy geroquery

# Front-end data: full build (needs `make data && make idmap` for the real corpus
# and the full NHANES DNAm cohort for the survival curves), or enrich the
# committed artifacts in place:
python3.13 -m geroquery.etl.build_frontend_data            # full build
python3.13 -m geroquery.etl.build_frontend_data --enrich   # layer derived fields

# Front end
cd frontend
npm ci
node validate_palette.js                # colour, both themes
npm run build && npx playwright test     # CI uses the matching browser
```

The survival curves and RMST populate only when `crosslayer.json` carries the
per-subject cohort, which the full ETL rebuilds; the committed artifact is the
full-cohort **summary** (hazard ratios, concordance), so the Mortality view shows
the hazard-ratio result today and gains the curves on the next full build. The
view degrades cleanly and the Playwright test skips when they are absent.

## Deferred, deliberately

**§3.6 proportional-hazards (Schoenfeld) diagnostics** are the one visual-plan
item not built. They gate whether RMST is the headline or a convenience, need the
full per-subject cohort to be meaningful, and are a diagnostic rather than a
headline view. Left as the clear next piece of survival work.

## Commit trail (author: Sophie-S-Z, no AI attribution per the standing request)

`7d14db1` influence + cross-species/structure joins + PI headline + palette validator ·
`73b522f` leave-one-out strip + tissue small-multiples ·
`67dc179` corpus-influence view ·
`43e2e72` cross-species Deming scatter ·
`47f37a1` weighted survival curves + RMST ·
`0cd6b23` overview landing ·
`d8fcb93` AlphaFold out-link · and this docs/demote pass.
