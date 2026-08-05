# GeroQuery front end

A static site. No server, no API, no runtime dependency on anything but a CDN.

```
Python (offline)                    browser
─────────────────────────────       ──────────────────────────
GeroStore ─┐
           ├─ random_effects ──► pooled.parquet    ─┐
           │                     contrasts.parquet ─┤
survival/ ─┴─ Cox + LR tests  ──► crosslayer.json  ─┼─► hyparquet ─► React + Observable Plot
                                  studies.json     ─┤
                                  meta.json        ─┘
```

**Every statistic is computed in Python.** The browser filters, sorts, bins and
draws; it never pools, never fits, never derives an interval. That line is
deliberate — a DerSimonian–Laird re-implementation in JavaScript would be a
second estimator for one quantity, and two estimators drift. The repo already
keeps `test_matrix_hedges_g_matches_scalar_implementation` around because a
vectorized copy of one estimator was worth pinning to its reference; across a
language boundary the risk is worse and the tooling is thinner.

## Running it

```bash
# 1. Export the data (needs a built store; run `make idmap` at least once first,
#    or 94% of genes will show as Ensembl accessions instead of symbols)
python -m geroquery.etl.build_frontend_data

# 2. Dev
npm install
npm run dev            # http://localhost:5173

# 3. Production build
npm run build          # -> dist/, 8.2 MB total
```

Or from the repo root: `make frontend-data` then `make frontend`.

## Deploying to Cloudflare Pages

Free, unlimited bandwidth, no credit card, does not sleep.

| Setting | Value |
|---|---|
| Build command | `npm run build` |
| Build output directory | `dist` |
| Root directory | `frontend` |
| Node version | 20 or later |

`public/_headers` is picked up automatically and sets caching plus
`Accept-Ranges`. Nothing else is needed — there is no server to configure and no
environment variable to set.

The data files must exist before the build. Either commit
`frontend/public/data/` or run the export in CI ahead of `npm run build`; the
build does not generate them.

**GitHub Pages** works identically (`dist/` is plain static files). If deploying
to a subpath, set `base` in `vite.config.js`; `DATA_BASE` in `src/lib/db.js`
already honours `import.meta.env.BASE_URL`.

## Why not DuckDB-WASM

The first cut used it, for range requests against Parquet row groups. Cloudflare
Pages' 25 MiB per-file limit is what exposed the mistake: the DuckDB wasm
binaries are **34 MB and 39 MB**, roughly seven times the size of the 6.8 MB
corpus they existed to query. Range requests only pay when the data dwarfs the
engine.

hyparquet is ~15 kB and reads the same files. Measured, on the built site:

| | DuckDB-WASM | hyparquet, ranged | hyparquet, whole-file |
|---|---|---|---|
| Engine payload | 73 MB (both bundles) | ~15 kB | ~15 kB |
| Requests for a gene | — | 8 | **2** |
| Bytes transferred | — | 7.58 MB | **7.06 MB** |
| Largest file | 39 MB (**over the limit**) | 5.1 MB | **5.1 MB** |

The ranged reader was also worse than plain fetches here: because every row
group is read anyway, its widening Range GETs overlap and re-fetch. One
`fetch` per file wins on requests, on bytes, and on not depending on the host
honouring `Range` at all.

## Layout

```
src/
  App.jsx                 shell, tabs, masthead, footer
  lib/db.js               the two tables + every query
  lib/state.js            URL-as-state; the query string *is* the app state
  lib/theme.js            light/dark, and resolved tokens for SVG contexts
  lib/format.js           signed effects, intervals, p-value floors
  components/
    GeneView.jsx          search rail + readout + contrast table
    ForestPlot.jsx        the signature chart
    Landscape.jsx         the whole corpus, binned by verdict
    Panel.jsx             32 contrasts, with the selection rule printed
    Mortality.jsx         the cross-layer survival result
```

## Design notes

**The zero line is the protagonist.** This product's claim is that an interval
can say "we cannot tell", so x = 0 is a real structural element: drawn heavier
than the axis and every other rule, in every chart. The mortality view's
reference is 1.0 rather than 0 — hazard ratios are multiplicative — but it does
the same job and is drawn the same weight.

**Type encodes epistemics.** Instrument Serif appears only on assertions (the
verdict, empty-state titles). IBM Plex Sans is the interface. IBM Plex Mono is
every measurement, with tabular figures. Three families because there are three
kinds of statement on the page.

**Direction colours are equal in weight.** A gene rising with age is not better
or worse than one falling, so neither gets the "good" colour. The null is a
neutral warm grey, never red: *we cannot tell* is an answer, not an error, and
colouring it as failure would undercut the one thing this tool is for.

**Filled versus hollow is the verdict grammar.** A filled dot means the interval
excludes the null; a hollow ring means it crosses. It reads the same way in the
gene rail, the forest plot, and the mortality chart.

## Accessibility and quality floor

Keyboard focus is visible everywhere (`:focus-visible`, 2px). `prefers-reduced-motion`
disables the one animation. Charts carry `role="img"` and a label. Layout is
responsive to 390 px — the rail collapses above the readout and charts scroll
inside their own container rather than the page. Light and dark are both
first-class, with an explicit toggle that beats the system preference in both
directions.

## Known gaps

- Charts are redrawn on theme change by keying the view subtree on the theme.
  It is correct and slightly wasteful.
- The corpus loads whole (6.8 MB) on first query. Fine on desktop and on a
  decent phone connection; a 2G user waits. Splitting `contrasts.parquet` per
  species would halve it and is the obvious next step.
- No test suite. Playwright against the built `dist` is the right harness and
  is not written yet — worth doing, given the repo's history of UI bugs that a
  green Python suite did not catch.
