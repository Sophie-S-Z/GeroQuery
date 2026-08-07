# 🧬 GeroQuery

**Type a gene. Get its measured relationship to aging across 32 published
studies — with an interval wide enough to admit when the data cannot support a
claim.**

[![CI](https://github.com/Sophie-S-Z/GeroQuery/actions/workflows/ci.yml/badge.svg)](https://github.com/Sophie-S-Z/GeroQuery/actions)
· License: MIT · Python ≥3.10

---

## What this is, in plain terms

Ask a simple question — *does gene X change as people get older?* — and there is
no good place to look. The evidence is scattered across hundreds of studies in
different databases, using different gene names, tissues, age definitions and
measurement platforms. Existing aging databases are **curated lists**: someone
read the papers and wrote down a verdict.

GeroQuery **re-derives the answer from the raw data**, then tells you how
confident you are allowed to be.

It downloads real published datasets, verifies every byte against a checksum,
computes a standardized effect size per gene per study, pools them with a
random-effects meta-analysis, and hands back **a number with an interval**. When
that interval crosses zero, it says so.

That last part is the whole point. A curated list can only say *yes* or *no*. An
interval can say **"we cannot tell"** — and for a lot of famous aging genes, that
turns out to be the honest answer.

### What makes it unusual

1. **Nothing is asserted that was not measured.** Every layer traces to a
   checksum-pinned upstream. There is exactly one synthetic table and it exists
   to test an estimator, not to fill a gap.
2. **It publishes its own null results.** The headline finding is that the
   most-cited aging gene in the literature does not replicate here.
3. **It says how likely it is to be wrong about the direction**, not just whether
   an interval clears zero.
4. **The clocks were checked against published values**, not merely run.
5. **A resilience layer** — critical slowing down, borrowed from dynamical
   systems — that no other aging portal exposes.

---

## Table of contents

- [The four things it can tell you](#the-four-things-it-can-tell-you)
- [Results](#results)
- [Quickstart](#quickstart) — **start here**
- [Running each piece](#running-each-piece)
- [Rebuilding the data](#rebuilding-the-data)
- [Testing](#testing)
- [The API](#the-api)
- [The MCP server](#the-mcp-server-using-geroquery-from-an-ai-agent)
- [Deploying the site](#deploying-the-site)
- [How to read the numbers](#how-to-read-the-numbers)
- [What is real data and what is not](#what-is-real-data-and-what-is-not)
- [Limitations](#limitations-read-these)
- [Repository layout](#repository-layout)
- [Troubleshooting](#troubleshooting)
- [Documentation map](#documentation-map)

---

## The four things it can tell you

| Layer | The question it answers | Backed by |
|---|---|---|
| **Gene signatures** | Does this gene change with age, and how sure can we be? | 31 checksum-pinned GEO DataSets → **485,905 effect sizes** over 46,091 genes |
| **Aging clocks** | How old does this sample look biologically? | 240 clocks, validated against two real 450K methylation series |
| **Resilience** | Is this population losing the ability to recover from perturbation? | Real NHANES, plus a planted-effect control |
| **Survival** | Does any of it predict death? | NHANES 1999–2002 DNAm subsample, n=2,517, 1,350 deaths, 20 y follow-up, **survey-weighted to 75.8 M US adults aged 50+** |

---

## Results

All published as found, including the nulls.

### 1 · CDKN2A/p16 does not replicate

The most-cited transcriptional marker of aging: pooled *g* = **+0.07**, 95% CI
[−0.20, +0.35], k = 14, I² = 14%. Fourteen contrasts — six up, eight down — and
they agree with each other that there is nothing to see.

The same pipeline *does* find p21/CDKN1A (+1.07), IGF1 (−0.53), the
metallothioneins and a coherent mitochondrial decline. **So the estimator works
and the p16 null is a result, not a bug** — a pair of tests pins both halves,
because that pair is what stops a null being dismissed as a broken estimator.
→ [`docs/RESULTS_GEO_SIGNATURES.md`](docs/RESULTS_GEO_SIGNATURES.md)

### 2 · Three published clocks add nothing over six routine blood tests

Survey-weighted Cox over 20 years: GrimAge2 carries HR **1.527** [1.410, 1.652]
per SD, and the trained-on hierarchy reproduces (mortality > phenotype >
chronological age). But **SkinBlood, Lin and Weidner carry no evidence** over a
six-marker blood panel. *The comparison is the product* — not "here is GrimAge's
hazard ratio", but "here is what it adds, and for three clocks the answer is
nothing." → [`docs/RESULTS_CROSSLAYER.md`](docs/RESULTS_CROSSLAYER.md)

### 3 · Half of critical slowing down replicates

On real NHANES, health-state **variance** rises with age in 20/20 analytic
configurations; marker **cross-correlation** in 0/20. `resilience_declines`
returns `False`, because it requires both. The synthetic fixture returns `True` —
which is exactly what it was built to do.
→ [`docs/RESULTS_NHANES_CSD.md`](docs/RESULTS_NHANES_CSD.md)

### 4 · The clocks are correct, not just non-crashing

Horvath vs the study authors' own published per-sample values: **r = 0.9984, MAE
1.39 years.** Cord blood: **0.23 years**. DNAmTL correlates −0.967 with age —
telomeres shorten, and the sign is the check.
→ [`docs/RESULTS_METHYLATION_CLOCKS.md`](docs/RESULTS_METHYLATION_CLOCKS.md)

---

## Quickstart

**Prerequisites**

| | Needed for | Notes |
|---|---|---|
| **Python 3.10+** | everything | 3.13/3.14 fine for the core |
| **Python 3.10–3.12** | the 240 real clocks only | `ecos` has no 3.14 wheel — see [Troubleshooting](#troubleshooting) |
| **Node.js 18+** | the static web front end | |
| **git** | cloning | |

### The fastest path — no network, no downloads

Every table has a committed **real** sample, so the whole thing builds and runs
offline in about a minute.

```bash
git clone https://github.com/Sophie-S-Z/GeroQuery.git
cd GeroQuery

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -e ".[dev,ui]"
make data-offline                  # builds the store from committed real samples
make test                          # 436 tests
```

> **The committed slices are real measurements, not generated stand-ins.** They
> are *samples*, so a number computed on them is not the full-panel number — and
> every adapter **returns** which mode produced it rather than only logging it.
> The web front end prints `curated slice` in its masthead when that is what it
> is showing.

### Then look at something

```bash
# The static web app — the best surface, no backend needed
make frontend
cd frontend && npx vite preview     # http://localhost:4173

# Or the REST API
python -m uvicorn geroquery.api.app:app --port 8010
# http://localhost:8010/docs
```

---

## Running each piece

Ports are stated explicitly because 8000 and 8501 are so often already taken.

### The static web front end (recommended)

React + Vite + Observable Plot. **No backend at all** — it reads two Parquet
files in the browser with `hyparquet`, so it deploys to any static host for free.

```bash
make frontend                       # exports data + builds into frontend/dist
cd frontend && npx vite preview     # http://localhost:4173
```

For live-reload development:

```bash
cd frontend
npm install
npm run dev                         # http://localhost:5173
```

**Five views:**

| View | What it shows |
|---|---|
| **Gene** | One gene's readout: pooled effect, interval, prediction interval, and a forest plot of every contrast behind it |
| **Landscape** | The whole corpus binned by effect size, plus the strongest replicated effects |
| **Certainty** | Every pooled estimate as effect against certainty — where the large effects are also the least trustworthy |
| **Panel** | All 32 contrasts with full provenance: accession, Series, platform, PMID, restrictions applied |
| **Mortality** | The survey-weighted survival result, weighted and unweighted side by side |

The **URL is the state**, so any result can be linked or cited:
`?view=gene&gene=CDKN2A&species=human`

### The REST API

```bash
python -m uvicorn geroquery.api.app:app --port 8010
```

Interactive docs at **http://localhost:8010/docs**.

### The Streamlit dashboard (legacy)

```bash
python -m streamlit run geroquery/ui/streamlit_app.py --server.port 8502
```

> ⚠️ **This surface is stale and slated for retirement.** It is correct but
> behind: it shows neither the prediction interval, the survey-weighted mortality
> result, nor the local false sign rate. Prefer the static front end.

### Docker

```bash
docker compose up                   # API on :8000, dashboard on :8501
```

---

## Rebuilding the data

Everything above runs on committed samples. To rebuild from the real upstreams:

```bash
make data          # ~679 MB, ~30 min. Fetches + SHA-256 verifies every artifact,
                   # then rebuilds every table and the store.
```

**Network is off by default.** `GEROQUERY_ALLOW_NETWORK` gates it, so tests and
CI are hermetic by construction rather than by convention. The `make` targets set
it themselves where they need it.

Individual stages, in dependency order:

| Command | What it does |
|---|---|
| `make fetch` | Download + SHA-256 verify every pinned upstream (55 artifacts) |
| `make idmap` | Rebuild the bulk gene identifier maps from mygene.info |
| `make signatures` | GEO aging panel + HAGR → signatures, studies, curated knowledge |
| `make crosslayer` | NHANES 1999–2002 DNAm clocks + health state + mortality, joined |
| `make data-offline` | Build the store from committed samples — no network |
| `make frontend-data` | Export the Parquet/JSON the browser app reads |

> ⚠️ **`make idmap` must run before `make frontend-data`.** Otherwise 94% of
> genes render as raw Ensembl accessions and the search rail cannot find them.
> The exporter raises if *no* map exists, but it cannot tell a stale map from a
> complete one.

---

## Testing

```bash
make ci                             # ruff + mypy + the full Python suite
make test                           # pytest with coverage
make frontend-test                  # export, build, and run 34 Playwright tests
npx playwright install chromium     # once, before the first frontend-test

pip install -e ".[dev,mcp]"         # the MCP transport tests need the SDK
pytest tests/test_mcp_transport.py

GEROQUERY_ALLOW_NETWORK=1 pytest -m live   # live upstream drift checks
```

**Current gate:** ruff · black · mypy clean on 65 files · **436 Python tests
passed, 1 skipped, 15 live deselected** · **34 Playwright tests** (desktop +
393px, against the built `dist`).

Notes that will save you time:

- **`.[dev]` is not enough.** The dashboard tests skip silently without `ui`.
  Use `.[dev,ui]`.
- **Live tests are deselected by default.** They are the only thing that catches
  upstream drift, so run them occasionally.
- **Playwright runs against `dist`, never the dev server** — two of this repo's
  silent UI failures were defects in the exported Parquet, which no component
  test can see.
- **Never run two Playwright suites at once.** They race on the preview port.

---

## The API

Base path `/v1`. Every error uses one envelope:
`{"error": {"code", "message", "detail"}}`.

| Endpoint | What it returns |
|---|---|
| `GET /v1/gene/{id}/signature` | Pooled aging signature + every per-study contrast. Filters: `species`, `tissue`, `sex`, `omic_layer` |
| `GET /v1/gene/{id}/card` | Signature + curated evidence + linked interventions, assembled |
| `POST /v1/geneset/signature` | Pooled signature for a gene set or pathway |
| `GET /v1/clocks` | 240 clocks with predicted-outcome and training-population metadata |
| `POST /v1/clock/apply` · `POST /v1/clock/compare` | Apply or compare clocks on a dataset id or an uploaded matrix |
| `GET /v1/intervention/{name}` | Intervention record + lifespan data, mammals first |
| `POST /v1/resilience/csd` · `POST /v1/resilience/recovery` | Critical slowing down; AR(1) recovery rate |
| `GET /v1/studies` · `/sources` · `/datasets` · `/version` · `/healthz` | Provenance and versioning |

Add `?format=json|csv|parquet` to any data endpoint.

```bash
curl "http://localhost:8010/v1/gene/CDKN2A/signature?species=human" | jq
```

---

## The MCP server: using GeroQuery from an AI agent

GeroQuery speaks the **Model Context Protocol**, so an agent can call it as a
tool. Every primary payload carries the interval and an honest verdict — an agent
that reads `no_evidence` gets the prediction interval too, so it knows what it
would see if it ran the experiment.

```bash
pip install -e ".[dev,mcp]"
geroquery-mcp                       # stdio transport
```

> ⚠️ **`mcp` is capped below 2.0 on purpose.** `mcp` 2.x requires
> starlette ≥ 1.0, which FastAPI 0.115 cannot use — an unbounded install upgrades
> starlette and **every API import then dies**. Lifting the cap means moving
> FastAPI forward in the same change.

---

## Deploying the site

The front end is fully static, so it hosts free anywhere.

```bash
make frontend                       # -> frontend/dist
```

Deploy `frontend/dist`. `frontend/public/data/` is **committed** even though it
is derived: it is the deployable artifact, and hosts like Cloudflare Pages build
from the repo and cannot run the Python ETL. `meta.json` carries the data version
that identifies it.

Full walkthrough for Cloudflare Pages and GitHub Pages:
[`docs/DEPLOY.md`](docs/DEPLOY.md).

---

## How to read the numbers

This section is the difference between using the tool and misreading it.

### Hedges' *g*

A standardized difference between old and young. Roughly: 0.2 small, 0.5 medium,
0.8 large. Positive means the gene is **higher** in the older group.

### The three intervals, because they answer three questions

| | What it says |
|---|---|
| **95% CI** (reported) | Where the *mean* effect across studies is. Hartung-Knapp-Sidik-Jonkman, widened to DerSimonian-Laird wherever HK comes out narrower |
| **95% CI (DL)** | What almost every other tool would report, kept for comparison |
| **95% prediction interval** | Where **the next study's** effect would fall — typically about twice as wide |

Confusing the first with the third is the most common misreading of a forest
plot. Median CI width here is 0.928; median PI width is 1.831.

> **Why not plain DerSimonian-Laird?** DL substitutes an *estimated* τ² into the
> weights and then proceeds as though it were known, which makes the interval too
> narrow — worst when *k* is small, and 24% of this corpus pools fewer than ten
> contrasts. Switching to Hartung-Knapp **retracted 1,164 claims and created
> none.**

### The verdict

`increases` / `decreases` / `no_evidence`, derived from the reported interval —
never from the point estimate, and never from a precision you cannot see on the
page.

### The local false sign rate (lfsr) — read this one

**The probability that the *direction* is wrong.** It is not a p-value. A p-value
answers "could this be exactly zero", which nobody believes of any gene; the lfsr
answers "if I tell you it rises with age, how often am I wrong" — the question
actually being asked.

It is a posterior under a prior fitted across the whole corpus, so **it carries
the multiplicity that a per-gene interval cannot**. That is why:

> **3,707 genes have an interval excluding zero, and only 916 have lfsr ≤ 0.05.**

An interval is a statement about *one* gene. A ranked table is a selection over
41,983 of them. Ranking that selection by a per-gene interval is the winner's
curse — and it was shipping here: the mouse table used to lead with `Kif21a`,
*g* = +1.154, CI [0.180, 2.127], and a **46% posterior chance of pointing the
wrong way**. Tables now rank on the **shrunken** effect and show the lfsr beside
it.

### Hazard ratios

Per standard deviation, from a survey-weighted Cox model. HR 1.527 means 52.7%
higher instantaneous risk of death per SD of clock age acceleration, for the
**population** the sample represents — not for the 2,517 people measured.

---

## What is real data and what is not

The first thing you should read, not a footnote.

| Layer | Status |
|---|---|
| **Gene aging signatures** | **Real.** 31 checksum-pinned GEO DataSets → 32 young-vs-old contrasts from 27 independent Series → **485,905 Hedges' *g* estimates** over 46,091 genes, 10 tissues, human + mouse. A 40,585-row slice is committed so tests run offline |
| **Aging clocks** | **Real, and validated.** 240 clocks from `biolearn` (63) and `pyaging` (173) plus the real Levine PhenoAge coefficients. **436 clock-dataset runs**; Horvath reproduces the authors' own published ages to r = 0.998, MAE 1.39 y |
| **DNA methylation** | **Real.** GSE64495 (450K, n=106, ages 2.3–73.7) and GSE30870 (newborns vs nonagenarians), checksum-pinned |
| **Clinical / resilience** | **Real.** NHANES 2017–2018, 4,895 complete cases aged 20–80, verified against pinned SHA-256 checksums |
| **Mortality** | **Real.** NHANES 1999–2002 DNAm subsample linked to the National Death Index: 2,517 subjects, 1,350 deaths, 20 y follow-up |
| **Curated knowledge & interventions** | **Real.** Five HAGR databases, checksum-pinned: GenAge, CellAge, LongevityMap (including nulls), DrugAge, GenDR |
| **Gene identifiers** | **Real.** Bundled canonical table with live batched mygene.info resolution behind it |
| **Tissue expression** | **Real.** GTEx Portal API v2, queried live |
| **`clinical_synthetic_csd`** | **Synthetic, deliberately.** 720 subjects with critical slowing down planted by construction. It is the only way to ask whether an estimator recovers an effect *known* to be there. Own dataset id, `SYNTHETIC` prefix, and a test asserting it never merges with NHANES |

**Six guarantees enforced in code, not by convention:**

1. Nothing enters the store unverified — **cache hits are re-verified, not
   trusted**.
2. Offline is hermetic: a cache miss with network off raises.
3. Licence gate before persistence; controlled sources refuse to be cached.
4. Full-vs-sample mode is *returned*, not logged.
5. Real and synthetic never merge.
6. The data version is a content digest, not a filename hash.

Full cache-vs-federate and licence table:
[`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md).

---

## Limitations (read these)

**Gene signatures**

- **Microarray-era only.** GEO stopped curating DataSets around 2016, so there is
  no RNA-seq. ARCHS4 is the highest-value next source.
- **Small groups** — 3–15 samples per study, median CI width 0.84. This panel
  detects large, consistent, cross-tissue effects and **cannot rule out moderate
  ones.** A null means "not detected by this instrument", never "absent".
- **Tissue coverage is lopsided** — half the human contrasts are skeletal muscle;
  there is one blood contrast (n = 9).
- **`omic_layer` reads `transcriptome` for every row.** The schema is multi-omic;
  the data is one layer deep.
- **Cross-study processing heterogeneity is uncorrected** (12 platforms, three
  decades of normalization practice).
- **32 contrasts come from 27 Series**, so a few share subjects. `series_id` is
  carried on every study row so this is checkable; the pooling does not yet
  account for it.
- Sex is `unspecified` for 24 of 32 contrasts; ~22% of probe rows drop at
  identifier resolution.

**Clocks and methylation**

- Two blood 450K series only, n = 106 and n = 40.
- Normalization is not the original authors'.
- `GPAge*` cannot run — GPy needs numpy < 2, biolearn's cvxpy needs ≥ 2.

**Clinical**

- NHANES 2017–2018 is **cross-sectional**, so age strata are a proxy for
  within-individual trajectories; no relaxation time is observed.
- Age is **topcoded at 80**, compressing exactly the stratum where critical
  slowing down should be strongest.
- No fasting, medication, or comorbidity adjustment.
- The CSD result is still a claim about a **sample**; the mortality result is
  survey-weighted and is a claim about a **population**.

**Curated**

- Mammalian rows only — 2,276 of 5,029 HAGR assertions are yeast, worm or fly and
  are parsed, counted and dropped.
- DrugAge records no gene targets, so **none are asserted**.

**Operational**

- HAGR pins go stale by design.
- GTEx open cannot be age-stratified — donor age is in the dbGaP-controlled tier.
- No wet-lab or clinical claims. These are research indicators with stated
  assumptions.

---

## Repository layout

```
geroquery/
  idmap/          any identifier -> canonical Ensembl; UBERON tissues; fractional age
  sources/        adapters + checksum-verifying fetch + cache-vs-federate licence gate
  harmonize/      Hedges' g, BH FDR, random-effects pooling (HKSJ + prediction
                  interval), corpus-wide adaptive shrinkage and the lfsr
  store/          partitioned Parquet (DuckDB) + SQLite metadata + content-hash versioning
  clocks/         240 clocks, predicted-outcome metadata, age acceleration, mortality risk
  resilience/     critical slowing down, AR(1) relaxation, control energy
  knowledge/      PubMed-verified citations + hallmark vocabulary
  survival/       Cox PH, survey-weighted (Binder) path, Harrell's C, Mahalanobis dysregulation
  api/            FastAPI, error envelope, pagination, format switch
  mcp/            Model Context Protocol tools over the API layer
  etl/            offline batch only, never on the request path; the living-evidence panel_diff
  ui/             Streamlit dashboard (legacy)
frontend/         the static React/Vite site — the primary surface
  public/data/    the committed, deployable Parquet/JSON export
  tests/          34 Playwright tests against the built dist
tests/            436 Python tests
evidence/         the living-evidence baseline and changelog
docs/             results, handoffs, roadmap, deployment
```

Strictly layered: **lower modules never import higher ones**, so every layer is
testable with the layer below mocked.

---

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `pip install ".[clocks]"` fails on Python 3.13/3.14 | `ecos` has no wheel there. Use a 3.10–3.12 interpreter for clocks — a second venv is the intended setup |
| Every API import dies with `Router.__init__() got an unexpected keyword argument 'on_startup'` | Something lifted the `mcp < 2` cap and upgraded starlette. Reinstall with the cap |
| 94% of genes show as `ENSG…` accessions | `make idmap` did not run before `make frontend-data` |
| Dashboard tests all skip | Installed `.[dev]` without `ui`. Use `.[dev,ui]` |
| A cache-miss raises instead of downloading | Network is off by default. `GEROQUERY_ALLOW_NETWORK=1` |
| biolearn seems to hang for ~2 minutes on first use | Its gallery initialisation. Not a hang |
| Playwright fails in `gene.spec.js` for no obvious reason | Two suites raced on the preview port. Re-run one at a time |
| A number differs from CI | You ran `make data` and have the full 485,905-row table; CI has the 40,585-row committed slice. Tests pin `prefer_full=False` — do the same in any new test |
| `npm run build` output looks fine but Playwright tests a stale page | Check for a build error; `tail` can hide it. Grep for `error` |

**Do not install GPy** — it breaks all 63 biolearn clocks.

---

## Documentation map

| Document | What it is |
|---|---|
| [`docs/HANDOFF.md`](docs/HANDOFF.md) | **The authoritative document.** Full state, every bug found, every gotcha |
| [`docs/VISUAL_AND_STATISTICAL_PLAN.md`](docs/VISUAL_AND_STATISTICAL_PLAN.md) | What to build next in charts and statistics, and what **not** to build |
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | Where this goes next, and why |
| [`docs/OVERVIEW.md`](docs/OVERVIEW.md) | Source table, architecture, use cases |
| [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) | Cache-vs-federate contract and licences |
| [`docs/DEPLOY.md`](docs/DEPLOY.md) | Putting the static site online, free |
| [`docs/RESULTS_*.md`](docs/) | The four scientific results, in full |
| [`evidence/README.md`](evidence/README.md) | The living-evidence loop and the three rules it must not break |
| `docs/HANDOFF_*.md` | Per-session records |

---

## License

Code: **MIT**.

Bundled data are redistributable slices of real upstream releases; every upstream
keeps its own licence, attributed per adapter and recorded in
`geroquery/sources/manifest.py`. HAGR data is free for non-commercial use with
attribution. See [`LICENSE`](LICENSE).
