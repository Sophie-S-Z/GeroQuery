# GeroQuery — roadmap

**Date:** 2026-08-05 · **Base:** `3b8cc47` · Supersedes `HANDOFF.md` §10.

> **Status, 2026-08-06.** §2.4's "survey weights on the CSD result" is done for
> the *mortality* result — `WTDN4YR` with `SDMVSTRA`/`SDMVPSU`, so it is now a
> claim about 75,754,006 US adults aged 50+ rather than 2,517 people. The
> 2017–2018 CSD result is still unweighted and is the last one that is not.
> §3.4's "the bug class that produced three shipped UI failures gets a real
> harness" is also done: 26 Playwright tests against the built `dist`, in CI.
> The research behind what to build next is in
> [`HANDOFF_2026-08-06.md`](HANDOFF_2026-08-06.md) §6, which supersedes §2.4 and
> §2.5 here.
>
> **Status, 2026-08-05.** Items 1–4 are **built**.
> §2.1 → [`RESULTS_CROSSLAYER.md`](RESULTS_CROSSLAYER.md) ·
> §2.2 → [`../evidence/`](../evidence/) + `.github/workflows/living-evidence.yml` ·
> §2.3 → `geroquery/mcp/` ·
> §3 → [`../frontend/`](../frontend/README.md), Vite + React + Observable Plot,
> static on Cloudflare Pages.
>
> **One §3 conclusion was wrong and is corrected in place.** §3.1 argues for
> DuckDB-WASM over Parquet row groups. Cloudflare's 25 MiB per-file limit
> exposed the flaw: the DuckDB wasm binaries are 34 MB and 39 MB — about seven
> times the 6.8 MB corpus they would query. Range requests only pay when the
> data dwarfs the engine, and here it is inverted. The shipped app uses
> **hyparquet** (~15 kB) and fetches both files whole: 2 requests, 7.06 MB,
> versus 73 MB of engine that would not have deployed at all. See
> `frontend/README.md` for the measured comparison.

Two questions drove this: what makes GeroQuery durably useful rather than a
finished artifact, and what replaces Streamlit. They have one answer in common,
so they are answered together in §4.

---

## 1. The diagnosis

GeroQuery today is a **snapshot with excellent provenance**. Every number traces
to a checksum-pinned upstream, the estimator is validated against planted
effects and published values, and the nulls are published. That is rare and it
is worth something.

It is also, right now, a thing that happened on 2026-08-04.

Three properties separate a research artifact from research infrastructure, and
GeroQuery has none of them yet:

| Property | Has it? | What it means |
|---|---|---|
| **A heartbeat** | No | The corpus is frozen at 31 GEO DataSets selected on 2026-08-03. Nothing re-derives. In twelve months the pins are stale and the panel is a historical document. |
| **Reachability** | Partly | An API on `localhost:8010` that a human must clone and build. Nothing else in the world can call it. |
| **A citable result** | No | Three real results, none published anywhere a citation could attach to. Tools that get cited get used for years; tools that don't, don't. |

Everything in §2 exists to install one of those three.

The strategic point: **the aging field's stated problem is exactly the one
GeroQuery is architected for.** The Biomarkers of Aging Consortium named its own
blockers as "limited reproducibility across cohorts, inconsistent dataset
structures, and insufficient validation against hard clinical end points."
GeroQuery is a checksum-pinned, schema-harmonized, provenance-first pipeline
with no hard-outcome validation. Two of three are already solved. The third is
now cheap to solve (§2.1), and solving it is the whole play.

---

## 2. What to build, in priority order

### 2.1 · Land the cross-layer result — **now unblocked**

`HANDOFF.md` §10.3 calls this "the reason this exists" and lists the obstacle:
clocks and resilience live in different cohorts. **That obstacle no longer
exists.** NCHS released NHANES DNA methylation epigenetic biomarkers as a
public-use file on 2024-07-31.

Verified live on 2026-08-05:

| File | URL | Bytes | Access |
|---|---|---|---|
| DNAm epigenetic biomarkers | `wwwn.cdc.gov/nchs/data/nhanes/dnam/dnmepi.sas7bdat` | 1,245,184 | **Public** |
| Linked mortality 1999–2000 | `ftp.cdc.gov/pub/HEALTH_STATISTICS/NCHS/datalinkage/linked_mortality/NHANES_1999_2000_MORT_2019_PUBLIC.dat` | 487,666 | Public |
| Linked mortality 2001–2002 | same directory, `..._2001_2002_...` | — | Public |
| CpG-level betas | `DNAm-CpG-Sites.csv` | — | **RDC proposal required** |

What is in `dnmepi`: per-participant Horvath, Hannum, SkinBlood, PhenoAge,
GrimAge, GrimAge2, DunedinPoAm and HorvathTelo, for **n ≈ 2,532 adults aged ≥50**
in NHANES 1999–2002. Join key is `SEQN` — the same key `sources/nhanes.py`
already uses. The linked mortality files carry follow-up through 2019-12-31 on
the same key.

That means **all three GeroQuery layers land on the same subjects for the first
time**, plus a hard outcome:

```
dnmepi (clocks)  ──SEQN──  1999-2002 labs (resilience)  ──SEQN──  mortality (outcome)
```

Three questions become answerable, and each is worth reporting whichever way it
comes out:

1. Does clock age acceleration predict mortality here? (Replication — checks the
   join and the pipeline. Expect yes; if no, something is wrong.)
2. Does the resilience/CSD signal predict mortality? (**Novel.** No aging portal
   reports this.)
3. Does resilience add anything *over* age acceleration? (**The claim.** This is
   the one that justifies the resilience layer existing at all.)

**Two honest obstacles, both surmountable, both must be stated in the result:**

- The current CSD panel is built on NHANES **2017–2018**. The DNAm subsample is
  **1999–2002**. The six markers (albumin, creatinine, glucose, CRP,
  lymphocyte %, RDW) all exist in the 1999–2002 lab files, but under different
  file and variable names — this is a real adapter, not a config change. Budget
  it as such, and treat the existing `MARKER_MAP` warning in §12 of the handoff
  as binding: the new cycle gets its own dataset id.
- **Age ≥50 truncates the range.** The CSD design measures variance across age
  strata; with a floor at 50 and NHANES topcoding at 85 for these cycles, you
  have roughly half the dynamic range the 2017–2018 analysis had. This will
  reduce power and must be reported as a limitation, not discovered later.

Note that using clocks NCHS already computed sidesteps the entire
normalization-mismatch limitation in `RESULTS_METHYLATION_CLOCKS.md` — these are
the survey's own batch-corrected values. It also means this work does not need
`biolearn`, `pyaging`, torch, or the `.venv312` interpreter at all.

**Effort:** the largest single item here. A 1999–2002 NHANES adapter, a
mortality-file parser (fixed-width, documented codebook), a survival model, and
a results document. **Value:** this is the citable result. Everything else in
this roadmap is distribution for it.

---

### 2.2 · Give it a heartbeat — the living evidence loop

Continuously-updated evidence synthesis is now an established pattern with its
own literature (*JMIR* 2026; *Mayo Clin Proc Digit Health* 2024). GeroQuery is
unusually well-placed for it because the panel is **rule-selected, not
hand-curated** — the selection query is already in `WORKING_NOTES.md`:

```
age[Subset Variable Type] AND gds[Filter] AND (Homo sapiens[Organism] OR Mus musculus[Organism])
```

189 records on 2026-08-03. A live test already pins that this stays ≥180. That
test is one line away from being a discovery mechanism.

**Build:** a scheduled GitHub Actions workflow (free, unlimited minutes on
public repos) that monthly re-runs the selection query, fetches anything new,
rebuilds the panel, and — the part that matters — **diffs the pooled estimates
against the previous run** and opens a PR containing the diff.

The output is not "the data was refreshed." It is a changelog of
**effect sizes that moved**:

```
CDKN2A   g +0.07 [-0.20,+0.35] k=14  →  g +0.11 [-0.14,+0.36] k=16   (still null)
FOXO3    g -0.31 [-0.62,+0.01] k=7   →  g -0.38 [-0.65,-0.11] k=9    (newly excludes zero)
```

Nothing else in aging publishes that. Curated databases publish *entries*;
this publishes *how the evidence changed*, with the interval, which is the only
honest way to say a claim strengthened.

Guardrails, because an auto-updating corpus can silently corrupt itself:

- The PR must never auto-merge. A human approves every panel change.
- The diff must flag any dataset that *left* the panel, not only additions —
  a GEO record being withdrawn is a signal.
- Pin the run's own manifest digest into the results doc, so a reader can tell
  which corpus version produced a number.
- Keep the rule. **Do not let the loop become a place where datasets get
  hand-picked** — that is the one thing keeping selection bias out (`HANDOFF.md`
  §12).

**Effort:** small — one workflow file, one diff script, one results template.
**Value:** highest ratio in this document. It is the difference between a
snapshot and a resource, and it costs a weekend.

---

### 2.3 · Make it agent-reachable — an MCP server

**MCPmed** (*Briefings in Bioinformatics*, 2026-01-07, article `bbag076`) is an
explicit call for bioinformatics web services to expose Model Context Protocol
endpoints, demonstrated on GEO, STRING and the UCSC Cell Browser. Open Targets
shipped an official MCP server in 2026. The pattern the paper recommends is
exactly the shape GeroQuery already has: a clean typed API layer, plus a thin
MCP endpoint over it.

`api/service.py` is 406 lines of already-typed domain methods with Pydantic
schemas above them. An MCP server over it is a wrapper, not a rewrite:

| Tool | Wraps |
|---|---|
| `gene_signature` | `GET /v1/gene/{id}/signature` |
| `geneset_signature` | `POST /v1/geneset/signature` |
| `intervention` | `GET /v1/intervention/{name}` |
| `studies` | `GET /v1/studies` |
| `gene_card` | `GET /v1/gene/{id}/card` |

The reason this is high-value and not a checkbox: **GeroQuery's differentiator
is precisely what an LLM cannot do without it.** Ask a model "does CDKN2A change
with age" and it answers from the literature — which is how CDKN2A became the
most-cited aging gene. Ask GeroQuery and you get *g* = +0.07, 95% CI
[−0.20, +0.35], k=14. The tool exists to contradict the consensus with a
measurement, and the consensus is exactly what a language model is made of.

Make each tool return the interval and `k` in its **primary** response payload,
never as an optional field. A tool that can return "we cannot tell" is worth
more to an agent than one that always returns a verdict.

**Effort:** small. `mcp` Python SDK, one module, ~200 lines.
**Value:** turns GeroQuery from a site someone has to find into a capability
that shows up inside other people's work. This is the distribution channel that
did not exist two years ago.

---

### 2.4 · Widen the evidence base

In rough order of value per unit of work:

- **ARCHS4** (`archs4py`, HDF5, ~2M uniformly processed RNA-seq samples from
  GEO/SRA). This removes the microarray-only limitation *and* the uncorrected
  cross-study processing heterogeneity in one move, because ARCHS4 samples went
  through one pipeline. It is lighter than recount3 — one HDF5 file, no R
  bridge. This is the single biggest quality upgrade available to the signature
  layer.
- **A second omic layer.** `omic_layer` reads `transcriptome` on all 485,905
  rows. The schema is multi-omic; the data is one layer, which makes the column
  a promise rather than a fact. Methylation is the natural second layer and the
  ingestion path already exists.
- **Nested random effects by series** — 32 contrasts from 27 Series share
  subjects, so the current pooling slightly overstates independence. This is a
  correctness fix, not a feature, and it changes published numbers.
- **Ortholog mapping for HAGR** — recovers 2,276 dropped assertions for the cost
  of a join.
- **Survey weights on the CSD result** — converts it from a claim about 4,895
  people to a claim about the US adult population.

---

### 2.5 · Things worth doing that do not fit a theme

- Expose `gene_report`, `list_curated_genes`, `references` over HTTP — the
  dashboard reaches past the API for these, which means the API is not actually
  the complete interface it is documented as.
- `?library=` filter on `/v1/clocks` (240 unfiltered entries is not a usable
  response).
- Finish the GSE30870 sweep (36 short; the runner resumes).
- Re-add a *sourced* gene→hallmark mapping.
- Pool multiple NHANES cycles past the age-80 topcode.

---

## 3. Hosting and UI — the options

### 3.1 · The constraint that decides everything

**GeroQuery's data is 7 MB of Parquet and its hot query path is pure OLAP.**

```
signatures/species=human/...  3.59 MB
signatures/species=mouse/...  3.59 MB
metadata.sqlite               2.93 MB
curated_knowledge.csv         2.09 MB
clinical parquet ×3           0.32 MB
```

Even the git-ignored full signature table is 60.7 MB as CSV — as Parquet it is
well under the 25 MiB per-file limit of every static host worth using.

A gene lookup is `SELECT ... WHERE gene = ?` over a partitioned columnar file.
That does not need a server. It needs a CDN and DuckDB-WASM, which reads Parquet
over HTTP **range requests** — the browser downloads only the row groups the
query touches, not the file.

Which means the honest answer to "what should replace Streamlit" is not another
app framework. It is: **most of GeroQuery should not have a backend at all.**

What genuinely needs Python:

| Capability | Needs a server? |
|---|---|
| Gene / geneset signature lookup, forest plots, panel browse, ranking, cross-species | **No** — SQL over static Parquet |
| DerSimonian–Laird pooling | **No** — ~30 lines of JS; Hedges' *g* is already precomputed per contrast |
| PhenoAge on an uploaded cohort | **No** — a linear combination of nine terms |
| CSD / critical slowing down | **No** — variance and Kendall τ across strata; SQL plus a little JS |
| The 240-clock library (biolearn, pyaging, torch) | **Yes** — but it is a *batch* job, so precompute it in CI and ship the results as Parquet |
| Arbitrary user methylation matrix through an exotic clock | **Yes** — the only irreducible case |

So: static site, CI as the compute layer, and one optional thin API for the last
row.

### 3.2 · Free hosts, compared (verified 2026-08-05)

| Host | Free tier | Sleeps? | Fits GeroQuery |
|---|---|---|---|
| **Cloudflare Pages** | **Unlimited bandwidth**, 500 builds/mo, 20k files, 25 MiB/file, custom domain, no card | **No** | **Best.** Range requests supported; largest asset is 3.6 MB |
| **GitHub Pages** | Static, 1 GB site / 100 GB mo soft limits, no card | No | Fine. Zero new accounts if the repo is already there |
| **Google Cloud Run** | 2M req/mo, 360k GB-s, 180k vCPU-s, always-free, commercial OK | Scales to zero | Best option *if* a Python API is needed. Requires a card |
| **Koyeb** | 1 web service, no card | Yes | Workable fallback for the API |
| **Render** | 512 MB / 0.1 CPU | Yes, cold starts | 512 MB is tight for duckdb + pandas + pyarrow |
| **Hugging Face Spaces** | **Static only.** Gradio/Docker now require PRO | — | **No longer viable** for a Python backend. Was the obvious answer; is not anymore |
| **Fly.io** | No free tier for new users | — | No |
| **Railway** | $1/mo credit, pauses when spent | — | No |
| **Streamlit Community Cloud** | 1 GB RAM, sleeps | Yes | The current ceiling |

The Hugging Face change is worth flagging explicitly because it inverts advice
that was correct until recently: **Docker and Gradio Spaces now require a paid
plan.** Only Static Spaces remain free.

### 3.3 · Frontend options, compared

**Recommended — Observable Framework** (static site generator, v1.13.4,
March 2026):

- Client-side SQL over Parquet via DuckDB is **built in**, not bolted on. This
  is the exact primitive GeroQuery needs and the reason it beats a from-scratch
  React app on time-to-quality.
- Observable Plot produces charts that look considerably better than
  Streamlit+Plotly defaults without design work — and the forest plot,
  effect-size distribution and tissue comparison in `HANDOFF.md` §10.2 are all
  ~15 lines of Plot each.
- **Data loaders run in any language, including Python.** The existing
  `etl/` scripts become build-time loaders unchanged. Python stays in the repo;
  it just stops being on the request path.
- Builds to plain static files → drop on Cloudflare Pages.
- *Honest risk:* Observable deprecated its own `deploy` command and its 2026
  investment is visibly going to Canvas and Notebooks 2.0. Framework is
  maintained, not accelerating. Mitigated by the fact that the output is plain
  static HTML with no runtime lock-in — if Framework stalls, the built site
  keeps working and porting is mechanical.

**Alternative — Vite + React + Observable Plot + `@duckdb/duckdb-wasm`:** more
control, better if a distinct visual identity matters, more work. `frontend/`
already has this stack scaffolded, so it is not a cold start. Pick this over
Framework if the design ambition is high; pick Framework if shipping is.

**If you want to stay in Python** (lower effort, lower ceiling):

| Framework | The case for | The case against |
|---|---|---|
| **marimo** (`export html-wasm`) | Reactive, no rerun-everything model, exports to a **static** WASM site — Streamlit's hosting problem disappears entirely. Already installed in `.venv312` | 20–50 MB payload, slow first load, Pyodide package limits (numpy/scipy/pandas yes; biolearn no) |
| **Shiny for Python** | Genuinely reactive dependency graph — recomputes only what changed. `shinylive` also compiles to static WASM | Smaller Python ecosystem; reactive model is a real learning curve |
| **Reflex** | Compiles to Next.js; a real web app in Python | Needs a Node build and a server; heavier than the problem |
| **Panel / Dash** | Mature, flexible | Still server-bound, still a free-tier ceiling |
| **NiceGUI** | Pleasant for app-like UIs | Server-bound |

**marimo is the correct choice if the goal is to escape Streamlit this week.**
Framework or React is the correct choice if the goal is the tool being good.

### 3.4 · Why static beats Streamlit on every axis named

Mapping directly onto the gaps in `HANDOFF.md` §10.1:

| Gap | Static resolution |
|---|---|
| No mobile / narrow viewport | Real CSS. Streamlit's reflow is not a responsive strategy |
| No loading / empty / error states | Ordinary component states; a first-run walkthrough is just a route |
| Forest plot is the only visualisation | Observable Plot; tissue comparison, effect-size distribution and per-study drill-down are each a short block |
| No keyboard nav, unmeasured contrast | Standard a11y tooling applies; axe and Lighthouse run in CI |
| Streamlit is a ceiling | Removed. No rerun-on-interaction, no 1 GB RAM cap, no sleep |
| **No shareable state** | **URL is the state.** `?gene=CDKN2A&species=human&tissue=muscle` — this alone is the strongest single argument, because a research tool you cannot link into is a research tool nobody cites |

And a quieter benefit that matters given this repo's history: **the bug class
that produced three shipped UI failures gets a real harness.** `AppTest` was the
only thing that would have caught them, and it is a weak instrument. A static
site is testable with Playwright against a real DOM, in CI, on every commit.

### 3.5 · The one real risk, and how to handle it

Moving pooling, PhenoAge and CSD to JavaScript means **two implementations of
the same estimator**, which will drift.

Do not accept that. Use the pattern this repo already has: freeze inputs and
expected outputs as golden files, and assert the JS matches the Python reference
to a tolerance in CI. `test_matrix_hedges_g_matches_scalar_implementation`
already pins a vectorized copy to a reference implementation for exactly this
reason — this is the same test, across a language boundary.

If a JS estimator cannot be pinned to the Python one, **do not port it.** Route
it to the API instead. A wrong interval is worse than a slow one.

---

## 4. Where the two questions converge

The static architecture is not only the better UI answer. It is what makes §2.2
possible.

A living evidence loop needs somewhere to publish each rebuild. With a server,
every refresh is a deploy against a free tier that sleeps, throttles, or changes
its terms — and the corpus is only as available as the cheapest host's current
policy. With static output, **a rebuild is a commit**: CI re-derives the panel,
writes Parquet, and the CDN serves it. The site cannot go down, cannot exceed a
quota, and costs nothing at any traffic level.

Same for the MCP server: a static Parquet corpus on a CDN is directly queryable
by anything — a browser, a notebook, an agent, `duckdb` on someone's laptop
reading the URL. The API becomes one consumer of the corpus rather than the only
door to it.

The three things that make GeroQuery infrastructure — a heartbeat, reachability,
a citable result — are all cheaper to build on static output than on a server.
That is the actual reason to move, and the UI improvement is a consequence.

---

## 5. Suggested sequence

| # | Work | Rough effort | Why here |
|---|---|---|---|
| 1 | **NHANES DNAm + linked mortality → the cross-layer result** (§2.1) | Large | The citable result. Everything else distributes it. Do it while the finding is still novel |
| 2 | **Living evidence loop** (§2.2) | Small | Best ratio in the document. Independent of everything else — can run in parallel |
| 3 | **MCP server** (§2.3) | Small | Wraps an API that already exists. Ship before the front end; it needs no UI |
| 4 | **Static front end** on Cloudflare Pages (§3) | Medium–large | Fixes every §10.1 gap and unlocks §2.2's publishing path |
| 5 | **Retire Streamlit** | Small | Only after 4 reaches parity. Keep it running until then |
| 6 | **ARCHS4** (§2.4) | Medium | Biggest data-quality upgrade, but it changes published numbers — do it after the loop exists to diff it |

Items 1–3 are independent and can proceed in any order or at once. Item 4
depends on nothing but benefits from 1 being done first, so the new front end
launches with the new result rather than being retrofitted around it.

---

## 6. Sources

- MCPmed — *Brief Bioinform* 2026-01-07, `bbag076` — https://pubmed.ncbi.nlm.nih.gov/41729821/
- Open Targets MCP server — https://blog.opentargets.org/official-open-targets-mcp/
- Biomarkers of Aging Challenge — https://www.agingconsortium.org/challenge · *Nature Aging* — https://www.nature.com/articles/s43587-026-01139-6
- NHANES DNA methylation — https://wwwn.cdc.gov/Nchs/Nhanes/DNAm/Default.aspx
- NCHS Public-Use Linked Mortality Files — https://cdc.gov/nchs/data-linkage/mortality-public.htm
- Living evidence synthesis — *JMIR* 2026 — https://www.jmir.org/2026/1/e76130 · *Mayo Clin Proc Digit Health* — https://pmc.ncbi.nlm.nih.gov/articles/PMC11975841/
- ARCHS4 — https://archs4.org/ · `archs4py` — https://github.com/MaayanLab/archs4py
- Observable Framework — https://observablehq.com/framework · DuckDB support — https://observablehq.com/framework/lib/duckdb
- DuckDB-WASM — https://github.com/duckdb/duckdb-wasm
- marimo WASM export — https://docs.marimo.io/guides/exporting/webassembly_html/
- Cloudflare Workers/Pages pricing — https://developers.cloudflare.com/workers/platform/pricing/
- Cloud Run pricing — https://cloud.google.com/run/pricing
- HF Spaces overview (paid-plan requirement for Gradio/Docker) — https://huggingface.co/docs/hub/en/spaces-overview
- Free-tier survey 2026 — https://render.com/articles/platforms-with-a-real-free-tier-for-developers-in-2026
