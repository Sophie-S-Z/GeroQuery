# GeroQuery — complete handoff

**Date:** 2026-08-05 · **Branch:** `main`
**Manifest:** 2026.3 (checksums verified 2026-08-05)
**Gate:** ruff · black clean · **389 tests + 15 live**

> **2026-08-05 session.** Three things landed, all from
> [`ROADMAP.md`](ROADMAP.md) §2. §10 below is superseded by that document.
>
> 1. **The cross-layer result** — §10.3's "reason this exists" is no longer
>    blocked and is now answered. NCHS released per-participant DNAm clocks for
>    NHANES 1999–2002 in July 2024, so clocks, health state, and death now land
>    on the same 2,517 people. → [`RESULTS_CROSSLAYER.md`](RESULTS_CROSSLAYER.md)
> 2. **A living evidence loop** — monthly rebuild that publishes *which pooled
>    estimates moved*, never auto-merges. → [`../evidence/`](../evidence/)
> 3. **An MCP server** — GeroQuery as a tool an agent can call, with the
>    interval in every primary payload. → `geroquery/mcp/`
>
> Also fixed: CI installed `[dev]` without `ui`, so the ten `AppTest` cases in
> §8 — the only ones that would have caught bugs #9–#11 — were being skipped in
> CI and had never actually run there.

This is the authoritative document. §1 explains what GeroQuery is in plain terms;
§2 recaps this session; everything after is reference.

---

## 1. What GeroQuery is, in plain terms

### The problem

Ask a simple question — *does gene X change as people get older?* — and there is
no good place to look. The evidence is scattered across hundreds of published
studies in different databases, using different gene names, different tissues,
different age definitions, and different measurement platforms. Existing aging
databases are curated lists: someone read the papers and wrote down a verdict.
That is useful, and it is not the same as the measurement.

### What this does

GeroQuery **re-derives the answer from the raw data**, then tells you how
confident you are allowed to be.

It downloads real published datasets, verifies every byte against a checksum,
computes a standardized effect size per gene per study, pools them with a
random-effects meta-analysis, and hands back **a number with a confidence
interval**. When that interval crosses zero, it says so.

That last part is the whole point. A curated list can only tell you *yes* or
*no*. An interval can tell you *we cannot tell* — and for a lot of famous aging
genes, that turns out to be the honest answer.

### The one-line version

> Type a gene, get its measured relationship to aging across 32 published
> contrasts, with an interval wide enough to admit when the data cannot support
> a claim.

### What makes it unusual

1. **Nothing is asserted that was not measured.** Every layer traces to a
   checksum-pinned upstream. There is exactly one synthetic table, and it exists
   to test the estimator, not to fill a gap.
2. **It publishes its own null results.** The headline finding is that the most
   cited aging gene in the literature does not replicate here (§4).
3. **It reports when it cannot answer.** Wide intervals are shown, not hidden.
4. **The clocks were checked against published values**, not just run.
5. **The resilience layer** — critical slowing down, borrowed from dynamical
   systems — is something no other aging portal exposes.

### Who it is for

A researcher checking whether a gene is worth pursuing; someone evaluating an
aging-biomarker claim; anyone who wants the evidence *behind* an aging database
entry rather than the entry.

---

## 2. What happened this session

Four pieces of work, in order.

### 2.1 Replaced every synthetic data layer

The repo was a well-architected, fully-tested system over generated data:
twelve genes with hand-assigned effect sizes, fabricated GEO accessions
(`GSE90001`–`GSE90228`), synthetic p-values, and a clinical cohort labelled
"NHANES" that was a random draw.

| Layer | Before | After |
|---|---|---|
| Gene signatures | 12 genes, invented numbers | **485,905 effect sizes over 46,091 genes**, 31 pinned GEO DataSets → 32 contrasts → 27 Series |
| Clinical | 720 synthetic subjects | **4,895 real NHANES 2017–2018 adults** |
| Curated knowledge | 12 hand-written genes | **2,720 assertions over 1,880 genes** (5 HAGR databases) |
| Interventions | 7 hand-written | **1,340 records** (DrugAge + GenDR) |
| Identifier table | 24 hand-picked genes | **2,560 records / 1,835 ortholog groups**, generated |
| Species lifespans | 9 hand-transcribed | parsed from pinned AnAge |

### 2.2 Integrated a competing revamp

While this ran, a separate session merged **PR #1** into `main` — an independent
attempt at the same goal by the opposite route. Its environment had biomedical
APIs blocked, so it *deleted* the fabricated data and replaced it with a
hand-encoded knowledge module. Rather than discard either side, both were
integrated (§7).

### 2.3 Ran the clocks on real methylation

236 clocks had been wired and **never run on real data**. Two 450K blood series
were pinned and swept: **436 clock-dataset runs**. Horvath's clock reproduces the
study authors' own published per-sample ages to **r = 0.998, MAE 1.39 years**,
and puts cord blood at **0.23 years**.

### 2.4 Rebuilt the dashboard, and found it had been lying

You found three bugs I shipped. All of them got past a green 321-test suite
**because nothing ran the UI**.

- **Every gene showed "isn't in the curated aging knowledge base yet."** I removed
  the hand-written per-gene table but ported the UI without rewiring it; it still
  branched on a key that is now always `None`, so the entire explorer below that
  line was unreachable.
- **The resilience tab could never load a cohort.** A patch rewrote the branch
  conditions but not the radio's options, so the selection never matched.
- **Seven of twenty-eight PMIDs pointed at unrelated papers** — while the file
  asserted it never invents identifiers.

All fixed, plus the requested redesign, plus **ten `AppTest` cases** that are the
only tests here that would have caught any of it.

---

## 3. Features, elements, and use cases

### 3.1 The three layers

| Layer | What it answers | Backed by |
|---|---|---|
| **Gene signatures** | Does this gene change with age, and how sure can we be? | 31 GEO DataSets → 485,905 effect sizes |
| **Aging clocks** | How old does this sample look biologically? | 240 clocks; validated on 2 real methylation series |
| **Resilience** | Is this population losing the ability to recover from perturbation? | Real NHANES + a planted-effect control |
| **Survival** *(new)* | Does any of it predict death? | NHANES 1999–2002 DNAm subsample, n=2,517, 1,350 deaths, 20 y follow-up |

The fourth row is what changed on 2026-08-05. Until then every result here
validated a measurement against another measurement; this is the first hard
outcome. → [`RESULTS_CROSSLAYER.md`](RESULTS_CROSSLAYER.md)

### 3.2 Everything you can do today

| # | Use case | How |
|---|---|---|
| 1 | Look up one gene's measured age effect, with CI and heterogeneity | Dashboard *Gene explorer*, or `GET /v1/gene/{id}/signature` |
| 2 | See the per-study evidence behind that pooled number | Forest plot in the explorer; `signatures[]` in the API |
| 3 | Compare human vs mouse for the same ortholog | Both species pooled separately and shown together |
| 4 | Filter by tissue, sex, or omic layer | Query parameters on the signature endpoint |
| 5 | Get a gene's curated evidence + linked interventions in one call | `GET /v1/gene/{id}/card` |
| 6 | Browse 1,880 genes carrying curated aging evidence | Dashboard dropdown; `list_curated_genes()` |
| 7 | Ask what a compound does to lifespan, per organism | `GET /v1/intervention/{name}` — mammals first |
| 8 | Pool a whole gene set / pathway | `POST /v1/geneset/signature` |
| 9 | Screen for replicated aging genes | 689 human genes at BH q<0.05 across ≥3 contrasts |
| 10 | **Falsify** a claimed aging marker | The CDKN2A result is the worked example |
| 11 | Compute biological age on real clinical data | Dashboard *Aging clock*; PhenoAge on 4,894 NHANES adults |
| 12 | Get a 10-year mortality risk from the same fit | `mortality_risk_10yr` on the clock result |
| 13 | Upload your own cohort CSV for clock or resilience | Both tabs accept uploads |
| 14 | Measure critical-slowing-down early-warning signals | Dashboard *Resilience*; `POST /v1/resilience/csd` |
| 15 | Verify an estimator against a known planted answer | Run resilience on the synthetic fixture and compare |
| 16 | Trace any number to its source | `/v1/studies` gives accession, Series, platform, PMID, and the exact restrictions applied |
| 17 | Check clock accuracy against published values | `docs/RESULTS_METHYLATION_CLOCKS.md` |
| 18 | Run entirely offline | Committed real slices; network off by default |
| 19 | Export as JSON, CSV, or Parquet | `?format=` on data endpoints |
| 20 | Read the citation layer | 28 PubMed-verified references |

### 3.3 API surface (`/v1`)

`GET /gene/{id}/signature` · `/gene/{id}/card` · `POST /geneset/signature` ·
`GET /clocks` · `POST /clock/apply` · `/clock/compare` ·
`GET /intervention/{name}` · `POST /resilience/csd` · `/resilience/recovery` ·
`GET /studies` · `/sources` · `/datasets` · `/version` · `/healthz`

Not yet over HTTP (dashboard calls them directly): `gene_report`,
`list_curated_genes`, `references`.

### 3.4 Architecture

```
ui  ->  api  ->  ( store, clocks, resilience, harmonize )  ->  ( sources, idmap, knowledge )
```

Strictly layered; lower never imports higher. 54 modules, 18 test files.

| Module | Responsibility |
|---|---|
| `idmap` | Any identifier → canonical Ensembl; UBERON tissues; fractional age |
| `sources` | Adapters + checksum-verifying fetch + cache-vs-federate licence gate |
| `harmonize` | Hedges' *g*, per-gene contrasts + BH FDR, DerSimonian–Laird pooling |
| `store` | Partitioned Parquet (DuckDB) + SQLite metadata + content-hash versioning |
| `clocks` | 240 clocks; predicted-outcome metadata; age acceleration; mortality risk |
| `resilience` | Critical slowing down, AR(1) relaxation, control energy |
| `knowledge` | PubMed-verified citations + hallmark vocabulary |
| `survival` | Cox PH (Breslow ties), Harrell's C, nested LR tests, Mahalanobis dysregulation |
| `api` | FastAPI, error envelope, pagination, format switch, versioned LRU |
| `mcp` | Model Context Protocol tools over the API layer; `tools.py` imports no SDK |
| `etl` | Offline batch only — never on the request path. Includes the living-evidence `panel_diff` |
| `ui` | Streamlit dashboard (`theme.py` + `streamlit_app.py`), React showcase |

---

## 4. The three scientific results

All published as found, including the nulls.

**1 · On real NHANES**, health-state variance rises with age in 20/20 analytic
configurations; marker cross-correlation in 0/20. `resilience_declines` returns
`False`. → [`RESULTS_NHANES_CSD.md`](RESULTS_NHANES_CSD.md)

**2 · CDKN2A/p16 does not replicate.** Pooled *g* = +0.07, 95% CI
[−0.20, +0.35], k=14, I²=14% — fourteen contrasts, six up and eight down, and
they agree with each other that there is nothing to see. The synthetic slice this
replaced had p16 planted at *g* = +1.20, and the old test asserted every direction
was "up". Meanwhile p21/CDKN1A (+1.07), IGF1 (−0.53), the metallothioneins and a
mitochondrial decline all replicate — so the estimator works.
→ [`RESULTS_GEO_SIGNATURES.md`](RESULTS_GEO_SIGNATURES.md)

**3 · The clocks are correct, not just non-crashing.** Horvath vs the authors'
own published values: **r = 0.9984, MAE 1.39 y**. Against chronological age:
3.61 y where the literature says ~3.6. Cord blood: **0.23 years**. DNAmTL
correlates −0.967 with age — telomeres shorten; the sign is the check.
→ [`RESULTS_METHYLATION_CLOCKS.md`](RESULTS_METHYLATION_CLOCKS.md)

---

## 5. Data inventory

### 5.1 Pinned artifacts — 44 files, ~679 MB

| Group | Files | Size | What |
|---|---|---|---|
| GEO expression | 31 | 313.0 MB | GDS SOFT — the signature panel |
| GEO methylation | 2 | 358.2 MB | 450K series — the clock validation panel |
| NHANES 2017–2018 | 4 | 7.2 MB | DEMO_J, BIOPRO_J, HSCRP_J, CBC_J |
| HAGR | 6 | 0.2 MB | GenAge ×2, CellAge, LongevityMap, DrugAge, GenDR |
| AnAge | 1 | 0.2 MB | Maximum lifespans |

### 5.2 Live-queried

mygene.info (identifier harmonization, ~71k resolved at build) · GTEx Portal v2
(tissue context) · biolearn (63 clocks) · pyaging (173 clocks)

### 5.3 Derived and committed

`genes.json` (0.88 MB) · `anage.json` · `signatures_curated.csv` (5.05 MB) ·
`curated_knowledge.csv` (2.09 MB) · `interventions.csv` · `studies.csv` ·
`nhanes_2017_2018_sample.csv`. Git-ignored because reproducible:
`signatures_full.csv` (60.7 MB), `clinical_nhanes_full.csv`.

### 5.4 The one synthetic table

`clinical_synthetic_csd` — 720 subjects with critical slowing down planted by
construction. **Not a stand-in for missing data.** It is the only way to ask
whether the estimator recovers an effect *known* to be there, which no real
dataset can answer. Own dataset id, `SYNTHETIC` prefix, and a test asserting it
never merges with NHANES.

### 5.5 Guarantees enforced in code

1. Nothing enters the store unverified; **cache hits are re-verified, not trusted**.
2. Offline is hermetic — a cache miss with network off raises.
3. Licence gate before persistence; controlled sources refuse.
4. Full-vs-sample mode is *returned*, not logged.
5. Real and synthetic never merge.
6. The data version is a content digest.

---

## 6. Every bug found

Each was exposed by running real data or the real UI. Each would have passed a stub.

| # | Bug | Consequence |
|---|---|---|
| 1 | `idmap/mygene.py` took `ensembl[0]`, which can be a homology cross-reference | Mouse genes got flycatcher and ferret Ensembl ids as their canonical join key. **Affected the live request path** |
| 2 | mygene resolves symbol `CDKN1A` to an alt-scaffold MHC copy; Entrez `1026` gives the primary | **p21 had twelve real contrasts and the API reported none** |
| 3 | Alias list truncated at 12, sorted alphabetically | `p16` fell off human CDKN2A; the resolver answered with the mouse ortholog |
| 4 | Download cache keyed on URL basename; HAGR serves two files as `/dataset.zip` | They shared a cache entry. The checksum caught it — that is what it is for |
| 5 | Data version hashed Parquet filenames and sizes; DuckDB writes partitions in parallel | Identical data could report different versions |
| 6 | Betas passed CpG-by-sample, and every gapped CpG deleted — including 3 Horvath coefficients | **236 clocks failed** with "features not present" |
| 7 | `/v1/intervention/{name}` returned `matches[0]` | Rapamycin answered with the *C. elegans* result, not mouse |
| 8 | A test relied on a git-ignored file being absent | Passed only if you had never run `make data` |
| 9 | Dashboard branched on `knowledge`, which is now always `None` | **Every gene showed "not in the knowledge base"**; the explorer was unreachable |
| 10 | Resilience radio options rewritten, branches not | **No cohort could ever load** |
| 11 | 7 of 28 PMIDs pointed at unrelated papers | A lamin B1 citation resolved to breast-cancer macrophages; PhenoAge to contraceptive implants |
| 12 | NHANES `XY_Estimation` is 1=female/2=male — the **inverse** of `RIAGENDR` | Both are valid-looking 1/2 integers, so the naive map gave 2.4% agreement instead of 97.6%: **DNAm sex was wrong for every subject**, silently |
| 13 | NHANES variable names are not stable across cycles (`LBXSCR`→`LBDSCR`, `LBXSAPSI`→`LBDSAPSI`) | Using one cycle's names for both does not raise — it drops 1,306 subjects and yields a half-size cohort that looks fine |
| 14 | MCP tools read `effect_size`/`intervention_type`; the model exposes `lifespan_effect_pct`/`itype` | Every rapamycin effect came back `None`, reading as "a drug with no measured effect" rather than as a mapping bug |
| 15 | CI installed `.[dev]`, not `.[dev,ui]` | `importorskip("streamlit")` skipped all ten `AppTest` cases. **The guard against bugs #9–#11 had never run in CI** |

Also fixed earlier: negative AR(1) coefficients reported as *high resilience*;
CSD accepting a positive slope alone (fires on noise ~half the time);
`resilience/service.py` had a live `NameError` on every `control_energy` call.

---

## 7. The divergence with PR #1, resolved

`origin/main` was independently revamped toward "no fabricated data" by encoding
cited knowledge by hand, because its environment had APIs blocked. This branch
ingested the real sources instead. Integration kept both sides' real work.

**Carried forward:** `clocks/phenoage.py` (real Levine coefficients) · the
Streamlit redesign · `knowledge/references.py` + hallmark vocabulary ·
`mortality_risk_10yr` (bare `except` narrowed) · `CanonicalGene.ortholog_group`.

**Not carried:** `knowledge/aging_knowledge.py` — its hand-written per-gene
evidence asserts what the panel now measures, and disagrees for CDKN2A. Also
`main`'s `api/app.py`, which drops three endpoints.

> ⚠️ Git merged five files **without conflict** because this branch had not
> touched the same lines — silently taking `main`'s versions, including the
> removal of `/gene/signature`, `/geneset/signature` and `/studies`.
> **"No conflict" is not "no decision."**

---

## 8. Quality state

```
ruff / black / mypy   clean (54 modules)
pytest                331 passed, 1 skipped (13 live deselected)
coverage              86%
```

| Well covered | Not covered, and why |
|---|---|
| `geo.py` 97% · `store.py` 95% · `methylation.py` 94% · `differential.py` 94% · `streamlit_app.py` 87% · `theme.py` 100% | `etl/build_idmap.py`, `build_fixtures.py`, `fetch_artifacts.py` **0%** — offline scripts run by `make`, never on a request path |

**18 test files.** The ones that encode judgement rather than mechanics:

- `test_p16_does_not_replicate` **and** `test_p21_does_replicate` — the pair is
  what stops a null being dismissed as a broken estimator.
- `test_matrix_hedges_g_matches_scalar_implementation` — pins the vectorized copy
  to the reference.
- `test_resilience_runs_and_separates_the_two_cohorts` — the estimator must find
  the planted effect *and* decline to claim it on real data.
- `test_live_every_pmid_resolves_to_the_paper_we_claim` — bug #11.
- `test_masthead_counts_are_computed_not_typed` — a stale corpus size on a page
  whose argument is provenance is a self-inflicted wound.

---

## 9. Running it

```bash
pip install -e ".[dev,ui]"          # add ",clocks" on Python 3.10-3.12

python -c "from geroquery.store import GeroStore; GeroStore().build()"
python -m uvicorn geroquery.api.app:app --port 8010
python -m streamlit run geroquery/ui/streamlit_app.py --server.port 8502
```

Dashboard **:8502** · API docs **:8010/docs**. Ports are explicit because 8000
and 8501 are often already taken.

Rebuild data: `make data` (~679 MB, ~30 min) · `make signatures` · `make idmap` ·
`make data-offline` (no network) · `make ci` · `GEROQUERY_ALLOW_NETWORK=1 pytest -m live`.

**Sanity numbers:** 31 datasets / 32 contrasts / 27 Series · 485,905 rows over
46,091 genes · 2,720 assertions over 1,880 genes · 1,340 interventions · NHANES
4,895 and 4,894 · 44 artifacts · rapamycin in mouse = +13.0% · Horvath vs
published MAE 1.39 y.

---

## 10. Next steps

### 10.1 The UI/UX is not finished — start here

The dashboard was rebuilt this session to fix three functional bugs and apply the
typography and dark-mode brief. **It is now correct, not yet excellent.** Known
gaps, roughly in order of value:

1. **No mobile or narrow-viewport work.** Streamlit's default reflow is the only
   thing keeping it usable below ~900px. The forest plot in particular does not
   degrade well.
2. **No loading, empty, or error states worth the name.** A gene with no data
   gets a plain panel; a slow query gets Streamlit's default spinner. First-run
   has no orientation at all — a new visitor sees CDKN1A and no explanation of
   why they should care.
3. **The forest plot is the only real visualisation.** Tissue comparison,
   effect-size distribution across the panel, and per-study drill-down are all
   absent. There is no way to *see* the 485,905 rows, only to query into them.
4. **No keyboard navigation or focus management**, and the contrast ratios were
   designed but never measured with a tool.
5. **Streamlit is a ceiling.** Tabs, radios and the rerun-on-every-interaction
   model fight a research tool that wants persistent state and deep linking.
   A React or Observable front end over the existing `/v1` API would remove that
   ceiling — the API is already complete enough to serve it, and `frontend/`
   exists as a stub.
6. **No shareable state.** You cannot link someone to a gene result.
7. **Motion, micro-interaction, and the small craft details** were out of scope
   this session and remain so.

### 10.2 Make the tool do more

The data is broader than the interface. Capabilities that exist in the store but
have no way in:

- **Compare a gene set / pathway visually** — `geneset_signature` is implemented
  and exposed over HTTP but has no UI.
- **Tissue-by-tissue view of one gene** — the data is per-tissue; the UI pools it.
- **Browse the panel itself** — 32 contrasts with full provenance, currently only
  reachable as JSON.
- **Rank and filter genes** — "show me every gene with q<0.05 that falls in
  muscle" is one DuckDB query and no interface.
- **Cross-species conservation view** — both species are computed; nothing
  displays the comparison as a comparison.
- **Clock comparison** — `/v1/clock/compare` exists; the dashboard runs one clock.
- **Intervention → gene → evidence path** — the joins exist; the journey does not.
- **Saved queries / export a result set** — no session persistence at all.

### 10.3 Correlate age acceleration against resilience — *the reason this exists*

Now unblocked. Clock age minus chronological age is computable on real samples;
resilience metrics are computable on real NHANES. **This cross-layer claim is
still the one headline claim with no evidence behind it in either direction.**

The obstacle: they live in different cohorts. In increasing order of effort —
find a GEO series carrying both methylation and clinical chemistry; pursue NHANES
DNA methylation subsamples; or apply for a cohort with both.

### 10.4 Widen the evidence base

- **recount3 or ARCHS4** — removes the microarray-only limitation *and* the
  uncorrected cross-study processing heterogeneity in one move. ARCHS4 is
  lighter: one HDF5, no R bridge.
- **A proteomic or methylation signature layer** — `omic_layer` currently reads
  `transcriptome` for every row. The schema is multi-omic; the data is one layer.
- **Longitudinal resilience** — `resilience/recovery.py` is implemented, tested,
  and has never run on real longitudinal data. NHANES III linked mortality is
  free; UK Biobank is ideal and the adapter exists.

### 10.5 Smaller, worth doing

- Finish the GSE30870 clock sweep (36 short; the runner resumes).
- Nested random effects by series — 32 contrasts from 27 Series share subjects.
- Expose `gene_report`, `list_curated_genes`, `references` over HTTP.
- Apply survey weights — makes the CSD result a population claim.
- Ortholog mapping for HAGR — recovers 2,276 dropped assertions.
- `?library=` filter on `/v1/clocks` (240 entries unfiltered).
- Re-add a *sourced* gene→hallmark mapping.
- Pool multiple NHANES cycles past the age-80 topcode.

---

## 11. Limitations

**Signatures:** microarray-era only, no RNA-seq · groups of 3–15, median 95% CI
width **0.84** — detects large consistent effects, cannot rule out moderate ones
· half the human contrasts are skeletal muscle · one blood contrast (n=9) ·
`omic_layer` is transcriptome everywhere · cross-study processing heterogeneity
uncorrected · 32 contrasts from 27 Series share some subjects · sex
`unspecified` for 24 of 32 · ~22% of probe rows dropped at resolution.

**Methylation/clocks:** two blood 450K series only · n=106 and n=40 · GSE30870
swept to 200/236 · `GPAge*` cannot run (GPy needs numpy<2, biolearn's cvxpy needs
≥2) · normalization is not the authors' · no mortality outcome.

**Clinical:** NHANES is cross-sectional — no relaxation time observed · age
topcoded at 80 · survey weights carried but not applied · no fasting/medication
adjustment.

**Curated:** mammalian rows only (2,276 dropped) · DrugAge has no gene links so
none are asserted · gene→hallmark annotation dropped.

**Operational:** HAGR pins go stale by design · GTEx open cannot be
age-stratified · React frontend not repointed.

---

## 12. Gotchas

- **`GEROQUERY_ALLOW_NETWORK` is off by default.**
- **Two interpreters on purpose** — 3.14 for the core, `.venv312` for clocks.
- **Do not install GPy.** It breaks all 63 biolearn clocks.
- **`signatures_full.csv` is gitignored, `signatures_curated.csv` is committed.**
  A developer who ran `make data` sees different counts from CI. **Tests pin
  `prefer_full=False`** — do the same in any new test.
- **Never widen `nhanes.MARKER_MAP`** — the resilience service infers biomarkers
  by exclusion, so a new column silently joins the health state. Extra markers go
  in `PHENOAGE_ONLY_MARKER_MAP` with their own dataset id.
- **Do not hand-edit the GEO panel** — it is rule-selected; that is what keeps
  selection bias out.
- **`get_registry()` memoizes** — reset `clocks.registry._REGISTRY` in tests.
- **HTTP 200 from Streamlit proves nothing.** It renders over websocket; the
  shell loads regardless. Use `AppTest` — that is how all three UI bugs hid.
- **Live tests are excluded by default** — the only thing catching upstream drift.
- **biolearn's gallery init takes ~2 minutes.** Not a hang.
- **`fetch.cache_path` is prefixed with the manifest key.** Do not "simplify" it.
- **NHANES URL trap:** `/Nchs/Nhanes/2017-2018/DEMO_J.XPT` returns HTML. Use
  `/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt`. A test pins it.
- **`XY_Estimation` is 1=female, 2=male.** The opposite of `RIAGENDR`. Never
  reuse `SEX_LABELS` for it; use `DNAM_SEX_LABELS`. Bug #12.
- **NHANES variable names differ between cycles.** Never assume a 2017-2018
  name works in 1999-2002. Check the file's actual columns first. Bug #13.
- **The two NHANES cohorts are not interchangeable.** 2017-2018 is age 20-80
  with hs-CRP; 1999-2002 is age 50-85 with standard CRP, and the age floor
  removes the CSD variance signal entirely (`RESULTS_CROSSLAYER.md` §5).
- **`geroquery/mcp/tools.py` must never import the MCP SDK.** That split is what
  lets the payload logic be tested without the extra installed.
- **Never let the living-evidence workflow merge itself,** and never hand-edit
  the GEO panel to accept a diff. A test asserts the first; only discipline
  enforces the second.
- **`.[dev]` is not enough to run the suite honestly** — the dashboard tests
  skip without `ui`. Use `.[dev,ui]`. Bug #15.

---

## 13. Commit history

```
8194c4a  chore: drop a stray empty root package-lock.json
996bed0  fix: rebuild the dashboard on measured data; correct seven wrong PMIDs
9b6f8f7  docs: single authoritative handoff
ef1dfd9  docs: correct the methylation sweep coverage
bb9a6da  docs: session handoff for the integration and methylation work
cd2ef49  merge: reconcile with the PR #1 revamp on main
10189b2  docs: methylation clock validation results
171b840  test: methylation parsing and the PhenoAge clock
6c91cf0  feat: integrate the main-branch revamp; ingest real methylation;
         source the last hand-written tables
4757514  docs: session handoff for the real-data ingestion work
e0c3301  feat: replace every synthetic data layer with checksum-pinned real sources
```

All authored by Sophie Zhang. `real-data-ingestion` retained at `e0c3301` as the
pre-merge state.

---

## 14. Document map

| Document | What it is |
|---|---|
| **`HANDOFF.md`** | This document. Authoritative |
| **`ROADMAP.md`** | Where this goes next, and why. Supersedes §10 |
| `RESULTS_CROSSLAYER.md` | Clocks × health state × mortality in one cohort |
| `../evidence/README.md` | The living-evidence loop and the three rules it must not break |
| `OVERVIEW.md` | Source table, architecture, use cases, limitations |
| `DATA_SOURCES.md` | Cache-vs-federate contract and licences |
| `RESULTS_NHANES_CSD.md` | Clinical resilience result (partial null) |
| `RESULTS_GEO_SIGNATURES.md` | Gene signature result (p16 does not replicate) |
| `RESULTS_METHYLATION_CLOCKS.md` | Clock validation result |
| `HANDOFF_2026-08-03.md`, `HANDOFF_2026-08-04.md` | Session records |
| `HANDOFF_TIER0_REVIEW.md`, `HANDOFF_TIER1_SIGNATURES.md`, `HANDOFF_PR1_REVAMP.md` | Historical |
| `WORKING_NOTES.md` | Verified URLs, variable maps, rebuild numbers |
