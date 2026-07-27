# 🧬 GeroQuery

**A search engine for the biology of aging.** Type a gene and instantly see what
the research says about how it changes with age — across species and data types —
plus a real biological-age clock and a dynamical-systems resilience tool.

License: MIT · Python ≥3.10

---

## What it does

Aging researchers normally cross-reference a dozen separate databases to answer
one simple question: *“does this gene change with age, and is anything known to
affect it?”* GeroQuery puts that in one place, gene-first, with plain-English
answers next to the evidence.

Three tools, one interface:

1. **Gene explorer** — type a gene (symbol, alias, Ensembl/Entrez/UniProt ID) and
   get: whether it goes **up or down with age** and how confident that is; a
   plain-English analysis of *what the data tells us*; the tissues and data types
   where it’s established; its **hallmarks of aging**; cross-species (human ↔
   mouse) conservation; which curated databases (GenAge, CellAge, LongevityMap,
   OpenGenes) flag it; linked lifespan **interventions**; answers to common
   questions; and every claim’s **PubMed citation**.
2. **Aging clock** — the real, published **PhenoAge** model (Levine et al., 2018).
   Feed it nine routine blood markers plus age and it returns a **biological age**
   in years and, from the same model, a **10-year mortality risk** — a clean way
   to see that “biological age” and “mortality risk” are two faces of one model.
3. **Resilience** — early-warning signals of aging from dynamical-systems theory.
   As a body loses resilience, its biomarkers get more variable and more
   synchronised; GeroQuery surfaces this as age-stratified **dispersion**,
   **co-fluctuation**, and a **Dynamic Network Biomarker** composite — each with
   a bootstrap confidence interval. It is *inspired by* critical-slowing-down
   theory but honestly labelled a cross-sectional proxy, not validated CSD.

## Data honesty (read this)

GeroQuery **only uses real data**. Concretely:

- **Gene aging biology** — every direction, tissue, hallmark, and analysis is
  curated from the published literature and carries a **verifiable PubMed
  citation**. Evidence strength is an honest ordinal label (`robust` /
  `established` / `emerging`), not a false-precision number.
- **Curated database flags** — genuine memberships in the HAGR databases
  (GenAge, CellAge, LongevityMap) and OpenGenes.
- **Interventions** — real lifespan-affecting compounds/diets with approximate
  reported rodent effects and primary citations (NIA ITP, DrugAge).
- **The clock** — the published PhenoAge coefficients, implemented faithfully
  with correct unit conversions.

The **one** simulated element is the bundled **example biomarker cohort**
(`example_cohort_simulated`): 720 synthetic subjects that exist only so the clock
and resilience tools can be tried without uploading data. It is **clearly
labelled `SIMULATED` everywhere** and is never presented as measurements of real
people — upload-your-own is the first-class path.

> Earlier versions of this project shipped fabricated GEO study accessions,
> effect sizes, and p-values to fake a “demo meta-analysis.” **Those have been
> removed entirely.** GeroQuery no longer invents evidence.

## Quickstart

```bash
pip install -e ".[dev,ui]"        # or: pip install -r requirements.txt
python -m geroquery.etl.build_fixtures   # (re)generate the bundled data slice
python -m streamlit run geroquery/ui/streamlit_app.py   # the dashboard
```

The dashboard opens in your browser. Start on **Gene explorer** and try
`CDKN2A` (p16), the classic aging gene.

Other entry points:

```bash
make dashboard      # streamlit dashboard
make api            # FastAPI service at http://localhost:8000/docs
make test           # 79 tests, ~0.4s
make ci             # ruff + black + mypy + tests
docker compose up   # API on :8000, dashboard on :8501
```

## API (`/v1`)

| Endpoint | What it returns |
|---|---|
| `GET /v1/gene/{id}/report` (alias `/card`) | the assembled, cited aging profile of a gene |
| `GET /v1/genes` | the curated gene set (for browse / autocomplete) |
| `POST /v1/geneset/summary` | direction-of-change summary across a set of genes |
| `GET /v1/clocks` · `POST /v1/clock/apply` · `POST /v1/clock/compare` | PhenoAge metadata & application (dataset ID or uploaded matrix) |
| `GET /v1/intervention/{name}` | intervention record + linked genes + citations |
| `POST /v1/clock/diagnostics` | clock applicability check (completeness, ranges, unit hints, bootstrap CI) |
| `POST /v1/resilience/csd` · `POST /v1/resilience/recovery` | resilience metrics (dispersion, co-fluctuation, DNB, recovery rate) |
| `GET /v1/references` · `GET /v1/sources` · `GET /v1/datasets` · `GET /v1/version` | provenance & versioning |

Every error uses one envelope: `{"error": {"code", "message", "detail"}}`.

## Architecture

GeroQuery is organised around deep modules, each hiding complexity behind a small,
stable interface and tested in isolation. Lower modules never import higher ones.

```
ui  ->  api  ->  ( knowledge, clocks, resilience, store, harmonize )  ->  ( sources, idmap )
```

| Module | Responsibility |
|---|---|
| **`idmap`** | resolve any gene ID → canonical; UBERON tissues; AnAge fractional age |
| **`knowledge`** | the curated, cited aging knowledge base (the scientific source of truth) |
| **`sources`** | uniform cached/federated source adapters with a licence gate |
| **`harmonize`** | Hedges’ g, random-effects meta-analysis, batch correction (real algorithms for when real per-study data is supplied) |
| **`store`** | SQLite metadata + Parquet phenotype datasets + versioning |
| **`clocks`** | the real PhenoAge clock; a seam for `pyaging`/`biolearn` DNAm clocks |
| **`resilience`** | age-stratified dispersion & co-fluctuation indicators + DNB composite (bootstrap CIs), AR(1)/OU recovery rate, documented fallback |
| **`api`** | FastAPI `/v1`, error envelope, `format=json\|csv\|parquet` |
| **`ui`** | the Streamlit dashboard |

## Limitations (read these)

- **PhenoAge** is a research instrument for populations, not a medical diagnosis;
  results on your own data depend on assay units matching the expected columns.
- The bundled example cohort is **simulated** and labelled as such; it exists to
  demonstrate the tools, not to make claims about real people.
- The resilience module’s cross-sectional analysis is an **age-stratified proxy**
  (clearly labelled `fallback_used=True`); true critical slowing down is a
  longitudinal phenomenon.
- Curated database memberships and intervention effects are real but
  approximate; follow the citations to the primary sources before publication.
- Live identifier resolution (mygene.info) and real DNA-methylation clocks
  (`pyaging`/`biolearn`) activate automatically where network/those libraries are
  available; the bundled demo is fully offline.

## Roadmap & strategy

Where GeroQuery competes, what is deliberately deferred, and why, is documented in
[`docs/ROADMAP.md`](docs/ROADMAP.md) (with the full external assessment in
[`docs/STRATEGY_2026.md`](docs/STRATEGY_2026.md)). Short version: it wins on a
methodological-rigor layer — cited, uncertainty-quantified, honestly-scoped
answers — not on breadth.

## License

Code: **MIT**. Curated data: redistributable, attributed to the upstream sources
(HAGR, OpenGenes, NIA ITP, and the cited literature). Controlled sources
(UK Biobank, protected GTEx) are federate/link-only and never re-hosted —
enforced in code by each adapter’s licence gate. See [`LICENSE`](LICENSE).
