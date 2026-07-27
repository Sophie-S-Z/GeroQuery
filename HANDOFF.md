# GeroQuery — Revamp Handoff

_Last updated: 2026-07-27 · Branch: `claude/geroquery-revamp-data-ui-459dby` · PR: [#1](https://github.com/Sophie-S-Z/GeroQuery/pull/1)_

This document is a complete record of the GeroQuery revamp so anyone picking this
up with fresh context can continue confidently. It covers the goal, every change
made, the file-by-file map, how to run and verify, the key design decisions, and
what could come next.

---

## 1. The goal

Revamp GeroQuery so that it:

1. **Uses only real data** — no fabricated elements or false evidence.
2. Has a **significantly upgraded UI/UX** (especially the UI): cleaner, sleeker,
   more polished, easier to use.
3. Gives **plain-English descriptions/analysis** of what the data tells us,
   alongside (not instead of) the statistics.
4. Ships **exceptional figures, captions, and answers** to as many pertinent
   questions as possible.

GeroQuery is a portfolio piece and a genuine research tool: a gene-first "search
engine for the biology of aging," plus a biological-age clock and a resilience
tool driven by biomarker data.

---

## 2. Starting state (what was wrong)

The repo was a mature, well-architected, fully-tested project (77 tests, 8 deep
modules). The problem was concentrated in the **data layer**:

- `geroquery/etl/build_fixtures.py` **generated fabricated evidence**: fake GEO
  accession numbers (`GSE90001`–`GSE90228`), synthetic effect sizes with random
  scatter, and synthetic p-values — then presented them through a real
  random-effects meta-analysis as if they were measured data (`signatures.csv`,
  `studies.csv`).
- The clinical cohort (`clinical_nhanes_slice.csv`) was **synthetic but
  mislabeled "NHANES."**
- The clocks (`clinical_phenoage_demo`, `clinical_mortality_demo`) were
  transparent linear toys with invented coefficients.

The genuinely real assets already present: gene identities (`idmap/data/genes.json`),
curated database memberships (GenAge/CellAge/LongevityMap/OpenGenes), intervention
lifespan effects, and the **directions** of aging change (established biology).

**Network constraint discovered:** outbound access to biomedical APIs
(mygene.info, NCBI) is policy-blocked in this environment (HTTP 403), so real data
had to be encoded from established, citable knowledge rather than fetched live.

---

## 3. What was done (summary)

| Area | Change |
|---|---|
| **Data honesty** | Deleted all fabricated data; rebuilt around cited curated knowledge |
| **Knowledge base** | New `geroquery/knowledge/` module — the scientific source of truth |
| **Clock** | Implemented the real published **Levine PhenoAge** model |
| **Backend** | Refactored `models`/`store`/`sources`/`service`/`api` to serve cited evidence |
| **UI** | Ground-up Streamlit redesign with custom theme, plain-English analysis, figures + captions |
| **Tests** | Rewrote affected tests; **79 pass**; ruff/black/mypy clean |
| **Docs** | Rewrote README, DATA_SOURCES note, end-to-end notebook |

---

## 4. Detailed changes

### 4.1 New: `geroquery/knowledge/` (the scientific core)

The single source of truth for what GeroQuery asserts about aging. **All real,
curated, cited — nothing simulated.**

- **`references.py`** — 28 real publications as `Reference` dataclasses (authors,
  year, journal, title, PMID). URLs point to the PubMed article when a confident
  PMID is known, else a PubMed title search (so we never assert a possibly-wrong
  identifier — an unverifiable PMID would itself be false evidence).
- **`aging_knowledge.py`** — for each of the 12 ortholog groups (CDKN2A, CDKN1A,
  TP53, LMNB1, GDF15, SIRT1, TERT, KL, FOXO3, IGF1, MTOR, IGF1R):
  - `direction_with_age` (up / down / context-dependent)
  - `confidence` (ordinal: robust / established / emerging — **not** a fake number)
  - `hallmarks` of aging (López-Otín framework, defined in `HALLMARKS`)
  - a plain-English `one_liner`, `role`, and `analysis` paragraph
  - `evidence` items (omic layer, direction, species, tissues, strength, finding,
    citation keys)
  - `faqs` — plain-English Q&A
  - Also `INTERVENTIONS` (6 real lifespan interventions with citations).
- **`__init__.py`** — public exports and accessors (`gene_knowledge`,
  `interventions_for_group`, `KNOWLEDGE`, `HALLMARKS`, `INTERVENTIONS`, `REFERENCES`).

> **To add a gene:** add a `GeneKnowledge` entry in `aging_knowledge.py` keyed by
> its ortholog group, add the gene to `idmap/data/genes.json` (with matching
> `ortholog_group`), add any new refs to `references.py`, and add curated flags in
> `build_fixtures.py::CURATED`. Then `python -m geroquery.etl.build_fixtures`.

### 4.2 New: `geroquery/clocks/phenoage.py` (real clock)

Faithful implementation of the published **Levine PhenoAge** model (Levine 2018;
methodology Liu 2018, NHANES-derived):

- Nine clinical markers + age (`REQUIRED_FEATURES`), taken in **conventional US
  units** and converted internally to the SI units the published coefficients
  require (albumin g/dL→g/L, creatinine mg/dL→µmol/L, glucose mg/dL→mmol/L,
  CRP mg/L→mg/dL then ln).
- `phenotypic_age()` returns biological age (years); `mortality_risk_10yr()`
  returns the Gompertz 10-year mortality risk from the same linear predictor.
- **The coefficients are the real published ones**, not invented.
- Validated: healthy 40yo → ≈28; unhealthy 70yo → ≈84 (mortality risk 0.6% vs 64%).

`clocks/registry.py` now registers a single real `phenoage` clock (wrapping the
above) and keeps the `pyaging`/`biolearn` seam for real DNAm clocks when installed.
`clocks/service.py` gained `mortality_risk_10yr` on `ClockResult` and reads `age`
from the matrix when a separate chronological-age array isn't passed.

### 4.3 Data files (`geroquery/sources/data/`)

- **Deleted:** `signatures.csv`, `studies.csv`, `clinical_nhanes_slice.csv`.
- **`curated_knowledge.csv`** — real GenAge/CellAge/LongevityMap/OpenGenes flags
  (50 rows), regenerated with source links.
- **`interventions.csv`** — 6 real interventions (rapamycin, caloric restriction,
  metformin, D+Q senolytics, 17α-estradiol, NMN) with approximate rodent effects
  and **primary-citation PubMed URLs**, generated from the knowledge base.
- **`example_cohort_simulated.csv`** — 720 synthetic subjects with the nine
  PhenoAge markers + age + sex. **Clearly labelled SIMULATED.** A shared latent
  "aging" factor whose variance grows with age drives all markers, so both
  biological-age acceleration and the resilience (critical-slowing-down) signals
  emerge realistically. `subject_id`s are `SIM####`.

`build_fixtures.py` regenerates all three deterministically (seeded) and deletes
any stale fabricated files.

### 4.4 Backend refactor

- **`models.py`** — added `ortholog_group` to `CanonicalGene`; removed the
  fabricated-signature types (`AgingSignature`, `MetaSignature`, `Study`, `GeneCard`).
- **`idmap/resolver.py`** — populates `ortholog_group` on resolved genes.
- **`sources/local_fixture.py`** — replaced `LocalSignatureSource` with
  `LocalEvidenceSource` (curated cached tier); kept `CuratedKnowledgeSource` and
  `InterventionSource`. `sources/__init__.py` updated accordingly.
- **`store/store.py`** — dropped the fabricated signatures Parquet + studies table
  and the `query_signatures`/`list_studies` API; keeps curated flags,
  interventions, the example-cohort Parquet dataset, and reproducible versioning
  (content hash of the bundled source files).
- **`api/service.py`** — new `gene_report()` composes resolver + knowledge +
  curated flags + interventions + references + hallmarks into the rich object the
  UI/API use. Added `list_curated_genes()`, `geneset_summary()`, `references()`.
  `gene_card()` is a back-compat alias. Removed the fake meta-analysis path.
- **`api/app.py`** — endpoints: `/v1/gene/{id}/report` (+ `/card` alias),
  `/v1/genes`, `/v1/geneset/summary`, `/v1/references` (replacing `/studies`);
  clocks/intervention/resilience/sources/datasets unchanged.

> The `harmonize` module (Hedges' g, random-effects meta-analysis, batch
> correction) is **kept as-is** — it's correct library code and its tests use
> in-test synthetic data. It's simply no longer fed fabricated numbers. When real
> per-study data is supplied, it's ready to use.

### 4.5 UI: `geroquery/ui/streamlit_app.py` (full rewrite) + `.streamlit/config.toml`

Custom CSS theme (Inter typography, cards, pills, hero, dark sidebar). Four tabs:

- **Gene explorer** — search + example chips; gene header; verdict hero
  (↑/↓/↔ glyph + direction/confidence pills + one-liner); "What the data tells us"
  analysis callout; **aging-signal-by-data-type figure** (bars point up/down,
  length = curated confidence, explicitly *not* a measured effect size);
  cross-species panel; hallmark chips; evidence rows with citations; curated DB
  badges; **intervention bars**; FAQ expanders; references; CSV download.
- **Aging clock** — PhenoAge explainer; simulated cohort (SIMULATED badge) or CSV
  upload (+ template); metric tiles; **PhenoAge-vs-chronological scatter colored by
  acceleration** with identity line; acceleration histogram; **mortality-vs-age
  scatter**; results download. Every figure has a caption.
- **Resilience** — plain-English CSD explainer; SIMULATED badge / upload; verdict
  card; **variance-vs-age** and **cross-correlation-vs-age** figures with captions;
  method/assumptions/limitations expander.
- **About & data** — what's real vs. what's simulated; sources & licences table;
  full reference list.

### 4.6 Tests / docs

- `tests/conftest.py` — `clinical_matrix` fixture now yields a realistic 9-marker
  PhenoAge cohort.
- Rewrote `test_clocks.py`, `test_store.py`, `test_sources.py`, `test_api.py`,
  `test_crosscutting.py` to assert the honest model. **`test_idmap.py`,
  `test_harmonize.py`, `test_resilience.py` unchanged.**
- `pyproject.toml` — added `[tool.ruff.lint.per-file-ignores]` allowing long
  content strings (E501) in the UI and knowledge modules.
- Rewrote `README.md`, added an honesty note to `docs/DATA_SOURCES.md`, rewrote
  `notebooks/01_end_to_end.ipynb` (executes cleanly).

---

## 5. How to run

```bash
pip install -e ".[dev,ui]"
python -m geroquery.etl.build_fixtures                 # regenerate bundled data
python -m streamlit run geroquery/ui/streamlit_app.py  # dashboard
```

Other: `make api` (FastAPI at :8000/docs), `make test`, `make ci`, `docker compose up`.

## 6. How to verify

```bash
ruff check geroquery tests      # clean
black --check geroquery tests   # clean
mypy geroquery                  # clean
python -m pytest                # 79 passed
```

The Streamlit app was verified rendering live (headless Chromium via Playwright:
`/opt/pw-browsers/chromium-1194/chrome-linux/chrome`) for all three functional
tabs, including clicking "Apply PhenoAge clock" and "Compute resilience."

---

## 7. Key design decision (important context)

**Why one simulated element remains.** The clock and resilience tools
fundamentally need a **per-subject biomarker matrix**. Outbound biomedical APIs
are blocked in this environment, so real de-identified cohort rows (e.g. NHANES)
could not be fetched. The chosen resolution:

- The **algorithm is real** (published PhenoAge coefficients; published CSD theory).
- The example cohort is a **transparent, clearly-labelled `SIMULATED` sandbox** for
  trying the tools — never presented as real people.
- **Upload-your-own CSV is the first-class path** (template provided in the app).

This is the honest distinction from the old behavior: fabricating GEO accessions
and p-values to fake a meta-analysis is false *evidence*; a labelled simulation
demonstrating a *real* algorithm is a standard worked example. This decision is
disclosed in the UI (About tab), README, and DATA_SOURCES.

---

## 8. Known limitations / future work

- **Curated gene set is 12 genes.** Searching outside it returns a clear "not
  curated yet" message rather than inventing an answer. To expand, follow the
  "To add a gene" note in §4.1.
- **Real cohort data:** if/when network access or a bundled open cohort is
  available, replace `example_cohort_simulated` with real de-identified NHANES
  rows and drop the SIMULATED badge. The clock/resilience code needs no changes —
  only the dataset.
- **Live resolution & DNAm clocks:** set `GEROQUERY_ALLOW_NETWORK=1` for
  mygene.info identifier resolution; install `pyaging`/`biolearn` to auto-register
  real epigenetic clocks (the registry seam already handles this).
- **Quantitative gene evidence:** currently ordinal confidence + direction. If a
  real, citable meta-analysis with reproducible effect sizes is encoded, the
  `harmonize` module can surface pooled numbers again — honestly this time.
- **PRD doc** (`docs/PRD_GeroQuery_Multi-Omic_Aging_Aggregator (1).md`) is a
  historical planning artifact describing the original (partly fabricated) design;
  left as-is. Update if it should reflect the new honest model.

---

## 9. Pointers

- Scientific truth: `geroquery/knowledge/aging_knowledge.py`, `references.py`
- Real clock: `geroquery/clocks/phenoage.py`
- Orchestration: `geroquery/api/service.py` (`gene_report`)
- UI: `geroquery/ui/streamlit_app.py`
- Data generation: `geroquery/etl/build_fixtures.py`
- PR: https://github.com/Sophie-S-Z/GeroQuery/pull/1
