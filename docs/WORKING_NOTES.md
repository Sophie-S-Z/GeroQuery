# GeroQuery — Tier 0 working notes

Live scratch state for the in-progress Tier 0 build. Delete once Tier 0 ships and
the content has moved into `README.md` and `docs/RESULTS_NHANES_CSD.md`.

---

## Decisions taken (owner, 2026-08-03)

| Question | Decision |
|---|---|
| Scope of this push | **Tier 0, done properly.** Real data end-to-end, real network paths, real clocks, all HIGH/MEDIUM defects fixed, README rewritten to the true real/synthetic boundary. Stop at the Tier 0 exit criteria. |
| If real NHANES shows weak/no CSD | **Publish the null; honesty is the asset.** Report the real result with CIs, keep the synthetic case beside it explicitly labelled as *method validation* (estimator recovers a planted effect), not as evidence about biology. |
| How real data lives in the repo | **Download-on-demand + pinned manifest.** Exact source URLs, release year, SHA-256. `make data` fetches and verifies. A small verified sample is committed so tests and CI run offline. |
| React frontend (`frontend/`) | **Keep**, and repoint at real data at the end of the work. Not cut. |

## Environment facts verified this session

- Network reachable for all Tier 0 upstreams: mygene.info, GTEx Portal API v2,
  CDC NHANES XPT, ClinicalTrials.gov v2.
- NHANES XPT parses with `pandas.read_sas(..., format="xport")` — **no `pyreadstat`
  dependency needed**.
- **Working NHANES URL pattern** (the `/Nchs/Nhanes/2017-2018/{FILE}.XPT` form is dead —
  returns HTML, not XPORT):
  `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/{FILE}.xpt`
- Not installed locally: `pyaging`, `biolearn`, `pyreadstat`, `statsmodels`.
- Local Python is 3.14.3; `pyproject.toml` targets `>=3.10` and CI matrixes 3.10–3.12.

## NHANES 2017–2018 variable map (verified against downloaded files)

| GeroQuery column | NHANES var | File | Units |
|---|---|---|---|
| `age` | `RIDAGEYR` | `DEMO_J` | years, **topcoded at 80** |
| `sex` | `RIAGENDR` | `DEMO_J` | 1 = male, 2 = female |
| `albumin` | `LBXSAL` | `BIOPRO_J` | g/dL |
| `creatinine` | `LBXSCR` | `BIOPRO_J` | mg/dL |
| `glucose` | `LBXSGL` | `BIOPRO_J` | mg/dL |
| `crp` | `LBXHSCRP` | `HSCRP_J` | mg/L (high-sensitivity) |
| `lymphocyte_pct` | `LBXLYPCT` | `CBC_J` | % |
| `rdw` | `LBXRDW` | `CBC_J` | % |

Join key is `SEQN`. These are exactly the six markers `clocks/registry.py:28`
(`CLINICAL_FEATURES`) already expects, so the clock path needs no schema change.

Complete cases, age >= 20: **n = 4,895** (354 subjects sit at the age-80 topcode).

---

## Headline finding — real NHANES CSD

**Promoted to [`docs/RESULTS_NHANES_CSD.md`](RESULTS_NHANES_CSD.md)**, regenerated
through the repo code path (`sources.nhanes` -> `ResilienceService.csd`) rather
than a standalone script. Summary: variance supported 20/20 configurations,
cross-correlation supported 0/20, `resilience_declines` -> False.

One correction to an earlier note in this file: the Kendall p-floor claim was
stated as "K=4 returns tau=+1.00, p=0.333". The measured K=4 cell is tau=+0.333,
p=0.75; the *floor* argument is what holds — with four strata the best achievable
two-sided Kendall p is 0.083, so K=4 can never reach p<0.05 even for a perfect
monotone increase. Verified empirically, recorded in the results doc.

---

## Tier 0 status: complete

| Item | State |
|---|---|
| Resilience defects (recovery / CSD / control) | done |
| Pinned manifest + checksum-verifying fetch layer | done — `sources/manifest.py`, `sources/fetch.py` |
| Real NHANES adapter | done — `sources/nhanes.py`, n=4,895; 600-row real sample committed |
| Synthetic generator renamed + separated | done — `build_clinical_synthetic`, dataset `clinical_synthetic_csd` |
| Real mygene.info batch resolution + cache | done — `idmap/mygene.py` |
| Real GTEx Portal v2 adapter | done — `sources/gtex.py` |
| Real clocks | done — `clocks/library.py` wraps biolearn |
| Store/API defects | done — connection reuse, glob-level partition pruning, cached `ensure_built`, narrowed exception, `InterventionNotFoundError` |
| CI | done — `workflow_dispatch` + weekly schedule; live job installs `[dev,clocks]` |
| Tests / README / React | done — 159 tests + 6 live, 82% coverage; README rewritten to the real/synthetic boundary; React shows the real NHANES resilience panel |

### Known gaps, deliberately left

1. **`pyaging` is not wrapped.** Needs torch + downloaded model artifacts; a
   wrapper CI cannot exercise would be untested code claiming real-clock support.
2. **GTEx cannot be age-stratified.** The open v2 API accepts
   `attributeSubset=ageBracket` but returns one undivided group — donor age is
   dbGaP-controlled. A live test pins this so we notice if it changes.
3. **The multi-omic signature slice is still synthetic.** This is the remaining
   Tier 1 item; the clinical/resilience path is fully real.
4. **biolearn could not be installed locally** (its `ecos` dependency has no
   wheel for Python 3.14). The wrapper was written against biolearn 0.9.1's
   actual source, and is exercised offline against a fake shaped like that API;
   the real library runs in the `smoke-live` CI job on Python 3.11.

---

## Tier 1.1 — real gene signatures (same day)

Complete. `signatures.csv` is gone; `signatures_full.csv` (git-ignored) and
`signatures_curated.csv` (committed) come from the GEO DataSets panel. Detail in
[`RESULTS_GEO_SIGNATURES.md`](RESULTS_GEO_SIGNATURES.md) and
[`OVERVIEW.md`](OVERVIEW.md).

### Facts verified this session

- **GDS full SOFT path**: `https://ftp.ncbi.nlm.nih.gov/geo/datasets/{SHARD}/{ACC}/soft/{ACC}_full.soft.gz`
  where `SHARD` is the accession with its last three digits replaced by `nnn`
  (`GDS707` -> `GDSnnn`, `GDS5226` -> `GDS5nnn`).
- **`acc.cgi?targ=self&form=text` does not work for GDS** — it redirects to the
  GDS browser. To triage headers cheaply, HTTP `Range`-request the first 128 KB of
  `{ACC}.soft.gz` and `zlib.decompressobj(MAX_WBITS|16)` it; the whole metadata
  block precedes `!dataset_table_begin`.
- **GEO DataSets query for the panel**:
  `age[Subset Variable Type] AND gds[Filter] AND (Homo sapiens[Organism] OR Mus musculus[Organism])`
  → 189 records on 2026-08-03. A live test pins that this stays ≥180.
- All six HAGR ZIPs download unauthenticated and are 8–49 KB.
- **DrugAge and GenDR are both served from a URL ending `/dataset.zip`.** This
  collided in the flat download cache; `fetch.cache_path` now prefixes the
  manifest key.

### Three real bugs found by running real data through the pipeline

1. **`idmap/mygene.py::_ensembl_gene` returned other species' identifiers.**
   mygene's `ensembl` field is sometimes a list containing homology
   cross-references; taking `[0]` gave mouse genes ids like `ENSFALG...`
   (collared flycatcher) and `ENSMPUG...` (ferret), which then became the
   canonical join key. Now filtered by the taxid's expected stable-id prefix,
   falling back to `ENTREZ:`. Affected the request path too, not just the ETL.
2. **The download cache collided on URL basename** (above).
3. **`/v1/intervention/{name}` returned `matches[0]`.** With real DrugAge there is
   one record per organism, so asking about rapamycin returned the *C. elegans*
   result (+19.8%) rather than the mouse one (+13.0%). Now returns all organisms,
   mammals first.

### Numbers to sanity-check a rebuild against

| | |
|---|---|
| Panel | 31 datasets, 32 contrasts, 27 GEO Series |
| Signature rows | 485,905 over 46,091 genes |
| Committed slice | 40,585 rows, 3,168 genes, both species |
| Curated assertions | 2,720 mammalian (5,029 parsed, 2,276 non-mammalian dropped) |
| Interventions | 1,340 |
| Rows dropped at id resolution | 141,156 of 627,155 (~22%) |
| Rapamycin, *Mus musculus* | +13.0% median over DrugAge's significant ITP rows |
