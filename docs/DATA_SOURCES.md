# Data-source contract

Cache-vs-federate is decided by **licence and size**, and enforced in code: each
adapter declares `license()` and `capabilities()`, and the store calls
`assert_cacheable()` before persisting anything. Controlled sources are
federate-only and non-redistributable — a test (`test_crosscutting.py`) verifies
they refuse caching.

**Status** distinguishes what is wired from what is contracted. A table that
lists a planned source identically to a live one is how a reader concludes the
system already ingests it.

| Source | Status | Adapter role | Access | Format | Cache vs Federate | Licence |
|---|---|---|---|---|---|---|
| **NCBI GEO DataSets** | **live** | gene aging signatures — 31 datasets, 32 young-vs-old contrasts, 485,905 effect sizes | FTP SOFT, SHA-256 pinned | SOFT (gz) | Cache harmonized extracts | Public, attribute |
| NHANES | **live** | clinical markers (clocks/resilience), 2017-2018, n=4,895 | CDC download, SHA-256 pinned | XPT | Cache | US public domain |
| mygene.info | **live** | gene ID harmonization | REST (batch) | JSON | Cache mappings | Public |
| GTEx (open) | **live** | tissue expression (median TPM + UBERON), **no age brackets** | Portal API v2 | JSON | Cache summaries; federate detail | Open summary |
| **HAGR — GenAge (human)** | **live** | 307 human ageing-associated genes | ZIP, SHA-256 pinned | CSV | Cache | HAGR, non-commercial + attribution |
| **HAGR — GenAge (models)** | **live** | 2,205 lifespan-manipulation genes (132 mammalian loaded) | ZIP, pinned | CSV | Cache | HAGR |
| **HAGR — CellAge** | **live** | 949 senescence regulators | ZIP, pinned | TSV | Cache | HAGR |
| **HAGR — LongevityMap** | **live** | 1,325 longevity-variant assertions, nulls included | ZIP, pinned | CSV | Cache | HAGR |
| **HAGR — DrugAge** | **live** | 3,423 lifespan experiments -> 1,334 compound x organism records | ZIP, pinned | CSV | Cache | HAGR |
| **HAGR — GenDR** | **live** | 214 dietary-restriction-essential genes | ZIP, pinned | CSV | Cache | HAGR |
| biolearn | **live** | 63 published clocks | Python package | — | n/a | Library terms |
| pyaging | **live** | 173 published clocks (metadata-registered) | Hugging Face | — | Lazy artifacts | Library terms |
| UBERON / Cell Ontology | **bundled** | tissue/cell harmonization | OBO/OLS extract | JSON | Cache | Open |
| AnAge | **bundled** | maximum lifespan -> fractional age | flat file | JSON | Cache | HAGR |
| recount3 | **stub** | uniform RNA-seq (human+mouse) — the highest-value next source | R/Bioconductor (ETL) | RSE | Federate / cache extracts | Public |
| ARCHS4 | **stub** | uniform RNA-seq counts, no R bridge needed | HDF5 | HDF5 | Cache slices | Public |
| GTEx (protected) | **stub** | raw genotype-linked; donor age lives here | dbGaP | — | **Federate/link only** | Controlled |
| CELLxGENE Census | **stub** | single-cell (stretch) | `cellxgene_census` | SOMA/AnnData | Federate | CC-BY |
| Human Protein Atlas | **planned** | protein aging views | programmatic/TSV | TSV/JSON | Cache | CC-BY-SA |
| Lehallier plasma proteome | **planned** | proteome aging slice | supp. tables | XLSX | Cache | Per-paper terms |
| ComputAgeBench | **planned** | methylation benchmark | Hugging Face | Parquet | Cache | Public |
| Open Genes | **planned** | curated gene knowledge | REST / SQLite | JSON/SQLite | Cache mirror | MPL-2.0 |
| **UK Biobank, dbGaP, ADNI, MIMIC, HRS-restricted** | **stub** | controlled cohorts | controlled-access | — | **Link/federate only, never re-host** | Controlled |

**41 pinned artifacts, ~320 MB, from 5 independent upstreams (NCBI, CDC, HAGR,
mygene.info, GTEx).** Full detail in [`OVERVIEW.md`](OVERVIEW.md).

## Build-time notes

- Database entry counts and sample sizes grow over time — re-verify before depending on any live source.
- Confirm current NIA ITP data embargo / acknowledgment terms before redistribution. DrugAge
  flags which of its experiments came from the NIA ITP; those rows carry `source = ITP`.
- HAGR releases are pinned by checksum, so a HAGR update fails the build loudly rather than
  changing results silently. That is intended; it means the pins need periodic maintenance.
- The GEO panel is not hand-picked. It is every GEO DataSet declaring an `age` subset variable
  (189 records) that survives the contrast rules in `sources/geo.py`. Change a rule and the
  panel changes; nothing was selected for what it showed.
- The R/Bioconductor bridge (recount3, GEOquery, gtexr) runs **only** inside offline, batch, versioned ETL — never on the API request path.

## Status vocabulary

- **live** — a real adapter that queries the upstream and returns real data.
- **bundled** — real curated content shipped in the repo; no live fetch yet.
- **stub** — `FederatedStub`: declares capabilities and licence so the store, API,
  and licence test can reason about it, but has no fetch body. Querying it errors.
- **planned** — contracted in this table only. No adapter exists.

## Correction: GTEx and age

An earlier version of this table listed GTEx (open) as providing "age brackets".
It does not. `/expression/geneExpression` accepts `attributeSubset=ageBracket` and
returns HTTP 200, but with a single record whose `subsetGroup` is `null` and the
full undivided sample list — donor age lives in the dbGaP-controlled tier. GTEx
open therefore contributes **tissue context**, not an aging signature. A live test
(`test_live_gtex_age_bracket_subsetting_still_does_not_work`) pins this so the
claim gets rechecked automatically rather than by memory.

## Checksum pinning

Every downloaded artifact is declared in `geroquery/sources/manifest.py` with its
exact URL, byte count, and SHA-256, and `geroquery/sources/fetch.py` refuses any
bytes that do not match. Upstreams re-publish files under the same path; without
a checksum that becomes an unexplained change in results months later with no way
to attribute it. Cache hits are re-verified rather than trusted, and a failed
download is removed rather than left where a later offline run would pick it up.

## Cache-path collision

`fetch.cache_path` prefixes each cache entry with its manifest key. URL basenames
are not unique across upstreams: HAGR serves both DrugAge and GenDR from a path
ending `/dataset.zip`, so keying on the basename alone made the two share one
cache entry and whichever downloaded second overwrote the first. The checksum
check turned that into a loud failure rather than a silent wrong answer — which
is what it is for — but the collision is fixed rather than relied upon.
