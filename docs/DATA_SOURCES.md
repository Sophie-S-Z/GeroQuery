# Data-source contract

Cache-vs-federate is decided by **licence and size**, and enforced in code: each
adapter declares `license()` and `capabilities()`, and the store calls
`assert_cacheable()` before persisting anything. Controlled sources are
federate-only and non-redistributable — a test (`test_crosscutting.py`) verifies
they refuse caching.

| Source | Adapter role | Access | Format | Cache vs Federate | Licence |
|---|---|---|---|---|---|
| NCBI GEO | transcriptome + methylation | GEOparse / E-utilities | SOFT / series matrix | Cache harmonized extracts | Public, attribute |
| recount3 | uniform RNA-seq (human+mouse) | R/Bioconductor (ETL) | RSE | Federate / cache extracts | Public |
| GTEx (open) | tissue expression, age brackets | Portal API v2 | JSON/GCT | Cache summaries; federate detail | Open summary |
| GTEx (protected) | raw genotype-linked | dbGaP | — | **Federate/link only** | Controlled |
| ARCHS4 | uniform RNA-seq counts | HDF5 | HDF5 | Cache slices | Public |
| CELLxGENE Census | single-cell (stretch) | `cellxgene_census` | SOMA/AnnData | Federate | CC-BY |
| Human Protein Atlas | protein aging views | programmatic/TSV | TSV/JSON | Cache | CC-BY-SA |
| Lehallier plasma proteome | proteome aging slice | supp. tables | XLSX | Cache | Per-paper terms |
| ComputAgeBench | methylation benchmark | Hugging Face | Parquet | Cache | Public |
| NHANES | clinical markers (clocks/resilience) | CDC download | XPT/CSV | Cache | Public |
| Open Genes | curated gene knowledge | REST / SQLite | JSON/SQLite | Cache mirror | MPL-2.0 |
| HAGR (GenAge/CellAge/LongevityMap/DrugAge/GenDR) | curated gene/drug/species | flat-file | CSV/TSV | Cache | Free, attribute |
| NIA ITP | mouse lifespan interventions | Mouse Phenome DB | CSV/portal | Cache | Public; embargo/ack rules |
| mygene.info | gene ID harmonization | REST (batch) | JSON | Cache mappings | Public |
| UBERON / Cell Ontology | tissue/cell harmonization | OBO/OLS | OBO/OWL | Cache | Open |
| **UK Biobank, dbGaP, ADNI, MIMIC, HRS-restricted** | controlled cohorts | controlled-access | — | **Link/federate only, never re-host** | Controlled |

## Build-time notes

- Database entry counts and sample sizes grow over time — re-verify before depending on any live source.
- Confirm current NIA ITP data embargo / acknowledgment terms before redistribution.
- The R/Bioconductor bridge (recount3, GEOquery, gtexr) runs **only** inside offline, batch, versioned ETL — never on the API request path.
