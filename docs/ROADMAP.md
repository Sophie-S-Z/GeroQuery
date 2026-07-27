# GeroQuery Roadmap

This roadmap distills the strategic assessment in [`STRATEGY_2026.md`](STRATEGY_2026.md) into
what is **done**, what is **deferred** (and why), and the **recommended future cuts**. It is the
honest-positioning document the strategy calls for: it states plainly which methodological gaps
GeroQuery fills today and which are aspirational, so nothing is over-claimed.

See [`../HANDOFF.md`](../HANDOFF.md) for the implementation history and
[`DATA_SOURCES.md`](DATA_SOURCES.md) for the data/federation contract.

## Where GeroQuery competes

GeroQuery cannot win on breadth — HAGR, Open Genes, Aging Atlas, recount3/ARCHS4, pyaging and
biolearn already cover curation, scale, and clocks. Its niche is a **methodological rigor layer**
those incumbents lack: cited, uncertainty-quantified, honestly-scoped answers about how a gene
relates to aging, plus biomarker tools that surface their own limitations. Naming stays
**GeroQuery** (the assessment endorses it).

## Done (this repo, offline-feasible)

| Item | Strategy ref | Status |
|---|---|---|
| Remove all synthetic magnitudes / fabricated evidence | Rec 1, Cut list | ✅ Done (prior revamp — see HANDOFF) |
| Real, cited curated gene knowledge (directions, hallmarks, plain-English, PubMed refs) | §Key Findings | ✅ Done |
| Real published clock (Levine PhenoAge) instead of toy clocks | §2 | ✅ Done |
| **Clock applicability / diagnostics** (completeness, out-of-range flags, unit-mismatch heuristics, bootstrap CI) | §2 top white space, Rec 1b | ✅ Done — `clocks/diagnostics.py`, `POST /v1/clock/diagnostics`, UI panel |
| **Uncertainty quantification** (bootstrap CIs on PhenoAge acceleration and on resilience trends) | Rec 3 | ✅ Done |
| **Honest resilience reframing** — "age-stratified dispersion & co-fluctuation indicator", contested-EWS caveats cited | §3, Rec 6 | ✅ Done |
| **DNB / l-DNB composite index** (pure Python, cross-sectional strata, with caveats) | §3, Rec 6 | ✅ Done — `resilience/dnb.py` |
| **Hartung–Knapp** option for random-effects meta-analysis (DL kept as default) | Rec 5 | ✅ Done — `harmonize/meta.py` |
| Hallmarks-of-aging mapping | §4 #9 | ✅ Done (in the knowledge base) |
| AR(1)/Ornstein–Uhlenbeck recovery-rate metric for longitudinal data | §3 | ✅ Present — `resilience/recovery.py` |

## Deferred — needs blocked data or heavy infrastructure (not built here)

These require outbound access to biomedical data services (**network-blocked in the current
environment**) and/or toolchains this repo intentionally does not carry (R, Next.js, Nextflow).
They are the highest-value next steps once real data ingestion is available.

| Item | Strategy ref | Why deferred |
|---|---|---|
| **Trajectory-aware gene profiles** (GAM + changepoint + wave classification per Shen 2024) | §2, Rec 1a | Needs real per-age omics (recount3/GTEx). The method is small (`pyGAM`+`ruptures`); ship it against real trajectories, not simulated data. Strongest single differentiator. |
| Real-data rebuild on Open Genes / HAGR / recount3 / GTEx-open / Tabula Muris Senis / ComputAgeBench | Rec 1 | Requires downloading those datasets (network-blocked here). |
| Clock diagnostics **benchmarked on ComputAgeBench** | §2 | Needs the 66-dataset HF Parquet download. |
| DNB **longitudinal validation** on Dog Aging Project / Mouse Phenome DB | §3, Rec 2 | Needs open longitudinal cohort download. |
| Cell-composition adjustment (EpiDISH / EPISCORE) | §2, Rec 2 | Needs reference panels + real bulk matrices. |
| Signature-reversal endpoint (L2S2 / LINCS) | §4 #7 | Needs LINCS access + caveat handling. |
| Cross-species Dutta–Sengupta dual normalization + orthology confidence (Ensembl Compara / DIOPT) | §2 | Needs Compara/DIOPT lookups. |
| Provenance/FAIR emission (PROV-O / RO-Crate, Zenodo DOIs) | §2 | Meaningful only once real datasets are ingested. |

## Recommended future cuts (documented; not executed by default)

The assessment recommends removing these; they are **left in place** for now because they are
non-harmful and a maintainer may still want them. Track them here and remove on decision:

- **Standalone React showcase (`frontend/`)** — the strategy calls it weakly differentiating and
  recommends a single Streamlit analyst surface (or, later, a proper API-first Next.js frontend for
  flagship views). It is isolated and does not affect the Python package. Keep only if you want a
  minimal API-consumer demo; otherwise remove.
- **Network control energy (`resilience/control.py` + `ResilienceService.control_energy`)** — a
  correct, unit-tested LTI minimum-control-energy implementation that is **exposed nowhere** (no
  API/UI/README). Marked "stretch" in the historical PRD. Safe to remove if you want to shrink the
  resilience surface to the two defensible metrics (dispersion/co-fluctuation + AR(1) recovery).
- **Universal-atlas ambition** — keep scope narrow and deep rather than chasing breadth.

## Positioning notes

- **Publication:** a NAR Database Issue is impractical for a solo, unhosted tool; target a
  Bioinformatics Applications Note or a bioRxiv preprint, and consider the Biomarkers of Aging
  Challenge as a community entry point.
- **Ecosystem:** rather than build a bespoke MCP server, contribute to the existing
  **longevity-genie** org, which already occupies the LLM-queryable-aging-data niche.
- **Hiring signal:** the durable value is a rigorously engineered, reproducible, real-data tool with
  a clear README + reproducible notebook — not breadth of features.
