<!--
Filed for reference: an external strategic assessment of GeroQuery (2026).
This is an advisory document, not a specification. The parts of it that are
feasible in the current offline environment have been actioned (see
docs/ROADMAP.md); data-/infrastructure-dependent items are tracked there as
deferred work. Filed verbatim below.
-->

# GeroQuery: A Technical Strategy to Elevate, Differentiate, and Ship a Credible Aging-Data Tool

## TL;DR
- **GeroQuery's aggregator layer sits in a crowded field** (HAGR, Open Genes, Aging Atlas, recount3, ARCHS4, pyaging, biolearn), and its "curated demonstration slice with synthetic magnitudes" would sink credibility with both scientists and hiring managers. The single most defensible differentiator is **not** the resilience module as currently specced (cross-sectional critical-slowing-down is a weak, largely unvalidated proxy) but a **trajectory-aware, uncertainty-quantified, cell-composition-adjusted meta-analysis service plus a "will this clock be valid on my data?" diagnostics checker** — none of the incumbents do these well.
- **Cut hard and go deep on one thing.** Drop the synthetic magnitudes, the network-control-energy view, the standalone React showcase, and the universal-atlas ambition. Build a small number of real, reproducible, well-sourced endpoints on genuinely open data (Open Genes, HAGR, recount3, GTEx-open, Tabula Muris Senis, ComputAgeBench) where every number carries a confidence interval and a provenance trail.
- **For the portfolio goal specifically:** longevity biotechs (Retro, NewLimit, BioAge, Insilico) hire for single-cell/multi-omics + ML + Nextflow + cloud + reproducibility. A solo toy tool is weakly differentiating; a rigorously engineered, real-data tool that demonstrates exactly those competencies — ideally contributed into the existing **longevity-genie** open ecosystem rather than competing with it — is what actually moves a hiring decision.

## Key Findings

### The aggregator white space is narrow; the methodological white space is wide
- **Curated gene/intervention databases are essentially filled.** HAGR (de Magalhães lab, University of Birmingham) bundles GenAge (>300 human + >2,000 model-organism genes), AnAge (>4,000 species), DrugAge (>500 compounds), CellAge (866 senescence genes), GenDR (214 dietary-restriction genes), and LongevityMap (3,144 variants / 884 genes). Open Genes (Rafikova et al., *Nucleic Acids Research* 2023; 2,402 genes, open-genes.org) is newer, ships an API and UI, and explicitly criticizes HAGR for making cross-experiment comparison "impossible." Aging Atlas (CNCB/NGDC, *NAR* 2021) spans transcriptomics/scRNA-seq/epigenomics/proteomics/pharmacogenomics but is manually curated. All are live and maintained.
- **The exact "LLM-queryable aging data" niche already exists.** The **longevity-genie** GitHub org publishes `gget-mcp`, `opengenes-mcp` (natural-language→SQL over Open Genes), `synergy-age-mcp`, `biothings-mcp`, `pharmacology-mcp`, `alphagenome-mcp`, and the `holy-bio-mcp` umbrella (50+ functions). A generic MCP interface over aging data is therefore **not** novel.
- **The unsolved scientific problems are the feature-generating core.** The biggest is that **aging is non-linear**: Shen et al. (*Nature Aging* 2024, 4:1619–1634, doi:10.1038/s43587-024-00692-2), profiling n=108 individuals aged 25–75 with median 1.7-year follow-up, found that "only 6.6% of the molecules and microbes changed with age in a linear fashion," while ~81.0% showed nonlinear patterns "with substantial dysregulation occurring at two major periods occurring at approximately 44 years and 60 years of chronological age." A single Hedges' g "effect with age" per gene is therefore scientifically wrong for the large majority of genes — and **no incumbent aging portal models trajectory shape.** This is the strongest white space.

### The clock ecosystem's known failures are directly monetizable as features
- pyaging (de Lima Camillo, *Bioinformatics* 2024; GPU, 30+→100+ clocks) and biolearn (Biomarkers of Aging Consortium; 39 biomarkers across 20,000+ individuals) are the incumbents. biolearn's own analysis shows chronological-age prediction does not correlate with mortality prediction (R=0.12, P=0.67) — i.e., the two clock classes measure different things.
- Higgins-Chen et al. (*Nature Aging* 2022, 2:644–661, doi:10.1038/s43587-022-00248-2) showed "technical noise produces deviations up to 9 years between replicates for six prominent epigenetic clocks," while retrained principal-component versions "show agreement between most replicates within 1.5 years" (clock ICCs 0.917–0.979; some age-acceleration ICCs as low as 0.755). Sugden et al. showed many individual CpGs are unreliable. Array attrition (450k→EPIC→EPICv2), missing-probe imputation, and cross-tissue/cross-population misapplication are all unhandled, and clock outputs ship without confidence intervals. ComputAgeBench (Kriukov et al., bioRxiv 2024 / KDD '25) — which "collected and harmonized 66 public datasets of blood DNA methylation, covering 19 such conditions... and tested 13 published clock models" — found second-generation clocks (GrimAge2, PhenoAge) generalize best while some clocks are implicitly disease-specific. A **"will this clock be valid on my data?" diagnostics service is strong, unoccupied white space.**

### The resilience/criticality module is the riskiest claim and needs reframing
- The DOSI work (Pyrkov, Fedichev et al., *Nat Commun* 2021, 12:2765) is real but used **longitudinal** complete-blood-count and wearable step data; its headline — "Extrapolation of this trend suggested that DOSI recovery time and variance would simultaneously diverge at a critical point of 120−150 years of age corresponding to a complete loss of resilience" — depends on that longitudinal structure. GeroQuery's cross-sectional, age-stratified proxy is a fundamentally weaker construct.
- Critical-slowing-down early-warning signals are contested: false positives are systematic (Jäger & Füllsack, *PLoS ONE* 2019), results are sensitive to detrending/window/bandwidth (Boettiger & Hastings 2012; Dakos et al. 2012/2015), and there is "no evidence for critical slowing down prior to human epileptic seizures" (arXiv:1908.08973). Cross-sectional variance-across-age-strata is **not** a validated CSD proxy.
- The strongest defensible upgrade is **Dynamic Network Biomarker (DNB) theory** (Chen Luonan et al., *Sci Rep* 2012, srep00342; single-sample l-DNB, *National Science Review* 2019), which detects pre-transition states from a composite of rising intra-module standard deviation, rising intra-module correlation, and falling module-to-background correlation, and ships a maintained `DNBr` R package. It is under-used in aging portals and is implementable on cross-sectional strata with explicit caveats.

## Details

### 1. Competitive / prior-art map

| Resource | What it does | What it does NOT do | Licence/API | Status |
|---|---|---|---|---|
| HAGR (GenAge/AnAge/DrugAge/CellAge/GenDR/LongevityMap) | Curated aging genes, species longevity, drugs, senescence, DR, variants | Cross-study quantitative meta-analysis; trajectory modeling; clocks | Free, terms-of-use; flat files | Maintained (Birmingham) |
| Open Genes | 2,402 curated genes with intervention/expression/association evidence; API+UI | Raw omics reprocessing; clocks; resilience | Open; REST API | Active (2023 NAR) |
| Aging Atlas (CNCB) | Multi-omics + scRNA-seq + pharmacogenomics browser | Programmatic meta-analysis; clocks; uncertainty | Free | Live; update cadence unclear |
| pyaging | 30+→100+ GPU clocks, multi-omics | Data aggregation; harmonization QC; applicability checks | Open (PyPI) | Active |
| biolearn | Harmonized clocks + unified datasets + evaluation | Trajectory; deconvolution; resilience | Open | Active (BoA Consortium) |
| ComputAgeBench | 66-dataset, 13-clock, 19-condition benchmark | Live service; user data; trajectory | Open (HuggingFace Parquet) | Active |
| recount3 / ARCHS4 / refine.bio | Uniformly processed RNA-seq at scale | Aging-specific harmonization; clocks; resilience | Open | Mature |
| L2S2 / CLUE / metaLINCS | LINCS L1000 signature reversal | Aging-specific signatures; longitudinal | Mixed | Active |
| longevity-genie (MCP org) | MCP servers over Open Genes, SynergyAge, gget, etc. | Harmonized meta-analysis; clocks-as-service; trajectory | Open | Active (2025) |

recount3 (Collado-Torres et al., *Genome Biology* 2021) comprises "over 750,000 publicly available human and mouse RNA sequencing (RNA-seq) samples uniformly processed by our new Monorail analysis pipeline" — 316,443 human and 416,803 mouse run accessions — so competing on scale is futile.

**Implication:** GeroQuery cannot win on breadth. It can win on a *methodological* layer the incumbents lack: trajectory-aware, composition-adjusted, uncertainty-quantified cross-study meta-analysis, plus clock applicability diagnostics.

### 2. Unsolved problems → concrete GeroQuery features

| Unsolved problem (source) | GeroQuery feature | Implementation |
|---|---|---|
| Aging is nonlinear; ~81% of markers non-linear, waves at 44/60 (Shen 2024) | Trajectory-aware gene profiles | Per-gene GAM/spline fits (`pyGAM`/`statsmodels`), changepoint detection (`ruptures`), wave classification; `GET /v1/gene/{id}/trajectory` returns shape class + changepoints + CI bands |
| Cell composition confounds bulk expression/DNAm (EpiDISH/EPISCORE; PMC10928575) | Composition-adjusted signatures as a first-class step | Deconvolve (EpiDISH/EPISCORE for DNAm; CIBERSORTx/Bisque for expression); report adjusted vs unadjusted effect; flag composition-driven genes |
| Clock replicate noise up to 9y; no CIs; misapplication (Higgins-Chen 2022; Sugden; ComputAgeBench) | **Clock diagnostics / applicability checker** | Report ICC-based reliability, missing-probe robustness, platform/tissue/ancestry-mismatch warnings, bootstrap CIs; `POST /v1/clock/diagnostics` |
| Cross-species normalization crude | Dutta–Sengupta + fractional-lifespan dual normalization | Offer both (Dutta & Sengupta 2016 piecewise: 150× to day 42, then 45×/30×/25×/20×); expose orthology confidence via Ensembl Compara / DIOPT with one-to-many flags |
| Intervention non-replication (ITP vs DrugAge) | Replication/confidence tier per intervention | Tier by number of independent labs, NIA-ITP status, effect size, Kaplan–Meier/log-rank availability |
| Signature-reversal pitfalls (L2S2; CMap reproducibility) | Reversal endpoint wrapping L2S2/LINCS with caveats | Submit up/down aging sets → candidate reversers; surface 978-landmark-gene imputation + cell-line + reproducibility caveats |
| Weak provenance/FAIR | PROV-O / RO-Crate emission, Bioschemas, dataset DOIs | Content-addressed `datasets.lock`, Zenodo DOIs, ontology-typed metadata (UBERON/EFO/MONDO/CL) |
| Statistical rigor gaps | Publication-bias + heterogeneity + concordance per result | Report I², Egger's test, trim-and-fill, cross-study sign concordance |

### 3. Making the resilience module defensible
- **Reframe honestly:** rename the cross-sectional metric an "age-stratified dispersion / co-fluctuation indicator," not "critical slowing down," and cite the EWS false-positive literature directly in both UI and docs.
- **Add DNB (l-DNB) as the flagship method:** composite index = (intra-module SD × intra-module correlation) / (module-to-background correlation), computable on cross-sectional strata; wrap `DNBr`. Cite Chen 2012 and Liu 2019.
- **Estimate an Ornstein–Uhlenbeck relaxation rate** wherever longitudinal data exists (recovery rate = OU mean-reversion coefficient); this is the DOSI-consistent, defensible metric.
- **Use real open longitudinal data** for genuine recovery-rate estimation: Dog Aging Project open data, Mouse Phenome Database / Nathan Shock longitudinal phenotyping, and open wearable datasets (PMData, LifeSnaps). Controlled cohorts (UK Biobank repeat assessments, HRS, ELSA, MESA, Framingham/dbGaP, BLSA, InCHIANTI, All of Us) require application and should be "federate/link-only."
- **What reviewers will demand for a preprint:** demonstration on ≥1 true longitudinal cohort; direct comparison of the cross-sectional proxy against a longitudinal ground truth; explicit false-positive/robustness analysis over window and detrending choices; and null models.

### 4. Feature ideas ranked by differentiation × effort

**High differentiation, moderate effort — build first:**
1. Trajectory-aware gene profiles (GAM + changepoint + wave classification) — directly operationalizes Shen 2024; nobody else does it.
2. Clock diagnostics/applicability service — operationalizes Higgins-Chen/Sugden/ComputAgeBench; strong white space.
3. Uncertainty quantification everywhere (bootstrap/conformal CIs on effects, clocks, resilience).
4. Reproducibility/confidence tier per gene–tissue result (I², sign concordance, Egger's, trim-and-fill).

**High differentiation, higher effort — build second:**
5. Cell-composition-adjusted signatures as a harmonization step.
6. DNB/l-DNB resilience module with longitudinal validation.
7. Signature-reversal endpoint wrapping L2S2/LINCS.

**Lower differentiation — borrow, defer, or skip:**
8. MCP server — valuable but longevity-genie already occupies it; contribute there rather than rebuild.
9. Hallmarks-of-aging mapping (López-Otín 2023) — cheap context layer.
10. Open Targets Genetics / GWAS Catalog longevity-evidence overlay — good context, moderate effort.
11. DuckDB-WASM in-browser demo — striking, but WASM is single-threaded by default and lacks feature parity (Parquet-compression and lazy-loading limitations reported by production users on Hacker News); restrict it to a small curated slice.

### 5. UI/UX
- **Framework:** Streamlit is fine for a fast analyst/demo surface but its full-rerun model and limited customization scale poorly and look template-y. For a tool that must impress engineers **and** serve biologists, build an **API-first FastAPI backend with a typed Next.js + React front end** for the flagship gene-card and trajectory views, keeping a single Streamlit page only as a quick analyst surface. Dash is the middle option if staying all-Python. Emulate Open Targets Platform and cBioPortal: fast gene search with disambiguation, deep-linkable permalinks, progressive disclosure of methodological caveats, and honest empty/loading/error states.
- **Charts:** forest plots (meta-analysis), volcano plots, trajectory plots with CI bands and marked changepoints/waves, survival/lifespan curves, clock-comparison scatter with an identity line, and — for resilience — a "recovery-rate vs age" plot with widening variance ribbons. Use viridis/cividis for magnitude and a colorblind-safe diverging map (blue–white–red) for signed effects, always labeling direction explicitly. Target WCAG AA.
- **Performance UX:** stream federated queries with server-sent events and progress indicators; cache aggressively per the licence gate; degrade gracefully when an upstream (GTEx-protected, GEO) is down.

### 6. Engineering
- DuckDB + Hive-partitioned Parquet fits gene×tissue×study×species well; prefer plain Hive partitioning over Iceberg/Delta at this scale. Serve DuckDB read-only with per-request connections and explicit memory limits; never share a mutable connection across web workers.
- Deployment: Hugging Face Spaces for the demo (aligns with the ComputAgeBench/HF ecosystem); Fly.io/Render for the API; DuckDB-WASM for a small curated client-side slice only.
- Testing/quality: property-based tests (Hypothesis), golden-file regression on effect sizes, contract tests for source adapters (pytest-recording/VCR), Pydantic v2 schemas, Pandera/Great Expectations data checks, and **Nextflow** for ETL (specifically, because Retro lists it as a requirement).

### 7. Positioning, naming, credibility
- **Naming collisions (checked):** "Resilio" is taken (Resilio Sync file-sync). "AgeAtlas" collides with the established Aging Atlas (*NAR* 2021) and generic "aging atlas" datasets. "Senescope" is taken (senescope.com senescence AI; also an IEEE tool). "OmniAge" is an existing company (omniage.com). "Senescan" is weakly distinctive amid SenePy/SenCID/SenoQuant. "GeroLake" is mostly clear (only geographic collisions). "GeroQuery" is mostly clear (only an OCR-typo overlap with the GeoQuery text-to-SQL dataset) but the "Gero-" prefix reads as affiliated with the Gero longevity company (gero.ai). **Recommendation: keep GeroQuery or pick GeroLake; avoid the rest.**
- **Publication venue:** the NAR Database Issue requires the resource to be freely available without login, live, and "updated or at least maintained in a fully functional form, ideally at the same URL for at least 5 years," with an explicit maintenance commitment from a senior author and host institution — impractical for a solo builder, since editors verify the site during review and reject stale resources. The realistic venue is a **Bioinformatics Applications Note or a bioRxiv preprint**; the Biomarkers of Aging Challenge (Ying et al., bioRxiv 2024) is a credible community entry point.
- **Hiring reality:** Retro's Bioinformatics Research Engineer spec lists Nextflow, scRNA-seq pipelines, Python testing/packaging, TensorFlow/PyTorch, cloud infrastructure, and database schema design; NewLimit lists a PhD-level comp-bio background, Scanpy/NumPy/SciPy/Pandas, AWS/GCP/Azure, and Docker/Singularity; Insilico emphasizes production Python plus foundation models and recruits through hackathons. Practitioner consensus (Teamblind engineers; bioinformatics career blogs) is blunt: recruiters rarely read code, a portfolio's payoff is as a credible interview story and proof of reproducible, real-data skill, and it only helps if it is non-trivial, well-documented (README + reproducible notebook), and tailored to the target company's stack. A solo toy is weakly differentiating; contributing to longevity-genie or shipping something with real users/citations carries far more weight.

## Recommendations
1. **First (highest differentiation per unit effort):** Kill the synthetic magnitudes. Rebuild on real open data (Open Genes, HAGR, recount3, GTEx-open, Tabula Muris Senis, ComputAgeBench). Ship (a) trajectory-aware gene profiles and (b) a clock diagnostics/applicability checker, both with bootstrap confidence intervals. These operationalize Shen 2024 and Higgins-Chen 2022 and fill genuine white space.
2. **Second:** Add cell-composition adjustment (EpiDISH/EPISCORE) and per-result reproducibility tiers (I², sign concordance, Egger's, trim-and-fill). Add the DNB/l-DNB resilience module with an honest longitudinal validation on Dog Aging Project or Mouse Phenome Database data.
3. **Third:** Signature-reversal endpoint wrapping L2S2; hallmarks and Open Targets/GWAS overlays; a DuckDB-WASM demo slice.
4. **Cut:** synthetic magnitudes; network control energy; the standalone React showcase; universal-atlas scope; rebuilding clocks; a bespoke MCP server (contribute to longevity-genie instead); streaming ingestion, auth/multi-tenancy, and GraphQL in the MVP.
5. **Meta-analytic method:** move from DerSimonian–Laird to **REML with the Hartung–Knapp adjustment** (better small-study coverage), keep Hedges' g but also report on the log-fold-change scale where appropriate, and always surface heterogeneity and publication-bias diagnostics.
6. **Benchmarks that change the plan:** if the cross-sectional resilience proxy fails validation against a longitudinal cohort, drop the "resilience" framing entirely and keep only DNB as an exploratory tool. If the clock-diagnostics service gains real users (GitHub stars, external citations), pursue a Bioinformatics Applications Note; if it doesn't, keep it as a preprint + portfolio artifact.

## Caveats
- The "aging in waves at 44/60" result is from a single 108-person cohort; treat it as motivating, not settled.
- Cross-sectional critical-slowing-down as an aging biomarker is contested and largely unvalidated; this report treats it strictly as a hypothesis-generating view.
- Hiring-signal quotes come largely from general software/data-science practitioners, not aging-biotech hiring managers specifically; treat them as directional.
- Some cited pages are secondary (news, aggregators); primary papers are named wherever possible, and the web_search budget was exhausted before the exemplar-portal-UX and REML/Hartung–Knapp searches could be corroborated with a second source, so those recommendations rest on established methodological consensus rather than a freshly fetched citation.
