# The field, and where GeroQuery fits (Nov 2026)

A survey of tools near GeroQuery's job, what it changed in this release, and
what stands between the project and commercial use. Sources are listed at the
end; claims that come only from secondary write-ups are marked as such.

## What GeroQuery does that the others do not

GeroQuery re-derives a gene's relationship to age from the published
expression data: 485,905 per-study effect sizes from 32 contrasts, pooled per
gene by random-effects meta-analysis. Each gene gets an interval, a prediction
interval for the next study, and a probability that its direction is wrong.
Everything else in the space either **asserts** (curated lists) or **shows one
dataset at a time** (expression browsers). None of them tells you how a curated
claim holds up when the raw data are pooled, or whether a new study would see
the same thing.

## The landscape

| Kind | Examples | What they give you | What they do not |
|---|---|---|---|
| Curated knowledge bases | **GenAge / HAGR** (307 human genes; record updated Sep 2026), **Open Genes** (2,402 genes, NAR 2023), **AgingHallmarksDB** (3,111 genes from 7 databases, Oct 2026), **Aging Atlas** (latest dataset added Jul 2026) | A list, a label, literature links | A measured effect, an interval, or a "we cannot tell" |
| Expression browsers | **voyAGEr** and **Age Effect Explorer** (GTEx tissues), **AgeMeta** (122 datasets, 2024), **Aging Signature** (mouse HSCs) | Per-tissue age trends, mostly from one consortium's data | Cross-study pooling with stated uncertainty |
| Single-cell atlases | **Tabula Muris Senis** (23 mouse tissues), **HuTAge** (GTEx + Tabula Sapiens), **CELLxGENE Census** | Cell-type resolution | Comparability across datasets without batch correction |
| Multi-species meta-analyses | A 2026 Nature study integrating >11,000 transcriptomes across four mammals into transcriptomic clocks *(seen only via a community summary; not verified against the paper)* | Conserved hallmarks, clocks | A per-gene, queryable answer |
| Clocks and biomarkers | **Biomarkers of Aging Consortium** open challenge (>550 clocks from 37 teams; established clocks err by ~5–8.5 years on the shared data), **biolearn**, organ-specific **plasma proteomic clocks** in UK Biobank (Nature Aging, 2025–26) | Biological-age estimates | Gene-level evidence |
| Agents and MCP | **BioMCP**, **BioContextAI**, **gget-mcp**, the **MCPmed** proposal (GEO, STRING, UCSC over MCP), **Biomni**, NIH **GeneAgent** | Natural-language access to many databases | Aging-specific, uncertainty-aware answers |
| Perturbation and "virtual cell" | **scPerturb**, the Replogle/Weissman Perturb-seq releases, and the several virtual-cell efforts training on them | What a cell does when one gene is switched off; predicted knockdown profiles | A link to aging, and usually a published comparison against the obvious baseline |
| Commercial | **Generation Lab** (SystemAge, $11M seed, Oct 2025), **Deep Longevity**, **BioAge** (Nasdaq), Longevity.Technology's **DLT** data platform | Diagnostics and market intelligence | Open, reproducible gene-level evidence |

Two field-wide problems recur in these sources and shape what was built:

1. **Cross-study heterogeneity is the main barrier to replication.**
   Meta-analyses of worm and mammal transcriptomes report that lab-to-lab
   variance survives batch correction. GeroQuery's prediction interval is built
   to state exactly that.
2. **Over-the-counter enrichment invites misuse.** AgingHallmarksDB warns that
   a pre-filtered list of significant genes breaks over-representation
   statistics. GeroQuery's list check therefore reports rates, not p-values.

## What this release adds, and why

| Change | The need it answers |
|---|---|
| **Gene list** view: paste up to 2,000 symbols, get each gene's verdict, interval, prediction interval and P(wrong sign), the list's rate beside the species rate, a CSV, and a link that reproduces the check | Researchers and students arrive with a *list* (a hit list, a pathway, a candidate panel), not one gene. Every curated base answers lists with a membership tick. This answers them with evidence. |
| **GenAge against the pooled evidence** (`curated.json`, computed in Python, test-first) | The question no other tool answers: how do curated aging genes fare in the raw data? |
| **Joining on Entrez, not symbol** | GenAge still lists four genes under retired HGNC symbols (ARNTL→BMAL1, ATP5O→ATP5PO, CTGF→CCN2, H2AFX→H2AX). A symbol join called them "not measured". |
| **Spreadsheet-date detection** in pasted lists | Spreadsheet software silently turns SEPT1/MARCH1 into dates, a long-documented corruption of published gene lists. A `1-Sep` is named as such, not filed as "not measured". |
| **Copy citation** on the gene page: the estimate as shown, the data version, and a permalink | A result that cannot be cited cannot be used in a paper or a thesis |
| **Contrasts CSV** per gene | Supplementary tables, teaching, re-analysis |
| **Links out** to AlphaFold, GenAge (only for genes it lists, under the name it lists them by), GTEx and NCBI Gene | GeroQuery answers one question. The links hand off the next ones without restating them. Open Genes is not linked: its pages answer 200 for any symbol, so a link cannot be verified to lead anywhere. |
| **⌘K / Ctrl+K** jumps to gene search from any view, beside the existing `/` | The command-key convention from Linear and Notion |
| **Perturbation layer** (CRISPRi screens, CC BY 4.0): coverage of which aging genes have ever been switched off, an alignment score adjusted for knockdown damage, cross-cell-line replication, and a published prediction benchmark | The panel is observational and cannot say what a gene *does*. Perturbation data can, in principle — and the layer reports, with the measurements behind it, how little of that promise this public data currently delivers |

The API (`POST /v1/geneset/signature`) and the MCP tool
(`geneset_aging_signature`) already answered list queries for agents. The new
view brings the same question to people who do not script.

## The result it surfaced

Of GenAge's 307 human aging genes, **all 307 are measured** in the panel now
that it includes RNA-seq (ten were missing from the microarray-only panel, among
them TERC, an RNA gene). **70 (23%) have a pooled direction with age**, against
13% of all 26,584 human genes. **Only 5 (2%) clear the prediction interval**,
against 2%. A curated aging gene is nearly twice as likely to resolve as a
random one, three in four still do not, and almost none is consistent enough
across studies that a new one would be expected to see it.

Read it carefully. GenAge admits genes for many reasons ("mammal", "model",
"cell", "putative", "upstream"), most of which make no claim about expression
changing with age. "No direction here" does not contradict GenAge's evidence;
it says expression change with age is not that evidence. The list view shows
each gene's GenAge evidence code beside its verdict so a reader can make that
distinction gene by gene.

## Commercial use

The code is MIT. The data are not uniformly licensed, so a commercial
deployment needs this audit first:

| Upstream | Pinned artifacts | Terms | Commercial consequence |
|---|---|---|---|
| NCBI GEO | 33 | US public domain; attribute. NCBI notes that submitters may assert their own rights in some records | Usable; keep attribution, and check any series whose submitter asserted rights |
| NHANES (NCHS) | 15 | US public domain | Usable |
| HAGR (GenAge, CellAge, LongevityMap, DrugAge, GenDR) | 7 | Free for **non-commercial** use, attribution required | **Must be licensed or removed** in a paid product. That covers the GenAge comparison, the overview finding built on it, and the GenAge links/evidence column |
| ARCHS4 (processed RNA-seq counts) | 2 files, pinned by release + size | The licence of the processed files is not stated; the underlying GEO/SRA samples are public | Confirm terms with the Ma'ayan Laboratory before a paid product ships the ARCHS4 contrasts |
| Perturb-seq (Replogle 2022, processed) | 2 | **CC BY 4.0** | Usable commercially; attribute Replogle et al. |
| Clocks (via biolearn / pyaging) | — | Per clock. Some widely used clocks are distributed by their developers under academic or non-commercial terms | Check each clock's terms before selling clock output |

A model consistent with those terms: the static site, API and MCP server stay
free for academic and student use. A commercial tier offers the same evidence
without the HAGR layer (or under a HAGR licence), a hosted API with an uptime
commitment, and private list checks. The browser check already sends nothing to
a server, which is the property a company's unpublished target list needs.

## What to build next, in order

1. ~~**Widen the panel.**~~ **Done (2026-11-14):** 244 rule-selected RNA-seq
   contrasts from ARCHS4 took the panel from 32 to 276 contrasts and narrowed
   the median interval width from 0.84 to 0.48; see `evidence/CHANGELOG.md`.
2. ~~**Account for shared subjects.**~~ **Done (2026-11-14):** contrasts are
   combined per GEO Series (rho = 1, the conservative bound) before pooling, so
   k counts independent Series. A three-level model that estimates the
   within-Series correlation instead of bounding it is the refinement.
3. **A three-level or mixed model for the perturbation score**, estimating the
   shared damage axis rather than regressing it out, and a third cell line.
   Two screens agreeing 54.6% of the time is too thin a base for the ranking to
   be worth publishing as more than leads; see `PERTURBATION.md`.
4. **A single-cell layer** from CELLxGENE Census / Tabula Muris Senis, kept as a
   separate pool per cell type (cross-dataset comparison there needs batch
   correction; it must not be pooled with bulk).
5. **Open Genes as a pinned snapshot**: its API returns JSON per symbol
   (`open-genes.com/api/gene/{symbol}`), so it can be checksum-pinned like HAGR,
   which would make its links verifiable.
6. The deferred **Schoenfeld proportional-hazards diagnostic** for the mortality
   layer, and the **Vite 5→8** upgrade.

## Sources

- Open Genes — NAR 2024: <https://academic.oup.com/nar/article/52/D1/D950/7259277>
- GenAge KG-Registry record: <https://kghub.org/kg-registry/resource/genage/genage.html>
- Human Ageing Genomic Resources, 2024 update: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10767973/>
- AgingHallmarksDB: <https://cb.imsc.res.in/aginghallmarksdb> · coverage: <https://bioengineer.org/new-database-maps-3111-genes-across-the-11-hallmarks-of-aging/>
- Aging Atlas: <https://ngdc.cncb.ac.cn/aging/index>
- voyAGEr (eLife): <https://elifesciences.org/articles/88623>
- Age Effect Explorer: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12889165/>
- AgeMeta: <https://link.springer.com/article/10.1134/S000629792402010X>
- Aging Signature (Groningen): <https://agingsignature.web.rug.nl/>
- Tabula Muris Senis (Nature 2020): <https://www.nature.com/articles/s41586-020-2496-1>
- HuTAge (bioRxiv): <https://www.biorxiv.org/content/10.1101/2024.08.07.606269.full.pdf>
- Multi-species transcriptomic hallmarks (community summary): <https://jadams5.github.io/aging_wiki/studies/tyshkovskiy-2026-universal-transcriptomic-hallmarks>
- Cross-tissue meta-analysis, Palmer et al. 2021: <https://www.aging-us.com/article/202648/text>
- C. elegans meta-analysis, inter-lab variance: <https://doaj.org/article/07979b06f16f479880fec58da412b7c9>
- Biomarkers of Aging challenge: <https://pmc.ncbi.nlm.nih.gov/articles/PMC13476092/> · <https://www.science.org/content/article/scientific-showdown-seeks-biological-clock-best-tracks-aging>
- Biolearn framework: <https://pmc.ncbi.nlm.nih.gov/articles/PMC13476091/>
- Organ-specific proteomic clocks (Nature Aging): <https://pubmed.ncbi.nlm.nih.gov/41299092/>
- MCPmed: <https://arxiv.org/html/2507.08055v1>
- BioContextAI knowledgebase MCP: <https://glama.ai/mcp/servers/@biocontext-ai/knowledgebase-mcp>
- gget-mcp: <https://archestra.ai/mcp-catalog/longevity-genie__gget-mcp>
- GeneAgent (NIH): <https://www.nextgov.com/artificial-intelligence/2025/07/nih-launches-new-ai-agent-analyzing-genomic-data/407034>
- Biomni Lab / GeneCards: <https://labcritics.com/?p=13978>
- Generation Lab seed round: <https://longevity.technology/news/generation-lab-raises-11m-to-advance-longevity-diagnostics/>
- Longevity.Technology DLT: <https://longevity.technology/news/longevity-technology-unveils-dlt-the-intelligence-platform-for-longevity/>
- HGNC multi-symbol checker: <https://www.genenames.org/tools/multi-symbol-checker/>
- Spreadsheet gene-name errors: Ziemann, Eren & El-Osta, *Genome Biology* 17:177 (2016)
