# The perturbation layer: what switching a gene off actually tells you

GeroQuery's panel is observational. It reads how a gene's expression differs
between young and old tissue, which is an association and says nothing about
what the gene *does*. Perturbation data is the obvious complement: knock the
gene down and watch what follows.

This layer adds that data, and reports what it supports — which is less than it
first appears. Everything below is computed in
`geroquery/etl/build_perturbation.py` and
`geroquery/etl/perturbation_benchmark.py` from checksum-pinned public files, and
is reproducible with `python -m geroquery.etl.build_perturbation`.

## The data

Replogle et al. (*Cell* 2022) knocked down one gene at a time with CRISPRi in
human cell lines and sequenced the transcriptome of the result. Two pseudobulk
screens are pinned, CC BY 4.0:

| Screen | Cell line | Knockdowns | Readout genes | With a measured effect |
|---|---|---|---|---|
| `replogle_rpe1` | RPE1, retinal pigment epithelium (not a cancer line) | 2,679 | 8,749 | 2,174 |
| `replogle_k562_gwps` | K562, chronic myeloid leukaemia (genome-wide) | 11,258 | 8,248 | 5,183 |

A knockdown is scored only when it targets a gene with an Ensembl id (so it can
be joined to the panel), is not a non-targeting control, and the screen's own
energy test rejects "nothing changed". The 3,105 RPE1 and 6,075 K562 rows that
fail are not scored: correlating a non-response with the aging direction would
produce a number with no content.

**Neither cell line is aged tissue.** Both are immortalized dividing cells, one
of them a leukaemia line. Nothing here measures aging; it measures what a cell
line does when a gene is removed.

## The positive control

Knocking a gene down must lower that gene.

| Screen | Median z of the target's own readout | Fraction that fall |
|---|---|---|
| RPE1 | −0.65 | 99.3% |
| K562 | −0.50 | 98.6% |

Against an off-target spread of about 0.29. This is pinned in the test suite,
and it matters more than usual here: most of this layer's results are negative,
and the control is what distinguishes "measured, and found little" from "loaded
the file wrong".

## 1. Which aging genes are even testable

Of the 3,503 human genes GeroQuery resolves a direction for:

| | Genes |
|---|---|
| Knocked down, and the cell measurably responded | **995** |
| Knocked down, and nothing measurable happened | 821 |
| Measured as a readout gene, never itself knocked down | 335 |
| Absent from both screens | 1,352 |

So for 1,352 resolved aging genes there is no perturbation evidence of this kind
at all, and for another 821 the experiment was done and returned nothing. This
is the layer's most durable contribution: a map of what can and cannot currently
be asked.

## 2. Toward or away from the aging signature

For each knockdown, its response across the readout genes is compared with the
panel's pooled aging direction over the same genes — a cosine, so a knockdown
whose downstream shifts mirror the age-related ones scores positive.

Two corrections are necessary, and both came from the data:

**The perturbed gene is excluded from its own score.** Its own readout is the
knockdown working, not a consequence of it. (Effect on the ranking: r = 0.9999.
It is excluded because leaving it in cannot be defended, not because it matters.)

**The score is adjusted for how hard the knockdown hit.** Raw alignment
correlates **+0.53 (RPE1) and +0.41 (K562)** with the size of the response, and
the average knockdown scores positive rather than zero. Knocking out anything
important makes a cell sick, and a sick cell's transcriptome resembles old
tissue, so the raw score ranks damage and would present that as biology. The
adjustment regresses alignment on log magnitude and standardizes the residual.

The adjusted score answers a deliberately narrow question: *among knockdowns
that disturbed the cell about this much, was this one unusually aging-like?*

### What the ranking is worth

| | RPE1 | K562 |
|---|---|---|
| Alignment spread, non-targeting controls | 0.072 | 0.072 |
| Alignment spread, real knockdowns | 0.108 | 0.126 |

Real knockdowns do spread wider than guides that target nothing, so the score
is reading something. But across the 1,973 genes scored in both screens:

- rank correlation between the two cell lines: **ρ = +0.14**
- they agree on the *sign* for **54.6%** of genes — a coin flip
- **no gene** reaches an FDR-significant alignment in either screen

Only the extremes reproduce. The 300 most extreme genes agree at ρ = +0.36, and
the top-50 lists overlap about four times above chance. The exporter therefore
reports a gene as replicating only when **both** screens rank it in their most
extreme 10% **and** agree on direction: **32 genes**, out of 1,973 scored twice.

Those 32 are the only part of this analysis worth treating as a lead, and a lead
is all they are: two cell lines, neither of them aged tissue, one statistic.

## 3. Do pooled aging claims hold under perturbation?

If genes that rise with age were driving the aging transcriptome, knocking them
down should move cells away from it — a negative association between a gene's
own pooled effect and its knockdown's alignment.

| Screen | Genes | Spearman ρ | p |
|---|---|---|---|
| RPE1 | 383 | −0.040 | 0.43 |
| K562 | 968 | −0.092 | 0.004 |

The sign is as predicted in both and the K562 association clears p = 0.05, but
ρ = −0.09 is negligible: a gene's aging direction tells you essentially nothing
about what its knockdown does. The honest reading is that these are two
different measurements, and the panel's claims are neither confirmed nor
contradicted by this data.

## 4. What a knockdown predictor is actually worth

The "virtual cell" demo trains on genes that were knocked down, holds one out,
predicts what its knockdown does, and reports how many of the biggest movers
went the right way. The numbers are impressive. The control is usually missing.

Most knockdowns of a consequential gene produce a similar response — the cell
stalls and stresses — so a model that **always predicts that one average
profile**, having learned nothing about which gene was hit, already gets most
directions right.

Scoring the top 25 movers, 5-fold cross-validation, ridge penalty chosen on an
inner split of the training genes only:

| Screen | Baseline: always predict the average knockdown | Trained model | Gain |
|---|---|---|---|
| RPE1 (n = 2,050) | **85.2%** | 85.9% | **+0.71 pp** (p = 2e−13) |
| K562 (n = 4,436) | **76.1%** | 76.7% | **+0.66 pp** (p = 5e−7) |

The gain is real and reproduces in both screens. It is also about **one
percentage point**, against a baseline that ignores the input entirely. Roughly
99% of such a demo's headline accuracy is the shared average response.

For strong-effect genes the baseline alone exceeds 95%, so "it got all 25 right
on a gene it had never seen" is the expected outcome of learning nothing.

Two cautions against over-reading this:

- It is scoped to this setup — these public pseudobulk screens, this feature
  choice (a gene's own column: how it responds when other genes are hit), ridge
  and nearest-neighbour models. It is not a claim about every virtual-cell
  effort, several of which use far richer data.
- The benchmark can detect a gain when one exists. That is pinned by a test
  that plants gene-specific signal and requires the benchmark to find it, so the
  small number here is a measurement rather than a dead pipeline.

**Nothing from this section is shipped into the product.** A predictor that
cannot clear "ignore the input" by a useful margin has no business serving
answers beside measured ones. The benchmark is kept because the comparison is
the useful artifact.

## How this layer is kept honest

- It is a separate surface. No modelled number is mixed into a pooled aging
  estimate, and no perturbation result changes a gene's verdict.
- Every figure is computed in Python and shipped; the browser renders them.
- The negative results above are pinned by tests, so a future change that
  quietly turns a coin flip into a finding fails the suite.

## Reproducing

```
python -m geroquery.etl.build_perturbation          # the layer and its artifact
python -m geroquery.etl.perturbation_benchmark      # section 4 (~3 minutes)
```

Both fetch the pinned files on first run and verify them against the checksums
in `geroquery/sources/manifest.py`.

## Attribution

Replogle JM, Saunders RA, Pogson AN, et al. Mapping information-rich
genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*.
2022;185(14):2559–2575.e28. Processed data CC BY 4.0.
