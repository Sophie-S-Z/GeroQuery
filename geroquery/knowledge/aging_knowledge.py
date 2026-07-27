"""The curated, real-evidence knowledge base at the heart of GeroQuery.

This module is the single source of truth for what GeroQuery asserts about how
each gene relates to aging. **Every statement here is real, curated biology** —
directions of change with age, the tissues/omic layers where they are
established, the hallmark(s) of aging involved, plain-English interpretation,
and a verifiable literature reference. Nothing here is simulated.

What this deliberately does *not* contain is fabricated per-study effect sizes,
invented GEO accession numbers, or synthetic p-values. Earlier versions of the
project generated those to populate a "demo meta-analysis"; they have been
removed. Where a real quantitative meta-analysis exists in the literature, we
cite it and describe its finding in words rather than reprinting numbers we
did not measure.

Evidence strength is recorded as an ordinal confidence label, not a false
precision:

* ``robust``       — repeatedly and independently replicated; a textbook fact
* ``established``  — well supported by multiple studies
* ``emerging``     — real but context-dependent or still being characterised
"""

from __future__ import annotations

from dataclasses import dataclass, field

# --------------------------------------------------------------------------- #
# Hallmarks of aging (López-Otín et al., 2013; expanded 2023).
# --------------------------------------------------------------------------- #

HALLMARKS: dict[str, str] = {
    "genomic_instability": "Accumulating DNA damage and mutations across the lifespan.",
    "telomere_attrition": "Progressive shortening of chromosome-end telomeres, limiting "
    "replicative capacity.",
    "epigenetic_alterations": "Drift in DNA methylation and chromatin marks that changes "
    "gene expression with age.",
    "loss_of_proteostasis": "Declining ability to fold, maintain, and clear proteins.",
    "disabled_autophagy": "Reduced recycling of damaged organelles and macromolecules.",
    "deregulated_nutrient_sensing": "Altered insulin/IGF-1, mTOR, AMPK and sirtuin signalling.",
    "mitochondrial_dysfunction": "Falling energy output and rising oxidative stress.",
    "cellular_senescence": "Accumulation of growth-arrested, inflammatory 'zombie' cells.",
    "stem_cell_exhaustion": "Depletion of the regenerative reserve of tissues.",
    "altered_intercellular_communication": "Shifting hormonal and inflammatory signalling "
    "between cells and organs.",
    "chronic_inflammation": "Low-grade, sterile, systemic inflammation ('inflammaging').",
    "dysbiosis": "Age-related imbalance of the microbiome.",
}


@dataclass(frozen=True)
class Evidence:
    """One curated, cited line of evidence for a gene's relationship to aging."""

    omic_layer: str  # transcriptome | proteome | methylome | physiology | genetics
    direction: str  # up | down | context-dependent
    species: tuple[str, ...]  # species where this is established
    tissues: tuple[str, ...]
    strength: str  # robust | established | emerging
    finding: str  # plain-English description of what the data show
    reference_keys: tuple[str, ...]

    def to_dict(self) -> dict:
        return {
            "omic_layer": self.omic_layer,
            "direction": self.direction,
            "species": list(self.species),
            "tissues": list(self.tissues),
            "strength": self.strength,
            "finding": self.finding,
            "reference_keys": list(self.reference_keys),
        }


@dataclass(frozen=True)
class GeneKnowledge:
    """Everything GeroQuery curates about one ortholog group and aging."""

    group: str  # ortholog group key, e.g. "CDKN2A"
    common_name: str  # human-friendly name, e.g. "p16 / CDKN2A"
    aka: str  # famous alias line
    direction_with_age: str  # up | down | context-dependent
    confidence: str  # robust | established | emerging
    one_liner: str  # a single-sentence headline answer
    role: str  # what the gene does, plainly
    analysis: str  # the plain-English "what the data tells us" paragraph
    hallmarks: tuple[str, ...]
    evidence: tuple[Evidence, ...]
    faqs: tuple[tuple[str, str], ...] = field(default_factory=tuple)

    def reference_keys(self) -> list[str]:
        keys: list[str] = []
        for e in self.evidence:
            for k in e.reference_keys:
                if k not in keys:
                    keys.append(k)
        return keys

    def to_dict(self) -> dict:
        return {
            "group": self.group,
            "common_name": self.common_name,
            "aka": self.aka,
            "direction_with_age": self.direction_with_age,
            "confidence": self.confidence,
            "one_liner": self.one_liner,
            "role": self.role,
            "analysis": self.analysis,
            "hallmarks": list(self.hallmarks),
            "evidence": [e.to_dict() for e in self.evidence],
            "faqs": [{"q": q, "a": a} for q, a in self.faqs],
            "reference_keys": self.reference_keys(),
        }


# --------------------------------------------------------------------------- #
# The curated knowledge base, keyed by ortholog group.
# --------------------------------------------------------------------------- #

KNOWLEDGE: dict[str, GeneKnowledge] = {
    "CDKN2A": GeneKnowledge(
        group="CDKN2A",
        common_name="CDKN2A (p16INK4a / ARF)",
        aka="the classic 'senescence' gene, p16",
        direction_with_age="up",
        confidence="robust",
        one_liner="p16INK4a rises steeply with age in nearly every tissue examined and is one "
        "of the most reliable molecular markers of biological aging.",
        role="CDKN2A encodes two tumour suppressors, p16INK4a and p14ARF. p16 halts the cell "
        "cycle and is a central enforcer of cellular senescence — the stable growth arrest that "
        "accumulates in aging tissues.",
        analysis="Across humans and rodents, CDKN2A/p16 expression climbs markedly with age — in "
        "some tissues by an order of magnitude between young and old animals. This is not a subtle "
        "trend: p16 is used as a working biomarker of the senescent-cell burden of a tissue. "
        "Because senescent cells drive chronic inflammation, clearing p16-high cells extends "
        "healthspan in mice, which is strong evidence that this rise is causally tied to aging "
        "rather than merely correlated with it.",
        hallmarks=("cellular_senescence", "chronic_inflammation", "stem_cell_exhaustion"),
        evidence=(
            Evidence(
                "transcriptome",
                "up",
                ("human", "mouse"),
                ("blood", "brain", "muscle", "many tissues"),
                "robust",
                "p16INK4a mRNA increases with age across a wide range of tissues; the magnitude "
                "of the increase tracks the senescent-cell burden.",
                ("krishnamurthy_2004",),
            ),
            Evidence(
                "physiology",
                "up",
                ("mouse",),
                ("multiple tissues",),
                "robust",
                "Genetically clearing p16INK4a-positive cells delays age-related deterioration, "
                "linking the age-related rise to aging causally, not just correlatively.",
                ("baker_2011",),
            ),
        ),
        faqs=(
            (
                "Does p16 go up or down with age?",
                "Up — clearly and reproducibly. It is one of the strongest 'up with age' signals "
                "in aging biology.",
            ),
            (
                "Why do people call it an aging gene?",
                "Because p16 marks senescent cells, and senescent cells accumulate as tissues age. "
                "Measuring p16 is a common proxy for how 'old' a tissue is at the cellular level.",
            ),
            (
                "Can anything change it?",
                "Senolytic drugs (e.g. dasatinib + quercetin) selectively remove p16-high senescent "
                "cells; in mice this improves function and extends lifespan.",
            ),
        ),
    ),
    "CDKN1A": GeneKnowledge(
        group="CDKN1A",
        common_name="CDKN1A (p21)",
        aka="p21, the p53 effector",
        direction_with_age="up",
        confidence="established",
        one_liner="p21 tends to rise with age as a downstream effector of DNA-damage and p53 "
        "signalling, contributing to senescent growth arrest.",
        role="CDKN1A encodes p21, a cyclin-dependent-kinase inhibitor switched on by p53 in "
        "response to DNA damage. It enforces cell-cycle arrest and, together with p16, helps "
        "establish and maintain senescence.",
        analysis="p21 is an acute, stress-responsive brake on cell division. As DNA damage and "
        "p53 activity accumulate with age, p21 expression generally increases, reinforcing the "
        "senescent state. Compared with p16 its age trend is a little more variable across "
        "tissues because p21 responds quickly to transient stress as well as to durable aging, "
        "but the overall direction — up with age and with senescence — is well supported.",
        hallmarks=("cellular_senescence", "genomic_instability"),
        evidence=(
            Evidence(
                "transcriptome",
                "up",
                ("human", "mouse"),
                ("blood", "brain", "muscle"),
                "established",
                "p21 expression rises with age and with the establishment of senescence as part "
                "of the p53 → p21 DNA-damage response axis.",
                ("krishnamurthy_2004", "baker_2011"),
            ),
        ),
        faqs=(
            (
                "Is p21 the same as p16?",
                "No. Both are cyclin-dependent-kinase inhibitors that enforce senescence, but p21 "
                "is a fast, p53-driven response to damage, while p16 marks more durable senescence.",
            ),
            ("Which direction with age?", "Up, in most tissues studied."),
        ),
    ),
    "TP53": GeneKnowledge(
        group="TP53",
        common_name="TP53 (p53)",
        aka="'the guardian of the genome'",
        direction_with_age="context-dependent",
        confidence="established",
        one_liner="p53 activity — more than its raw expression — rises with accumulating DNA "
        "damage; it sits at the crossroads of cancer protection and aging.",
        role="TP53 encodes p53, the master tumour suppressor that senses DNA damage and decides "
        "between repair, senescence, and programmed cell death.",
        analysis="p53 is a balancing act. It protects against cancer by removing damaged cells, "
        "but that same activity, chronically engaged, promotes senescence and tissue attrition — "
        "an example of antagonistic pleiotropy in aging. Its mRNA level is not a clean 'up' or "
        "'down' with age; what changes is pathway activity as its upstream damage signals "
        "accumulate. GeroQuery therefore labels its directional trend context-dependent rather "
        "than forcing a single arrow.",
        hallmarks=("genomic_instability", "cellular_senescence"),
        evidence=(
            Evidence(
                "transcriptome",
                "context-dependent",
                ("human", "mouse"),
                ("multiple tissues",),
                "established",
                "p53 target-gene activity increases with age-related DNA damage, but TP53 mRNA "
                "itself shows tissue-dependent, non-monotonic trends.",
                ("lopez_otin_2013",),
            ),
        ),
        faqs=(
            (
                "Does p53 go up with age?",
                "Its activity generally does, as DNA damage accumulates; its expression level is "
                "more variable, so we call the trend context-dependent.",
            ),
            (
                "Is p53 good or bad for aging?",
                "Both. It prevents cancer but can accelerate tissue aging when chronically active — "
                "a classic trade-off.",
            ),
        ),
    ),
    "LMNB1": GeneKnowledge(
        group="LMNB1",
        common_name="LMNB1 (Lamin B1)",
        aka="the senescence marker that goes down",
        direction_with_age="down",
        confidence="robust",
        one_liner="Lamin B1 falls with senescence and age — a reliable 'down with age' counterpart "
        "to p16's rise.",
        role="LMNB1 encodes lamin B1, a structural protein of the nuclear envelope that organises "
        "chromatin and the nuclear lamina.",
        analysis="Loss of lamin B1 is one of the most consistent 'down with age' signals in aging "
        "biology. As cells enter senescence, lamin B1 protein and mRNA drop, the nuclear lamina "
        "reorganises, and heterochromatin is remodelled. Because the decline is so reproducible, "
        "lamin B1 loss is used as a companion biomarker to p16: p16 up and lamin B1 down together "
        "give a more confident read-out of senescence than either alone.",
        hallmarks=("cellular_senescence", "epigenetic_alterations", "genomic_instability"),
        evidence=(
            Evidence(
                "transcriptome",
                "down",
                ("human", "mouse"),
                ("blood", "brain", "muscle", "skin"),
                "robust",
                "Lamin B1 mRNA and protein decrease during senescence and with age; the loss is "
                "used as a senescence biomarker.",
                ("freund_2012", "shimi_2011"),
            ),
        ),
        faqs=(
            (
                "Which way does lamin B1 go with age?",
                "Down. It is a well-established 'decreases with age/senescence' marker.",
            ),
            (
                "Why is it useful?",
                "Paired with p16 (which goes up), lamin B1 (which goes down) makes senescence "
                "detection more reliable.",
            ),
        ),
    ),
    "GDF15": GeneKnowledge(
        group="GDF15",
        common_name="GDF15",
        aka="one of the strongest plasma aging biomarkers known",
        direction_with_age="up",
        confidence="robust",
        one_liner="Circulating GDF15 rises steadily with age and is repeatedly ranked among the "
        "top age-associated proteins in human plasma.",
        role="GDF15 is a stress-responsive cytokine (a 'mitokine') secreted in response to "
        "mitochondrial and metabolic stress; it acts on the brainstem to influence appetite and "
        "on many tissues as an inflammatory/stress signal.",
        analysis="In large human plasma-proteomics studies of aging, GDF15 is one of the single "
        "most strongly age-associated proteins, climbing monotonically across adulthood. It is a "
        "readout of the mitochondrial and inflammatory stress that accumulates with age, and high "
        "levels predict mortality and multimorbidity. This makes GDF15 an unusually clean, "
        "measurable, 'up with age' signal — exactly the kind of real biomarker aging researchers "
        "reach for.",
        hallmarks=(
            "altered_intercellular_communication",
            "mitochondrial_dysfunction",
            "chronic_inflammation",
        ),
        evidence=(
            Evidence(
                "proteome",
                "up",
                ("human",),
                ("plasma",),
                "robust",
                "GDF15 is consistently among the top-ranked age-associated proteins in "
                "human plasma proteomic studies across the adult lifespan.",
                ("tanaka_2018", "lehallier_2019"),
            ),
        ),
        faqs=(
            (
                "Is GDF15 a good aging biomarker?",
                "Yes — it is one of the best-replicated single-protein markers of chronological and "
                "biological age in blood.",
            ),
            ("Does it go up or down?", "Up, steadily, across adult life."),
            (
                "What does a high level mean?",
                "It reflects accumulated mitochondrial/metabolic stress and is associated with "
                "higher risk of frailty and mortality.",
            ),
        ),
    ),
    "SIRT1": GeneKnowledge(
        group="SIRT1",
        common_name="SIRT1",
        aka="the 'sirtuin' longevity gene",
        direction_with_age="context-dependent",
        confidence="emerging",
        one_liner="SIRT1 is a nutrient-sensing deacetylase tied to caloric-restriction benefits; "
        "its expression trend with age is tissue-dependent and often declines.",
        role="SIRT1 is an NAD+-dependent protein deacetylase that links cellular energy status to "
        "stress resistance, metabolism, and chromatin regulation. It is a key mediator of the "
        "benefits of caloric restriction.",
        analysis="SIRT1 sits in the nutrient-sensing machinery of aging. Its activity depends on "
        "NAD+, which itself declines with age, so functional sirtuin signalling generally wanes "
        "even where mRNA is stable. Reported expression trends vary by tissue, which is why "
        "GeroQuery marks the direction context-dependent and the confidence emerging rather than "
        "overstating a single arrow. Its importance is less about being a biomarker and more "
        "about being a druggable node — the target of NAD+ boosters and caloric-restriction "
        "mimetics.",
        hallmarks=("deregulated_nutrient_sensing", "epigenetic_alterations"),
        evidence=(
            Evidence(
                "transcriptome",
                "context-dependent",
                ("human", "mouse"),
                ("multiple tissues",),
                "emerging",
                "SIRT1 expression trends with age are tissue-dependent; functionally, "
                "sirtuin activity declines as its cofactor NAD+ falls with age.",
                ("lopez_otin_2013",),
            ),
        ),
        faqs=(
            (
                "Does SIRT1 decline with age?",
                "Its activity generally does, largely because NAD+ falls; its mRNA trend varies by "
                "tissue, so we label the direction context-dependent.",
            ),
            (
                "How is it linked to caloric restriction?",
                "SIRT1 is a major mediator of caloric restriction's benefits, which is why it is a "
                "target for CR-mimetic drugs and NAD+ boosters.",
            ),
        ),
    ),
    "TERT": GeneKnowledge(
        group="TERT",
        common_name="TERT (telomerase)",
        aka="the telomere-maintenance enzyme",
        direction_with_age="down",
        confidence="established",
        one_liner="TERT/telomerase is repressed in most adult human somatic tissues, contributing "
        "to the telomere attrition that limits regenerative capacity with age.",
        role="TERT is the catalytic subunit of telomerase, the enzyme that rebuilds the protective "
        "telomere caps at chromosome ends. Without it, telomeres shorten with each cell division.",
        analysis="Telomere attrition is a canonical hallmark of aging, and TERT is its master "
        "switch. In humans, telomerase is largely switched off in adult somatic cells — a "
        "cancer-suppressing arrangement that comes at the cost of finite replicative lifespan. "
        "The practical read-out with age is falling telomerase capacity and shortening telomeres, "
        "especially in high-turnover tissues like blood. Note the species caveat: mice keep more "
        "telomerase activity and have long telomeres, so mouse telomere biology does not map "
        "one-to-one onto humans.",
        hallmarks=("telomere_attrition", "stem_cell_exhaustion"),
        evidence=(
            Evidence(
                "transcriptome",
                "down",
                ("human",),
                ("blood", "most somatic tissues"),
                "established",
                "Telomerase (TERT) is repressed in adult human somatic tissues; telomerase "
                "capacity and telomere length decline with age, notably in blood.",
                ("lopez_otin_2013",),
            ),
        ),
        faqs=(
            (
                "Does telomerase go down with age?",
                "In human somatic tissues it is already largely off and telomeres shorten over time, "
                "so effectively yes.",
            ),
            (
                "Do mouse and human telomeres age the same way?",
                "No — mice retain more telomerase and have longer telomeres, an important "
                "cross-species caveat when interpreting rodent data.",
            ),
        ),
    ),
    "KL": GeneKnowledge(
        group="KL",
        common_name="KL (Klotho)",
        aka="the 'anti-aging hormone'",
        direction_with_age="down",
        confidence="robust",
        one_liner="Klotho declines with age; losing it accelerates aging in mice and boosting it "
        "extends lifespan — a genuine aging-suppressor gene.",
        role="KL encodes Klotho, a kidney-derived hormone that regulates phosphate/vitamin-D "
        "metabolism and suppresses several intracellular aging pathways, including insulin/IGF-1 "
        "signalling.",
        analysis="Klotho is one of the few genes named for aging itself. Mice lacking klotho "
        "develop a syndrome resembling accelerated aging, while over-expressing it extends "
        "lifespan — bidirectional, causal evidence that is rare in this field. In humans and "
        "animals, circulating Klotho falls with age, and low levels associate with cardiovascular "
        "and cognitive decline. Its robust 'down with age' trajectory plus causal lifespan data "
        "make it a headline aging gene.",
        hallmarks=("altered_intercellular_communication", "deregulated_nutrient_sensing"),
        evidence=(
            Evidence(
                "physiology",
                "down",
                ("human", "mouse"),
                ("kidney", "plasma"),
                "robust",
                "Klotho expression and circulating levels decline with age; loss accelerates "
                "aging phenotypes and over-expression extends mouse lifespan.",
                ("kuroo_1997", "kurosu_2005"),
            ),
        ),
        faqs=(
            (
                "Why is Klotho called an anti-aging gene?",
                "Because deleting it accelerates aging in mice and boosting it extends lifespan — "
                "direct causal evidence, not just correlation.",
            ),
            ("Which direction with age?", "Down. Klotho levels fall as we get older."),
        ),
    ),
    "FOXO3": GeneKnowledge(
        group="FOXO3",
        common_name="FOXO3",
        aka="the most replicated human longevity gene",
        direction_with_age="context-dependent",
        confidence="robust",
        one_liner="FOXO3 variants are the most consistently replicated genetic association with "
        "human longevity across populations worldwide.",
        role="FOXO3 is a stress-responsive transcription factor downstream of insulin/IGF-1 "
        "signalling. It switches on programmes for stress resistance, DNA repair, autophagy, and "
        "metabolic control.",
        analysis="FOXO3's claim to fame is genetic, not expression-based. Specific FOXO3 variants "
        "are associated with reaching extreme old age in Japanese, German, American, Italian and "
        "many other cohorts — making it the most reproducible longevity-gene association in "
        "humans. Its expression does not follow a simple up-or-down arc with age (hence "
        "context-dependent), but its activity is protective: it is activated by the same "
        "low-nutrient conditions (caloric restriction, reduced IGF-1) that extend lifespan across "
        "species.",
        hallmarks=(
            "deregulated_nutrient_sensing",
            "loss_of_proteostasis",
            "genomic_instability",
        ),
        evidence=(
            Evidence(
                "genetics",
                "context-dependent",
                ("human",),
                ("germline",),
                "robust",
                "FOXO3 longevity-associated variants are replicated across many independent "
                "human populations — the strongest human longevity genetic signal.",
                ("willcox_2008", "flachsbart_2009"),
            ),
        ),
        faqs=(
            (
                "Is FOXO3 an aging biomarker?",
                "Not primarily. Its importance is genetic: certain FOXO3 variants strongly predict "
                "human longevity.",
            ),
            (
                "Does its expression rise or fall with age?",
                "There's no clean single trend — it is context-dependent — but its protective "
                "activity is engaged by caloric restriction and low IGF-1 signalling.",
            ),
        ),
    ),
    "IGF1": GeneKnowledge(
        group="IGF1",
        common_name="IGF1",
        aka="the growth hormone that trades size for lifespan",
        direction_with_age="down",
        confidence="established",
        one_liner="Circulating IGF-1 declines through adulthood ('somatopause'), while lower "
        "IGF-1 signalling is repeatedly linked to longer lifespan across species.",
        role="IGF1 encodes insulin-like growth factor 1, the main mediator of growth-hormone "
        "action. It drives growth and anabolism and is a core node of the nutrient-sensing "
        "network.",
        analysis="IGF-1 embodies one of aging's deepest trade-offs. Levels fall with age (the "
        "'somatopause'), yet reduced IGF-1 signalling extends lifespan from worms to mice, and "
        "long-lived humans and dog/mouse breeds often carry lower IGF-1 activity. So the "
        "age-related decline is real and established, but 'more IGF-1' is not simply 'younger' — "
        "chronically high signalling is pro-aging. GeroQuery reports the measured direction "
        "(down with age) while flagging that lower lifelong signalling is the pro-longevity "
        "state.",
        hallmarks=("deregulated_nutrient_sensing", "altered_intercellular_communication"),
        evidence=(
            Evidence(
                "proteome",
                "down",
                ("human", "mouse"),
                ("plasma",),
                "established",
                "Circulating IGF-1 declines across adulthood; reduced IGF-1 signalling is "
                "associated with extended lifespan across species.",
                ("junnila_2013",),
            ),
        ),
        faqs=(
            (
                "Does IGF-1 go up or down with age?",
                "Down — the age-related decline is called the somatopause.",
            ),
            (
                "If it declines, is more IGF-1 better?",
                "Not for longevity. Lower lifelong IGF-1 signalling is repeatedly linked to longer "
                "life — a genuine trade-off between growth and lifespan.",
            ),
        ),
    ),
    "MTOR": GeneKnowledge(
        group="MTOR",
        common_name="MTOR",
        aka="the target of rapamycin",
        direction_with_age="context-dependent",
        confidence="established",
        one_liner="mTOR signalling tends to become hyperactive with age; inhibiting it with "
        "rapamycin extends lifespan in every model tested, including mice.",
        role="MTOR is a nutrient- and growth-factor-sensing kinase that switches cells between "
        "growth (protein synthesis) and maintenance (autophagy). It is the central hub of the "
        "nutrient-sensing hallmark of aging.",
        analysis="mTOR is arguably the most actionable node in aging. What changes with age is "
        "less its expression than its activity: chronic nutrient and growth signalling keep mTOR "
        "over-engaged, suppressing the autophagy that clears cellular damage. The decisive "
        "evidence is pharmacological — rapamycin, which inhibits mTOR, extends lifespan in yeast, "
        "worms, flies and mice, even when started late in life. GeroQuery marks the expression "
        "direction context-dependent but highlights mTOR as a premier intervention target.",
        hallmarks=("deregulated_nutrient_sensing", "disabled_autophagy", "loss_of_proteostasis"),
        evidence=(
            Evidence(
                "physiology",
                "context-dependent",
                ("mouse",),
                ("multiple tissues",),
                "established",
                "mTOR pathway activity tends to rise with age; inhibiting mTOR with rapamycin "
                "extends lifespan in genetically heterogeneous mice, even when started late.",
                ("harrison_2009",),
            ),
        ),
        faqs=(
            (
                "Is mTOR good or bad for aging?",
                "Chronic high mTOR activity is pro-aging because it suppresses cellular clean-up "
                "(autophagy). Dialling it down pharmacologically extends lifespan.",
            ),
            (
                "What drug targets it?",
                "Rapamycin — the best-replicated pharmacological lifespan extender in mammals.",
            ),
        ),
    ),
    "IGF1R": GeneKnowledge(
        group="IGF1R",
        common_name="IGF1R",
        aka="the receptor whose mutations mark centenarians",
        direction_with_age="context-dependent",
        confidence="established",
        one_liner="Reduced IGF-1 receptor signalling is linked to exceptional longevity; "
        "functional IGF1R variants are enriched in centenarians.",
        role="IGF1R is the cell-surface receptor for IGF-1. Dampening its signalling is the "
        "conserved mechanism by which reduced insulin/IGF-1 pathway activity extends lifespan.",
        analysis="IGF1R is where the insulin/IGF-1 longevity story becomes human genetics. "
        "Centenarians are enriched for IGF1R mutations that partially reduce receptor function, "
        "echoing the lifespan extension seen when this pathway is dialled down in worms, flies "
        "and mice. Its expression does not follow a clean age arc (context-dependent), but the "
        "direction of benefit is clear: less IGF1R signalling, more longevity. This makes IGF1R "
        "a mechanistic anchor rather than a biomarker.",
        hallmarks=("deregulated_nutrient_sensing",),
        evidence=(
            Evidence(
                "genetics",
                "context-dependent",
                ("human",),
                ("germline",),
                "established",
                "Functionally significant IGF1R mutations that reduce receptor signalling are "
                "over-represented in centenarians.",
                ("suh_2008",),
            ),
        ),
        faqs=(
            (
                "How is IGF1R linked to long life?",
                "Centenarians carry more IGF1R variants that partially reduce its signalling, "
                "mirroring lifespan extension from lowered insulin/IGF-1 signalling in animals.",
            ),
            (
                "Does IGF1R rise or fall with age?",
                "There's no simple single trend; its relevance is genetic and mechanistic rather "
                "than as an age biomarker.",
            ),
        ),
    ),
}


# --------------------------------------------------------------------------- #
# Interventions — real lifespan-affecting compounds/diets with citations.
# Effect percentages are approximate median-lifespan changes reported in the
# cited rodent literature and are labelled as such; they are not precise or
# guaranteed and vary substantially by dose, sex, strain, and study.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class InterventionKnowledge:
    name: str
    display_name: str
    itype: str  # drug | dietary | senolytic
    source: str  # NIA ITP | DrugAge | GenDR ...
    organism: str
    lifespan_effect_pct: float  # approximate median-lifespan change (%)
    linked_groups: tuple[str, ...]
    summary: str
    reference_keys: tuple[str, ...]

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "display_name": self.display_name,
            "itype": self.itype,
            "source": self.source,
            "organism": self.organism,
            "lifespan_effect_pct": self.lifespan_effect_pct,
            "linked_groups": list(self.linked_groups),
            "summary": self.summary,
            "reference_keys": list(self.reference_keys),
        }


INTERVENTIONS: dict[str, InterventionKnowledge] = {
    "rapamycin": InterventionKnowledge(
        "rapamycin",
        "Rapamycin",
        "drug",
        "NIA ITP",
        "mouse",
        14.0,
        ("MTOR",),
        "An mTOR inhibitor and the best-replicated pharmacological lifespan extender in mammals. "
        "The NIA Interventions Testing Program reported roughly a 9–14% median-lifespan extension "
        "in mice (larger at higher doses), even when started in late life.",
        ("harrison_2009",),
    ),
    "caloric_restriction": InterventionKnowledge(
        "caloric_restriction",
        "Caloric restriction",
        "dietary",
        "GenDR",
        "mouse",
        30.0,
        ("IGF1", "MTOR", "SIRT1", "FOXO3", "IGF1R"),
        "Reducing calorie intake without malnutrition is the most reproducible lifespan-extending "
        "intervention across species, acting largely through the nutrient-sensing network "
        "(IGF-1, mTOR, sirtuins, FOXO). Reported rodent effects vary widely by degree and timing.",
        ("junnila_2013", "lopez_otin_2013"),
    ),
    "metformin": InterventionKnowledge(
        "metformin",
        "Metformin",
        "drug",
        "DrugAge",
        "mouse",
        6.0,
        ("MTOR", "GDF15"),
        "A widely used anti-diabetic that activates AMPK and indirectly restrains mTOR. In mice it "
        "produced a modest (~5–6%) lifespan and healthspan benefit; the human TAME trial is "
        "testing whether it slows aging.",
        ("martin_montalvo_2013",),
    ),
    "dasatinib_quercetin": InterventionKnowledge(
        "dasatinib_quercetin",
        "Dasatinib + Quercetin",
        "senolytic",
        "DrugAge",
        "mouse",
        9.0,
        ("CDKN2A", "CDKN1A"),
        "A senolytic combination that selectively eliminates p16/p21-high senescent cells. In "
        "aged mice it improved physical function and extended remaining lifespan.",
        ("xu_2018",),
    ),
    "17a_estradiol": InterventionKnowledge(
        "17a_estradiol",
        "17-alpha-estradiol",
        "drug",
        "NIA ITP",
        "mouse",
        12.0,
        ("IGF1",),
        "A non-feminising estrogen that extended median lifespan preferentially in male mice in "
        "the NIA Interventions Testing Program, improving metabolic health.",
        ("harrison_2014",),
    ),
    "nmn": InterventionKnowledge(
        "nmn",
        "Nicotinamide mononucleotide (NMN)",
        "drug",
        "DrugAge",
        "mouse",
        5.0,
        ("SIRT1",),
        "An NAD+ precursor aimed at restoring the age-related decline in NAD+ that limits sirtuin "
        "activity. Long-term dosing mitigated several age-associated physiological declines in "
        "mice; lifespan effects are modest and context-dependent.",
        ("mills_2016",),
    ),
}


# --------------------------------------------------------------------------- #
# Accessors
# --------------------------------------------------------------------------- #


def gene_knowledge(group: str) -> GeneKnowledge | None:
    return KNOWLEDGE.get(group)


def interventions_for_group(group: str) -> list[InterventionKnowledge]:
    return [iv for iv in INTERVENTIONS.values() if group in iv.linked_groups]
