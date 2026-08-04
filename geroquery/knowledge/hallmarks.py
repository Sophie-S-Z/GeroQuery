"""The López-Otín hallmarks of aging — a published vocabulary, not data.

Used to group findings for a reader. The taxonomy is from López-Otín et al.,
*Cell* 2013 ("The hallmarks of aging") and its 2023 update, which added
disabled autophagy, chronic inflammation, and dysbiosis to the original nine.
See :data:`geroquery.knowledge.REFERENCES` for the citations.
"""

from __future__ import annotations

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
