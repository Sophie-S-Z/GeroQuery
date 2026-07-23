"""M-knowledge — GeroQuery's curated, real-evidence aging knowledge base.

This module is the scientific source of truth. It contains only real, curated,
and cited biology — no simulated effect sizes, invented accessions, or synthetic
statistics. Higher layers (the service, the API, the dashboard) compose it with
gene-ID resolution to answer questions about how a gene relates to aging.
"""

from __future__ import annotations

from .aging_knowledge import (
    HALLMARKS,
    INTERVENTIONS,
    KNOWLEDGE,
    Evidence,
    GeneKnowledge,
    InterventionKnowledge,
    gene_knowledge,
    interventions_for_group,
)
from .references import REFERENCES, Reference, ref

__all__ = [
    "HALLMARKS",
    "INTERVENTIONS",
    "KNOWLEDGE",
    "Evidence",
    "GeneKnowledge",
    "InterventionKnowledge",
    "gene_knowledge",
    "interventions_for_group",
    "REFERENCES",
    "Reference",
    "ref",
]
