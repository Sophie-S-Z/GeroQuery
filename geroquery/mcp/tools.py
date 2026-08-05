"""Tool implementations for the MCP server — plain functions, no SDK import.

Split from :mod:`geroquery.mcp.server` for the same reason ``api/service.py`` is
split from ``api/app.py``: the protocol binding is trivial and the payload
shaping is not, so the part worth testing should not require the transport to be
installed.

**The design rule for every tool here: the interval and the study count go in
the primary response, never in an optional field.**

That is not house style, it is the entire point. A language model asked "does
CDKN2A change with age" answers from the literature — which is how CDKN2A became
the most-cited aging gene despite pooling to *g* = +0.07, 95% CI
[-0.20, +0.35] over fourteen contrasts here. The value of putting GeroQuery
behind a tool call is that the model gets a measurement instead of a consensus.
A tool that returns a point estimate and hides the interval throws that away and
launders a null into a number.

So each payload carries, unconditionally:

* ``verdict`` — one of ``increases``, ``decreases``, ``no_evidence``,
  ``not_measured``. ``no_evidence`` is a real answer, not an error.
* ``interval`` — the CI, always.
* ``n_studies`` / ``n_contrasts`` — how much evidence is behind it.
* ``caveats`` — the limitations that apply to *this* answer.
"""

from __future__ import annotations

from typing import Any

from ..api.service import GeroService, get_service
from ..exceptions import GeroQueryError

# Applies to every signature answer. Kept here rather than in prose docs because
# a tool response is often the only thing a caller ever reads.
SIGNATURE_CAVEATS: tuple[str, ...] = (
    "Microarray-era GEO DataSets only; no RNA-seq.",
    "Groups of 3-15 subjects; median 95% CI width 0.84 — detects large consistent "
    "effects, cannot rule out moderate ones.",
    "Half the human contrasts are skeletal muscle; one blood contrast (n=9).",
    "32 contrasts come from 27 GEO Series, so some share subjects.",
    "Cross-study processing heterogeneity is not corrected.",
)

# A pooled estimate whose interval contains zero is reported as no evidence,
# regardless of the point estimate's sign.
NULL_VERDICT = "no_evidence"


def _verdict(meta: dict) -> str:
    """Turn a pooled estimate into a word, honestly.

    The rule is the interval, not the point estimate and not the p-value. An
    effect of +0.07 with an interval spanning [-0.20, +0.35] is not a small
    increase; it is an absence of evidence, and saying so is the only reason
    this tool is worth calling.
    """
    low, high = meta.get("ci_low"), meta.get("ci_high")
    if low is None or high is None:
        return NULL_VERDICT
    if low > 0:
        return "increases"
    if high < 0:
        return "decreases"
    return NULL_VERDICT


def _summarize(meta: dict) -> dict:
    return {
        "species": meta.get("species"),
        "omic_layer": meta.get("omic_layer"),
        "verdict": _verdict(meta),
        "pooled_effect_hedges_g": meta.get("pooled_effect"),
        "interval": {
            "ci_low": meta.get("ci_low"),
            "ci_high": meta.get("ci_high"),
            "level": 0.95,
        },
        "n_studies": meta.get("n_studies"),
        "heterogeneity_i2": meta.get("heterogeneity_i2"),
        "p_value": meta.get("p_value"),
    }


def gene_aging_signature(
    gene: str,
    species: str | None = None,
    tissue: str | None = None,
    service: GeroService | None = None,
) -> dict[str, Any]:
    """Measured age association for one gene, re-derived from raw GEO data.

    Args:
        gene: symbol, alias, Ensembl or Entrez id. ``p16``, ``CDKN2A`` and
            ``ENSG00000147889`` all resolve.
        species: ``human`` or ``mouse``. Omit to get both.
        tissue: UBERON term or common name, e.g. ``skeletal muscle``.
    """
    svc = service or get_service()
    try:
        raw = svc.gene_signature(gene, species=species, tissue=tissue)
    except GeroQueryError as exc:
        return {
            "query": gene,
            "verdict": "not_measured",
            "error": {"code": exc.code, "message": exc.message},
            "caveats": list(SIGNATURE_CAVEATS),
        }

    metas = [_summarize(m) for m in raw.get("meta_signatures", [])]
    contrasts = raw.get("signatures", [])
    return {
        "query": gene,
        "gene": raw.get("gene"),
        # If nothing was measured, say that rather than returning an empty list
        # a caller might read as "no effect".
        "verdict": metas[0]["verdict"] if metas else "not_measured",
        "pooled": metas,
        "n_contrasts": len(contrasts),
        "contrasts": [
            {
                "study": c.get("study_id"),
                "tissue": c.get("tissue"),
                "species": c.get("species"),
                "sex": c.get("sex"),
                "age_range": c.get("age_range"),
                "effect_size_hedges_g": c.get("effect_size"),
                "standard_error": c.get("standard_error"),
                "direction": c.get("direction"),
                "q_value": c.get("q_value"),
            }
            for c in contrasts
        ],
        "caveats": list(SIGNATURE_CAVEATS),
        "data_version": svc.version(),
    }


def geneset_aging_signature(
    genes: list[str],
    species: str | None = None,
    service: GeroService | None = None,
) -> dict[str, Any]:
    """Pool a gene set or pathway. Unresolved genes are named, with the reason."""
    svc = service or get_service()
    raw = svc.geneset_signature(genes, species=species)
    per_gene = [
        {
            "gene": entry["gene"].get("symbol"),
            "pooled": [_summarize(m) for m in entry.get("meta_signatures", [])],
        }
        for entry in raw.get("per_gene", [])
    ]
    measured = [g for g in per_gene if g["pooled"]]
    return {
        "n_requested": len(genes),
        "n_resolved": len(raw.get("resolved", [])),
        "n_with_measurements": len(measured),
        "aggregate_pooled_effect": raw.get("aggregate_pooled_effect"),
        # An aggregate over genes is a summary of summaries and has no honest
        # interval; saying so beats inventing one.
        "aggregate_interval": None,
        "aggregate_note": (
            "The aggregate is the mean of per-gene pooled effects and carries no "
            "confidence interval. Use the per-gene intervals to judge the set."
        ),
        "per_gene": per_gene,
        "unresolved": raw.get("unresolved", []),
        "unresolved_detail": raw.get("unresolved_detail", {}),
        "caveats": list(SIGNATURE_CAVEATS),
        "data_version": svc.version(),
    }


def intervention_effect(name: str, service: GeroService | None = None) -> dict[str, Any]:
    """Lifespan effect of a compound or intervention, per organism.

    Every organism is returned, mammals first. DrugAge holds one record per
    organism, so answering with a single number means answering a question about
    a drug with a result about nematodes.
    """
    svc = service or get_service()
    try:
        raw = svc.intervention(name)
    except GeroQueryError as exc:
        return {
            "query": name,
            "found": False,
            "error": {"code": exc.code, "message": exc.message},
        }
    return {
        "query": name,
        "found": True,
        "name": raw.get("name"),
        "n_records": raw.get("n"),
        # Field names are taken from Intervention.to_dict() and pinned by
        # test_intervention_fields_are_populated_not_silently_none. Guessing
        # them produces a payload of nulls that reads as "no effect measured"
        # rather than as a mapping bug — the same failure mode that made the
        # dashboard report every gene as absent from the knowledge base.
        "by_organism": [
            {
                "organism": iv.get("organism"),
                "median_lifespan_change_percent": iv.get("lifespan_effect_pct"),
                "intervention_type": iv.get("itype"),
                "source": iv.get("source"),
                "url": iv.get("url"),
            }
            for iv in raw.get("interventions", [])
        ],
        "caveats": [
            "DrugAge records model-organism lifespan experiments. A mammalian "
            "result is not a human result.",
            "DrugAge carries no gene links, so no gene associations are asserted.",
        ],
        "data_version": svc.version(),
    }


def list_studies(service: GeroService | None = None) -> dict[str, Any]:
    """Every contrast in the panel, with full provenance.

    The panel is rule-selected, not hand-picked: every GEO DataSet declaring an
    ``age`` subset variable that survives the contrast rules. That is what keeps
    selection bias out, and it is why this tool exists — a caller can audit the
    corpus rather than trust it.
    """
    svc = service or get_service()
    studies = svc.studies()
    return {
        "n_studies": len(studies),
        "selection_rule": (
            "GEO DataSets query: age[Subset Variable Type] AND gds[Filter] AND "
            "(Homo sapiens[Organism] OR Mus musculus[Organism]); then an adult "
            "young-vs-old split restricted to the control arm of every other "
            "subset variable, with at least three samples per group."
        ),
        "studies": studies,
        "data_version": svc.version(),
    }


def data_provenance(service: GeroService | None = None) -> dict[str, Any]:
    """Which upstreams back these answers, under which licences, at which version."""
    svc = service or get_service()
    from ..sources.manifest import MANIFEST, MANIFEST_VERSION, VERIFIED_ON

    return {
        "manifest_version": MANIFEST_VERSION,
        "checksums_verified_on": VERIFIED_ON,
        "n_pinned_artifacts": len(MANIFEST),
        "guarantee": (
            "Every pinned artifact is SHA-256 verified on fetch, and cache hits are "
            "re-verified rather than trusted. Nothing enters the store unverified."
        ),
        "sources": svc.sources(),
        "datasets": svc.datasets(),
        "data_version": svc.version(),
    }


# Tool name -> (callable, one-line description). The server binds these; the
# tests call them directly.
TOOLS: dict[str, tuple[Any, str]] = {
    "gene_aging_signature": (
        gene_aging_signature,
        "Measured age association for one gene, re-derived from raw GEO data, with "
        "a confidence interval. Returns 'no_evidence' when the interval crosses zero.",
    ),
    "geneset_aging_signature": (
        geneset_aging_signature,
        "Pool a gene set or pathway; per-gene intervals plus an aggregate.",
    ),
    "intervention_effect": (
        intervention_effect,
        "Lifespan effect of a compound, per organism, mammals first.",
    ),
    "list_studies": (
        list_studies,
        "Every contrast in the panel with full provenance and the selection rule.",
    ),
    "data_provenance": (
        data_provenance,
        "Pinned upstreams, licences, checksum guarantee, and the data version.",
    ),
}
