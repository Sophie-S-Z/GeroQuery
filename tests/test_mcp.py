"""Tests for the MCP tool surface.

Two things are worth testing here and one is not. Worth testing: that the
payloads carry the interval and the honest verdict, and that every field a tool
claims to expose is actually populated. Not worth testing: the SDK's transport,
which is why ``tools.py`` imports no SDK and these tests need no ``mcp``
install.

The populated-field tests exist because the failure mode of a wrong key is
``None``, and a payload of nulls reads to a caller as "measured, no effect"
rather than as a bug. That is exactly how the dashboard came to report every
gene as absent from the knowledge base.
"""

from __future__ import annotations

import pytest

from geroquery.mcp import tools
from geroquery.mcp.server import SERVER_INSTRUCTIONS, TOOL_SCHEMAS

# --- the verdict rule --------------------------------------------------------


@pytest.mark.parametrize(
    ("ci_low", "ci_high", "expected"),
    [
        (0.55, 1.59, "increases"),
        (-1.20, -0.30, "decreases"),
        # The CDKN2A case: a positive point estimate whose interval spans zero.
        # Reporting this as a small increase is the single most damaging thing
        # this tool could do.
        (-0.20, 0.35, "no_evidence"),
        (-0.01, 2.00, "no_evidence"),
        (0.0, 1.0, "no_evidence"),
        (None, None, "no_evidence"),
    ],
)
def test_verdict_follows_the_interval_not_the_point_estimate(ci_low, ci_high, expected):
    assert tools._verdict({"ci_low": ci_low, "ci_high": ci_high}) == expected


def test_cdkn2a_returns_no_evidence_with_its_interval_intact():
    """The worked example. p16 is the most-cited gene in aging and does not
    replicate here; the tool has to say so, and show the interval that says it."""
    result = tools.gene_aging_signature("CDKN2A", species="human")
    assert result["verdict"] == "no_evidence"

    pooled = result["pooled"][0]
    assert pooled["interval"]["ci_low"] < 0 < pooled["interval"]["ci_high"]
    assert pooled["n_studies"] >= 10
    assert result["n_contrasts"] >= 10


def test_a_gene_that_does_replicate_is_not_reported_as_null():
    """The pair matters: without this, 'no_evidence' everywhere would look like
    a working tool rather than a broken estimator."""
    result = tools.gene_aging_signature("CDKN1A", species="human")
    assert result["verdict"] == "increases"
    assert result["pooled"][0]["interval"]["ci_low"] > 0


def test_an_unknown_gene_is_not_measured_rather_than_no_evidence():
    """'We never looked' and 'we looked and cannot tell' are different answers."""
    result = tools.gene_aging_signature("NOT_A_REAL_GENE_XYZZY")
    assert result["verdict"] == "not_measured"
    assert "error" in result


# --- populated, not silently None -------------------------------------------


def test_intervention_fields_are_populated_not_silently_none():
    """Pins the Intervention.to_dict() key names.

    The first version of this tool read ``effect_size`` and ``intervention_type``;
    the model exposes ``lifespan_effect_pct`` and ``itype``. Every value came
    back None and the payload looked like a drug with no measured effect.
    """
    result = tools.intervention_effect("rapamycin")
    assert result["found"] is True

    by_organism = result["by_organism"]
    assert by_organism, "rapamycin must return at least one organism"
    # Mammals lead, so a caller asking about a drug does not get a nematode.
    assert by_organism[0]["organism"] == "Mus musculus"

    for row in by_organism:
        assert row["organism"]
        assert row["intervention_type"] is not None, "itype mapped to the wrong key"
    assert any(
        r["median_lifespan_change_percent"] is not None for r in by_organism
    ), "lifespan_effect_pct mapped to the wrong key: every effect is None"

    mouse = by_organism[0]
    assert mouse["median_lifespan_change_percent"] == pytest.approx(13.0, abs=0.1)


def test_contrast_fields_are_populated_not_silently_none():
    """Same check for the per-study rows."""
    contrasts = tools.gene_aging_signature("CDKN1A", species="human")["contrasts"]
    assert contrasts
    for row in contrasts:
        assert row["study"], "study_id mapped to the wrong key"
        assert row["effect_size_hedges_g"] is not None
        assert row["species"]
    assert any(r["standard_error"] is not None for r in contrasts)


def test_unknown_intervention_is_reported_as_not_found_not_as_a_crash():
    result = tools.intervention_effect("definitely-not-a-compound")
    assert result["found"] is False
    assert "error" in result


# --- payload contract --------------------------------------------------------


def test_every_gene_answer_carries_caveats_and_a_data_version():
    """A tool response is often the only thing a caller reads, so the
    limitations travel with the number rather than living in a doc."""
    result = tools.gene_aging_signature("CDKN2A")
    assert result["caveats"]
    assert any("muscle" in c for c in result["caveats"])
    assert any("CI" in c or "interval" in c for c in result["caveats"])
    assert result["data_version"]["data_version"]


def test_geneset_aggregate_declines_to_invent_an_interval():
    """A mean of pooled effects has no honest CI. Returning one would be the
    most quietly wrong number this package could emit."""
    result = tools.geneset_aging_signature(["CDKN1A", "CDKN2A", "IGF1"], species="human")
    assert result["aggregate_interval"] is None
    assert "no confidence interval" in result["aggregate_note"]
    assert result["n_requested"] == 3
    for entry in result["per_gene"]:
        for pooled in entry["pooled"]:
            assert "interval" in pooled


def test_studies_tool_publishes_the_selection_rule():
    """The panel is rule-selected; a caller should be able to audit the corpus
    rather than trust it."""
    result = tools.list_studies()
    assert result["n_studies"] > 0
    assert "Subset Variable Type" in result["selection_rule"]


def test_provenance_reports_the_checksum_guarantee():
    result = tools.data_provenance()
    assert result["n_pinned_artifacts"] > 0
    assert "SHA-256" in result["guarantee"]
    assert result["manifest_version"]


# --- server binding ----------------------------------------------------------


def test_every_tool_has_a_schema_and_a_description():
    assert sorted(TOOL_SCHEMAS) == sorted(tools.TOOLS)
    for name, (fn, description) in tools.TOOLS.items():
        assert callable(fn)
        assert len(description) > 30, f"{name} needs a description an agent can act on"
        assert TOOL_SCHEMAS[name]["type"] == "object"


def test_server_instructions_tell_a_model_not_to_override_a_null():
    """The instructions are the only defence against a model answering from the
    literature consensus the measurement contradicts."""
    assert "no_evidence" in SERVER_INSTRUCTIONS
    assert "not_measured" in SERVER_INSTRUCTIONS
    assert "background knowledge" in SERVER_INSTRUCTIONS
