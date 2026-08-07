"""The MCP binding, exercised over an actual client-server session.

``tests/test_mcp.py`` calls the tool functions directly, which is the right place
for payload logic and is deliberately SDK-free. It cannot see anything that goes
wrong *between* a client and those functions: a schema that does not match the
signature it documents, a tool registered under a name nothing dispatches, a
payload that is not JSON-serializable, an error that reaches the client as a
crash rather than as a tool error. Those are the failures an agent would actually
hit, and until now nothing ran them.

The session is in-memory, so this is a real protocol round-trip with no
subprocess and no stdio.
"""

from __future__ import annotations

import json

import pytest

anyio = pytest.importorskip("anyio")
memory = pytest.importorskip("mcp.shared.memory")

from geroquery.mcp.server import SERVER_NAME, TOOL_SCHEMAS, build_server  # noqa: E402
from geroquery.mcp.tools import TOOLS  # noqa: E402


def _session():
    return memory.create_connected_server_and_client_session(build_server(), raise_exceptions=True)


def _run(coroutine_factory):
    return anyio.run(coroutine_factory)


def test_the_server_initializes_and_hands_over_its_instructions():
    async def go():
        async with _session() as client:
            return await client.initialize()

    init = _run(go)
    assert init.serverInfo.name == SERVER_NAME
    # The instructions are the only thing telling a model not to override a
    # measured null with what it remembers from the literature. If they do not
    # survive the handshake, the tool is a lookup with no guardrail.
    assert "no_evidence" in (init.instructions or "")
    assert "background knowledge" in (init.instructions or "")


def test_every_registered_tool_is_advertised_with_a_schema_that_matches_it():
    async def go():
        async with _session() as client:
            await client.initialize()
            return await client.list_tools()

    listed = _run(go)
    names = {tool.name for tool in listed.tools}
    assert names == set(TOOLS)
    assert names == set(TOOL_SCHEMAS)

    import inspect

    for tool in listed.tools:
        assert tool.description
        properties = set(tool.inputSchema.get("properties", {}))
        parameters = set(inspect.signature(TOOLS[tool.name][0]).parameters) - {"service"}
        # A property the function does not accept becomes a TypeError at call
        # time; a parameter the schema omits is one an agent can never reach.
        assert properties <= parameters, f"{tool.name}: schema advertises unknown args"
        required = set(tool.inputSchema.get("required", []))
        assert required <= properties


def test_a_gene_query_round_trips_with_its_interval_intact():
    async def go():
        async with _session() as client:
            await client.initialize()
            return await client.call_tool("gene_aging_signature", {"gene": "CDKN2A"})

    result = _run(go)
    assert not result.isError
    payload = json.loads(result.content[0].text)

    # The product's whole claim, checked at the far end of the wire rather than
    # in the function that builds it.
    assert payload["verdict"] in {"increases", "decreases", "no_evidence", "not_measured"}
    human = next(p for p in payload["pooled"] if p["species"] == "human")
    assert human["interval"]["ci_low"] < human["interval"]["ci_high"]
    assert human["n_studies"] >= 3
    # And the prediction interval, which an agent needs precisely because it
    # will otherwise read the confidence interval as a claim about the next study.
    assert human["prediction_interval"]["pi_low"] < human["prediction_interval"]["pi_high"]


def test_a_gene_that_was_never_measured_says_so_rather_than_erroring():
    async def go():
        async with _session() as client:
            await client.initialize()
            return await client.call_tool("gene_aging_signature", {"gene": "NOTAGENE"})

    result = _run(go)
    assert not result.isError
    assert json.loads(result.content[0].text)["verdict"] == "not_measured"


def test_a_missing_required_argument_comes_back_as_a_tool_error():
    """Not as a crashed session. An agent can recover from the first."""

    async def go():
        async with _session() as client:
            await client.initialize()
            return await client.call_tool("gene_aging_signature", {})

    result = _run(go)
    assert result.isError
    assert "gene" in result.content[0].text


def test_every_tool_returns_json_over_the_wire():
    """Serializability is a transport property, not a payload one.

    ``tools.py`` returns dicts, and a dict holding a numpy float or a date is
    perfectly valid Python that fails only when something tries to encode it.
    """
    arguments = {
        "gene_aging_signature": {"gene": "CDKN1A"},
        "geneset_aging_signature": {"genes": ["CDKN1A", "IGF1"]},
        "intervention_effect": {"name": "rapamycin"},
        "list_studies": {},
        "data_provenance": {},
    }
    assert set(arguments) == set(TOOLS)

    async def go():
        async with _session() as client:
            await client.initialize()
            out = {}
            for name, args in arguments.items():
                result = await client.call_tool(name, args)
                out[name] = (result.isError, result.content[0].text)
            return out

    for name, (is_error, text) in _run(go).items():
        assert not is_error, f"{name} errored: {text[:200]}"
        assert isinstance(json.loads(text), dict), name
