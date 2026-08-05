"""MCP server — GeroQuery as a tool an agent can call.

Run:  ``python -m geroquery.mcp.server``            (stdio, for a local client)
Install: ``pip install -e ".[mcp]"``

Why this exists. MCPmed (*Brief Bioinform* 2026;27(1):bbag076) argues that
bioinformatics web services should expose a machine-actionable layer rather than
only a human-facing one, and demonstrates it on GEO, STRING and the UCSC Cell
Browser; Open Targets shipped an official server on the same reasoning. The
pattern it recommends is a typed API plus a thin MCP endpoint over it, which is
the shape this package already had — ``api/service.py`` is 400 lines of typed
domain methods, so this module is a binding, not a rewrite.

The substantive argument is narrower than "agents are useful". GeroQuery's
differentiator is precisely the thing a language model cannot do unaided: ask a
model whether CDKN2A changes with age and it answers from the literature, which
is how CDKN2A became the most-cited gene in aging despite pooling to
*g* = +0.07, 95% CI [-0.20, +0.35] across fourteen contrasts here. The tool
returns the measurement and the interval. That is the whole product.

Thin by construction: routing and schema only. Every payload is shaped in
:mod:`geroquery.mcp.tools`, which imports no SDK and is therefore testable
without one.
"""

from __future__ import annotations

from typing import Any

from .tools import TOOLS

SERVER_NAME = "geroquery"
SERVER_INSTRUCTIONS = """\
GeroQuery re-derives aging evidence from raw published data rather than serving \
a curated verdict. Every gene answer carries a confidence interval and the \
number of contrasts behind it.

Read `verdict` before `pooled_effect_hedges_g`. A verdict of `no_evidence` means \
the confidence interval crosses zero — it is a real answer and should be \
reported as "the data cannot support a claim", not as a small effect or as a \
missing value. `not_measured` is different and means the gene is not in the \
panel at all.

Do not substitute background knowledge for a `no_evidence` verdict. Several \
famous aging genes do not replicate here, and reporting the literature consensus \
over the measurement defeats the purpose of the query.
"""

# JSON Schemas, one per tool. Written out rather than generated so the
# descriptions an agent actually reads are visible and editable here.
TOOL_SCHEMAS: dict[str, dict[str, Any]] = {
    "gene_aging_signature": {
        "type": "object",
        "properties": {
            "gene": {
                "type": "string",
                "description": "Symbol, alias, Ensembl or Entrez id. 'p16', 'CDKN2A' "
                "and 'ENSG00000147889' all resolve.",
            },
            "species": {
                "type": "string",
                "enum": ["human", "mouse"],
                "description": "Omit to get both species pooled separately.",
            },
            "tissue": {
                "type": "string",
                "description": "UBERON term or common name, e.g. 'skeletal muscle'.",
            },
        },
        "required": ["gene"],
    },
    "geneset_aging_signature": {
        "type": "object",
        "properties": {
            "genes": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Gene symbols or identifiers.",
            },
            "species": {"type": "string", "enum": ["human", "mouse"]},
        },
        "required": ["genes"],
    },
    "intervention_effect": {
        "type": "object",
        "properties": {
            "name": {
                "type": "string",
                "description": "Compound or intervention name, e.g. 'rapamycin'.",
            }
        },
        "required": ["name"],
    },
    "list_studies": {"type": "object", "properties": {}},
    "data_provenance": {"type": "object", "properties": {}},
}


def build_server() -> Any:
    """Construct the MCP server, binding every tool in :data:`TOOLS`.

    Imported lazily so that ``geroquery.mcp.tools`` — the part with logic worth
    testing — stays importable without the SDK installed.
    """
    try:
        from mcp.server import Server
        from mcp.types import TextContent, Tool
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "The MCP server needs the 'mcp' package. Install it with:\n"
            '    pip install -e ".[mcp]"'
        ) from exc

    import json

    server = Server(SERVER_NAME, instructions=SERVER_INSTRUCTIONS)

    @server.list_tools()
    async def _list_tools() -> list[Tool]:
        return [
            Tool(
                name=name,
                description=description,
                inputSchema=TOOL_SCHEMAS[name],
            )
            for name, (_fn, description) in TOOLS.items()
        ]

    @server.call_tool()
    async def _call_tool(name: str, arguments: dict[str, Any]) -> list[TextContent]:
        entry = TOOLS.get(name)
        if entry is None:
            raise ValueError(f"Unknown tool {name!r}. Known: {sorted(TOOLS)}")
        fn, _description = entry
        payload = fn(**(arguments or {}))
        return [TextContent(type="text", text=json.dumps(payload, default=str))]

    return server


def main() -> int:  # pragma: no cover - transport loop
    import anyio
    from mcp.server.stdio import stdio_server

    server = build_server()

    async def _run() -> None:
        async with stdio_server() as (read, write):
            await server.run(read, write, server.create_initialization_options())

    anyio.run(_run)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
