"""Federated / controlled source adapters (interface-complete stubs).

They declare real capabilities and licences so the store, API, and licence test
can reason about them; the fetch bodies are where production wiring lands. The
controlled adapters are federate-only and non-redistributable — `assert_cacheable`
will refuse to persist them, which the licence test verifies.
"""

from __future__ import annotations

from .base import FederatedStub, License

# Open but large / long-tail: federate live, extracts may be cached.
RECOUNT3 = FederatedStub(
    name="recount3",
    _omics=("transcriptome",),
    _species=("human", "mouse"),
    _federated=True,
    _cacheable=True,
    _license=License("public", True, "recount3 / Bioconductor"),
    _note="Uniform RNA-seq. Wired via the R/Bioconductor ETL bridge, off the request path.",
)

CELLXGENE_CENSUS = FederatedStub(
    name="cellxgene-census",
    _omics=("single_cell",),
    _species=("human", "mouse"),
    _federated=True,
    _cacheable=True,
    _license=License("CC-BY", True, "CZ CELLxGENE Census"),
    _note="Single-cell; federated via the cellxgene_census API (stretch layer).",
)

# gtex-open is a real adapter now — see sources/gtex.py. It is registered through
# all_adapters() directly rather than as a stub here.

# Controlled-access: federate/link only, NEVER cached or re-hosted.
UK_BIOBANK = FederatedStub(
    name="uk-biobank",
    _omics=("clinical", "proteome"),
    _species=("human",),
    _federated=True,
    _cacheable=False,
    _license=License("controlled-access", redistributable=False, attribution="UK Biobank"),
    _note="Controlled access. Linked/federated only; redistribution is prohibited.",
)

GTEX_PROTECTED = FederatedStub(
    name="gtex-protected",
    _omics=("transcriptome",),
    _species=("human",),
    _federated=True,
    _cacheable=False,
    _license=License(
        "dbGaP-controlled", redistributable=False, attribution="GTEx protected / dbGaP"
    ),
    _note=(
        "Protected raw GTEx via dbGaP. Federate/link only. This is where donor age "
        "lives — the open API in sources/gtex.py cannot stratify by age."
    ),
)

FEDERATED_ADAPTERS = [RECOUNT3, CELLXGENE_CENSUS, UK_BIOBANK, GTEX_PROTECTED]
