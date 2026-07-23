"""M2 sources — source adapters (public surface)."""

from .base import Capabilities, FederatedStub, License, SourceAdapter
from .federated import FEDERATED_ADAPTERS
from .local_fixture import (
    CuratedKnowledgeSource,
    InterventionSource,
    LocalSignatureSource,
)


def all_adapters() -> list[SourceAdapter]:
    """Every registered adapter, cached + federated, for /v1/sources."""
    return [
        LocalSignatureSource(),
        CuratedKnowledgeSource(),
        InterventionSource(),
        *FEDERATED_ADAPTERS,
    ]


__all__ = [
    "SourceAdapter",
    "Capabilities",
    "License",
    "FederatedStub",
    "LocalSignatureSource",
    "CuratedKnowledgeSource",
    "InterventionSource",
    "FEDERATED_ADAPTERS",
    "all_adapters",
]
