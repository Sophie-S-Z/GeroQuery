"""M2 sources — source adapters (public surface)."""

from . import nhanes
from .base import Capabilities, FederatedStub, License, SourceAdapter
from .federated import FEDERATED_ADAPTERS
from .fetch import fetch_all, fetch_artifact, is_cached
from .geo import GeoDataSetsSource
from .gtex import GtexOpenSource
from .hagr import HagrSource
from .local_fixture import (
    CuratedKnowledgeSource,
    InterventionSource,
    LocalSignatureSource,
)
from .manifest import MANIFEST, MANIFEST_VERSION, RemoteArtifact, get_artifact
from .nhanes import NhanesClinicalSource


def all_adapters() -> list[SourceAdapter]:
    """Every registered adapter, cached + federated, for /v1/sources."""
    return [
        LocalSignatureSource(),
        GeoDataSetsSource(),
        HagrSource(),
        CuratedKnowledgeSource(),
        InterventionSource(),
        NhanesClinicalSource(),
        GtexOpenSource(),
        *FEDERATED_ADAPTERS,
    ]


__all__ = [
    "SourceAdapter",
    "Capabilities",
    "License",
    "FederatedStub",
    "LocalSignatureSource",
    "GeoDataSetsSource",
    "HagrSource",
    "CuratedKnowledgeSource",
    "InterventionSource",
    "NhanesClinicalSource",
    "GtexOpenSource",
    "nhanes",
    "FEDERATED_ADAPTERS",
    "all_adapters",
    "MANIFEST",
    "MANIFEST_VERSION",
    "RemoteArtifact",
    "get_artifact",
    "fetch_artifact",
    "fetch_all",
    "is_cached",
]
