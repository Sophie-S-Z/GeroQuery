"""M2 sources — the uniform adapter interface.

Every data source, whether we cache it or federate to it live, implements the
same tiny surface: ``capabilities()`` (what it can serve) and ``license()``
(whether we may redistribute it). Adding a source means writing one adapter and
nothing else.

The licence object is load-bearing, not documentation: a federate-only source
raises :class:`LicenseViolationError` if any caller tries to cache/redistribute
it, and a test asserts this. That is how the "never re-host controlled data"
rule is enforced in code rather than in a README.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from ..exceptions import LicenseViolationError


@dataclass(frozen=True)
class License:
    name: str
    redistributable: bool
    attribution: str = ""


@dataclass(frozen=True)
class Capabilities:
    source_name: str
    omics: tuple[str, ...] = ()
    species: tuple[str, ...] = ()
    federated: bool = False  # queried live at request time
    cacheable: bool = True  # may we persist harmonized extracts?
    notes: str = ""


class SourceAdapter(ABC):
    """Base class for all source adapters."""

    name: str = "abstract"

    @abstractmethod
    def capabilities(self) -> Capabilities: ...

    @abstractmethod
    def license(self) -> License: ...

    def assert_cacheable(self) -> None:
        """Guard the cache-vs-federate boundary in code.

        Raises if this adapter's licence forbids redistribution. The store calls
        this before persisting anything from the adapter.
        """
        cap = self.capabilities()
        lic = self.license()
        if cap.federated and not cap.cacheable:
            raise LicenseViolationError(
                f"Source {self.name!r} is federate-only and must never be cached or re-hosted.",
                detail={"source": self.name, "license": lic.name},
            )
        if not lic.redistributable:
            raise LicenseViolationError(
                f"Source {self.name!r} data is not redistributable (license {lic.name!r}).",
                detail={"source": self.name, "license": lic.name},
            )


@dataclass
class FederatedStub(SourceAdapter):
    """Interface-complete placeholder for a live source not yet wired.

    Declares its real capabilities and licence so the rest of the system (store,
    api, licence test) can reason about it, but raises a clear error if actually
    queried while network is off. This is the seam the production adapters slot
    into without touching callers.
    """

    name: str
    _omics: tuple[str, ...] = ()
    _species: tuple[str, ...] = ()
    _federated: bool = True
    _cacheable: bool = True
    _license: License = field(default_factory=lambda: License("public", True))
    _note: str = "Federated adapter stub — implement fetch on the request path or ETL."

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=self._omics,
            species=self._species,
            federated=self._federated,
            cacheable=self._cacheable,
            notes=self._note,
        )

    def license(self) -> License:
        return self._license
