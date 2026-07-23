"""Central exception hierarchy.

One place so the API layer (M7) can map domain errors to HTTP status codes and
the consistent ``{error: {code, message, detail}}`` envelope without importing
every module's internals.
"""

from __future__ import annotations

from typing import Any


class GeroQueryError(Exception):
    """Base for every domain error. Carries a stable machine-readable code."""

    code = "geroquery_error"
    http_status = 400

    def __init__(self, message: str, detail: Any | None = None):
        super().__init__(message)
        self.message = message
        self.detail = detail


class GeneNotFoundError(GeroQueryError):
    code = "gene_not_found"
    http_status = 404


class TissueNotFoundError(GeroQueryError):
    code = "tissue_not_found"
    http_status = 404


class UnknownSpeciesError(GeroQueryError):
    code = "unknown_species"
    http_status = 422


class ClockNotFoundError(GeroQueryError):
    code = "clock_not_found"
    http_status = 404


class ClockInputError(GeroQueryError):
    """Uploaded matrix is missing features a clock requires, etc."""

    code = "clock_input_error"
    http_status = 422


class ResilienceInputError(GeroQueryError):
    code = "resilience_input_error"
    http_status = 422


class SourceError(GeroQueryError):
    """An upstream adapter failed (network, parse, rate limit)."""

    code = "source_error"
    http_status = 502


class LicenseViolationError(GeroQueryError):
    """Attempted to cache/redistribute a federate-only source."""

    code = "license_violation"
    http_status = 403
