"""M7 api — FastAPI service layer (public surface)."""

from .app import create_app
from .service import GeroService, get_service

__all__ = ["create_app", "GeroService", "get_service"]
