"""M8 mcp — Model Context Protocol surface over the existing API layer.

``tools`` holds the payload shaping and imports no SDK, so it is testable
without one. ``server`` is the transport binding.
"""

from .tools import TOOLS

__all__ = ["TOOLS"]
