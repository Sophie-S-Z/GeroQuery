"""Runtime configuration and canonical filesystem locations.

Kept deliberately tiny: paths to bundled reference data plus a couple of
network toggles. Everything is overridable by environment variable so the same
code runs identically in dev, CI, and a container.
"""

from __future__ import annotations

import os
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent

# Bundled, redistributable reference data shipped inside the wheel.
IDMAP_DATA = PACKAGE_ROOT / "idmap" / "data"
SOURCES_DATA = PACKAGE_ROOT / "sources" / "data"


def _flag(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


# Live network access to upstreams (mygene.info, GEO, GTEx, ...). Off by
# default so tests and CI are hermetic and never hammer external APIs.
ALLOW_NETWORK: bool = _flag("GEROQUERY_ALLOW_NETWORK", False)

# Where DuckDB/Parquet and the SQLite metadata DB live for a hosted deployment.
# Defaults to a per-process temp build assembled from bundled fixtures.
DATA_HOME: Path = Path(os.environ.get("GEROQUERY_DATA_HOME", PACKAGE_ROOT / "_build"))

# On-disk cache directory for federated responses.
CACHE_HOME: Path = Path(os.environ.get("GEROQUERY_CACHE_HOME", PACKAGE_ROOT / "_cache"))

# Default HTTP timeout (seconds) for any federated request.
HTTP_TIMEOUT: float = float(os.environ.get("GEROQUERY_HTTP_TIMEOUT", "20"))
