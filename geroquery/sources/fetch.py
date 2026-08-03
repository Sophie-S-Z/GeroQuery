"""Checksum-verifying fetch layer for pinned remote artifacts.

Contract, in one line: :func:`fetch_artifact` returns a path to bytes whose SHA-256 equals
the value pinned in :mod:`geroquery.sources.manifest`, or it raises. There is no
path through this module that hands back unverified data.

Behaviour:
  * A cache hit is re-verified, not trusted. A truncated or edited cache entry is
    discarded and re-fetched rather than silently used.
  * Downloads land on a temp file and are only moved into place after the digest
    matches, so an interrupted download can never leave a valid-looking cache
    entry behind.
  * With ``GEROQUERY_ALLOW_NETWORK`` off (the default), a cache miss raises
    :class:`NetworkDisabledError` instead of reaching out. Tests and CI are
    hermetic by construction, not by convention.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Iterable
from pathlib import Path

from .. import config
from ..exceptions import ChecksumMismatchError, NetworkDisabledError, SourceError
from .manifest import MANIFEST, RemoteArtifact, get_artifact

# Read in chunks so a multi-GB artifact never has to fit in memory.
_CHUNK = 1 << 20  # 1 MiB


def cache_dir() -> Path:
    """Directory holding verified downloads."""
    return Path(config.CACHE_HOME) / "downloads"


def cache_path(artifact: RemoteArtifact, directory: Path | None = None) -> Path:
    """Where this artifact's verified bytes live on disk.

    Prefixed with the manifest key because URL basenames are not unique across
    upstreams: HAGR serves both DrugAge and GenDR from a path ending
    ``/dataset.zip``. Keying on the basename alone made the two share one cache
    entry, so whichever downloaded second overwrote the first. The checksum check
    turned that into a loud failure rather than a silent wrong answer, which is
    the point of it — but the collision is fixed here rather than relied upon.
    """
    return (directory or cache_dir()) / f"{artifact.key}__{artifact.filename}"


def sha256_file(path: Path) -> str:
    """Streaming SHA-256 of a file on disk."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(_CHUNK):
            h.update(chunk)
    return h.hexdigest()


def verify(path: Path, artifact: RemoteArtifact) -> None:
    """Raise unless ``path`` matches the artifact's pinned size and digest."""
    size = path.stat().st_size
    if size != artifact.n_bytes:
        raise ChecksumMismatchError(
            f"{artifact.key}: expected {artifact.n_bytes} bytes, got {size}.",
            detail={
                "artifact": artifact.key,
                "expected_bytes": artifact.n_bytes,
                "got_bytes": size,
            },
        )
    digest = sha256_file(path)
    if digest != artifact.sha256:
        raise ChecksumMismatchError(
            f"{artifact.key}: SHA-256 mismatch. Upstream file changed, or the download is corrupt.",
            detail={
                "artifact": artifact.key,
                "url": artifact.url,
                "expected_sha256": artifact.sha256,
                "got_sha256": digest,
            },
        )


def is_cached(artifact: RemoteArtifact, directory: Path | None = None) -> bool:
    """True only if a cache entry exists *and* still verifies."""
    path = cache_path(artifact, directory)
    if not path.exists():
        return False
    try:
        verify(path, artifact)
    except ChecksumMismatchError:
        return False
    return True


def _download(artifact: RemoteArtifact, dest: Path, timeout: float) -> None:
    """Stream the artifact to ``dest``. Caller verifies; this only transports."""
    import httpx

    try:
        with httpx.stream(
            "GET",
            artifact.url,
            timeout=timeout,
            follow_redirects=True,
            headers={"User-Agent": "GeroQuery/1.0 (+https://github.com/praneel/geroquery)"},
        ) as resp:
            resp.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in resp.iter_bytes(_CHUNK):
                    fh.write(chunk)
    except httpx.HTTPError as exc:
        dest.unlink(missing_ok=True)
        raise SourceError(
            f"Failed to download {artifact.key} from {artifact.url}: {exc}",
            detail={"artifact": artifact.key, "url": artifact.url},
        ) from exc


def fetch_artifact(
    artifact: RemoteArtifact | str,
    *,
    directory: Path | None = None,
    allow_network: bool | None = None,
    force: bool = False,
    timeout: float | None = None,
) -> Path:
    """Return a path to verified bytes for ``artifact``, downloading if needed.

    Args:
        artifact: a :class:`RemoteArtifact` or its manifest key.
        directory: cache directory override (defaults to ``CACHE_HOME/downloads``).
        allow_network: override ``GEROQUERY_ALLOW_NETWORK`` for this call.
        force: re-download even if a valid cache entry exists.
        timeout: per-request timeout in seconds.

    Raises:
        NetworkDisabledError: cache miss while network access is disabled.
        ChecksumMismatchError: downloaded bytes do not match the pinned digest.
        SourceError: transport failure.
    """
    if isinstance(artifact, str):
        artifact = get_artifact(artifact)
    allow_network = config.ALLOW_NETWORK if allow_network is None else allow_network
    timeout = config.HTTP_TIMEOUT if timeout is None else timeout

    target = cache_path(artifact, directory)
    if not force and is_cached(artifact, directory):
        return target

    if not allow_network:
        raise NetworkDisabledError(
            f"{artifact.key} is not in the local cache and network access is disabled. "
            f"Run `make data` (or set GEROQUERY_ALLOW_NETWORK=1) to fetch it.",
            detail={"artifact": artifact.key, "url": artifact.url, "cache_path": str(target)},
        )

    target.parent.mkdir(parents=True, exist_ok=True)
    # Download to a sibling temp file; only a verified file gets the real name.
    tmp = target.with_suffix(target.suffix + ".part")
    _download(artifact, tmp, timeout)
    try:
        verify(tmp, artifact)
    except ChecksumMismatchError:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, target)
    return target


def fetch_all(
    artifacts: Iterable[RemoteArtifact | str] | None = None,
    *,
    directory: Path | None = None,
    allow_network: bool | None = None,
    force: bool = False,
) -> dict[str, Path]:
    """Fetch and verify several artifacts, keyed by manifest key."""
    items = list(artifacts) if artifacts is not None else list(MANIFEST.values())
    out: dict[str, Path] = {}
    for item in items:
        art = get_artifact(item) if isinstance(item, str) else item
        out[art.key] = fetch_artifact(
            art, directory=directory, allow_network=allow_network, force=force
        )
    return out
