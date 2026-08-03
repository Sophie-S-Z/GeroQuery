"""CLI to download and SHA-256 verify every pinned upstream artifact.

Run: ``python -m geroquery.etl.fetch_artifacts [KEY ...]``

Lives in the ETL package rather than next to the fetch layer because
``geroquery.sources`` re-exports from ``geroquery.sources.fetch``, which makes
``python -m geroquery.sources.fetch`` emit a "found in sys.modules" RuntimeWarning.
"""

from __future__ import annotations

import argparse

from ..exceptions import SourceError
from ..sources.fetch import fetch_all
from ..sources.manifest import MANIFEST, MANIFEST_VERSION, VERIFIED_ON


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch and verify pinned data artifacts.")
    parser.add_argument("keys", nargs="*", help="manifest keys (default: all)")
    parser.add_argument("--force", action="store_true", help="re-download even if cached")
    parser.add_argument(
        "--offline", action="store_true", help="fail instead of downloading on a cache miss"
    )
    args = parser.parse_args(argv)

    keys = args.keys or sorted(MANIFEST)
    print(f"manifest {MANIFEST_VERSION} (checksums verified {VERIFIED_ON})")
    try:
        paths = fetch_all(
            keys,
            allow_network=not args.offline,
            force=args.force,
        )
    except SourceError as exc:
        # A failed checksum or a disabled network is an expected operator-facing
        # condition, not a crash. Print the actionable message, skip the traceback.
        print(f"error: {exc.message}")
        if exc.detail:
            print(f"detail: {exc.detail}")
        return 1

    for key in keys:
        print(f"{key}: verified -> {paths[key]}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
