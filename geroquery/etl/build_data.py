"""Real-data ETL: fetch pinned upstream artifacts and materialize local tables.

Run: ``python -m geroquery.etl.build_data``

Distinct from :mod:`geroquery.etl.build_fixtures`, which generates the *synthetic*
method-validation slice. This module only ever writes data that came off a real
upstream and passed its pinned SHA-256.

Outputs land in ``geroquery/sources/data/``:

``clinical_nhanes_full.csv``
    All 4,895 complete NHANES 2017-2018 cases aged 20+. Git-ignored — it is
    reproducible from the manifest, so committing it would just be a second,
    unverifiable copy.

``nhanes_2017_2018_sample.csv`` (with ``--write-sample``)
    Seeded 600-row subset of the same real rows. Committed, so tests and CI run
    offline. Regenerate only when the manifest changes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..config import SOURCES_DATA
from ..sources.manifest import MANIFEST_VERSION, VERIFIED_ON
from ..sources.nhanes import RELEASE, full_path, load_full, write_sample


def build_nhanes(
    data_dir: Path | None = None,
    *,
    allow_network: bool | None = None,
    write_offline_sample: bool = False,
) -> dict[str, object]:
    """Fetch, verify, harmonize, and write the real NHANES clinical table."""
    data_dir = data_dir or SOURCES_DATA
    frame = load_full(allow_network=allow_network)

    out = full_path(data_dir)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)

    result: dict[str, object] = {
        "release": RELEASE,
        "manifest_version": MANIFEST_VERSION,
        "checksums_verified_on": VERIFIED_ON,
        "subjects": len(frame),
        "age_min": float(frame["age"].min()),
        "age_max": float(frame["age"].max()),
        "path": str(out),
    }
    if write_offline_sample:
        result["sample_path"] = str(write_sample(frame, data_dir))
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--write-sample",
        action="store_true",
        help="also regenerate the committed offline sample (only when the manifest changes)",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="use only already-cached artifacts; fail rather than download",
    )
    args = parser.parse_args(argv)

    summary = build_nhanes(
        allow_network=False if args.offline else True,
        write_offline_sample=args.write_sample,
    )
    for key, value in summary.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
