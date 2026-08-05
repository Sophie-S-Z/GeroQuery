"""Real-data ETL: build the NHANES 1999-2002 cross-layer cohort.

Run: ``python -m geroquery.etl.build_crosslayer [--write-sample]``

Joins three pinned upstreams on ``SEQN`` — NCHS-computed DNAm clocks, the
1999-2002 clinical panel, and 2019 linked mortality — into the one table in this
repo where a clock, a health state, and an outcome are observed on the same
person.

Outputs land in ``geroquery/sources/data/``:

``crosslayer_nhanes_dnam_full.csv``
    All 2,517 complete cases. Git-ignored: reproducible from the manifest, so a
    committed copy would just be a second unverifiable one.

``nhanes_dnam_1999_2002_sample.csv`` (with ``--write-sample``)
    Seeded 600-row subset of the same real rows. Committed, so tests run
    offline. Estimates computed on it are *not* the reported result — the
    adapter returns its mode so a caller can tell which produced a number.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..config import SOURCES_DATA
from ..sources.manifest import MANIFEST_VERSION, VERIFIED_ON
from ..sources.nhanes_dnam import (
    AGE_LIKE_CLOCKS,
    CLOCKS,
    RELEASE,
    full_path,
    load_full,
    sex_discordance,
    write_sample,
)


def build_crosslayer(
    data_dir: Path | None = None,
    *,
    allow_network: bool | None = None,
    write_offline_sample: bool = False,
) -> dict[str, object]:
    """Fetch, verify, join, and write the cross-layer cohort."""
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
        "deaths": int(frame["died"].sum()),
        "n_clocks": len(CLOCKS),
        "n_age_like_clocks": len(AGE_LIKE_CLOCKS),
        "age_min": float(frame["age"].min()),
        "age_max": float(frame["age"].max()),
        "followup_years_median": round(float(frame["followup_years"].median()), 2),
        "followup_years_max": round(float(frame["followup_years"].max()), 2),
        # Reported every build: a rising rate is the cheapest signal that
        # samples were mixed up upstream.
        "dnam_sex_qc": sex_discordance(frame),
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

    summary = build_crosslayer(
        allow_network=False if args.offline else True,
        write_offline_sample=args.write_sample,
    )
    for key, value in summary.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
