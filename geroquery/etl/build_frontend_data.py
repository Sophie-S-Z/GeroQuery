"""Export the static assets the browser front end queries.

Run: ``python -m geroquery.etl.build_frontend_data``

The front end has no backend. It is a static bundle on a CDN that reads these
two Parquet files in the browser with hyparquet, a ~15 kB pure-JS reader.

The first cut of this used DuckDB-WASM. That was wrong, and the Cloudflare
Pages 25 MiB per-file limit is what surfaced it: the DuckDB wasm binaries are
34 MB and 39 MB, which is roughly seven times the size of the data they were
there to query. Range requests only pay when the data dwarfs the engine. Here
it is the other way round, so the whole corpus is simply fetched once — 6.6 MB,
less than a fifth of the engine alone — and every query after that is local and
instant.

**Every statistic in the output is computed here, in Python.** The browser does
no meta-analysis: it reads pooled effects that this module already pooled, using
the same :func:`geroquery.harmonize.random_effects` the API serves. That is a
deliberate architectural line. Re-implementing DerSimonian-Laird in JavaScript
would create two estimators for one quantity, and they would drift — the repo
already keeps `test_matrix_hedges_g_matches_scalar_implementation` around
because a vectorized copy of one estimator was worth pinning to its reference.
Across a language boundary the risk is worse and the tooling is thinner, so the
line is: **JavaScript renders numbers, it never derives them.**

Outputs land in ``frontend/public/data/``:

``pooled.parquet``
    One row per gene per species: the pooled effect, its interval, k, I-squared.
    Sorted by symbol so the search path prunes.

``contrasts.parquet``
    Every per-study effect size, sorted by gene_id so a symbol lookup scans a
    contiguous run.

Both are written with **snappy** rather than zstd: the browser reader
(hyparquet) has snappy built in, and zstd would mean shipping a second
decompressor to save a few hundred kilobytes on a file that is already under
five megabytes.

``studies.json`` · ``crosslayer.json`` · ``meta.json``
    Small enough to fetch whole.

``meta.json`` records whether the store was built from the full signature table
or the committed curated slice, and the front end displays it. A number derived
from the 40,585-row slice must never be presented as the 485,905-row result —
same contract as ``NhanesClinicalSource.clinical_frame`` returning its mode.
"""

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import pandas as pd

from .. import __version__
from ..harmonize import random_effects
from ..sources.manifest import MANIFEST, MANIFEST_VERSION, VERIFIED_ON
from ..store import GeroStore

# Repo-root-relative output directory. Vite serves `public/` at the site root,
# so these land at /data/*.
DEFAULT_OUT = Path(__file__).resolve().parents[2] / "frontend" / "public" / "data"

# Row count of the committed curated slice. Used only to label the build mode;
# a developer who has run `make data` gets the full table and CI does not.
CURATED_SLICE_ROWS = 40_585

# Genes pooled from fewer contrasts than this get no pooled row. With one or two
# contrasts the between-study variance is not estimable, so the interval would
# be a shape rather than a summary.
MIN_CONTRASTS = 3


def _pool(frame: pd.DataFrame) -> pd.DataFrame:
    """Pool every (gene, species) group with enough contrasts.

    Uses the same estimator the API serves, so the site and the API cannot
    disagree about what a gene's effect is.
    """
    rows = []
    grouped = frame.dropna(subset=["effect_size", "standard_error"]).groupby(
        ["gene_id", "species"], observed=True
    )
    for (gene_id, species), group in grouped:
        if len(group) < MIN_CONTRASTS:
            continue
        pooled = random_effects(group["effect_size"].tolist(), group["standard_error"].tolist())
        rows.append(
            {
                "gene_id": gene_id,
                "species": species,
                "g": round(pooled.pooled_effect, 4),
                "se": round(pooled.standard_error, 4),
                "ci_low": round(pooled.ci_low, 4),
                "ci_high": round(pooled.ci_high, 4),
                "p_value": pooled.p_value,
                "i2": round(pooled.i2, 1),
                "tau2": round(pooled.tau2, 4),
                "k": int(pooled.n_studies),
                # The verdict is computed here, once, from the interval — not in
                # the browser from a threshold someone might change.
                "verdict": (
                    "increases"
                    if pooled.ci_low > 0
                    else "decreases" if pooled.ci_high < 0 else "no_evidence"
                ),
                "n_tissues": int(group["tissue"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def _symbols(gene_ids: list[str]) -> dict[str, str]:
    """Map canonical Ensembl ids to display symbols, so search is by name.

    Reads the consolidated bulk id maps that ``make idmap`` already writes, and
    inverts them. Those maps are keyed by the identifier that was *queried*
    (Entrez id or symbol) but every resolved record carries its Ensembl
    canonical id, so the mapping needed here is already on disk — no network,
    no second source of truth.

    The obvious implementation, looping ``resolver.resolve_gene`` over 46,091
    ids, resolves **5.6%** of them: the bundled ``genes.json`` is a curated
    2,560-record table for the request path, not a genome-wide symbol index.
    That left the search rail unable to find 94% of the corpus and the
    landscape table showing raw accessions. Worth stating because the failure
    was silent — every id still got a label, the label was just its own id.
    """
    from ..idmap.bulk import map_path

    wanted = set(gene_ids)
    out: dict[str, str] = {}
    for species in ("human", "mouse"):
        path = map_path(species)
        if not path.exists():
            continue
        with open(path, encoding="utf-8") as handle:
            mapping = json.load(handle)
        for record in mapping.values():
            if not record:
                continue
            canonical, symbol = record.get("canonical_id"), record.get("symbol")
            if canonical in wanted and symbol and canonical not in out:
                out[canonical] = symbol

    if not out:
        # Better to say so than to ship a site where every gene is an accession.
        raise FileNotFoundError(
            "No bulk id maps found under CACHE_HOME/idmap. Run `make idmap` "
            "before building the front-end data, or the search rail will be "
            "unusable."
        )
    return out


def _crosslayer_payload() -> dict:
    """The mortality result, precomputed. Absent if the cohort was never built.

    Not fatal: the cross-layer table is git-ignored and reproducible, so a clean
    checkout can build the site without it. The front end renders a stated
    "not built" panel rather than an empty chart.
    """
    from ..sources import nhanes_dnam as nd

    try:
        frame, mode = nd.NhanesDnamSource().crosslayer_frame(prefer_full=True)
    except Exception as exc:  # noqa: BLE001
        return {"available": False, "reason": str(exc)}

    from ..survival import crosslayer_analysis, mahalanobis_dysregulation

    dysregulation = mahalanobis_dysregulation(frame, list(nd.MARKERS))
    predicts = {c.column: c.predicts for c in nd.CLOCKS}
    notes = {c.column: c.note for c in nd.CLOCKS}

    clocks = []
    for column in nd.AGE_LIKE_CLOCKS:
        result = crosslayer_analysis(
            frame,
            nd.age_acceleration(frame, column),
            dysregulation.values,
            clock_name=column,
            predicts=predicts.get(column, ""),
        )
        joint = {row["covariate"]: row for row in result.joint_model.summary_rows()}
        alone = {row["covariate"]: row for row in result.clock_model.summary_rows()}
        clocks.append(
            {
                "clock": column,
                "predicts": predicts.get(column, ""),
                "note": notes.get(column, ""),
                "hr_alone": alone["clock_acceleration"]["hazard_ratio"],
                "hr": joint["clock_acceleration"]["hazard_ratio"],
                "ci_low": joint["clock_acceleration"]["ci_low"],
                "ci_high": joint["clock_acceleration"]["ci_high"],
                "excludes_null": joint["clock_acceleration"]["excludes_null"],
                "dysregulation_hr": joint["dysregulation"]["hazard_ratio"],
                "dysregulation_ci_low": joint["dysregulation"]["ci_low"],
                "dysregulation_ci_high": joint["dysregulation"]["ci_high"],
                "c_baseline": round(result.baseline.concordance, 4),
                "c_clock": round(result.clock_model.concordance, 4),
                "c_dysregulation": round(result.dysregulation_model.concordance, 4),
                "c_joint": round(result.joint_model.concordance, 4),
                "dysregulation_adds_p": result.dysregulation_adds_over_clock["p_value"],
                "clock_adds_p": result.clock_adds_over_dysregulation["p_value"],
            }
        )

    clocks.sort(key=lambda row: -row["hr"])
    return {
        "available": True,
        "mode": mode,
        "n_subjects": int(len(frame)),
        "n_deaths": int(frame["died"].sum()),
        "followup_years_median": round(float(frame["followup_years"].median()), 2),
        "age_min": float(frame["age"].min()),
        "age_max": float(frame["age"].max()),
        "sex_qc": nd.sex_discordance(frame),
        "dysregulation": dysregulation.to_dict(),
        "clocks": clocks,
        "caveats": list(nd.CAVEATS),
    }


def build(out_dir: Path | None = None, *, store: GeroStore | None = None) -> dict:
    out_dir = out_dir or DEFAULT_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    store = (store or GeroStore()).ensure_built()

    contrasts = pd.read_parquet(store.sig_dir)
    for column in ("species", "omic_layer"):
        contrasts[column] = contrasts[column].astype(str)

    pooled = _pool(contrasts)
    symbols = _symbols(sorted(pooled["gene_id"].unique()))
    pooled["symbol"] = pooled["gene_id"].map(symbols).fillna(pooled["gene_id"])
    # Sorted by symbol: the search path scans this file, and a sorted column
    # lets DuckDB's row-group statistics skip almost all of it.
    pooled = pooled.sort_values(["symbol", "species"]).reset_index(drop=True)
    pooled.to_parquet(out_dir / "pooled.parquet", index=False, compression="snappy")

    # Sorted by gene_id: a forest plot reads one contiguous run of row groups
    # instead of touching every one.
    contrasts = contrasts.sort_values(["gene_id", "study_id"]).reset_index(drop=True)
    contrasts["symbol"] = contrasts["gene_id"].map(symbols).fillna(contrasts["gene_id"])
    contrasts.to_parquet(out_dir / "contrasts.parquet", index=False, compression="snappy")

    studies = [s.to_dict() for s in store.list_studies()]
    (out_dir / "studies.json").write_text(json.dumps(studies), encoding="utf-8")

    crosslayer = _crosslayer_payload()
    (out_dir / "crosslayer.json").write_text(json.dumps(crosslayer), encoding="utf-8")

    n_rows = int(len(contrasts))
    meta = {
        "built_on": date.today().isoformat(),
        "code_version": __version__,
        "data_version": store.version(),
        "manifest_version": MANIFEST_VERSION,
        "checksums_verified_on": VERIFIED_ON,
        "n_pinned_artifacts": len(MANIFEST),
        "n_contrasts_rows": n_rows,
        "n_genes": int(contrasts["gene_id"].nunique()),
        "n_pooled": int(len(pooled)),
        "n_studies": len(studies),
        "n_series": len({s.get("series_id") for s in studies if s.get("series_id")}),
        "n_tissues": int(contrasts["tissue"].nunique()),
        "species": sorted(contrasts["species"].unique()),
        # The mode contract, carried onto the page. A masthead that says
        # "485,905" when the store held the 40,585-row slice would be exactly
        # the self-inflicted wound test_masthead_counts_are_computed_not_typed
        # exists to prevent.
        "signature_mode": "full" if n_rows > CURATED_SLICE_ROWS else "curated_slice",
        "min_contrasts_for_pooling": MIN_CONTRASTS,
        "crosslayer_available": crosslayer.get("available", False),
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")

    sizes = {p.name: p.stat().st_size for p in sorted(out_dir.glob("*")) if p.is_file()}
    return {**meta, "out_dir": str(out_dir), "bytes": sizes}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, default=None, help="output directory")
    args = parser.parse_args(argv)

    summary = build(args.out)
    for key, value in summary.items():
        if key == "bytes":
            for name, size in value.items():
                print(f"  {name}: {size / 1_048_576:.2f} MB")
        else:
            print(f"{key}: {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
