"""Real-data ETL for gene aging signatures and curated knowledge.

Run: ``GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_signatures``

This is the module that removed the last generated layer from GeroQuery. It
fetches the checksum-pinned GEO DataSets panel, derives a young-vs-old contrast
per dataset under the rules in :mod:`geroquery.sources.geo`, estimates a Hedges'
*g* per gene, resolves every probe to a canonical gene id, and writes the tables
the store builds from. It also loads the five HAGR curated databases.

Outputs in ``geroquery/sources/data/``:

``signatures_full.csv``
    Every gene of every contrast. Git-ignored: it is reproducible from the
    manifest, and a committed copy would be a second, unverifiable one.

``signatures_curated.csv``
    The rows whose gene carries a HAGR assertion, plus each such gene's
    cross-species ortholog. Committed, so tests and CI run offline. This is a
    *gene-restricted slice of the same real estimates*, not a different or
    simulated calculation — but it is a slice, so a number computed on it is not
    the full-panel number.

``studies.csv``, ``curated_knowledge.csv``, ``interventions.csv``
    Small enough to commit whole.

Runtime is dominated by parsing ~300 MB of gzipped SOFT, not by the network.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ..config import SOURCES_DATA
from ..harmonize.differential import contrast_effects
from ..idmap.bulk import resolve_bulk
from ..sources import geo, hagr
from ..sources.fetch import fetch_artifact
from ..sources.manifest import (
    GEO_AGING_PANEL,
    GEO_LICENSE,
    HAGR,
    MANIFEST_VERSION,
    VERIFIED_ON,
)

# Every dataset in the panel is an expression array.
OMIC_LAYER = "transcriptome"

FULL_SIGNATURES = "signatures_full.csv"
CURATED_SIGNATURES = "signatures_curated.csv"
STUDIES = "studies.csv"
CURATED_KNOWLEDGE = "curated_knowledge.csv"
INTERVENTIONS = "interventions.csv"

SIGNATURE_COLUMNS = [
    "gene_id",
    "study_id",
    "omic_layer",
    "species",
    "tissue",
    "sex",
    "age_range",
    "effect_size",
    "direction",
    "p_value",
    "q_value",
    "standard_error",
    "source",
]

# Effect sizes are reported to four decimals. More would imply a precision these
# sample sizes do not support; fewer would collapse distinct small effects.
ROUND_TO = 4


def _log(message: str, *, quiet: bool = False) -> None:
    if not quiet:
        print(message, flush=True)


def contrast_frames(
    *, allow_network: bool | None = None, accessions: list[str] | None = None, quiet: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Estimate every contrast in the panel.

    Returns ``(signatures, studies, skipped)``. ``signatures`` is keyed by
    ``gene_key`` (Entrez id, or ``SYM:<symbol>`` where the platform annotation
    gave no Entrez id) — resolution to canonical ids happens once, afterwards,
    across the whole panel.
    """
    keys = accessions or sorted(GEO_AGING_PANEL)
    signature_parts: list[pd.DataFrame] = []
    study_rows: list[dict] = []
    skipped: list[str] = []

    for index, accession in enumerate(keys, 1):
        artifact = GEO_AGING_PANEL[accession]
        path = fetch_artifact(artifact, allow_network=allow_network)
        header, table = geo.read_soft(path)
        contrasts, reasons = geo.build_contrasts(header)
        skipped.extend(reasons)
        _log(
            f"[{index:2d}/{len(keys)}] {accession} {header.organism} "
            f"{len(table):>6d} probes -> {len(contrasts)} contrast(s)",
            quiet=quiet,
        )

        gene_keys, _identifiers = geo.probe_annotation(table)
        for contrast in contrasts:
            effects = contrast_effects(
                table,
                list(contrast.young_samples),
                list(contrast.old_samples),
                genes=gene_keys,
            )
            if effects.empty:
                skipped.append(f"{contrast.study_id}: no gene passed probe annotation")
                continue
            frame = effects.reset_index(names="gene_key")
            frame["study_id"] = contrast.study_id
            frame["omic_layer"] = OMIC_LAYER
            frame["species"] = contrast.species
            frame["tissue"] = contrast.tissue
            frame["sex"] = contrast.sex
            frame["age_range"] = contrast.age_range
            frame["source"] = "GEO"
            signature_parts.append(frame)

            study_rows.append(
                {
                    "study_id": contrast.study_id,
                    "source": "GEO",
                    "omic_layer": OMIC_LAYER,
                    "species": contrast.species,
                    "tissue": contrast.tissue,
                    "sample_size": len(contrast.young_samples) + len(contrast.old_samples),
                    "processing_method": (
                        f"GDS curated matrix ({header.value_type}); "
                        f"MaxMean probe collapse; Hedges' g"
                    ),
                    "license": GEO_LICENSE,
                    "url": header.url,
                    "version": artifact.release,
                    "series_id": header.series_id,
                    "platform": header.platform,
                    "pubmed_id": header.pubmed_id,
                    "n_young": len(contrast.young_samples),
                    "n_old": len(contrast.old_samples),
                    "age_range": contrast.age_range,
                    "sex": contrast.sex,
                    "restrictions": " ; ".join(contrast.restrictions),
                }
            )

    signatures = (
        pd.concat(signature_parts, ignore_index=True) if signature_parts else pd.DataFrame()
    )
    return signatures, pd.DataFrame(study_rows), skipped


def resolve_gene_ids(
    signatures: pd.DataFrame, *, allow_network: bool | None = None, quiet: bool = False
) -> tuple[pd.DataFrame, dict[str, int], dict[str, str]]:
    """Attach a canonical ``gene_id`` to every signature row.

    Also returns ``{gene_id: symbol}``, which is what lets the committed offline
    slice span both species: HAGR curates mostly human genes, and matching on
    upper-cased symbol is how this project already pairs orthologs.

    Probes whose identifier resolves to nothing keep an ``ENTREZ:<id>`` id when
    the platform gave an Entrez id — that is a stable identifier, just not an
    Ensembl one — and are dropped when all we had was an unrecognised symbol.
    Keeping a made-up id would put rows in the store that can never join to
    curated knowledge and cannot be traced to a gene.
    """
    stats: dict[str, int] = {}
    resolved: dict[tuple[str, str], str] = {}
    symbols: dict[str, str] = {}

    for species, group in signatures.groupby("species", sort=True):
        keys = sorted(set(group["gene_key"]))
        queries = [k[len("SYM:") :] if k.startswith("SYM:") else k for k in keys]
        mapping = resolve_bulk(
            queries, str(species), allow_network=allow_network, progress=not quiet
        )
        hits = 0
        for key, query in zip(keys, queries, strict=True):
            record = mapping.get(query)
            if record is not None:
                resolved[(str(species), key)] = record["canonical_id"]
                symbols.setdefault(record["canonical_id"], record["symbol"])
                hits += 1
            elif not key.startswith("SYM:"):
                resolved[(str(species), key)] = f"ENTREZ:{key}"
        stats[f"{species}_keys"] = len(keys)
        stats[f"{species}_resolved"] = hits

    out = signatures.copy()
    out["gene_id"] = [
        resolved.get((species, key))
        for species, key in zip(out["species"], out["gene_key"], strict=True)
    ]
    stats["rows_before"] = len(out)
    out = out[out["gene_id"].notna()]
    stats["rows_dropped_unresolved"] = stats["rows_before"] - len(out)

    # One canonical id can collect several probe keys (deprecated Entrez ids that
    # now point at the same gene). Keep the most significant row per gene per
    # study rather than emitting the gene twice with different effect sizes.
    out = out.sort_values("p_value", kind="stable").drop_duplicates(["study_id", "gene_id"])
    stats["rows_dropped_duplicate_gene"] = (
        len(signatures) - stats["rows_dropped_unresolved"] - len(out)
    )
    for column in ("effect_size", "standard_error", "p_value", "q_value"):
        out[column] = out[column].round(ROUND_TO)
    ordered = out[SIGNATURE_COLUMNS].sort_values(["gene_id", "study_id"]).reset_index(drop=True)
    return ordered, stats, symbols


def build_curated(
    *, allow_network: bool | None = None, quiet: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Load HAGR into ``(curated_knowledge, interventions, stats)``."""
    paths = {key: fetch_artifact(HAGR[key], allow_network=allow_network) for key in HAGR}
    assertions = hagr.load_assertions(paths)
    interventions = hagr.load_interventions(paths)

    stats = {"assertions_parsed": len(assertions), "interventions_parsed": len(interventions)}

    by_species: dict[str, list[str]] = {}
    for item in assertions:
        if item.species:
            by_species.setdefault(item.species, []).append(item.entrez or item.symbol)
    mapping: dict[tuple[str, str], str] = {}
    for species, queries in sorted(by_species.items()):
        resolved = resolve_bulk(queries, species, allow_network=allow_network, progress=not quiet)
        for query, record in resolved.items():
            if record is not None:
                mapping[(species, query)] = record["canonical_id"]

    knowledge_rows = []
    symbol_to_gene: dict[str, str] = {}
    for item in assertions:
        if not item.species:
            continue
        gene_id = mapping.get((item.species, item.entrez or item.symbol))
        if gene_id is None:
            continue
        symbol_to_gene.setdefault(item.symbol.upper(), gene_id)
        knowledge_rows.append(
            {
                "gene_id": gene_id,
                "database": item.database,
                "assertion": item.assertion,
                "url": item.url,
            }
        )
    knowledge = pd.DataFrame(knowledge_rows).drop_duplicates()
    stats["curated_symbols"] = len({a.symbol.upper() for a in assertions if a.species})
    stats["assertions_loaded"] = len(knowledge)
    stats["assertions_dropped_non_mammalian"] = sum(1 for a in assertions if not a.species)

    intervention_rows = [
        {
            "intervention_id": item.intervention_id,
            "name": item.name,
            "itype": item.itype,
            "source": item.source,
            "organism": item.organism,
            "lifespan_effect_pct": item.lifespan_effect_pct,
            "linked_gene_ids": "|".join(
                sorted({g for s in item.gene_symbols if (g := symbol_to_gene.get(s.upper()))})
            ),
            "url": item.url,
            "n_experiments": item.n_experiments,
        }
        for item in interventions
    ]
    return knowledge, pd.DataFrame(intervention_rows), stats


def build_all(
    data_dir: Path | None = None,
    *,
    allow_network: bool | None = None,
    quiet: bool = False,
) -> dict[str, object]:
    """Run the whole signature + curated-knowledge ETL and write every table."""
    data_dir = data_dir or SOURCES_DATA
    data_dir.mkdir(parents=True, exist_ok=True)

    _log("Estimating GEO age contrasts...", quiet=quiet)
    raw, studies, skipped = contrast_frames(allow_network=allow_network, quiet=quiet)
    if raw.empty:
        raise RuntimeError("The GEO panel produced no signatures; refusing to write empty tables.")

    _log("Resolving gene identifiers...", quiet=quiet)
    signatures, id_stats, symbols = resolve_gene_ids(raw, allow_network=allow_network, quiet=quiet)

    _log("Loading HAGR curated knowledge...", quiet=quiet)
    knowledge, interventions, hagr_stats = build_curated(allow_network=allow_network, quiet=quiet)

    signatures.to_csv(data_dir / FULL_SIGNATURES, index=False)
    studies.to_csv(data_dir / STUDIES, index=False)
    knowledge.to_csv(data_dir / CURATED_KNOWLEDGE, index=False)
    interventions.to_csv(data_dir / INTERVENTIONS, index=False)

    # The committed slice: every gene carrying a HAGR assertion, plus the
    # cross-species ortholog of each. Without the ortholog step the slice would
    # be almost entirely human — HAGR curates human genes — and the mouse half
    # of the panel would be invisible to anyone who has not run `make data`.
    curated_genes = set(knowledge["gene_id"])
    curated_symbols = {s.upper() for gene, s in symbols.items() if gene in curated_genes}
    keep = signatures["gene_id"].isin(curated_genes) | signatures["gene_id"].map(
        lambda g: symbols.get(g, "").upper() in curated_symbols
    )
    sample = signatures[keep]
    sample.to_csv(data_dir / CURATED_SIGNATURES, index=False)

    return {
        "manifest_version": MANIFEST_VERSION,
        "checksums_verified_on": VERIFIED_ON,
        "datasets": len(GEO_AGING_PANEL),
        "contrasts": len(studies),
        "series": int(studies["series_id"].nunique()),
        "signature_rows": len(signatures),
        "genes": int(signatures["gene_id"].nunique()),
        "curated_slice_rows": len(sample),
        "curated_slice_genes": int(sample["gene_id"].nunique()),
        "curated_slice_species": int(sample["species"].nunique()),
        "curated_assertions": len(knowledge),
        "interventions": len(interventions),
        "skipped": len(skipped),
        **id_stats,
        **hagr_stats,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--offline",
        action="store_true",
        help="use only already-verified cached artifacts; fail rather than download",
    )
    parser.add_argument("--quiet", action="store_true", help="suppress per-dataset progress")
    args = parser.parse_args(argv)

    summary = build_all(allow_network=not args.offline, quiet=args.quiet)
    print()
    for key, value in summary.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
