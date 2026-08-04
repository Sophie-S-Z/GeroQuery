"""Generate the bundled identifier tables from real sources.

Run: ``GEROQUERY_ALLOW_NETWORK=1 python -m geroquery.etl.build_idmap``

These two files were the last hand-written tables in GeroQuery. Both were
*accurate* — the identifiers and lifespans in them are real — but neither was
sourced, so nothing could say when they went stale or where a number came from.
That is the same objection this project raised against the synthetic signature
layer, and it applies to a 24-row table exactly as much as to a 200,000-row one.

``idmap/data/genes.json``
    Was 24 hand-picked "well-studied aging genes". Now: every gene carrying a
    HAGR curated assertion, plus its cross-species ortholog, resolved through
    mygene.info. That is ~1,880 genes rather than 24, which is also a real
    functional gain — the offline resolver could previously name 24 of the 1,880
    genes in its own curated table, and silently returned nothing for the rest.

``idmap/data/anage.json``
    Was 9 hand-transcribed maximum lifespans. Now: parsed from the pinned AnAge
    release, which is the database those numbers were copied out of.

Both are committed, because the resolver must work with no network. The
difference is that they are now *derived* artifacts with a recorded provenance
and a re-run command, not typed-in constants.
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

from .. import config
from ..idmap.bulk import resolve_bulk
from ..sources import hagr
from ..sources.fetch import fetch_artifact
from ..sources.manifest import ANAGE, HAGR, MANIFEST_VERSION

GENES_FILENAME = "genes.json"
ANAGE_FILENAME = "anage.json"

# Species the rest of GeroQuery models, plus the model organisms that appear in
# HAGR's intervention tables — a fractional-lifespan figure for a DrugAge result
# in C. elegans is only meaningful if the worm's maximum lifespan is known.
ANAGE_SPECIES: dict[str, tuple[str, tuple[str, ...]]] = {
    "human": ("Homo sapiens", ("human", "homo sapiens")),
    "mouse": ("Mus musculus", ("mouse", "mus musculus", "house mouse")),
    "rat": ("Rattus norvegicus", ("rat", "rattus norvegicus", "norway rat")),
    "worm": ("Caenorhabditis elegans", ("worm", "c. elegans", "caenorhabditis elegans")),
    "fly": ("Drosophila melanogaster", ("fly", "drosophila", "drosophila melanogaster")),
    "yeast": ("Saccharomyces cerevisiae", ("yeast", "saccharomyces cerevisiae")),
    "zebrafish": ("Danio rerio", ("zebrafish", "danio rerio")),
    "killifish": ("Nothobranchius furzeri", ("killifish", "nothobranchius furzeri")),
    "naked_mole_rat": ("Heterocephalus glaber", ("naked mole rat", "heterocephalus glaber")),
}

ANAGE_MEMBER = "anage_data.txt"


def build_anage(*, allow_network: bool | None = None) -> dict:
    """Parse maximum longevity for the modelled species out of the AnAge release."""
    path = fetch_artifact(ANAGE, allow_network=allow_network)
    with zipfile.ZipFile(path) as archive:
        name = next((n for n in archive.namelist() if n.endswith(ANAGE_MEMBER)), None)
        if name is None:
            raise RuntimeError(f"{ANAGE_MEMBER} not in {path.name}: {archive.namelist()}")
        text = archive.read(name).decode("utf-8", "replace")

    rows = [line.split("\t") for line in text.splitlines() if line.strip()]
    header = rows[0]
    genus_i, species_i = header.index("Genus"), header.index("Species")
    longevity_i = header.index("Maximum longevity (yrs)")

    by_latin: dict[str, float] = {}
    for row in rows[1:]:
        if len(row) <= longevity_i:
            continue
        latin = f"{row[genus_i].strip()} {row[species_i].strip()}"
        try:
            by_latin[latin] = float(row[longevity_i])
        except ValueError:
            continue

    species = []
    missing = []
    for key, (latin, common) in ANAGE_SPECIES.items():
        value = by_latin.get(latin)
        if value is None:
            missing.append(latin)
            continue
        species.append(
            {
                "key": key,
                "latin": latin,
                "common": list(common),
                "max_lifespan_years": value,
            }
        )
    return {
        "_provenance": (
            f"Maximum longevity (years) parsed from the pinned AnAge release "
            f"({ANAGE.url}, SHA-256 {ANAGE.sha256[:16]}...). Regenerate with "
            f"`python -m geroquery.etl.build_idmap`. Manifest {MANIFEST_VERSION}."
        ),
        "_missing_from_anage": missing,
        "species": species,
    }


# Well-known synonyms that mygene does not carry. Vocabulary, not evidence: each
# is a name the literature actually uses for the gene, and losing it makes the
# resolver answer "not found" to a query a biologist would consider obvious.
# Kept deliberately short — this is not a place to accumulate guesses.
ALIAS_SUPPLEMENT: dict[str, tuple[str, ...]] = {
    "KL": ("klotho", "alpha-klotho"),
    "CDKN2A": ("p16ink4a", "p16-ink4a"),
    "CDKN1A": ("p21cip1", "p21waf1"),
    "TERT": ("telomerase", "htert"),
    "MTOR": ("mtor", "frap1"),
    "IGF1": ("igf-1", "somatomedin c"),
}


def build_genes(*, allow_network: bool | None = None, quiet: bool = False) -> dict:
    """Resolve every HAGR-curated gene and its ortholog into the bundled table.

    **Resolved by Entrez id wherever HAGR supplies one**, and only by symbol
    otherwise. That is not a preference, it is a correctness requirement: mygene
    resolves the symbol ``CDKN1A`` to ``ENSG00000205643`` — a copy of the gene on
    an alternative MHC haplotype scaffold — while the Entrez id ``1026`` resolves
    to the primary ``ENSG00000124762``. The signature table is keyed on Entrez,
    because that is what GEO platform annotations supply, so a symbol-built gene
    table silently fails to join: p21 had twelve real contrasts and the API
    reported none.
    """
    paths = {key: fetch_artifact(HAGR[key], allow_network=allow_network) for key in HAGR}
    assertions = hagr.load_assertions(paths)

    # Entrez first, symbol only as a fallback (LongevityMap indexes by symbol).
    by_species: dict[str, set[str]] = {"human": set(), "mouse": set()}
    symbol_only: set[str] = set()
    for item in assertions:
        if not item.species:
            continue
        if item.entrez and item.entrez.isdigit():
            by_species[item.species].add(item.entrez)
        else:
            symbol_only.add(item.symbol)

    # HAGR curates mostly human genes. Querying every symbol against mouse too is
    # what makes the offline table cross-species, which the ortholog view needs.
    records: dict[str, dict] = {}
    for species in ("human", "mouse"):
        queries = sorted(by_species[species] | symbol_only)
        resolved = resolve_bulk(queries, species, allow_network=allow_network, progress=not quiet)
        for query, record in resolved.items():
            if record is None:
                continue
            # An Entrez-keyed hit wins over a symbol-keyed one for the same gene.
            existing = records.get(record["canonical_id"])
            if existing is None or (query.isdigit() and not existing.get("_from_entrez")):
                record = {**record, "_from_entrez": query.isdigit()}
                records[record["canonical_id"]] = record

    genes = [
        {
            "canonical_id": r["canonical_id"],
            "symbol": r["symbol"],
            "species": r["species"],
            "entrez": r.get("entrez"),
            "ensembl": r.get("ensembl"),
            "uniprot": r.get("uniprot"),
            "name": r.get("name"),
            "ortholog_group": r["ortholog_group"],
            # Not truncated. Slicing the sorted list to a fixed length looks
            # harmless and is not: aliases sort alphabetically, so a cap drops
            # the tail, and "p16" is late in CDKN2A's alias list. The offline
            # resolver then answered "p16" with the mouse ortholog, because the
            # human record no longer claimed the alias.
            "aliases": sorted(
                {a for a in (r.get("aliases") or [])}
                | set(ALIAS_SUPPLEMENT.get(r["symbol"].upper(), ()))
            ),
        }
        for r in sorted(records.values(), key=lambda x: (x["ortholog_group"], x["species"]))
    ]
    return {
        "_provenance": (
            "Every gene carrying a HAGR curated assertion (GenAge, CellAge, "
            "LongevityMap, GenDR), plus its cross-species ortholog, resolved through "
            "mygene.info. Regenerate with `python -m geroquery.etl.build_idmap`. "
            f"Manifest {MANIFEST_VERSION}."
        ),
        # Retired Ensembl ids that must still resolve. Kept by hand because they
        # are a property of *our* history — ids we have published — not of any
        # upstream release, so no upstream can supply them.
        "deprecated_ids": {
            "ENSG00000147888": "ENSG00000147889",
            "ENSG00000267174": "ENSG00000141510",
        },
        "genes": genes,
    }


def build_all(
    out_dir: Path | None = None, *, allow_network: bool | None = None, quiet: bool = False
) -> dict[str, int]:
    out_dir = out_dir or config.IDMAP_DATA
    out_dir.mkdir(parents=True, exist_ok=True)

    anage = build_anage(allow_network=allow_network)
    genes = build_genes(allow_network=allow_network, quiet=quiet)

    (out_dir / ANAGE_FILENAME).write_text(json.dumps(anage, indent=1), encoding="utf-8")
    (out_dir / GENES_FILENAME).write_text(json.dumps(genes, indent=1), encoding="utf-8")
    return {
        "species": len(anage["species"]),
        "species_missing_from_anage": len(anage["_missing_from_anage"]),
        "genes": len(genes["genes"]),
        "ortholog_groups": len({g["ortholog_group"] for g in genes["genes"]}),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--offline", action="store_true", help="use only cached artifacts")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    for key, value in build_all(allow_network=not args.offline, quiet=args.quiet).items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
