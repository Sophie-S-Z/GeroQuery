"""Deterministically build the bundled, honest data slice.

Run: ``python -m geroquery.etl.build_fixtures``

What this writes is intentionally, verifiably real:

* ``curated_knowledge.csv`` — real membership of each gene in the HAGR/OpenGenes
  curated databases (GenAge, CellAge, LongevityMap, OpenGenes). These are
  genuine curated facts, linked back to the source portal.
* ``interventions.csv`` — real lifespan-affecting interventions (rapamycin,
  caloric restriction, metformin, senolytics, …) with approximate reported
  rodent effect sizes and a link to the primary publication. Generated from the
  curated knowledge base so numbers and citations stay in sync.
* ``example_cohort_simulated.csv`` — a **clearly-labelled simulated** cohort of
  nine clinical-chemistry markers (the exact inputs of the real PhenoAge clock)
  plus age and sex. It exists only so the clock and resilience tools can be
  tried without uploading data; it is a transparent worked example, never
  presented as measurements of real people. Upload-your-own is the first-class
  path.

There are **no fabricated GEO accessions, effect sizes, or p-values** here — the
earlier "demo signature slice" that invented those has been removed. Everything
GeroQuery asserts about a gene's relationship to aging now comes from the cited
knowledge base in :mod:`geroquery.knowledge`.

Everything is seeded, so re-running produces byte-identical CSVs.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ..config import SOURCES_DATA
from ..idmap import get_resolver
from ..knowledge import INTERVENTIONS, REFERENCES

SEED = 20240117

# Real HAGR / OpenGenes curated-database memberships, keyed by ortholog group.
# ortholog_group -> [(database, assertion)]
CURATED: dict[str, list[tuple[str, str]]] = {
    "CDKN2A": [("GenAge", "human ageing-associated gene"), ("CellAge", "senescence inducer")],
    "CDKN1A": [("CellAge", "senescence inducer"), ("GenAge", "ageing-associated")],
    "TP53": [("GenAge", "human ageing-associated gene"), ("CellAge", "senescence regulator")],
    "SIRT1": [("GenAge", "ageing-associated"), ("OpenGenes", "longevity-associated")],
    "FOXO3": [
        ("LongevityMap", "human longevity variant"),
        ("GenAge", "ageing-associated"),
        ("OpenGenes", "longevity-associated"),
    ],
    "TERT": [("GenAge", "ageing-associated"), ("LongevityMap", "human longevity variant")],
    "KL": [
        ("GenAge", "ageing-associated"),
        ("LongevityMap", "klotho longevity variant"),
        ("OpenGenes", "lifespan-extending overexpression"),
    ],
    "MTOR": [("GenAge", "ageing-associated"), ("OpenGenes", "lifespan-extending knockdown")],
    "IGF1": [("GenAge", "ageing-associated"), ("OpenGenes", "insulin/IGF-1 signalling")],
    "IGF1R": [
        ("GenAge", "ageing-associated"),
        ("OpenGenes", "lifespan-extending knockdown"),
        ("LongevityMap", "human longevity variant"),
    ],
    "GDF15": [("OpenGenes", "aging biomarker")],
    "LMNB1": [("CellAge", "senescence marker (down)")],
}

# Reference marker profile for the simulated example cohort, in conventional US
# clinical units. (base_young, base_old, worsening_loading). The loading sign
# encodes the direction a marker moves as health declines; a shared aging factor
# with age-growing variance drives all markers, reproducing the critical-slowing
# -down signature the resilience tool looks for.
_MARKERS = {
    #                 young   old    loading
    "albumin_gdl": (4.55, 3.95, -0.35),
    "creatinine_mgdl": (0.85, 1.12, 0.12),
    "glucose_mgdl": (90.0, 110.0, 9.0),
    "crp_mgl": (1.2, 3.4, 1.4),
    "lymphocyte_pct": (34.0, 23.0, -3.2),
    "mcv_fl": (89.0, 92.8, 1.6),
    "rdw_pct": (12.9, 14.7, 0.7),
    "alp_ul": (66.0, 84.0, 6.0),
    "wbc_1000ul": (5.9, 6.9, 0.7),
}


def build_curated(resolver) -> list[dict]:
    rows = []
    for group, flags in CURATED.items():
        for gene in resolver.genes_in_group(group):
            for db, assertion in flags:
                rows.append(
                    {
                        "gene_id": gene.canonical_id,
                        "database": db,
                        "assertion": assertion,
                        "url": f"https://genomics.senescence.info/genes/search?search={gene.symbol}",
                    }
                )
    return rows


def build_interventions(resolver) -> list[dict]:
    rows = []
    for i, iv in enumerate(INTERVENTIONS.values(), 1):
        linked: list[str] = []
        for group in iv.linked_groups:
            linked.extend(g.canonical_id for g in resolver.genes_in_group(group))
        primary_ref = REFERENCES[iv.reference_keys[0]] if iv.reference_keys else None
        rows.append(
            {
                "intervention_id": f"IV{i:03d}",
                "name": iv.name,
                "itype": iv.itype,
                "source": iv.source,
                "organism": iv.organism,
                "lifespan_effect_pct": iv.lifespan_effect_pct,
                "linked_gene_ids": "|".join(sorted(set(linked))),
                "url": primary_ref.url if primary_ref else "",
            }
        )
    return rows


def build_example_cohort(rng: np.random.Generator) -> list[dict]:
    """A SIMULATED example cohort carrying the nine real PhenoAge markers.

    Not real patient data. A shared latent 'aging' factor whose variance grows
    with age drives all markers, so biological-age acceleration and the
    resilience (critical-slowing-down) signature both emerge realistically while
    the cohort remains an honest, transparent worked example.
    """
    rows = []
    n = 720
    for i in range(n):
        age = float(rng.integers(20, 86))
        aging = (age - 20.0) / 65.0  # 0 (young) .. 1 (old)
        common = rng.normal(0.0, 0.35 + 1.6 * aging)  # variance grows with age
        row = {"subject_id": f"SIM{i:04d}", "age": age, "sex": rng.choice(["male", "female"])}
        for marker, (young, old, loading) in _MARKERS.items():
            trend = young + (old - young) * aging
            idio = rng.normal(0.0, abs(loading) * 0.6)
            val = trend + loading * common + idio
            # Keep values physiologically non-negative.
            row[marker] = round(max(val, 0.01), 4)
        rows.append(row)
    return rows


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def build_all(out_dir: Path | None = None) -> dict[str, int]:
    out_dir = out_dir or SOURCES_DATA
    rng = np.random.default_rng(SEED)
    resolver = get_resolver()

    curated_rows = build_curated(resolver)
    iv_rows = build_interventions(resolver)
    cohort_rows = build_example_cohort(rng)

    _write_csv(out_dir / "curated_knowledge.csv", curated_rows)
    _write_csv(out_dir / "interventions.csv", iv_rows)
    _write_csv(out_dir / "example_cohort_simulated.csv", cohort_rows)

    # Remove any stale fabricated files from earlier versions of the project.
    for stale in ("signatures.csv", "studies.csv", "clinical_nhanes_slice.csv"):
        p = out_dir / stale
        if p.exists():
            p.unlink()

    return {
        "curated": len(curated_rows),
        "interventions": len(iv_rows),
        "example_cohort_subjects": len(cohort_rows),
    }


if __name__ == "__main__":
    counts = build_all()
    for k, v in counts.items():
        print(f"{k}: {v}")
