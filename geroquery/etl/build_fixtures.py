"""Deterministically build the bundled harmonized data slice.

Run: ``python -m geroquery.etl.build_fixtures``

Everything is seeded, so re-running produces byte-identical CSVs — a property the
reproducibility test relies on. The *directions* of effect encode established
aging biology (p16/CDKN2A up, LMNB1 down, GDF15 up in plasma, klotho/IGF1 down);
magnitudes and per-study scatter are synthetic and clearly labelled as a
demonstration slice. The production ETL replaces the generators, not the schema.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from ..config import SOURCES_DATA
from ..idmap import get_resolver

SEED = 20240117

# ortholog_group -> per-omic "true" standardized age effect (+ up with age).
TRUE_EFFECTS: dict[str, dict[str, float]] = {
    "CDKN2A": {"transcriptome": 1.20, "methylome": 0.60},
    "CDKN1A": {"transcriptome": 0.90},
    "LMNB1": {"transcriptome": -1.00, "methylome": -0.40},
    "GDF15": {"transcriptome": 0.70, "proteome": 1.10},
    "SIRT1": {"transcriptome": -0.50, "methylome": -0.30},
    "TERT": {"transcriptome": -0.80},
    "KL": {"transcriptome": -0.90, "proteome": -0.60},
    "FOXO3": {"transcriptome": -0.20},
    "IGF1": {"transcriptome": -0.60, "proteome": -0.70},
    "TP53": {"transcriptome": 0.30},
    "MTOR": {"transcriptome": 0.25},
    "IGF1R": {"transcriptome": -0.30},
}

TISSUES_BY_OMIC = {
    "transcriptome": ["blood", "brain", "muscle"],
    "methylome": ["blood"],
    "proteome": ["plasma"],
}
# proteome/methylome are human-only in this slice; transcriptome covers both.
SPECIES_BY_OMIC = {
    "transcriptome": ["human", "mouse"],
    "methylome": ["human"],
    "proteome": ["human"],
}

CURATED = {
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

# name -> (type, source, organism, lifespan_effect_pct, [ortholog_groups])
INTERVENTIONS = {
    "rapamycin": ("drug", "ITP", "mouse", 14.0, ["MTOR"]),
    "caloric_restriction": ("dietary", "GenDR", "mouse", 30.0, ["IGF1", "MTOR", "SIRT1", "FOXO3"]),
    "metformin": ("drug", "DrugAge", "mouse", 6.0, ["MTOR", "GDF15"]),
    "resveratrol": ("drug", "DrugAge", "mouse", 4.0, ["SIRT1"]),
    "dasatinib_quercetin": ("drug", "DrugAge", "mouse", 9.0, ["CDKN2A", "CDKN1A"]),
    "nmn": ("drug", "DrugAge", "mouse", 5.0, ["SIRT1"]),
    "17a_estradiol": ("drug", "ITP", "mouse", 12.0, ["IGF1"]),
}


def _se_from_g(g: float, n1: int, n2: int) -> float:
    return math.sqrt((n1 + n2) / (n1 * n2) + (g * g) / (2.0 * (n1 + n2)))


def _p_from_z(g: float, se: float) -> float:
    from scipy import stats

    z = g / se if se > 0 else 0.0
    return float(2.0 * stats.norm.sf(abs(z)))


def build_signatures_and_studies(rng: np.random.Generator, resolver):
    sig_rows, study_rows = [], []
    study_counter = 90000
    for group, omics in TRUE_EFFECTS.items():
        records = resolver.genes_in_group(group)
        for omic, true_eff in omics.items():
            for species in SPECIES_BY_OMIC[omic]:
                gene = next((r for r in records if r.species == species), None)
                if gene is None:
                    continue
                mag = true_eff * (0.9 if species == "mouse" else 1.0)
                for tissue in TISSUES_BY_OMIC[omic]:
                    n_studies = 3 if omic == "transcriptome" else 2
                    for _ in range(n_studies):
                        n1 = int(rng.integers(8, 30))
                        n2 = int(rng.integers(8, 30))
                        g = float(mag + rng.normal(0, 0.15))
                        se = _se_from_g(g, n1, n2)
                        p = _p_from_z(g, se)
                        study_counter += 1
                        study_id = f"GEO:GSE{study_counter}"
                        sig_rows.append(
                            {
                                "gene_id": gene.canonical_id,
                                "study_id": study_id,
                                "omic_layer": omic,
                                "species": species,
                                "tissue": tissue,
                                "sex": rng.choice(["male", "female", "both"]),
                                "age_range": rng.choice(["young_vs_old", "20-80", "3m_vs_24m"]),
                                "effect_size": round(g, 4),
                                "direction": "up" if g >= 0 else "down",
                                "p_value": round(p, 6),
                                "q_value": round(min(1.0, p * 1.5), 6),
                                "standard_error": round(se, 4),
                                "source": "GEO",
                            }
                        )
                        study_rows.append(
                            {
                                "study_id": study_id,
                                "source": "GEO",
                                "omic_layer": omic,
                                "species": species,
                                "tissue": tissue,
                                "sample_size": n1 + n2,
                                "processing_method": "harmonized_demo_slice",
                                "license": "public-attribute",
                                "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE{study_counter}",
                                "version": "2024.1",
                            }
                        )
    return sig_rows, study_rows


def build_curated(resolver):
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


def build_interventions(resolver):
    rows = []
    for i, (name, (itype, source, organism, eff, groups)) in enumerate(INTERVENTIONS.items(), 1):
        linked = []
        for group in groups:
            linked.extend(g.canonical_id for g in resolver.genes_in_group(group))
        rows.append(
            {
                "intervention_id": f"IV{i:03d}",
                "name": name,
                "itype": itype,
                "source": source,
                "organism": organism,
                "lifespan_effect_pct": eff,
                "linked_gene_ids": "|".join(sorted(set(linked))),
                "url": f"https://genomics.senescence.info/drugs/search?search={name}",
            }
        )
    return rows


def build_clinical(rng: np.random.Generator):
    """A clinical/NHANES-style slice engineered with critical slowing down.

    A shared 'aging' latent factor whose variance grows with age drives all
    biomarkers, so both the variance of the health state and the cross-correlation
    among biomarkers rise with age — the cross-sectional CSD signature the
    resilience module reproduces.
    """
    markers = ["albumin", "creatinine", "glucose", "crp", "lymphocyte_pct", "rdw"]
    base = {
        "albumin": 4.3,
        "creatinine": 0.9,
        "glucose": 95.0,
        "crp": 1.5,
        "lymphocyte_pct": 32.0,
        "rdw": 13.0,
    }
    loading = {
        "albumin": -0.4,
        "creatinine": 0.3,
        "glucose": 5.0,
        "crp": 0.8,
        "lymphocyte_pct": -2.5,
        "rdw": 0.6,
    }
    rows = []
    n = 720
    for i in range(n):
        age = float(rng.integers(20, 86))
        aging = (age - 20.0) / 65.0
        common = rng.normal(0.0, 0.4 + 1.8 * aging)  # variance grows with age
        row = {"subject_id": f"S{i:04d}", "age": age, "sex": rng.choice(["male", "female"])}
        for m in markers:
            idio = rng.normal(0.0, 0.5)  # fixed idiosyncratic noise
            row[m] = round(base[m] + loading[m] * common + idio, 4)
        rows.append(row)
    return markers, rows


def _write_csv(path: Path, rows: list[dict]):
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

    sig_rows, study_rows = build_signatures_and_studies(rng, resolver)
    curated_rows = build_curated(resolver)
    iv_rows = build_interventions(resolver)
    _markers, clinical_rows = build_clinical(rng)

    _write_csv(out_dir / "signatures.csv", sig_rows)
    _write_csv(out_dir / "studies.csv", study_rows)
    _write_csv(out_dir / "curated_knowledge.csv", curated_rows)
    _write_csv(out_dir / "interventions.csv", iv_rows)
    _write_csv(out_dir / "clinical_nhanes_slice.csv", clinical_rows)

    return {
        "signatures": len(sig_rows),
        "studies": len(study_rows),
        "curated": len(curated_rows),
        "interventions": len(iv_rows),
        "clinical_subjects": len(clinical_rows),
    }


if __name__ == "__main__":
    counts = build_all()
    for k, v in counts.items():
        print(f"{k}: {v}")
