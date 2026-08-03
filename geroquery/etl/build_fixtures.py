"""Build the one deliberately synthetic table GeroQuery ships.

Run: ``python -m geroquery.etl.build_fixtures``

This module used to generate the gene signature table, the curated-knowledge
table, and the intervention table as well. All three are now real:

* signatures — :mod:`geroquery.etl.build_signatures`, from the checksum-pinned
  GEO DataSets aging panel.
* curated knowledge and interventions — the same module, from the HAGR releases
  (GenAge, CellAge, LongevityMap, DrugAge, GenDR).

What remains is a *method-validation fixture* with critical slowing down planted
in it. It is not a stand-in for missing data and never becomes one: it exists to
answer "does the estimator recover an effect that is known to be there?", which
no real dataset can answer, because in real data you do not know the answer.

It is written to its own dataset id (``clinical_synthetic_csd``), described with
a ``SYNTHETIC`` prefix, and a test asserts it never merges with the real NHANES
table. Output is seeded, so re-running produces a byte-identical CSV.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ..config import SOURCES_DATA

SEED = 20240117

SYNTHETIC_CLINICAL = "clinical_synthetic_csd.csv"

MARKERS = ("albumin", "creatinine", "glucose", "crp", "lymphocyte_pct", "rdw")

# Population means, on the same units as the NHANES columns they mirror.
BASELINE = {
    "albumin": 4.3,
    "creatinine": 0.9,
    "glucose": 95.0,
    "crp": 1.5,
    "lymphocyte_pct": 32.0,
    "rdw": 13.0,
}

# How strongly each marker loads on the shared latent factor, and in which
# direction. Albumin and lymphocyte percent fall as health declines; the rest rise.
LOADING = {
    "albumin": -0.4,
    "creatinine": 0.3,
    "glucose": 5.0,
    "crp": 0.8,
    "lymphocyte_pct": -2.5,
    "rdw": 0.6,
}

# The planted effect: the latent factor's SD grows linearly with normalized age.
# This is what makes both the variance and the cross-correlation signal rise —
# by construction, not by biology.
LATENT_SD_INTERCEPT = 0.4
LATENT_SD_SLOPE = 1.8

IDIOSYNCRATIC_SD = 0.5

N_SUBJECTS = 720
AGE_MIN, AGE_MAX = 20, 86


def build_clinical_synthetic(rng: np.random.Generator) -> tuple[tuple[str, ...], list[dict]]:
    """A synthetic clinical slice with critical slowing down planted in it.

    **This is a method-validation fixture, not data about humans.** A shared
    'aging' latent factor whose variance is explicitly engineered to grow with
    age drives every biomarker, so the health-state variance and the marker
    cross-correlation both rise with age *by construction*.

    The real clinical data lives in :mod:`geroquery.sources.nhanes`; on it, the
    variance signal replicates and the cross-correlation signal does not. See
    ``docs/RESULTS_NHANES_CSD.md``.
    """
    rows = []
    for i in range(N_SUBJECTS):
        age = float(rng.integers(AGE_MIN, AGE_MAX))
        aging = (age - AGE_MIN) / (AGE_MAX - 1 - AGE_MIN)
        common = rng.normal(0.0, LATENT_SD_INTERCEPT + LATENT_SD_SLOPE * aging)
        row: dict = {"subject_id": f"S{i:04d}", "age": age, "sex": rng.choice(["male", "female"])}
        for marker in MARKERS:
            idiosyncratic = rng.normal(0.0, IDIOSYNCRATIC_SD)
            row[marker] = round(BASELINE[marker] + LOADING[marker] * common + idiosyncratic, 4)
        rows.append(row)
    return MARKERS, rows


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_all(out_dir: Path | None = None) -> dict[str, int]:
    out_dir = out_dir or SOURCES_DATA
    rng = np.random.default_rng(SEED)
    _markers, rows = build_clinical_synthetic(rng)
    _write_csv(out_dir / SYNTHETIC_CLINICAL, rows)
    return {"synthetic_clinical_subjects": len(rows)}


if __name__ == "__main__":
    for key, value in build_all().items():
        print(f"{key}: {value}")
