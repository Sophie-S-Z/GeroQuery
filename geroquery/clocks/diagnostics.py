"""Clock applicability diagnostics — "will this clock be valid on my data?".

Published aging clocks are routinely misapplied: markers in the wrong units,
values outside physiological range, or missing inputs silently coerced. Clock
outputs also ship without uncertainty. This module checks an uploaded matrix
against the real PhenoAge input contract *before* trusting the number, and
attaches a bootstrap confidence interval — directly operationalising the
clock-reliability literature (Higgins-Chen et al., 2022; ComputAgeBench).

It only uses the data supplied; nothing is fabricated.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import phenoage

# Plausible per-sample physiological ranges, conventional US units.
# (soft_low, soft_high) flag out-of-typical values; wildly outside suggests a unit error.
REFERENCE_RANGES: dict[str, tuple[float, float]] = {
    "albumin_gdl": (3.0, 5.5),
    "creatinine_mgdl": (0.4, 1.5),
    "glucose_mgdl": (50.0, 250.0),
    "crp_mgl": (0.0, 50.0),
    "lymphocyte_pct": (5.0, 60.0),
    "mcv_fl": (70.0, 115.0),
    "rdw_pct": (10.0, 25.0),
    "alp_ul": (20.0, 300.0),
    "wbc_1000ul": (2.0, 30.0),
    "age": (18.0, 110.0),
}

# Median-based unit-mismatch heuristics: (predicate on median, message).
_UNIT_HINTS: dict[str, tuple] = {
    "albumin_gdl": (
        lambda m: m > 20,
        "median looks like g/L — PhenoAge expects g/dL (divide by 10).",
    ),
    "glucose_mgdl": (lambda m: m < 30, "median looks like mmol/L — PhenoAge expects mg/dL (×18)."),
    "creatinine_mgdl": (
        lambda m: m > 20,
        "median looks like µmol/L — PhenoAge expects mg/dL (÷88.4).",
    ),
    "crp_mgl": (
        lambda m: m > 60,
        "median very high — check units (expects mg/L high-sensitivity CRP).",
    ),
}


def phenoage_diagnostics(df: pd.DataFrame, n_boot: int = 500, seed: int = 0) -> dict:
    """Applicability report + bootstrap-CI'd PhenoAge for an uploaded matrix."""
    required = list(phenoage.REQUIRED_FEATURES)
    provided = list(df.columns)
    missing = [c for c in required if c not in provided]

    report: dict = {
        "clock_id": "phenoage",
        "n_samples": int(len(df)),
        "required_features": required,
        "missing_features": missing,
        "applicable": not missing,
        "warnings": [],
        "range_flags": {},
        "unit_warnings": [],
    }
    if missing:
        report["warnings"].append(
            f"Missing {len(missing)} required feature(s): {', '.join(missing)}. "
            "The clock cannot be applied until these are provided."
        )
        return report

    # Per-marker out-of-range and unit checks.
    for col in required:
        series = pd.to_numeric(df[col], errors="coerce")
        n_nan = int(series.isna().sum())
        lo, hi = REFERENCE_RANGES[col]
        out = int(((series < lo) | (series > hi)).sum())
        if n_nan:
            report["warnings"].append(f"{col}: {n_nan} non-numeric/missing value(s).")
        if out:
            report["range_flags"][col] = {
                "out_of_range": out,
                "range": [lo, hi],
                "observed_median": float(np.nanmedian(series)),
            }
        if col in _UNIT_HINTS:
            pred, msg = _UNIT_HINTS[col]
            med = float(np.nanmedian(series))
            if np.isfinite(med) and pred(med):
                report["unit_warnings"].append(f"{col}: {msg} (median={med:.2f})")

    if report["range_flags"]:
        n_cols = len(report["range_flags"])
        report["warnings"].append(
            f"{n_cols} marker(s) have values outside typical physiological range — "
            "verify units and outliers before trusting the estimate."
        )

    # PhenoAge with a bootstrap CI on the cohort mean (needs clean numeric input).
    try:
        clean = df[required].apply(pd.to_numeric, errors="coerce").dropna()
        if len(clean):
            ages = phenoage.phenotypic_age(clean)
            report["mean_phenoage"] = float(np.mean(ages))
            if len(ages) >= 10:
                rng = np.random.default_rng(seed)
                boot = ages[rng.integers(0, len(ages), (n_boot, len(ages)))].mean(axis=1)
                report["mean_phenoage_ci"] = [
                    float(np.percentile(boot, 2.5)),
                    float(np.percentile(boot, 97.5)),
                ]
    except Exception as exc:  # noqa: BLE001 — diagnostics must never crash
        report["warnings"].append(f"Could not compute PhenoAge: {exc}")

    if report["applicable"] and not report["warnings"] and not report["unit_warnings"]:
        report["warnings"].append("No applicability issues detected — inputs look well-formed.")
    return report
