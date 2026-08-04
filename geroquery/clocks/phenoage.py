"""The real Levine PhenoAge clinical aging clock.

PhenoAge is a published, peer-reviewed biological-age clock built from nine
routine clinical-chemistry markers plus chronological age (Levine et al., 2018,
*Aging*; methodology from Liu et al., 2018, *PLoS Medicine*, derived on NHANES).
Unlike the transparent teaching clocks it replaces, **the coefficients below are
the real published coefficients** — GeroQuery does not invent them.

The model works in two steps:

1. A Gompertz proportional-hazards mortality model combines the markers into a
   predicted 10-year mortality risk.
2. That risk is mapped back onto the age scale of the reference population, so
   the output ("phenotypic age") is expressed in years and is directly
   comparable to chronological age.

Because the same model yields both a mortality risk and an age, PhenoAge is a
clean, honest way to show that "how old you look biologically" and "your
short-term mortality risk" are two facets of one calibrated model.

Inputs are supplied in conventional US clinical units and converted internally
to the SI units the published coefficients require.

The model is ported unchanged from the coefficients as published. What is new
here is that it now has real data to run on: :mod:`geroquery.sources.nhanes`
carries all nine markers out of the already-pinned NHANES XPORT files, so
``phenoage`` is exercised against 4,894 real US adults rather than a fixture.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Required input columns, in conventional US clinical units.
REQUIRED_FEATURES: tuple[str, ...] = (
    "albumin_gdl",  # g/dL
    "creatinine_mgdl",  # mg/dL
    "glucose_mgdl",  # mg/dL
    "crp_mgl",  # mg/L (high-sensitivity CRP)
    "lymphocyte_pct",  # %
    "mcv_fl",  # fL
    "rdw_pct",  # %
    "alp_ul",  # U/L (alkaline phosphatase)
    "wbc_1000ul",  # 10^3 cells/uL
    "age",  # chronological age, years
)

# Published PhenoAge linear-predictor coefficients (SI units).
_INTERCEPT = -19.9067
_COEF = {
    "albumin_gL": -0.0336,
    "creatinine_umolL": 0.0095,
    "glucose_mmolL": 0.1953,
    "ln_crp_mgdl": 0.0954,
    "lymphocyte_pct": -0.0120,
    "mcv_fl": 0.0268,
    "rdw_pct": 0.3306,
    "alp_ul": 0.0019,
    "wbc_1000ul": 0.0554,
    "age": 0.0804,
}

# Gompertz mortality constants from the published derivation.
_GAMMA = 0.0076927
_MONTHS = 120.0  # 10-year horizon
# Age-scale mapping constants.
_A = 141.50225
_B = -0.00553
_C = 0.090165


# GeroQuery's clinical frame names its columns after the marker, not the unit,
# because the resilience module and the NHANES adapter share that vocabulary.
# Accepting both spellings here means neither side has to rename to talk to the
# other, and a caller cannot accidentally pair `albumin` (g/dL) with a
# coefficient expecting g/L.
COLUMN_ALIASES: dict[str, str] = {
    "albumin": "albumin_gdl",
    "creatinine": "creatinine_mgdl",
    "glucose": "glucose_mgdl",
    "crp": "crp_mgl",
    "rdw": "rdw_pct",
    "mcv": "mcv_fl",
    "wbc": "wbc_1000ul",
    "alp": "alp_ul",
}


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename GeroQuery marker columns to the names the coefficients expect.

    Columns already in PhenoAge's own naming are left alone, so a frame that
    mixes the two conventions still resolves.
    """
    rename = {
        source: target
        for source, target in COLUMN_ALIASES.items()
        if source in df.columns and target not in df.columns
    }
    return df.rename(columns=rename) if rename else df


def _linear_predictor(df: pd.DataFrame) -> np.ndarray:
    albumin_gL = df["albumin_gdl"].to_numpy(float) * 10.0
    creatinine_umolL = df["creatinine_mgdl"].to_numpy(float) * 88.4
    glucose_mmolL = df["glucose_mgdl"].to_numpy(float) / 18.0
    # CRP: mg/L -> mg/dL, then natural log (guard against non-positive values).
    crp_mgdl = np.clip(df["crp_mgl"].to_numpy(float) / 10.0, 1e-4, None)
    ln_crp = np.log(crp_mgdl)

    xb = (
        _INTERCEPT
        + _COEF["albumin_gL"] * albumin_gL
        + _COEF["creatinine_umolL"] * creatinine_umolL
        + _COEF["glucose_mmolL"] * glucose_mmolL
        + _COEF["ln_crp_mgdl"] * ln_crp
        + _COEF["lymphocyte_pct"] * df["lymphocyte_pct"].to_numpy(float)
        + _COEF["mcv_fl"] * df["mcv_fl"].to_numpy(float)
        + _COEF["rdw_pct"] * df["rdw_pct"].to_numpy(float)
        + _COEF["alp_ul"] * df["alp_ul"].to_numpy(float)
        + _COEF["wbc_1000ul"] * df["wbc_1000ul"].to_numpy(float)
        + _COEF["age"] * df["age"].to_numpy(float)
    )
    return xb


def mortality_risk_10yr(df: pd.DataFrame) -> np.ndarray:
    """Predicted 10-year mortality risk (0–1) from the PhenoAge mortality model."""
    xb = _linear_predictor(normalize_columns(df))
    return 1.0 - np.exp(-np.exp(xb) * (np.exp(_GAMMA * _MONTHS) - 1.0) / _GAMMA)


def phenotypic_age(df: pd.DataFrame) -> np.ndarray:
    """PhenoAge (biological age, in years)."""
    risk = mortality_risk_10yr(normalize_columns(df))
    # Numerically guard the double-log mapping.
    risk = np.clip(risk, 1e-8, 1 - 1e-8)
    return _A + np.log(_B * np.log(1.0 - risk)) / _C
