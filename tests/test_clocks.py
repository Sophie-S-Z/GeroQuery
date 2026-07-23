"""M5 clocks — the real PhenoAge model: validation, application, age acceleration."""

import numpy as np
import pandas as pd
import pytest

from geroquery.clocks import ClockService
from geroquery.clocks.phenoage import REQUIRED_FEATURES
from geroquery.exceptions import ClockInputError, ClockNotFoundError


@pytest.fixture(scope="module")
def svc():
    return ClockService()


def test_every_clock_declares_outcome_and_features(svc):
    for c in svc.list_clocks():
        assert c.predicted_outcome in {"chronological_age", "mortality", "pace_of_aging"}
        assert c.training_population
        assert len(c.required_features) > 0


def test_phenoage_is_registered(svc):
    ids = {c.clock_id for c in svc.list_clocks()}
    assert "phenoage" in ids


def test_phenoage_tracks_age_and_reports_acceleration(svc, clinical_matrix):
    feats, true_age = clinical_matrix
    res = svc.apply_clock("phenoage", feats, chronological_age=true_age)
    pred = np.array(res.predictions)
    # Biological age should correlate strongly with chronological age.
    assert np.corrcoef(pred, true_age)[0, 1] > 0.6
    assert res.age_acceleration is not None
    assert res.mean_age_acceleration is not None
    # PhenoAge is built on a mortality model, so it also reports a risk.
    assert res.mortality_risk_10yr is not None
    assert all(0.0 <= r <= 1.0 for r in res.mortality_risk_10yr)


def test_phenoage_reads_age_from_matrix_when_not_passed(svc, clinical_matrix):
    feats, _ = clinical_matrix
    res = svc.apply_clock("phenoage", feats)  # age column present in matrix
    assert res.age_acceleration is not None


def test_unhealthy_profile_is_biologically_older(svc):
    healthy = {
        "albumin_gdl": 4.7,
        "creatinine_mgdl": 0.8,
        "glucose_mgdl": 85,
        "crp_mgl": 0.5,
        "lymphocyte_pct": 36,
        "mcv_fl": 88,
        "rdw_pct": 12.5,
        "alp_ul": 58,
        "wbc_1000ul": 5.2,
        "age": 50,
    }
    unhealthy = {
        "albumin_gdl": 3.6,
        "creatinine_mgdl": 1.4,
        "glucose_mgdl": 140,
        "crp_mgl": 8.0,
        "lymphocyte_pct": 16,
        "mcv_fl": 96,
        "rdw_pct": 16.0,
        "alp_ul": 130,
        "wbc_1000ul": 9.5,
        "age": 50,
    }
    df = pd.DataFrame([healthy, unhealthy])
    res = svc.apply_clock("phenoage", df, chronological_age=[50, 50])
    assert res.predictions[1] > res.predictions[0]  # same age, worse markers → older
    assert res.mortality_risk_10yr[1] > res.mortality_risk_10yr[0]


def test_missing_feature_raises_precise_error(svc, clinical_matrix):
    feats, _ = clinical_matrix
    with pytest.raises(ClockInputError) as exc:
        svc.apply_clock("phenoage", feats.drop(columns=["crp_mgl"]))
    assert "crp_mgl" in exc.value.detail["missing_features"]


def test_nan_input_rejected(svc, clinical_matrix):
    feats, _ = clinical_matrix
    bad = feats.copy()
    bad.iloc[0, 0] = np.nan
    with pytest.raises(ClockInputError):
        svc.apply_clock("phenoage", bad)


def test_compare_runs_multiple_clocks(svc, clinical_matrix):
    feats, true_age = clinical_matrix
    results = svc.compare_clocks(["phenoage"], feats, true_age)
    assert len(results) == 1
    assert results[0].clock_id == "phenoage"


def test_unknown_clock_raises(svc, clinical_matrix):
    feats, _ = clinical_matrix
    with pytest.raises(ClockNotFoundError):
        svc.apply_clock("no_such_clock", feats)


def test_empty_matrix_rejected(svc):
    with pytest.raises(ClockInputError):
        svc.apply_clock("phenoage", pd.DataFrame())


def test_required_features_are_the_nine_markers_plus_age():
    assert len(REQUIRED_FEATURES) == 10
    assert "age" in REQUIRED_FEATURES
