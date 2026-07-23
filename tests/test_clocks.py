"""M5 clocks — validation, application, age acceleration, metadata."""

import numpy as np
import pandas as pd
import pytest

from geroquery.clocks import ClockService
from geroquery.exceptions import ClockInputError, ClockNotFoundError


@pytest.fixture(scope="module")
def svc():
    return ClockService()


def test_every_clock_declares_outcome_and_features(svc):
    for c in svc.list_clocks():
        assert c.predicted_outcome in {"chronological_age", "mortality", "pace_of_aging"}
        assert c.training_population
        assert len(c.required_features) > 0


def test_clock_recovers_known_ages(svc, clinical_matrix):
    feats, true_age = clinical_matrix
    res = svc.apply_clock("clinical_phenoage_demo", feats, chronological_age=true_age)
    mae = np.mean(np.abs(np.array(res.predictions) - true_age))
    assert mae < 1.5  # within noise
    assert res.mean_age_acceleration == pytest.approx(0.0, abs=0.5)
    assert res.age_acceleration is not None


def test_missing_feature_raises_precise_error(svc, clinical_matrix):
    feats, _ = clinical_matrix
    with pytest.raises(ClockInputError) as exc:
        svc.apply_clock("clinical_phenoage_demo", feats.drop(columns=["crp"]))
    assert "crp" in exc.value.detail["missing_features"]


def test_nan_input_rejected(svc, clinical_matrix):
    feats, _ = clinical_matrix
    bad = feats.copy()
    bad.iloc[0, 0] = np.nan
    with pytest.raises(ClockInputError):
        svc.apply_clock("clinical_phenoage_demo", bad)


def test_mortality_clock_has_no_age_acceleration(svc, clinical_matrix):
    feats, true_age = clinical_matrix
    res = svc.apply_clock("clinical_mortality_demo", feats, chronological_age=true_age)
    assert res.predicted_outcome == "mortality"
    assert res.age_acceleration is None  # meaningless for a non-age clock


def test_compare_runs_multiple_clocks(svc, clinical_matrix):
    feats, true_age = clinical_matrix
    results = svc.compare_clocks(
        ["clinical_phenoage_demo", "clinical_mortality_demo"], feats, true_age
    )
    assert len(results) == 2
    assert {r.clock_id for r in results} == {"clinical_phenoage_demo", "clinical_mortality_demo"}


def test_unknown_clock_raises(svc, clinical_matrix):
    feats, _ = clinical_matrix
    with pytest.raises(ClockNotFoundError):
        svc.apply_clock("no_such_clock", feats)


def test_empty_matrix_rejected(svc):
    with pytest.raises(ClockInputError):
        svc.apply_clock("clinical_phenoage_demo", pd.DataFrame())
