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
    """Every clock must state what it predicts, what it was trained on, and what
    it needs as input.

    Deliberately not a closed vocabulary. When biolearn is installed the registry
    holds 60+ real models predicting BMI, smoking status, cholesterol and cell
    proportions alongside age — an allow-list of {chronological_age, mortality,
    pace_of_aging} would force those into a bucket they do not belong in, which
    is the exact mislabelling this metadata exists to prevent.
    """
    for c in svc.list_clocks():
        assert c.predicted_outcome, c.clock_id
        assert c.predicted_outcome == c.predicted_outcome.strip().lower()
        assert c.training_population, c.clock_id
        assert c.units, c.clock_id
        # required_features is legitimately empty for the pyaging tier: its
        # model artifacts download lazily, so the feature list is not knowable
        # at registration time. Everything else must declare its inputs.
        if c.library != "pyaging":
            assert len(c.required_features) > 0, c.clock_id


def test_non_age_predictors_are_not_labelled_as_age_clocks(svc):
    """A predictor of smoking, BMI, cholesterol or sex is not an aging clock.
    Only runs where the library tier is loaded; the reference tier has none."""
    NOT_AGE = ("smoking", "bmi", "cholesterol", "sex", "bodyfat", "alcohol", "education")
    checked = 0
    for c in svc.list_clocks():
        if any(token in c.clock_id.lower() for token in NOT_AGE):
            checked += 1
            assert c.predicted_outcome != "chronological_age", c.clock_id
    # No assertion on `checked` — it is legitimately 0 without biolearn.


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
