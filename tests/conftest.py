"""Shared fixtures. A single freshly-built store/service is reused across the
suite; store-specific tests build isolated stores in tmp dirs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from geroquery.api.app import create_app
from geroquery.api.service import GeroService
from geroquery.clocks.phenoage import REQUIRED_FEATURES
from geroquery.store import GeroStore


@pytest.fixture(scope="session")
def store(tmp_path_factory) -> GeroStore:
    home = tmp_path_factory.mktemp("store_home")
    return GeroStore(data_home=home).build()


@pytest.fixture(scope="session")
def service(store) -> GeroService:
    return GeroService(store=store)


@pytest.fixture(scope="session")
def client(service) -> TestClient:
    return TestClient(create_app(service))


@pytest.fixture
def clinical_matrix():
    """A realistic cohort of the nine PhenoAge clinical markers plus age, where
    the markers worsen with age. Biological age from the real PhenoAge model must
    therefore track chronological age."""
    rng = np.random.default_rng(7)
    n = 80
    ages = rng.integers(25, 85, n).astype(float)
    aging = (ages - 25) / 60.0
    df = pd.DataFrame(
        {
            "albumin_gdl": 4.6 - 0.5 * aging + rng.normal(0, 0.1, n),
            "creatinine_mgdl": 0.8 + 0.3 * aging + rng.normal(0, 0.05, n),
            "glucose_mgdl": 88 + 20 * aging + rng.normal(0, 3, n),
            "crp_mgl": np.clip(0.8 + 2.5 * aging + rng.normal(0, 0.3, n), 0.05, None),
            "lymphocyte_pct": 35 - 10 * aging + rng.normal(0, 1, n),
            "mcv_fl": 89 + 3 * aging + rng.normal(0, 1, n),
            "rdw_pct": 12.8 + 1.6 * aging + rng.normal(0, 0.2, n),
            "alp_ul": 64 + 18 * aging + rng.normal(0, 4, n),
            "wbc_1000ul": 5.8 + 1.0 * aging + rng.normal(0, 0.4, n),
            "age": ages,
        }
    )
    assert list(df.columns) == list(REQUIRED_FEATURES)
    return df, ages
