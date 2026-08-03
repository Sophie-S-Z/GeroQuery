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
from geroquery.clocks.registry import CLINICAL_FEATURES, get_registry
from geroquery.store import GeroStore


@pytest.fixture(scope="session")
def store(tmp_path_factory) -> GeroStore:
    """Built in sample mode on purpose.

    Both ``prefer_full_*=False`` flags pin the store to the committed offline
    slices: the 600-row NHANES sample and the curated-gene signature slice.
    Without them the suite would assert against 4,895 clinical rows and the full
    GEO panel on a machine that has run ``make data``, and against the slices in
    CI — so the same test would pass or fail depending on the developer's
    download cache.
    """
    home = tmp_path_factory.mktemp("store_home")
    return GeroStore(data_home=home).build(prefer_full_clinical=False, prefer_full_signatures=False)


@pytest.fixture(scope="session")
def service(store) -> GeroService:
    return GeroService(store=store)


@pytest.fixture(scope="session")
def client(service) -> TestClient:
    return TestClient(create_app(service))


@pytest.fixture
def clinical_matrix():
    """A clinical feature matrix generated from the reference clock's own model,
    so predictions must recover the planted ages within noise tolerance."""
    clock = get_registry().get("clinical_phenoage_demo")
    rng = np.random.default_rng(7)
    n = 80
    feats = pd.DataFrame({f: rng.normal(1.0, 0.4, n) for f in CLINICAL_FEATURES})
    w = np.array([clock.weights[f] for f in CLINICAL_FEATURES])
    true_age = clock.intercept + feats.to_numpy() @ w + rng.normal(0, 1.0, n)
    return feats, true_age
