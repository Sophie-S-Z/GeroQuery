"""M7 api — endpoint behavior, envelopes, format switch, OpenAPI."""

import io

import pandas as pd

_CLINICAL_ROW = {
    "albumin_gdl": 4.4,
    "creatinine_mgdl": 0.9,
    "glucose_mgdl": 95,
    "crp_mgl": 1.2,
    "lymphocyte_pct": 30,
    "mcv_fl": 90,
    "rdw_pct": 13.2,
    "alp_ul": 72,
    "wbc_1000ul": 6.2,
    "age": 55,
}


def test_healthz_and_version(client):
    assert client.get("/healthz").json()["status"] == "ok"
    v = client.get("/v1/version").json()
    assert "code_version" in v and "data_version" in v


def test_gene_report_shape_and_cross_species(client):
    body = client.get("/v1/gene/CDKN2A/report").json()
    assert body["gene"]["symbol"] == "CDKN2A"
    assert body["knowledge"]["direction_with_age"] == "up"
    species = {o["species"] for o in body["orthologs"]}
    assert {"human", "mouse"} <= species  # conservation view
    assert body["references"]


def test_gene_card_alias(client):
    card = client.get("/v1/gene/klotho/card").json()
    assert card["gene"]["symbol"] == "KL"
    assert len(card["curated_flags"]) >= 1
    assert card["knowledge"]["direction_with_age"] == "down"


def test_genes_list(client):
    genes = client.get("/v1/genes").json()["genes"]
    symbols = {g["symbol"] for g in genes}
    assert {"CDKN2A", "LMNB1", "GDF15", "FOXO3"} <= symbols
    assert all("direction_with_age" in g for g in genes)


def test_geneset_summary(client):
    r = client.post(
        "/v1/geneset/summary",
        json={"genes": ["CDKN2A", "LMNB1", "KL"], "species": "human"},
    ).json()
    assert set(r["resolved"]) == {"CDKN2A", "LMNB1", "KL"}
    assert r["direction_counts"]["up"] >= 1
    assert r["direction_counts"]["down"] >= 1


def test_clocks_list_carry_outcome_metadata(client):
    clocks = client.get("/v1/clocks").json()["clocks"]
    assert clocks and all("predicted_outcome" in c for c in clocks)
    assert any(c["clock_id"] == "phenoage" for c in clocks)


def test_clock_apply_on_dataset(client):
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "phenoage", "dataset_id": "example_cohort_simulated"},
    ).json()
    assert r["n_samples"] == 720
    assert r["predicted_outcome"] == "chronological_age"
    assert "mortality_risk_10yr" in r


def test_clock_apply_on_uploaded_matrix(client):
    records = [_CLINICAL_ROW, {**_CLINICAL_ROW, "age": 70, "crp_mgl": 4.0, "glucose_mgdl": 120}]
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "phenoage", "matrix": {"records": records}},
    )
    assert r.status_code == 200
    assert r.json()["n_samples"] == 2


def test_clock_apply_reports_mean_acceleration_ci(client):
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "phenoage", "dataset_id": "example_cohort_simulated"},
    ).json()
    ci = r["mean_age_acceleration_ci"]
    assert len(ci) == 2 and ci[0] <= r["mean_age_acceleration"] <= ci[1]


def test_clock_diagnostics_on_dataset(client):
    r = client.post(
        "/v1/clock/diagnostics",
        json={"dataset_id": "example_cohort_simulated"},
    ).json()
    assert r["clock_id"] == "phenoage"
    assert r["applicable"] is True
    assert r["missing_features"] == []
    assert len(r["mean_phenoage_ci"]) == 2


def test_clock_diagnostics_flags_missing_and_units(client):
    # Missing features -> not applicable.
    miss = client.post(
        "/v1/clock/diagnostics",
        json={"matrix": {"records": [{"albumin_gdl": 4.0}]}},
    ).json()
    assert miss["applicable"] is False and miss["missing_features"]
    # Albumin supplied in g/L (×10) -> unit warning.
    bad_units = {**_CLINICAL_ROW, "albumin_gdl": 44.0}
    ru = client.post("/v1/clock/diagnostics", json={"matrix": {"records": [bad_units]}}).json()
    assert any("g/L" in w for w in ru["unit_warnings"])


def test_clock_apply_missing_feature_envelope(client):
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "phenoage", "matrix": {"records": [{"albumin_gdl": 4.0}]}},
    )
    assert r.status_code == 422
    err = r.json()["error"]
    assert err["code"] == "clock_input_error"
    assert "missing_features" in err["detail"]


def test_intervention_and_not_found(client):
    ok = client.get("/v1/intervention/rapamycin").json()
    assert ok["intervention"]["name"] == "rapamycin"
    assert ok["intervention"]["references"]
    missing = client.get("/v1/intervention/notadrug")
    assert missing.status_code == 404
    assert missing.json()["error"]["code"] == "intervention_not_found"


def test_resilience_csd_endpoint(client):
    r = client.post(
        "/v1/resilience/csd",
        json={"dataset_id": "example_cohort_simulated", "n_strata": 6},
    ).json()
    assert r["resilience_declines"] is True
    assert r["fallback_used"] is True
    assert r["assumptions"]


def test_resilience_recovery_endpoint(client):
    series = [3.0, 1.5, 0.8, 0.4, 0.2, 0.12, 0.07, 0.04]
    r = client.post("/v1/resilience/recovery", json={"series": series}).json()
    assert r["recovery_rate"] > 0


def test_references_pagination_and_format(client):
    j = client.get("/v1/references", params={"limit": 5}).json()
    assert j["n"] == 5
    csv = client.get("/v1/references", params={"format": "csv", "limit": 5})
    assert csv.headers["content-type"].startswith("text/csv")
    df = pd.read_csv(io.StringIO(csv.text))
    assert "url" in df.columns


def test_sources_and_datasets(client):
    names = {s["name"] for s in client.get("/v1/sources").json()["sources"]}
    assert {"local-curated", "uk-biobank"} <= names
    ds = {d["dataset_id"] for d in client.get("/v1/datasets").json()["datasets"]}
    assert "example_cohort_simulated" in ds


def test_error_envelope_structure(client):
    r = client.get("/v1/gene/___nope___/card")
    assert r.status_code == 404
    err = r.json()["error"]
    assert set(err.keys()) == {"code", "message", "detail"}


def test_openapi_schema_valid(client):
    schema = client.get("/openapi.json").json()
    assert schema["info"]["title"] == "GeroQuery API"
    assert "/v1/gene/{gene_id}/report" in schema["paths"]
