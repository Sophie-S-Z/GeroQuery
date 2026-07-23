"""M7 api — endpoint behavior, envelopes, pagination, format switch, OpenAPI."""

import io

import pandas as pd


def test_healthz_and_version(client):
    assert client.get("/healthz").json()["status"] == "ok"
    v = client.get("/v1/version").json()
    assert "code_version" in v and "data_version" in v


def test_gene_signature_shape_and_cross_species(client):
    r = client.get("/v1/gene/CDKN2A/signature", params={"omic_layer": "transcriptome"})
    assert r.status_code == 200
    body = r.json()
    assert body["gene"]["symbol"] == "CDKN2A"
    species = {m["species"] for m in body["meta_signatures"]}
    assert {"human", "mouse"} <= species  # conservation view
    assert all(m["direction"] == "up" for m in body["meta_signatures"])


def test_gene_signature_pagination(client):
    full = client.get("/v1/gene/CDKN2A/signature").json()["signatures"]
    page = client.get("/v1/gene/CDKN2A/signature", params={"limit": 3, "offset": 2}).json()
    assert page["signatures"] == full[2:5]


def test_gene_signature_format_csv_and_parquet(client):
    csv = client.get("/v1/gene/CDKN2A/signature", params={"format": "csv"})
    assert csv.headers["content-type"].startswith("text/csv")
    assert "effect_size" in csv.text

    pq = client.get("/v1/gene/CDKN2A/signature", params={"format": "parquet"})
    df = pd.read_parquet(io.BytesIO(pq.content))
    assert "gene_id" in df.columns and len(df) > 0


def test_gene_card(client):
    card = client.get("/v1/gene/klotho/card").json()
    assert card["gene"]["symbol"] == "KL"
    assert len(card["curated_flags"]) >= 1


def test_geneset_signature(client):
    r = client.post(
        "/v1/geneset/signature",
        json={
            "genes": ["CDKN2A", "LMNB1", "SIRT1"],
            "species": "human",
            "omic_layer": "transcriptome",
        },
    ).json()
    assert set(r["resolved"]) == {"CDKN2A", "LMNB1", "SIRT1"}
    assert r["aggregate_pooled_effect"] is not None


def test_clocks_list_carry_outcome_metadata(client):
    clocks = client.get("/v1/clocks").json()["clocks"]
    assert all("predicted_outcome" in c for c in clocks)


def test_clock_apply_on_dataset(client):
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "clinical_phenoage_demo", "dataset_id": "clinical_nhanes_slice"},
    ).json()
    assert r["n_samples"] == 720
    assert r["predicted_outcome"] == "chronological_age"


def test_clock_apply_on_uploaded_matrix(client):
    records = [
        {
            "albumin": 4.2,
            "creatinine": 0.9,
            "glucose": 95,
            "crp": 1.2,
            "lymphocyte_pct": 30,
            "rdw": 13,
        },
        {
            "albumin": 3.8,
            "creatinine": 1.1,
            "glucose": 110,
            "crp": 3.0,
            "lymphocyte_pct": 22,
            "rdw": 15,
        },
    ]
    r = client.post(
        "/v1/clock/apply",
        json={
            "clock_id": "clinical_phenoage_demo",
            "matrix": {"records": records},
            "chronological_age": [40, 70],
        },
    )
    assert r.status_code == 200
    assert r.json()["n_samples"] == 2


def test_clock_apply_missing_feature_envelope(client):
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "clinical_phenoage_demo", "matrix": {"records": [{"albumin": 4.0}]}},
    )
    assert r.status_code == 422
    err = r.json()["error"]
    assert err["code"] == "clock_input_error"
    assert "missing_features" in err["detail"]


def test_clock_compare(client):
    r = client.post(
        "/v1/clock/compare",
        json={
            "clock_ids": ["clinical_phenoage_demo", "clinical_mortality_demo"],
            "dataset_id": "clinical_nhanes_slice",
        },
    ).json()
    assert len(r["results"]) == 2


def test_intervention_and_not_found(client):
    ok = client.get("/v1/intervention/rapamycin").json()
    assert ok["intervention"]["name"] == "rapamycin"
    missing = client.get("/v1/intervention/notadrug")
    assert missing.status_code == 404
    assert missing.json()["error"]["code"] == "intervention_not_found"


def test_resilience_csd_endpoint(client):
    r = client.post(
        "/v1/resilience/csd", json={"dataset_id": "clinical_nhanes_slice", "n_strata": 6}
    ).json()
    assert r["resilience_declines"] is True
    assert r["fallback_used"] is True
    assert r["assumptions"]


def test_resilience_recovery_endpoint(client):
    series = [3.0, 1.5, 0.8, 0.4, 0.2, 0.12, 0.07, 0.04]
    r = client.post("/v1/resilience/recovery", json={"series": series}).json()
    assert r["recovery_rate"] > 0


def test_studies_pagination_and_format(client):
    j = client.get("/v1/studies", params={"limit": 5}).json()
    assert j["n"] == 5
    csv = client.get("/v1/studies", params={"format": "csv", "limit": 5})
    assert csv.headers["content-type"].startswith("text/csv")


def test_sources_and_datasets(client):
    names = {s["name"] for s in client.get("/v1/sources").json()["sources"]}
    assert {"local-harmonized", "uk-biobank"} <= names
    ds = {d["dataset_id"] for d in client.get("/v1/datasets").json()["datasets"]}
    assert "clinical_nhanes_slice" in ds


def test_error_envelope_structure(client):
    r = client.get("/v1/gene/___nope___/card")
    assert r.status_code == 404
    err = r.json()["error"]
    assert set(err.keys()) == {"code", "message", "detail"}


def test_openapi_schema_valid(client):
    schema = client.get("/openapi.json").json()
    assert schema["info"]["title"] == "GeroQuery API"
    assert "/v1/gene/{gene_id}/signature" in schema["paths"]
