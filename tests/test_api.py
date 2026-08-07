"""M7 api — endpoint behavior, envelopes, pagination, format switch, OpenAPI."""

import io

import pandas as pd
import pytest

from geroquery.exceptions import GeroQueryError


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
    for meta in body["meta_signatures"]:
        assert meta["ci_low"] <= meta["pooled_effect"] <= meta["ci_high"]
        assert meta["n_studies"] >= 2
        assert 0.0 <= meta["heterogeneity_i2"] <= 100.0


def test_p16_does_not_replicate_across_the_real_panel(client):
    """Pins the actual finding rather than the expected one.

    CDKN2A/p16 is the most-cited transcriptional marker of aging, and across
    these 32 real GEO contrasts its pooled effect is indistinguishable from
    zero. The synthetic slice this replaced had p16 planted at g=+1.20, so the
    old version of this test asserted every direction was "up" and passed. See
    docs/RESULTS_GEO_SIGNATURES.md for why arrays are a poor instrument for this
    particular transcript.
    """
    body = client.get("/v1/gene/CDKN2A/signature").json()
    human = [m for m in body["meta_signatures"] if m["species"] == "human"]
    assert human, "CDKN2A should still be on the human panel"
    pooled = human[0]
    assert (
        pooled["ci_low"] < 0 < pooled["ci_high"]
    ), "CDKN2A now has a CI excluding zero; the panel or the pipeline changed."


def test_p21_does_replicate_across_the_real_panel(client):
    """The counterweight to the test above: the pipeline is capable of finding a
    real effect, so a null result is a result and not a broken estimator."""
    body = client.get("/v1/gene/CDKN1A/signature").json()
    human = [m for m in body["meta_signatures"] if m["species"] == "human"]
    assert human
    pooled = human[0]
    assert pooled["pooled_effect"] > 0.5
    assert pooled["ci_low"] > 0


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


def test_geneset_reports_why_a_gene_was_unresolved(client):
    """An unknown gene is reported with a reason, not silently dropped."""
    r = client.post(
        "/v1/geneset/signature",
        json={"genes": ["CDKN2A", "___nope___"], "species": "human"},
    ).json()
    assert r["resolved"] == ["CDKN2A"]
    assert r["unresolved"] == ["___nope___"]
    assert r["unresolved_detail"]["___nope___"]


def test_geneset_does_not_swallow_unexpected_resolver_failures(service, monkeypatch):
    """Only gene/species lookup failures count as 'unresolved'.

    A bare `except Exception` here reported a dead id-mapping backend as "these
    gene names are unknown" — the caller got a plausible-looking empty result and
    no indication anything was broken.
    """

    def _boom(self, query, species=None):
        raise RuntimeError("id-mapping backend is down")

    monkeypatch.setattr("geroquery.idmap.resolver.GeneResolver.resolve_gene", _boom)
    with pytest.raises(RuntimeError, match="id-mapping backend is down"):
        service.geneset_signature(["CDKN2A"], species="human")


def test_clocks_list_carry_outcome_metadata(client):
    clocks = client.get("/v1/clocks").json()["clocks"]
    assert all("predicted_outcome" in c for c in clocks)


def test_clock_apply_on_dataset(client):
    r = client.post(
        "/v1/clock/apply",
        json={"clock_id": "clinical_phenoage_demo", "dataset_id": "clinical_nhanes_slice"},
    ).json()
    assert r["n_samples"] == 600  # committed NHANES sample; see the `store` fixture
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


def test_intervention_returns_every_organism_mammals_first(client):
    """DrugAge holds one record per organism. Returning only the first meant
    answering "how much does rapamycin extend lifespan?" with a nematode."""
    ok = client.get("/v1/intervention/rapamycin").json()
    assert ok["name"].lower() == "rapamycin"
    organisms = [i["organism"] for i in ok["interventions"]]
    assert "Mus musculus" in organisms and "Caenorhabditis elegans" in organisms
    assert organisms.index("Mus musculus") < organisms.index("Caenorhabditis elegans")
    mouse = next(i for i in ok["interventions"] if i["organism"] == "Mus musculus")
    assert mouse["lifespan_effect_pct"] == 13.0  # median of DrugAge's significant ITP results


def test_unknown_intervention_is_a_404(client):
    missing = client.get("/v1/intervention/notadrug")
    assert missing.status_code == 404
    assert missing.json()["error"]["code"] == "intervention_not_found"


def test_resilience_csd_endpoint(client):
    r = client.post(
        "/v1/resilience/csd", json={"dataset_id": "clinical_nhanes_slice", "n_strata": 6}
    ).json()
    assert r["fallback_used"] is True  # cross-sectional strata, not a time series
    assert r["assumptions"]
    assert r["n_samples"] == 600
    for block in ("variance_evidence", "crosscorr_evidence"):
        assert set(r[block]) >= {"slope", "ci_low", "ci_high", "supported", "verdict"}


def test_csd_excludes_survey_weight_from_biomarkers(client):
    """WTMEC2YR is a design variable. Auto-inferring biomarkers must not pick it
    up, or the sampling design ends up inside the measured health state."""
    r = client.post(
        "/v1/resilience/csd", json={"dataset_id": "clinical_nhanes_slice", "n_strata": 6}
    ).json()
    assert "survey_weight" not in r["biomarker_cols"]
    assert set(r["biomarker_cols"]) == {
        "albumin",
        "creatinine",
        "glucose",
        "crp",
        "lymphocyte_pct",
        "rdw",
    }


def test_csd_recovers_planted_effect_in_synthetic_fixture(client):
    """Method validation: on data with critical slowing down built in, the
    estimator must find it. This is the only claim the synthetic slice supports."""
    r = client.post(
        "/v1/resilience/csd", json={"dataset_id": "clinical_synthetic_csd", "n_strata": 6}
    ).json()
    assert r["variance_evidence"]["supported"] is True
    assert r["crosscorr_evidence"]["supported"] is True
    assert r["resilience_declines"] is True


def test_csd_does_not_claim_evidence_on_underpowered_real_sample(client):
    """The 600-row NHANES sample is underpowered: the variance slope is positive
    but its bootstrap CI straddles zero, so the gate must decline to claim.

    This is the discriminating test. Before the CSD fix, a positive slope alone
    was enough to return `resilience_declines: True`, so this dataset and the
    planted-effect fixture above would have been indistinguishable. Note this is
    absence of evidence at n=600, not evidence of absence — on the full 4,895-case
    cohort the variance signal *is* supported. See docs/RESULTS_NHANES_CSD.md.
    """
    r = client.post(
        "/v1/resilience/csd", json={"dataset_id": "clinical_nhanes_slice", "n_strata": 6}
    ).json()
    assert r["variance_evidence"]["slope"] > 0
    assert r["variance_evidence"]["ci_low"] < 0 < r["variance_evidence"]["ci_high"]
    assert r["variance_evidence"]["supported"] is False
    assert r["resilience_declines"] is False


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


# ---- one verdict, three surfaces (bug #20 / #22) ---------------------------


def test_every_api_verdict_is_readable_off_its_own_shipped_interval(service):
    """The invariant the API used to break for four genes.

    `MetaSignature.verdict` came from `PooledEffect.verdict` judged at full
    precision, while `ci_low`/`ci_high` shipped rounded to 4 dp and without the
    -0.0 collapse. The result was rows like `[-1.195, -0.0]` labelled
    "decreases", contradicting both the published site and their own interval.
    A caller can only ever see the shipped numbers, so those are what the
    verdict has to be derived from.
    """
    seen = 0
    for gene in ("CDKN2A", "TP53", "IL6", "SIRT1", "FOXO3"):
        try:
            payload = service.gene_signature(gene)
        except GeroQueryError:
            continue
        for meta in payload["meta_signatures"]:
            seen += 1
            low, high = meta["ci_low"], meta["ci_high"]
            expected = "increases" if low > 0 else "decreases" if high < 0 else "no_evidence"
            assert (
                meta["verdict"] == expected
            ), f"{gene}: verdict {meta['verdict']!r} contradicts [{low}, {high}]"
            # -0.0 prints as "-0.000" and is `not < 0` in JavaScript.
            assert str(low) != "-0.0" and str(high) != "-0.0"
    assert seen, "fixture resolved no pooled signatures — the assertion never ran"


def test_the_api_and_the_static_exporter_agree_gene_for_gene(service):
    """Same estimator, same rounding, same verdict rule — so same answer.

    These are two independent code paths over the same store, and they
    disagreed for four genes until the rounding and the verdict rule were
    unified in `harmonize.meta`.
    """
    from geroquery.harmonize.meta import report_bound

    for gene in ("CDKN2A", "TP53", "IL6"):
        try:
            payload = service.gene_signature(gene)
        except GeroQueryError:
            continue
        for meta in payload["meta_signatures"]:
            # What the exporter writes for the same pooled estimate.
            assert meta["ci_low"] == report_bound(meta["ci_low"])
            assert meta["ci_high"] == report_bound(meta["ci_high"])
