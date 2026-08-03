"""M8 ui — Streamlit dashboard.

Runs the GeroQuery service in-process (so the demo is a single one-click deploy
on Hugging Face Spaces / Streamlit Community Cloud with no separate backend),
while still exercising exactly the same service layer the REST API exposes.

Launch:  streamlit run geroquery/ui/streamlit_app.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from geroquery.api.service import GeroService

st.set_page_config(page_title="GeroQuery", page_icon="🧬", layout="wide")


@st.cache_resource
def get_service() -> GeroService:
    return GeroService()


svc = get_service()

st.title("🧬 GeroQuery")
st.caption(
    "Harmonized, multi-omic, cross-species, clock- and resilience-aware aging data — "
    f"data version `{svc.version()['data_version']}`"
)

tab_gene, tab_clock, tab_resilience, tab_about = st.tabs(
    ["Gene explorer", "Aging clocks", "Resilience", "About"]
)

# ---- Gene explorer -------------------------------------------------------
with tab_gene:
    col_q, col_s = st.columns([3, 1])
    query = col_q.text_input("Gene (symbol / Ensembl / Entrez / alias)", value="CDKN2A")
    species = col_s.selectbox("Species", ["both", "human", "mouse"], index=0)
    species_arg = None if species == "both" else species

    if not query.strip():
        st.info("Type a gene to begin — try CDKN2A (p16), LMNB1, KLOTHO, GDF15.")
    else:
        try:
            card = svc.gene_card(query, species_arg)
        except Exception as exc:  # noqa: BLE001 - surface any resolution error to the UI
            st.error(f"Could not resolve “{query}”: {exc}")
            st.stop()

        g = card["gene"]
        st.subheader(f"{g['symbol']} — {g.get('name', '')}")
        st.write(
            f"Canonical `{g['canonical_id']}` · Entrez `{g.get('entrez', '—')}` · "
            f"UniProt `{g.get('uniprot', '—')}` · species **{g['species']}**"
        )

        metas = pd.DataFrame(card["meta_signatures"])
        sigs = pd.DataFrame(card["signatures"])

        if metas.empty:
            st.warning("No GEO contrast in this panel covers this gene.")
        else:
            left, right = st.columns(2)
            with left:
                st.markdown("**Pooled aging effect (random-effects meta-analysis)**")
                st.dataframe(
                    metas[
                        [
                            "omic_layer",
                            "species",
                            "pooled_effect",
                            "ci_low",
                            "ci_high",
                            "heterogeneity_i2",
                            "n_studies",
                            "direction",
                        ]
                    ],
                    use_container_width=True,
                    hide_index=True,
                )
            with right:
                st.markdown("**Effect by omic layer & species**")
                fig = px.bar(
                    metas,
                    x="omic_layer",
                    y="pooled_effect",
                    color="species",
                    barmode="group",
                    error_y=metas["ci_high"] - metas["pooled_effect"],
                    labels={"pooled_effect": "pooled effect (Hedges g)"},
                )
                fig.add_hline(y=0, line_dash="dot", opacity=0.4)
                st.plotly_chart(fig, use_container_width=True)

            if not sigs.empty and {"effect_size", "p_value"} <= set(sigs.columns):
                st.markdown("**Per-study volcano**")
                sigs = sigs.copy()
                sigs["neg_log10_p"] = -np.log10(sigs["p_value"].clip(lower=1e-300))
                vol = px.scatter(
                    sigs,
                    x="effect_size",
                    y="neg_log10_p",
                    color="omic_layer",
                    symbol="species",
                    hover_data=["study_id", "tissue"],
                    labels={"effect_size": "effect size (Hedges g)", "neg_log10_p": "-log10 p"},
                )
                vol.add_vline(x=0, line_dash="dot", opacity=0.4)
                st.plotly_chart(vol, use_container_width=True)

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Curated knowledge**")
            flags = pd.DataFrame(card["curated_flags"])
            st.dataframe(
                flags if not flags.empty else pd.DataFrame({"info": ["none in slice"]}),
                use_container_width=True,
                hide_index=True,
            )
        with c2:
            st.markdown("**Linked interventions**")
            ivs = pd.DataFrame(card["interventions"])
            if ivs.empty:
                st.write(
                    "No linked interventions. DrugAge records no gene targets, so only "
                    "GenDR dietary-restriction genes carry intervention links."
                )
            else:
                st.dataframe(
                    ivs[["name", "itype", "source", "lifespan_effect_pct"]],
                    use_container_width=True,
                    hide_index=True,
                )

        if not sigs.empty:
            st.download_button(
                "⬇ Download signatures (CSV)",
                sigs.to_csv(index=False),
                file_name=f"{g['symbol']}_signatures.csv",
                mime="text/csv",
            )

# ---- Aging clocks --------------------------------------------------------
with tab_clock:
    clocks = pd.DataFrame(svc.list_clocks())
    st.markdown(
        "**Registered clocks** — note the predicted-outcome column: "
        "chronological-age accuracy does *not* imply mortality prediction."
    )
    st.dataframe(
        clocks[["clock_id", "predicted_outcome", "input_type", "training_population"]],
        use_container_width=True,
        hide_index=True,
    )

    clock_id = st.selectbox("Clock", clocks["clock_id"].tolist())
    dataset = st.selectbox("Dataset", [d["dataset_id"] for d in svc.datasets()])
    if st.button("Apply clock"):
        try:
            df = svc.store.get_dataset(dataset)
            chrono = df["age"].tolist() if "age" in df.columns else None
            result = svc.apply_clock(clock_id, df, chrono)
            st.success(
                f"Applied {clock_id} to {result['n_samples']} samples "
                f"({result['predicted_outcome']})."
            )
            out = pd.DataFrame({"predicted": result["predictions"]})
            if chrono is not None and result.get("age_acceleration"):
                out["chronological_age"] = chrono
                out["age_acceleration"] = result["age_acceleration"]
                fig = px.scatter(
                    out,
                    x="chronological_age",
                    y="predicted",
                    labels={"predicted": f"predicted ({result['units']})"},
                )
                lo, hi = min(chrono), max(chrono)
                fig.add_trace(
                    go.Scatter(
                        x=[lo, hi], y=[lo, hi], mode="lines", line=dict(dash="dot"), name="identity"
                    )
                )
                st.plotly_chart(fig, use_container_width=True)
                st.metric("Mean age acceleration (yrs)", round(result["mean_age_acceleration"], 2))
            st.dataframe(out.head(50), use_container_width=True, hide_index=True)
        except Exception as exc:  # noqa: BLE001
            st.error(str(exc))

# ---- Resilience ----------------------------------------------------------
with tab_resilience:
    st.markdown(
        "**Critical slowing down** — early-warning signals of resilience loss. "
        "Rising variance *and* rising cross-correlation with age indicate the system "
        "is approaching a tipping point."
    )
    dataset = st.selectbox("Dataset", [d["dataset_id"] for d in svc.datasets()], key="res_ds")
    n_strata = st.slider("Age strata", 3, 12, 6)
    if st.button("Compute resilience"):
        try:
            res = svc.resilience_csd(dataset_id=dataset, n_strata=n_strata)
            verdict = "declines with age ✓" if res["resilience_declines"] else "no clear decline"
            st.success(f"Resilience {verdict}. Method: `{res['method']}`")
            frame = pd.DataFrame(
                {
                    "age": res["strata_midpoints"],
                    "variance": res["variance"],
                    "cross_correlation": res["cross_correlation"],
                }
            )
            f1 = px.line(frame, x="age", y="variance", markers=True, title="Variance vs age")
            f2 = px.line(
                frame,
                x="age",
                y="cross_correlation",
                markers=True,
                title="Cross-correlation vs age",
            )
            st.plotly_chart(f1, use_container_width=True)
            st.plotly_chart(f2, use_container_width=True)
            with st.expander("Assumptions & limits"):
                for a in res["assumptions"]:
                    st.write("• " + a)
        except Exception as exc:  # noqa: BLE001
            st.error(str(exc))

with tab_about:
    st.markdown("""
GeroQuery unifies transcriptome + methylome + proteome + curated knowledge +
interventions behind one harmonized, cross-species API — plus a dynamical-systems
**resilience** layer no other aging portal exposes.

This dashboard calls the same service layer as the REST API (`/v1`). See the
README for architecture, the eight deep modules, and how to run the API and tests.
        """)
