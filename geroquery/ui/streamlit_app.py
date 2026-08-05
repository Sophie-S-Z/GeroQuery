"""GeroQuery — research dashboard.

A gene-first reader for the biology of aging, over measured data.

The central design decision: **this surface shows intervals, not verdicts.** The
version of this dashboard it replaces was built around a hand-written table that
asserted a direction and a confidence level per gene, so it could render an arrow
and a label. The data underneath is now 485,905 random-effects estimates pooled
from real GEO contrasts, and the honest primitive for that is a forest plot with
a confidence interval — including when the interval crosses zero and the answer
is *we cannot tell*.

Everything here is measured or ingested from a checksum-pinned upstream. The one
exception is the method-validation fixture, badged wherever it appears.

Launch:  python -m streamlit run geroquery/ui/streamlit_app.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from geroquery.api.service import GeroService
from geroquery.exceptions import GeroQueryError
from geroquery.ui.theme import PALETTES, css, plotly_layout

st.set_page_config(
    page_title="GeroQuery",
    page_icon="◐",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_resource
def get_service() -> GeroService:
    return GeroService()


svc = get_service()

# ---- theme ---------------------------------------------------------------

if "theme" not in st.session_state:
    st.session_state.theme = "dark"

with st.sidebar:
    st.markdown("## Appearance")
    dark = st.toggle(
        "Dark",
        value=st.session_state.theme == "dark",
        help="Default is dark: this page is read for long stretches, usually indoors.",
    )
    st.session_state.theme = "dark" if dark else "light"

P = PALETTES[st.session_state.theme]
st.markdown(css(P), unsafe_allow_html=True)


# ---- small helpers -------------------------------------------------------


def panel(body: str) -> None:
    st.markdown(f"<div class='gq-panel'>{body}</div>", unsafe_allow_html=True)


def badge(kind: str, text: str) -> str:
    return f"<span class='gq-badge gq-{kind}'>{text}</span>"


def read_interval(meta: dict) -> tuple[str, str, str]:
    """Turn a pooled estimate into a claim, a colour, and the reason for it.

    This is the only place the dashboard decides what a number means, and it
    decides it from the interval rather than the point estimate.
    """
    lo, hi, g = meta["ci_low"], meta["ci_high"], meta["pooled_effect"]
    if lo > 0:
        return ("rises with age", P.up, "the interval lies entirely above zero")
    if hi < 0:
        return ("falls with age", P.down, "the interval lies entirely below zero")
    direction = "higher" if g > 0 else "lower"
    return (
        "no detectable change",
        P.null,
        f"the interval spans zero, so a {direction} value in older samples is not "
        f"distinguishable from chance in this panel",
    )


def forest(rows: list[dict], pooled: dict, title: str) -> go.Figure:
    """Per-study estimates with their intervals, and the pooled result beneath.

    A forest plot rather than a bar chart because the width of each interval is
    the point: a study with three samples per group and one with fifteen should
    not look equally certain.
    """
    fig = go.Figure()
    labels, xs, los, his = [], [], [], []
    for r in rows:
        se = r.get("standard_error") or 0.0
        labels.append(f"{r['study_id'].replace('GEO:', '')} · {r.get('tissue') or '—'}")
        xs.append(r["effect_size"])
        los.append(r["effect_size"] - 1.96 * se)
        his.append(r["effect_size"] + 1.96 * se)

    for i, (lo, hi) in enumerate(zip(los, his, strict=True)):
        fig.add_shape(type="line", x0=lo, x1=hi, y0=i, y1=i, line=dict(color=P.line, width=1.4))
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=list(range(len(xs))),
            mode="markers",
            marker=dict(size=7, color=P.muted),
            hovertemplate="%{text}<br>g = %{x:.3f}<extra></extra>",
            text=labels,
        )
    )

    y0 = -1.6
    _claim, colour, _why = read_interval(pooled)
    fig.add_shape(
        type="line",
        x0=pooled["ci_low"],
        x1=pooled["ci_high"],
        y0=y0,
        y1=y0,
        line=dict(color=colour, width=3),
    )
    fig.add_trace(
        go.Scatter(
            x=[pooled["pooled_effect"]],
            y=[y0],
            mode="markers",
            marker=dict(size=13, color=colour, symbol="diamond"),
            hovertemplate=(
                f"pooled g = {pooled['pooled_effect']:.3f}<br>"
                f"95% CI [{pooled['ci_low']:.2f}, {pooled['ci_high']:.2f}]<extra></extra>"
            ),
        )
    )
    fig.add_vline(x=0, line=dict(color=P.muted, width=1, dash="dot"))

    layout = plotly_layout(P, height=max(230, 30 * len(xs) + 130))
    layout["title"] = dict(
        text=title, font=dict(family="Courier Prime, monospace", size=11, color=P.muted), x=0
    )
    layout["yaxis"] = dict(
        tickmode="array",
        tickvals=[*range(len(labels)), y0],
        ticktext=[*labels, "<b>POOLED</b>"],
        gridcolor="rgba(0,0,0,0)",
        zeroline=False,
        linecolor="rgba(0,0,0,0)",
        tickfont=dict(family="Courier Prime, monospace", size=9.5, color=P.muted),
    )
    layout["xaxis"] = dict(
        title=dict(
            text="Hedges' g   (positive = higher in older samples)",
            font=dict(size=10.5, color=P.muted),
        ),
        gridcolor=P.plot_grid,
        zeroline=False,
        linecolor=P.line,
    )
    fig.update_layout(**layout)
    return fig


# ---- masthead ------------------------------------------------------------


@st.cache_data(show_spinner=False)
def corpus_summary() -> dict:
    """Counted from the store, never typed into the page.

    A masthead that hardcodes its own numbers goes stale the first time the panel
    changes, and a stale count on a page whose whole argument is provenance is a
    self-inflicted wound.
    """
    signatures = svc.store.query_signatures()
    studies = svc.studies()
    return {
        "rows": len(signatures),
        "genes": len({s.gene_id for s in signatures}),
        "contrasts": len(studies),
        "datasets": len({s["study_id"].split(":")[1] for s in studies if ":" in s["study_id"]}),
        "series": len({s.get("series_id") for s in studies if s.get("series_id")}),
    }


version = svc.version()
corpus = corpus_summary()
st.markdown(
    f"""<div class='gq-mast'>
<h1>GeroQuery</h1>
<p class='lede'>Ask what a gene actually does with age, and get an interval wide
enough to say when the answer is <em>we cannot tell</em>.</p>
<div class='meta'>{corpus["rows"]:,} effect sizes &middot; {corpus["genes"]:,} genes
&middot; {corpus["contrasts"]} contrasts from {corpus["datasets"]} GEO DataSets
&middot; NHANES 2017&ndash;2018 &middot; data {version["data_version"]}</div>
</div>""",
    unsafe_allow_html=True,
)

tab_gene, tab_clock, tab_res, tab_about = st.tabs(
    ["Gene explorer", "Aging clock", "Resilience", "About & data"]
)


# ---- gene explorer -------------------------------------------------------


def render_gene() -> None:
    curated = svc.list_curated_genes(limit=400)
    symbols = sorted({g["symbol"] for g in curated})

    left, right = st.columns([3, 2], gap="large")
    with left:
        typed = st.text_input(
            "Gene",
            value=st.session_state.get("gene_q", "CDKN1A"),
            placeholder="symbol, alias, Ensembl or Entrez id",
            label_visibility="collapsed",
        )
    with right:
        picked = st.selectbox(
            "Browse curated", ["Browse the curated set…", *symbols], label_visibility="collapsed"
        )

    query = typed
    if picked != "Browse the curated set…" and picked != st.session_state.get("last_pick"):
        st.session_state.last_pick = picked
        st.session_state.gene_q = picked
        query = picked

    if not query.strip():
        panel(
            "<div class='gq-note'>Type a gene symbol, or pick one of the "
            f"{len(curated):,} genes carrying curated aging evidence.</div>"
        )
        return

    try:
        report = svc.gene_report(query.strip())
    except GeroQueryError as exc:
        panel(
            f"<h4>Not resolved</h4><div class='gq-note'>{exc.message}</div>"
            "<div class='gq-note' style='margin-top:.6rem'>Try a HGNC symbol "
            "(<span class='gq-measure'>CDKN1A</span>), an alias "
            "(<span class='gq-measure'>p21</span>), or an Ensembl id.</div>"
        )
        return

    gene = report["gene"]
    metas = report["meta_signatures"]
    sigs = report["signatures"]

    st.markdown(
        f"<div class='gq-panel'>"
        f"<div style='display:flex;justify-content:space-between;align-items:baseline;"
        f"flex-wrap:wrap;gap:.8rem'>"
        f"<div style='font-size:1.9rem;font-weight:500'>{gene['symbol']}"
        f"<span style='color:{P.muted};font-size:1.05rem;font-weight:400'> &middot; "
        f"{gene.get('name') or 'no description'}</span></div>"
        f"<div>{badge('tag', gene['species'])}</div></div>"
        f"<div class='gq-note gq-measure' style='margin-top:.7rem;font-size:.8rem'>"
        f"{gene['canonical_id']}"
        f"{' · entrez ' + gene['entrez'] if gene.get('entrez') else ''}"
        f"{' · uniprot ' + gene['uniprot'] if gene.get('uniprot') else ''}</div>"
        f"<div class='gq-note' style='margin-top:.45rem;font-size:.85rem'>aliases: "
        f"<span class='gq-measure'>{', '.join(gene.get('aliases', [])[:8]) or '—'}</span>"
        f"</div></div>",
        unsafe_allow_html=True,
    )

    if not metas:
        panel(
            "<h4>No pooled estimate</h4>"
            "<div class='gq-note'>This gene resolves, but no contrast in the "
            "31-dataset GEO panel carries a probe for it. That is a coverage gap in a "
            "microarray-era panel, not a statement about the gene.</div>"
        )
    else:
        for meta in metas:
            claim, colour, why = read_interval(meta)
            rows = [
                s
                for s in sigs
                if s["species"] == meta["species"] and s["omic_layer"] == meta["omic_layer"]
            ]
            st.markdown(
                f"<div class='gq-panel'>"
                f"<h4>{meta['species']} &middot; {meta['omic_layer']} &middot; "
                f"{meta['n_studies']} contrasts</h4>"
                f"<div class='gq-verdict'>"
                f"<div class='claim' style='color:{colour}'>{gene['symbol']} {claim}</div>"
                f"<div class='because'>Pooled Hedges' <em>g</em> "
                f"<span class='gq-interval'>{meta['pooled_effect']:+.3f}"
                f"<span class='ci'>95% CI [{meta['ci_low']:+.2f}, {meta['ci_high']:+.2f}]"
                f"</span></span> &mdash; {why}.</div></div></div>",
                unsafe_allow_html=True,
            )
            c1, c2, c3 = st.columns(3)
            c1.metric("Contrasts", meta["n_studies"])
            c2.metric("Heterogeneity I²", f"{meta['heterogeneity_i2']:.0f}%")
            c3.metric("p (pooled)", f"{meta['p_value']:.2e}")
            if rows:
                st.plotly_chart(
                    forest(rows, meta, f"{meta['species']} — per-contrast estimates"),
                    width="stretch",
                    config={"displayModeBar": False},
                )

        st.markdown(
            "<div class='gq-caveat'>Each horizontal line is one published contrast "
            "and its 95% interval; the diamond is the random-effects pool. Groups here "
            "are 3–15 samples, so the median interval is 0.84 wide — this panel detects "
            "large, consistent effects and cannot rule out moderate ones.</div>",
            unsafe_allow_html=True,
        )

    flags = report.get("curated_flags", [])
    if flags:
        rows_html = "".join(
            f"<tr><td>{f['database']}</td><td>{f['assertion']}</td>"
            f"<td><a href='{f.get('url', '')}' target='_blank'>source</a></td></tr>"
            for f in flags[:14]
        )
        panel(
            f"<h4>Curated evidence &middot; {len(flags)} assertions</h4>"
            f"<table class='gq-table'><thead><tr><th>Database</th><th>Assertion</th>"
            f"<th></th></tr></thead><tbody>{rows_html}</tbody></table>"
        )

    interventions = report.get("interventions", [])
    if interventions:
        rows_html = "".join(
            f"<tr><td>{i['name']}</td><td>{i.get('organism') or '—'}</td>"
            f"<td class='num'>{i.get('lifespan_effect_pct', '—')}</td>"
            f"<td>{i['source']}</td></tr>"
            for i in interventions[:10]
        )
        panel(
            f"<h4>Linked interventions</h4>"
            f"<table class='gq-table'><thead><tr><th>Intervention</th><th>Organism</th>"
            f"<th class='num'>Median lifespan Δ%</th><th>Source</th></tr></thead>"
            f"<tbody>{rows_html}</tbody></table>"
            "<div class='gq-note' style='margin-top:.8rem'>Only GenDR records "
            "gene-level dependence. DrugAge lists no gene targets, so no drug is linked "
            "to a gene here — that edge was never in the data.</div>"
        )


# ---- aging clock ---------------------------------------------------------

PHENOAGE_INPUTS = "albumin, creatinine, glucose, crp, lymphocyte_pct, mcv, rdw, alp, wbc"


def render_clock() -> None:
    st.markdown(
        "<div class='gq-panel'><h4>PhenoAge &middot; Levine et al. 2018</h4>"
        "<p>Nine routine blood markers and chronological age, combined by a published "
        "Gompertz mortality model and mapped back onto the age scale. The same fit "
        "yields a biological age <em>and</em> a calibrated 10-year mortality risk, which "
        "is the cleanest way to show those are two readouts of one model rather than "
        "interchangeable outputs.</p></div>",
        unsafe_allow_html=True,
    )

    source = st.radio(
        "Cohort", ["Real NHANES cohort", "Upload a CSV"], horizontal=True, key="clock_src"
    )

    frame = None
    if source == "Real NHANES cohort":
        frame = svc.store.get_dataset("clinical_nhanes_phenoage")
        st.markdown(
            f"{badge('real', 'measured')} <span class='gq-note'>{len(frame):,} US adults, "
            "NHANES 2017–2018, carrying all nine markers PhenoAge requires. Fetched from "
            "CDC and verified against a pinned SHA-256.</span>",
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            f"<div class='gq-note'>A CSV with <span class='gq-measure'>age</span> plus "
            f"<span class='gq-measure'>{PHENOAGE_INPUTS}</span>, in conventional US "
            "clinical units.</div>",
            unsafe_allow_html=True,
        )
        upload = st.file_uploader("Cohort CSV", type=["csv"], label_visibility="collapsed")
        if upload is not None:
            frame = pd.read_csv(upload)

    if frame is None or not st.button("Run PhenoAge"):
        return

    try:
        result = svc.apply_clock("phenoage", frame, chronological_age=frame["age"].tolist())
    except (GeroQueryError, KeyError) as exc:
        message = exc.message if isinstance(exc, GeroQueryError) else f"missing column {exc}"
        panel(f"<h4>Could not run</h4><div class='gq-note'>{message}</div>")
        return

    predicted = np.asarray(result["predictions"], dtype=float)
    actual = frame["age"].to_numpy(dtype=float)
    accel = np.asarray(result["age_acceleration"], dtype=float)

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", f"{result['n_samples']:,}")
    c2.metric("r with chronological age", f"{np.corrcoef(predicted, actual)[0, 1]:.3f}")
    c3.metric("MAE", f"{np.abs(accel).mean():.2f} y")
    c4.metric("Age acceleration SD", f"{accel.std():.2f} y")

    scatter = go.Figure()
    lo = float(min(actual.min(), predicted.min()))
    hi = float(max(actual.max(), predicted.max()))
    scatter.add_trace(
        go.Scatter(
            x=[lo, hi],
            y=[lo, hi],
            mode="lines",
            line=dict(color=P.muted, width=1, dash="dot"),
            hoverinfo="skip",
        )
    )
    scatter.add_trace(
        go.Scatter(
            x=actual,
            y=predicted,
            mode="markers",
            marker=dict(size=4.5, color=P.accent, opacity=0.45),
            hovertemplate="chronological %{x:.0f}<br>PhenoAge %{y:.1f}<extra></extra>",
        )
    )
    layout = plotly_layout(P, height=380)
    layout["xaxis"] = dict(
        title=dict(text="chronological age", font=dict(size=11)),
        gridcolor=P.plot_grid,
        linecolor=P.line,
    )
    layout["yaxis"] = dict(
        title=dict(text="PhenoAge", font=dict(size=11)), gridcolor=P.plot_grid, linecolor=P.line
    )
    scatter.update_layout(**layout)
    st.plotly_chart(scatter, width="stretch", config={"displayModeBar": False})

    st.markdown(
        "<div class='gq-caveat'>Points above the dotted line are biologically older than "
        "their years. The spread around it is the signal PhenoAge exists to measure — a "
        "clock sitting exactly on the line would be reporting the birth certificate back "
        "to you.</div>",
        unsafe_allow_html=True,
    )

    risk = result.get("mortality_risk_10yr")
    if risk:
        hist = go.Figure(
            go.Histogram(
                x=np.asarray(risk, dtype=float) * 100,
                nbinsx=48,
                marker=dict(color=P.accent, opacity=0.75),
            )
        )
        layout = plotly_layout(P, height=260)
        layout["xaxis"] = dict(
            title=dict(text="predicted 10-year mortality risk (%)", font=dict(size=11)),
            gridcolor=P.plot_grid,
            linecolor=P.line,
        )
        layout["yaxis"] = dict(
            title=dict(text="participants", font=dict(size=11)),
            gridcolor=P.plot_grid,
            linecolor=P.line,
        )
        hist.update_layout(**layout)
        panel("<h4>Ten-year mortality risk, from the same fit</h4>")
        st.plotly_chart(hist, width="stretch", config={"displayModeBar": False})


# ---- resilience ----------------------------------------------------------

COHORTS = {
    "Real NHANES cohort": ("clinical_nhanes_slice", "real"),
    "Method-validation fixture": ("clinical_synthetic_csd", "sim"),
}


def render_resilience() -> None:
    st.markdown(
        "<div class='gq-panel'><h4>Critical slowing down</h4>"
        "<p>A system near breakdown gets slower to recover from small knocks: its "
        "signals wobble more, and they start moving together. GeroQuery looks for both "
        "across age strata — rising <em>variance</em> and rising "
        "<em>cross-correlation</em>. Claiming a resilience decline requires both, "
        "because a single rising indicator fires on noise about half the time.</p></div>",
        unsafe_allow_html=True,
    )

    choice = st.radio("Cohort", [*COHORTS, "Upload a CSV"], horizontal=True, key="res_src")
    n_strata = st.slider("Age strata", 3, 12, 6)

    dataset_id, frame = None, None
    if choice in COHORTS:
        dataset_id, kind = COHORTS[choice]
        frame = svc.store.get_dataset(dataset_id)
        if kind == "real":
            st.markdown(
                f"{badge('real', 'measured')} <span class='gq-note'>{len(frame):,} US "
                "adults, NHANES 2017–2018, six routine blood markers.</span>",
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                f"{badge('sim', 'constructed')} <span class='gq-note'>720 subjects with "
                "critical slowing down <em>planted by construction</em>. Not evidence "
                "about people — it exists so the estimator can be checked against a known "
                "answer. Run both and compare.</span>",
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            "<div class='gq-note'>A cross-sectional CSV with an "
            "<span class='gq-measure'>age</span> column and numeric biomarker columns."
            "</div>",
            unsafe_allow_html=True,
        )
        upload = st.file_uploader("Cohort CSV", type=["csv"], label_visibility="collapsed")
        if upload is not None:
            frame = pd.read_csv(upload)

    if frame is None or not st.button("Compute early-warning signals"):
        return

    try:
        if dataset_id:
            res = svc.resilience_csd(dataset_id=dataset_id, n_strata=n_strata)
        else:
            cols = [c for c in frame.columns if c not in ("subject_id", "age", "sex")]
            res = svc.resilience_csd(data=frame, biomarker_cols=cols, n_strata=n_strata)
    except GeroQueryError as exc:
        panel(f"<h4>Could not compute</h4><div class='gq-note'>{exc.message}</div>")
        return

    declines = res["resilience_declines"]
    colour = P.up if declines else P.real
    headline = "Resilience declines with age" if declines else "No resilience decline detected"
    st.markdown(
        f"<div class='gq-panel'><div class='gq-verdict'>"
        f"<div class='claim' style='color:{colour}'>{headline}</div>"
        f"<div class='because'>{res['verdict']}</div></div></div>",
        unsafe_allow_html=True,
    )

    mid = res["strata_midpoints"]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=mid,
            y=res["variance"],
            mode="lines+markers",
            name="variance",
            line=dict(color=P.up, width=2),
            marker=dict(size=7),
            hovertemplate="age %{x:.0f}<br>variance %{y:.4f}<extra></extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=mid,
            y=res["cross_correlation"],
            mode="lines+markers",
            name="cross-correlation",
            line=dict(color=P.down, width=2, dash="dash"),
            marker=dict(size=7),
            hovertemplate="age %{x:.0f}<br>cross-corr %{y:.4f}<extra></extra>",
        )
    )
    layout = plotly_layout(P, height=340)
    layout["showlegend"] = True
    layout["legend"] = dict(
        orientation="h", y=1.12, x=0, font=dict(family="Courier Prime, monospace", size=10)
    )
    layout["xaxis"] = dict(
        title=dict(text="age stratum midpoint", font=dict(size=11)),
        gridcolor=P.plot_grid,
        linecolor=P.line,
    )
    layout["yaxis"] = dict(
        title=dict(text="indicator", font=dict(size=11)), gridcolor=P.plot_grid, linecolor=P.line
    )
    fig.update_layout(**layout)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

    rows_html = "".join(
        f"<tr><td>{name}</td><td class='num'>{res[f'{key}_trend_slope']:+.5f}</td>"
        f"<td class='num'>{res[f'{key}_trend_tau']:+.3f}</td>"
        f"<td class='num'>{res[f'{key}_trend_p']:.3f}</td></tr>"
        for name, key in (("variance", "variance"), ("cross-correlation", "crosscorr"))
    )
    panel(
        f"<h4>Trend across strata</h4>"
        f"<table class='gq-table'><thead><tr><th>Indicator</th><th class='num'>slope</th>"
        f"<th class='num'>Kendall τ</th><th class='num'>p</th></tr></thead>"
        f"<tbody>{rows_html}</tbody></table>"
        f"<div class='gq-note' style='margin-top:.9rem'>n = {res['n_samples']:,} &middot; "
        f"{res['n_bootstrap']} bootstrap resamples &middot; age-detrended within stratum: "
        f"{str(res['detrended']).lower()}</div>"
    )

    if res.get("fallback_used"):
        st.markdown(
            "<div class='gq-caveat'>These are age-stratified approximations, not "
            "within-person trajectories. NHANES is cross-sectional, so no relaxation time "
            "is observed — which is precisely what would distinguish critical slowing down "
            "from ordinary accumulating heterogeneity between people.</div>",
            unsafe_allow_html=True,
        )


# ---- about ---------------------------------------------------------------


def render_about() -> None:
    st.markdown(
        "<div class='gq-panel'><h4>What is measured</h4>"
        "<p>Every number on this page comes from a checksum-pinned upstream. Nothing is "
        "generated to fill a gap.</p></div>",
        unsafe_allow_html=True,
    )

    sources = svc.sources()
    rows_html = "".join(
        f"<tr><td>{s.get('source_name', '—')}</td>"
        f"<td>{'federated' if s.get('federated') else 'cached'}</td>"
        f"<td>{s.get('license') or '—'}</td></tr>"
        for s in sources
    )
    panel(
        f"<h4>Adapters &middot; {len(sources)}</h4>"
        f"<table class='gq-table'><thead><tr><th>Source</th><th>Mode</th><th>Licence</th>"
        f"</tr></thead><tbody>{rows_html}</tbody></table>"
    )

    studies = svc.studies()
    series = {s.get("series_id") for s in studies if s.get("series_id")}
    panel(
        f"<h4>Evidence base</h4>"
        f"<table class='gq-table'><tbody>"
        f"<tr><td>GEO contrasts</td><td class='num'>{len(studies)}</td></tr>"
        f"<tr><td>Independent GEO Series</td><td class='num'>{len(series)}</td></tr>"
        f"<tr><td>Curated genes</td><td class='num'>{len(svc.list_curated_genes()):,}</td></tr>"
        f"<tr><td>Clocks registered</td><td class='num'>{len(svc.list_clocks())}</td></tr>"
        f"</tbody></table>"
        f"<div class='gq-note' style='margin-top:.9rem'>{len(studies)} contrasts come from "
        f"{len(series)} Series: GEO splits some experiments across two array halves, so a "
        "few share subjects. The pooling assumes independence and does not yet account "
        "for that.</div>"
    )

    st.markdown(
        f"<div class='gq-panel'><h4>The one constructed table</h4>"
        f"<p>{badge('sim', 'constructed')} <span style='margin-left:.5rem'>"
        "<span class='gq-measure'>clinical_synthetic_csd</span> — 720 subjects with "
        "critical slowing down planted by construction. It is not a stand-in for missing "
        "data: it is the only way to ask whether the estimator recovers an effect that is "
        "<em>known</em> to be there, which no real dataset can answer. It has its own "
        "dataset id and a test asserting it never merges with the NHANES table.</span>"
        "</p></div>",
        unsafe_allow_html=True,
    )

    references = svc.references()
    rows_html = "".join(
        f"<tr><td>{r['authors']} ({r['year']})</td><td>{r['title']}</td>"
        f"<td><a href='{r['url']}' target='_blank'>"
        f"{'PMID ' + r['pmid'] if r.get('pmid') else 'search'}</a></td></tr>"
        for r in references
    )
    panel(
        f"<h4>References &middot; {len(references)}</h4>"
        f"<table class='gq-table'><tbody>{rows_html}</tbody></table>"
        "<div class='gq-note' style='margin-top:.9rem'>Every PMID here is checked against "
        "PubMed by a live test that compares the recorded title with the one the API "
        "returns. Seven were wrong when that check was first run.</div>"
    )


with tab_gene:
    render_gene()
with tab_clock:
    render_clock()
with tab_res:
    render_resilience()
with tab_about:
    render_about()
