"""GeroQuery — Streamlit dashboard.

A gene-first explorer for the biology of aging, plus a real biological-age clock
and a dynamical-systems resilience tool. Runs the GeroQuery service in-process,
so the whole app is a single one-click deploy with no separate backend.

Every number shown here is measured or ingested from a checksum-pinned upstream:
gene evidence from the GEO DataSets aging panel, curated flags and interventions
from the HAGR releases, the clinical cohorts from NHANES 2017-2018. The one
exception is the method-validation fixture, labelled wherever it appears, which
exists so an estimator can be checked against a known answer.

Launch:  python -m streamlit run geroquery/ui/streamlit_app.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from geroquery.api.service import GeroService
from geroquery.clocks.phenoage import REQUIRED_FEATURES

st.set_page_config(
    page_title="GeroQuery · the aging-gene search engine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --------------------------------------------------------------------------- #
# Palette & theme
# --------------------------------------------------------------------------- #

INK = "#1f2430"
MUTED = "#5b6472"
LINE = "#e6e8ef"
ACCENT = "#4f46e5"
UP = "#e0554f"  # increases with age
DOWN = "#2f7bd6"  # decreases with age
CONTEXT = "#8a94a6"  # context-dependent
CONF = {"robust": "#0f9d8f", "established": "#4f46e5", "emerging": "#94a3b8"}

DIR_META = {
    "up": ("↑", UP, "Increases with age"),
    "down": ("↓", DOWN, "Decreases with age"),
    "context-dependent": ("↔", CONTEXT, "Context-dependent"),
    "unknown": ("·", MUTED, "Not curated"),
}


def _plotly_theme(fig: go.Figure, height: int = 340, legend: bool = True) -> go.Figure:
    fig.update_layout(
        template="simple_white",
        font=dict(family="Inter, ui-sans-serif, system-ui, sans-serif", size=13, color=INK),
        margin=dict(l=10, r=16, t=32, b=10),
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        hoverlabel=dict(font_size=12, font_family="Inter, sans-serif"),
        showlegend=legend,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, font=dict(size=11)),
    )
    fig.update_xaxes(
        gridcolor=LINE, zeroline=False, linecolor=LINE, ticks="outside", tickcolor=LINE
    )
    fig.update_yaxes(
        gridcolor=LINE, zeroline=False, linecolor=LINE, ticks="outside", tickcolor=LINE
    )
    return fig


CSS = """
<style>
:root { --ink:#1f2430; --muted:#5b6472; --line:#e6e8ef; --accent:#4f46e5; --soft:#f6f7fb; }
html, body, [class*="css"] { font-family: 'Inter', ui-sans-serif, system-ui, -apple-system, sans-serif; }
.block-container { padding-top: 1.4rem; max-width: 1180px; }
h1,h2,h3,h4 { color: var(--ink); letter-spacing: -0.01em; }
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }

/* Hero */
.gq-hero { background: linear-gradient(120deg,#eef0ff 0%, #f6f7fb 55%, #ffffff 100%);
  border: 1px solid var(--line); border-radius: 18px; padding: 22px 26px; margin-bottom: 4px; }
.gq-hero h1 { font-size: 1.9rem; margin: 0 0 2px 0; }
.gq-hero p { color: var(--muted); margin: 0; font-size: 0.98rem; }
.gq-brandmark { font-size: 1.15rem; font-weight: 700; color: var(--accent); }

/* Cards */
.gq-card { background:#fff; border:1px solid var(--line); border-radius:14px; padding:18px 20px;
  margin-bottom:14px; box-shadow: 0 1px 2px rgba(20,24,40,0.03); }
.gq-card h4 { margin:0 0 10px 0; font-size:0.82rem; text-transform:uppercase; letter-spacing:0.06em;
  color: var(--muted); font-weight:700; }

/* Verdict block */
.gq-verdict { display:flex; align-items:center; gap:16px; }
.gq-glyph { font-size:2.6rem; line-height:1; font-weight:700; width:56px; text-align:center; }
.gq-verdict-main { flex:1; }
.gq-verdict-main .lead { font-size:1.05rem; color:var(--ink); margin:2px 0 0 0; }

/* Pills / badges */
.gq-pill { display:inline-block; padding:3px 11px; border-radius:999px; font-size:0.76rem;
  font-weight:600; margin:2px 6px 2px 0; border:1px solid transparent; }
.gq-chip { display:inline-block; padding:5px 12px; border-radius:10px; font-size:0.82rem;
  background:var(--soft); border:1px solid var(--line); color:var(--ink); margin:3px 6px 3px 0; }
.gq-db { display:inline-block; padding:6px 12px; border-radius:10px; font-size:0.82rem; font-weight:600;
  background:#eef0ff; border:1px solid #dfe2ff; color:#3a37b3; margin:3px 6px 3px 0; }

/* Analysis callout */
.gq-analysis { background:#fbfbfe; border-left:3px solid var(--accent); border-radius:0 10px 10px 0;
  padding:14px 18px; color:var(--ink); font-size:0.98rem; line-height:1.55; }

/* Evidence rows */
.gq-ev { border:1px solid var(--line); border-radius:12px; padding:12px 16px; margin-bottom:10px; background:#fff; }
.gq-ev .finding { color:var(--ink); font-size:0.94rem; margin:6px 0 4px 0; line-height:1.5; }
.gq-ev .meta { color:var(--muted); font-size:0.8rem; }

/* metric tiles */
.gq-tiles { display:flex; gap:12px; flex-wrap:wrap; }
.gq-tile { flex:1; min-width:150px; background:#fff; border:1px solid var(--line); border-radius:14px;
  padding:14px 16px; }
.gq-tile .v { font-size:1.5rem; font-weight:700; color:var(--ink); }
.gq-tile .l { font-size:0.78rem; color:var(--muted); text-transform:uppercase; letter-spacing:0.05em; }

/* captions */
.gq-cap { color:var(--muted); font-size:0.83rem; line-height:1.5; margin:2px 2px 14px 2px; }
.gq-note { color:var(--muted); font-size:0.86rem; }

/* Same shape as .gq-sim so the pair reads as a matched set: a viewer should be
   able to tell at a glance which panel is measured and which is constructed. */
.gq-real { display:inline-block; background:#ecfdf5; color:#047857;
  border:1px solid #a7f3d0;
  border-radius:8px; padding:2px 9px; font-size:0.74rem; font-weight:700;
  letter-spacing:0.03em; }
.gq-sim { display:inline-block; background:#fff7ed; color:#b45309; border:1px solid #fed7aa;
  border-radius:8px; padding:2px 9px; font-size:0.74rem; font-weight:700; letter-spacing:0.03em; }

section[data-testid="stSidebar"] { background:#0f1222; }
section[data-testid="stSidebar"] * { color:#e7e9f5 !important; }
section[data-testid="stSidebar"] .gq-side-mut { color:#9aa0c0 !important; font-size:0.82rem; }
section[data-testid="stSidebar"] code { background:#1b1f38 !important; color:#c7cbf0 !important;
  padding:1px 6px; border-radius:5px; font-size:0.78rem; }
.stTabs [data-baseweb="tab-list"] { gap: 4px; }
.stTabs [data-baseweb="tab"] { font-weight:600; }
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)


@st.cache_resource
def get_service() -> GeroService:
    return GeroService()


svc = get_service()


# --------------------------------------------------------------------------- #
# Small HTML helpers
# --------------------------------------------------------------------------- #


def pill(text: str, color: str) -> str:
    return (
        f"<span class='gq-pill' style='background:{color}1a;color:{color};"
        f"border-color:{color}44'>{text}</span>"
    )


def dir_pill(direction: str) -> str:
    glyph, color, label = DIR_META.get(direction, DIR_META["unknown"])
    return pill(f"{glyph} {label}", color)


def conf_pill(conf: str) -> str:
    color = CONF.get(conf, MUTED)
    return pill(conf.capitalize() + " evidence", color)


def caption(text: str) -> None:
    st.markdown(f"<div class='gq-cap'>{text}</div>", unsafe_allow_html=True)


# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #

with st.sidebar:
    st.markdown("<div class='gq-brandmark'>🧬 GeroQuery</div>", unsafe_allow_html=True)
    st.markdown(
        "<div class='gq-side-mut'>A search engine for the biology of aging. Type a gene, "
        "see what the research says about how it changes with age — plus a real biological-age "
        "clock and a resilience tool.</div>",
        unsafe_allow_html=True,
    )
    st.markdown("---")
    st.markdown("**Curated genes**", unsafe_allow_html=True)
    genes = svc.list_curated_genes()
    st.markdown(
        "<div class='gq-side-mut'>"
        + "  ".join(f"<code>{g['symbol']}</code>" for g in genes)
        + "</div>",
        unsafe_allow_html=True,
    )
    st.markdown("---")
    st.markdown(
        "<div class='gq-side-mut'><b>Data honesty.</b> Gene biology below is real, curated, and "
        "cited. The clock is the published PhenoAge model. The example biomarker cohort is "
        "<b>simulated</b> and labelled as such.</div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        f"<div class='gq-side-mut'>Data version <code>{svc.version()['data_version']}</code></div>",
        unsafe_allow_html=True,
    )


# --------------------------------------------------------------------------- #
# Hero
# --------------------------------------------------------------------------- #

st.markdown(
    "<div class='gq-hero'><h1>GeroQuery</h1>"
    "<p>Type a gene and instantly see what the research says about its relationship to aging — "
    "across species and data types — then score biological age and resilience from biomarkers.</p>"
    "</div>",
    unsafe_allow_html=True,
)

tab_gene, tab_clock, tab_res, tab_about = st.tabs(
    ["🔎  Gene explorer", "⏱  Aging clock", "📉  Resilience", "ℹ️  About & data"]
)


# ======================================================================= #
# GENE EXPLORER
# ======================================================================= #


def render_evidence_figure(evidence: list[dict]) -> None:
    """A qualitative 'aging signal map': one row per line of evidence, pointing
    up (increases) or down (decreases), length = curated confidence (NOT a
    measured effect size)."""
    strength_len = {"robust": 1.0, "established": 0.68, "emerging": 0.42}
    rows = list(reversed(evidence))  # first item on top
    labels, xs, colors, texts = [], [], [], []
    for e in rows:
        d = e["direction"]
        _, color, _ = DIR_META.get(d, DIR_META["unknown"])
        length = strength_len.get(e["strength"], 0.4)
        val = 0.0 if d == "context-dependent" else (length if d == "up" else -length)
        labels.append(f"{e['omic_layer']} · {'/'.join(e['species'])}")
        xs.append(val)
        colors.append(color)
        texts.append(e["strength"])
    fig = go.Figure()
    for lab, x, color, txt in zip(labels, xs, colors, texts, strict=True):
        if x == 0:  # context-dependent marker
            fig.add_trace(
                go.Scatter(
                    x=[0],
                    y=[lab],
                    mode="markers",
                    marker=dict(symbol="diamond", size=13, color=color),
                    hovertext=[f"{txt} · context-dependent"],
                    hoverinfo="text",
                    showlegend=False,
                )
            )
        else:
            fig.add_trace(
                go.Bar(
                    x=[x],
                    y=[lab],
                    orientation="h",
                    marker=dict(color=color, line=dict(width=0)),
                    width=0.5,
                    hovertext=[f"{txt} evidence"],
                    hoverinfo="text",
                    showlegend=False,
                )
            )
    fig.add_vline(x=0, line=dict(color=MUTED, width=1))
    fig.add_annotation(
        x=-0.82,
        y=1.12,
        xref="x",
        yref="paper",
        text="↓ decreases with age",
        showarrow=False,
        font=dict(color=DOWN, size=11),
    )
    fig.add_annotation(
        x=0.82,
        y=1.12,
        xref="x",
        yref="paper",
        text="increases with age ↑",
        showarrow=False,
        font=dict(color=UP, size=11),
    )
    fig.update_xaxes(range=[-1.15, 1.15], showticklabels=False, title="")
    _plotly_theme(fig, height=90 + 46 * len(rows), legend=False)
    st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False})


def render_gene_explorer():
    st.markdown("### Search a gene")
    c1, c2 = st.columns([4, 1])
    query = c1.text_input(
        "Gene symbol, alias, or ID",
        value=st.session_state.get("gene_query", "CDKN2A"),
        label_visibility="collapsed",
        placeholder="e.g. CDKN2A, p16, KLOTHO, GDF15, FOXO3, TP53…",
    )
    species_sel = c2.selectbox("Species", ["both", "human", "mouse"], label_visibility="collapsed")

    chips = ["CDKN2A", "GDF15", "KL", "LMNB1", "FOXO3", "TERT", "SIRT1", "MTOR"]
    cols = st.columns(len(chips))
    for col, name in zip(cols, chips, strict=True):
        if col.button(name, use_container_width=True, key=f"chip_{name}"):
            st.session_state["gene_query"] = name
            st.rerun()

    if not query.strip():
        st.info("Type a gene to begin — try **CDKN2A** (p16), the classic aging gene.")
        return

    species_arg = None if species_sel == "both" else species_sel
    try:
        report = svc.gene_report(query, species_arg)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Couldn't resolve “{query}”. {exc}")
        st.caption(
            "This demo ships a curated set of well-studied aging genes — see the sidebar "
            "for the full list. Live identifier resolution is available when network access "
            "is enabled."
        )
        return

    gene = report["gene"]
    know = report["knowledge"]

    # --- Header + verdict ---
    st.markdown(
        f"<div class='gq-card'>"
        f"<div style='display:flex;justify-content:space-between;align-items:baseline;flex-wrap:wrap'>"
        f"<div><span style='font-size:1.5rem;font-weight:700'>{gene['symbol']}</span> "
        f"<span style='color:{MUTED}'>· {gene.get('name','')}</span></div>"
        f"<div class='gq-note'>{know['aka'] if know else ''}</div></div>"
        f"<div class='gq-note' style='margin-top:6px'>Canonical <code>{gene['canonical_id']}</code>"
        f" · Entrez <code>{gene.get('entrez','—')}</code> · UniProt <code>{gene.get('uniprot','—')}</code>"
        f" · aliases {', '.join(gene.get('aliases',[])[:6]) or '—'}</div>"
        f"</div>",
        unsafe_allow_html=True,
    )

    if not know:
        st.warning("This gene resolves, but it isn't in the curated aging knowledge base yet.")
        return

    # Verdict hero
    glyph, color, label = DIR_META.get(know["direction_with_age"], DIR_META["unknown"])
    st.markdown(
        f"<div class='gq-card'><div class='gq-verdict'>"
        f"<div class='gq-glyph' style='color:{color}'>{glyph}</div>"
        f"<div class='gq-verdict-main'>"
        f"<div>{dir_pill(know['direction_with_age'])}{conf_pill(know['confidence'])}</div>"
        f"<div class='lead'>{know['one_liner']}</div></div></div></div>",
        unsafe_allow_html=True,
    )

    # Plain-English analysis
    st.markdown(
        "<div class='gq-card'><h4>What the data tells us</h4>"
        f"<div class='gq-analysis'>{know['analysis']}</div>"
        f"<div class='gq-note' style='margin-top:10px'><b>Role.</b> {know['role']}</div>"
        "</div>",
        unsafe_allow_html=True,
    )

    # Two columns: evidence figure + cross-species / hallmarks
    left, right = st.columns([3, 2])
    with left:
        st.markdown(
            "<div class='gq-card'><h4>Aging signal by data type</h4>", unsafe_allow_html=True
        )
        render_evidence_figure(know["evidence"])
        st.markdown("</div>", unsafe_allow_html=True)
        caption(
            "Each bar is one curated, cited line of evidence. Direction shows whether the gene "
            "goes <b>up</b> (right, red) or <b>down</b> (left, blue) with age; a diamond marks a "
            "context-dependent gene. Bar length encodes <b>curated confidence</b> "
            "(robust &gt; established &gt; emerging) — it is deliberately <b>not</b> a measured "
            "effect size, because GeroQuery does not invent numbers it did not measure."
        )
    with right:
        # Cross-species conservation
        orths = report["orthologs"]
        human = next((o for o in orths if o["species"] == "human"), None)
        mouse = next((o for o in orths if o["species"] == "mouse"), None)
        cons = "conserved" if human and mouse else "single species"
        st.markdown(
            "<div class='gq-card'><h4>Cross-species</h4>"
            + (f"<div class='gq-chip'>🧑 human · <b>{human['symbol']}</b></div>" if human else "")
            + (f"<div class='gq-chip'>🐭 mouse · <b>{mouse['symbol']}</b></div>" if mouse else "")
            + f"<div class='gq-note' style='margin-top:8px'>The aging relationship is "
            f"<b>{cons}</b> across the orthologs curated here.</div></div>",
            unsafe_allow_html=True,
        )
        # Hallmarks
        hm = report["hallmarks"]
        st.markdown(
            "<div class='gq-card'><h4>Hallmarks of aging</h4>"
            + "".join(
                f"<span class='gq-chip' title=\"{h['description']}\">{h['name']}</span>" for h in hm
            )
            + "</div>",
            unsafe_allow_html=True,
        )

    # Evidence detail
    st.markdown("<div class='gq-card'><h4>Evidence &amp; citations</h4>", unsafe_allow_html=True)
    refmap = {r["key"]: r for r in report["references"]}
    for e in know["evidence"]:
        _, ecolor, elabel = DIR_META.get(e["direction"], DIR_META["unknown"])
        cites = " · ".join(
            f"<a href='{refmap[k]['url']}' target='_blank'>{refmap[k]['authors']} {refmap[k]['year']}</a>"
            for k in e["reference_keys"]
            if k in refmap
        )
        st.markdown(
            f"<div class='gq-ev'>"
            f"<div>{pill(elabel, ecolor)}{pill(e['strength'].capitalize(), CONF.get(e['strength'], MUTED))}"
            f"<span class='gq-chip'>{e['omic_layer']}</span>"
            f"<span class='gq-chip'>{', '.join(e['tissues'])}</span></div>"
            f"<div class='finding'>{e['finding']}</div>"
            f"<div class='meta'>Species: {', '.join(e['species'])} &nbsp;·&nbsp; {cites}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)

    # Curated databases + interventions
    cdb, civ = st.columns(2)
    with cdb:
        st.markdown(
            "<div class='gq-card'><h4>Curated in aging databases</h4>", unsafe_allow_html=True
        )
        flags = report["curated_flags"]
        if flags:
            for f in flags:
                st.markdown(
                    f"<a href='{f.get('url','#')}' target='_blank' class='gq-db'>{f['database']}</a>"
                    f"<span class='gq-note'> — {f['assertion']}</span><br>",
                    unsafe_allow_html=True,
                )
        else:
            st.markdown(
                "<div class='gq-note'>Not flagged in the bundled curated set.</div>",
                unsafe_allow_html=True,
            )
        st.markdown("</div>", unsafe_allow_html=True)
        caption(
            "Membership in the HAGR databases (GenAge, CellAge, LongevityMap) and OpenGenes — "
            "real, human-curated catalogs of aging- and longevity-associated genes."
        )
    with civ:
        st.markdown("<div class='gq-card'><h4>Linked interventions</h4>", unsafe_allow_html=True)
        ivs = report["interventions"]
        if ivs:
            render_intervention_bars(ivs)
        else:
            st.markdown(
                "<div class='gq-note'>No linked lifespan interventions in the bundled set.</div>",
                unsafe_allow_html=True,
            )
        st.markdown("</div>", unsafe_allow_html=True)
        if ivs:
            caption(
                "Approximate <i>median-lifespan</i> changes reported in the cited rodent "
                "studies (NIA ITP / DrugAge). Effects vary widely by dose, sex, and strain — "
                "these are directional, not guarantees."
            )

    # FAQ
    st.markdown("<div class='gq-card'><h4>Common questions</h4></div>", unsafe_allow_html=True)
    for qa in know["faqs"]:
        with st.expander(qa["q"]):
            st.markdown(qa["a"])

    # References + download
    st.markdown("<div class='gq-card'><h4>References</h4>", unsafe_allow_html=True)
    for r in report["references"]:
        st.markdown(
            f"<div class='gq-note' style='margin-bottom:4px'>{r['citation']} "
            f"<a href='{r['url']}' target='_blank'>↗ PubMed</a></div>",
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)

    ev_df = pd.DataFrame(
        [
            {
                "gene": gene["symbol"],
                "direction_with_age": e["direction"],
                "omic_layer": e["omic_layer"],
                "species": "/".join(e["species"]),
                "tissues": "/".join(e["tissues"]),
                "strength": e["strength"],
                "finding": e["finding"],
                "references": "; ".join(e["reference_keys"]),
            }
            for e in know["evidence"]
        ]
    )
    st.download_button(
        "⬇ Download this gene's curated evidence (CSV)",
        ev_df.to_csv(index=False),
        file_name=f"{gene['symbol']}_aging_evidence.csv",
        mime="text/csv",
    )


def render_intervention_bars(ivs: list[dict]) -> None:
    ivs = sorted(ivs, key=lambda d: d["lifespan_effect_pct"])
    fig = go.Figure(
        go.Bar(
            x=[i["lifespan_effect_pct"] for i in ivs],
            y=[i["display_name"] for i in ivs],
            orientation="h",
            marker=dict(color=ACCENT, line=dict(width=0)),
            text=[f"+{i['lifespan_effect_pct']:.0f}%" for i in ivs],
            textposition="outside",
            hovertext=[i["summary"] for i in ivs],
            hoverinfo="text",
        )
    )
    fig.update_xaxes(
        title="≈ median lifespan change (%)",
        range=[0, max(i["lifespan_effect_pct"] for i in ivs) * 1.25],
    )
    _plotly_theme(fig, height=70 + 44 * len(ivs), legend=False)
    st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False})


with tab_gene:
    render_gene_explorer()


# ======================================================================= #
# AGING CLOCK
# ======================================================================= #

CLOCK_TEMPLATE_COLS = list(REQUIRED_FEATURES)
_UNIT_HELP = {
    "albumin_gdl": "Albumin (g/dL)",
    "creatinine_mgdl": "Creatinine (mg/dL)",
    "glucose_mgdl": "Glucose (mg/dL)",
    "crp_mgl": "C-reactive protein (mg/L)",
    "lymphocyte_pct": "Lymphocytes (%)",
    "mcv_fl": "Mean cell volume (fL)",
    "rdw_pct": "Red cell distribution width (%)",
    "alp_ul": "Alkaline phosphatase (U/L)",
    "wbc_1000ul": "White blood cells (10³/µL)",
    "age": "Chronological age (years)",
}


def render_clock():
    st.markdown("### PhenoAge · a real biological-age clock")
    st.markdown(
        "<div class='gq-card'><div class='gq-analysis'>"
        "A biological-age clock asks: <b>given your blood chemistry, how old does your body look?</b> "
        "GeroQuery ships <b>PhenoAge</b> (Levine et al., 2018) — a real, published clock built from "
        "nine routine blood markers plus age. It returns a <b>phenotypic age</b> in years and, from "
        "the very same model, a <b>10-year mortality risk</b> — a clean way to see that "
        "‘biological age’ and ‘mortality risk’ are two faces of one calibrated model."
        "</div></div>",
        unsafe_allow_html=True,
    )

    clocks = svc.list_clocks()
    ci = clocks[0]
    st.markdown(
        f"<div class='gq-note'>Predicts <b>{ci['predicted_outcome'].replace('_',' ')}</b> · "
        f"units <b>{ci['units']}</b> · trained on <b>{ci['training_population']}</b>.</div>",
        unsafe_allow_html=True,
    )

    src = st.radio("Data source", ["Real NHANES cohort", "Upload my own CSV"], horizontal=True)
    df = None
    if src == "Real NHANES cohort":
        df = svc.store.get_dataset("clinical_nhanes_phenoage")
        st.markdown(
            "<span class='gq-real'>REAL DATA</span> "
            f"<span class='gq-note'>{len(df):,} US adults from NHANES 2017-2018 "
            "carrying all nine markers the published PhenoAge model needs. Fetched "
            "from CDC, verified against a pinned SHA-256.</span>",
            unsafe_allow_html=True,
        )
    else:
        st.caption(
            "Upload a CSV with these columns (conventional US clinical units): "
            + ", ".join(f"`{c}`" for c in CLOCK_TEMPLATE_COLS)
        )
        template = pd.DataFrame([{c: "" for c in CLOCK_TEMPLATE_COLS}])
        st.download_button(
            "⬇ Download blank template",
            template.to_csv(index=False),
            file_name="phenoage_template.csv",
            mime="text/csv",
        )
        up = st.file_uploader("Upload biomarker CSV", type=["csv"])
        if up is not None:
            df = pd.read_csv(up)

    if df is None:
        st.info("Pick the simulated cohort or upload a CSV, then run the clock.")
        return

    if not st.button("Apply PhenoAge clock", type="primary"):
        return

    try:
        res = svc.apply_clock("phenoage", df, df["age"].tolist() if "age" in df.columns else None)
    except Exception as exc:  # noqa: BLE001
        st.error(str(exc))
        return

    pred = np.array(res["predictions"])
    chrono = np.array(df["age"]) if "age" in df.columns else None
    risk = np.array(res.get("mortality_risk_10yr", []))

    # Tiles
    if res.get("mean_age_acceleration") is not None:
        faster = int(np.sum(np.array(res["age_acceleration"]) > 0))
        pct_faster = 100 * faster / len(pred)
        st.markdown(
            "<div class='gq-tiles'>"
            f"<div class='gq-tile'><div class='v'>{res['n_samples']}</div><div class='l'>subjects</div></div>"
            f"<div class='gq-tile'><div class='v'>{np.mean(pred):.1f} yr</div><div class='l'>mean biological age</div></div>"
            f"<div class='gq-tile'><div class='v'>{res['mean_age_acceleration']:+.1f} yr</div><div class='l'>mean age acceleration</div></div>"
            f"<div class='gq-tile'><div class='v'>{pct_faster:.0f}%</div><div class='l'>aging faster than birthday</div></div>"
            "</div>",
            unsafe_allow_html=True,
        )
        st.markdown("<div style='height:12px'></div>", unsafe_allow_html=True)

    col1, col2 = st.columns(2)
    if chrono is not None:
        with col1:
            accel = np.array(res["age_acceleration"])
            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=chrono,
                    y=pred,
                    mode="markers",
                    marker=dict(
                        size=7,
                        color=accel,
                        colorscale=[[0, DOWN], [0.5, "#e9ecf5"], [1, UP]],
                        cmid=0,
                        showscale=True,
                        colorbar=dict(title="accel (yr)", thickness=10),
                        line=dict(width=0),
                    ),
                    hovertemplate="chrono %{x:.0f} → PhenoAge %{y:.1f}<extra></extra>",
                )
            )
            lo, hi = float(np.min(chrono)), float(np.max(chrono))
            fig.add_trace(
                go.Scatter(
                    x=[lo, hi],
                    y=[lo, hi],
                    mode="lines",
                    line=dict(dash="dot", color=MUTED),
                    name="same age",
                )
            )
            fig.update_xaxes(title="chronological age (years)")
            fig.update_yaxes(title="biological age — PhenoAge (years)")
            _plotly_theme(fig, height=380, legend=False)
            st.plotly_chart(fig, use_container_width=True)
            caption(
                "Each point is a subject. Points <b>above</b> the dotted line are biologically "
                "older than their birthday (red, positive acceleration); points <b>below</b> are "
                "younger (blue). A tight cloud around the line means the clock and the calendar "
                "mostly agree."
            )
        with col2:
            fig2 = go.Figure(
                go.Histogram(
                    x=res["age_acceleration"],
                    nbinsx=30,
                    marker=dict(color=ACCENT, line=dict(width=0)),
                )
            )
            fig2.add_vline(x=0, line=dict(color=MUTED, dash="dot"))
            fig2.update_xaxes(title="age acceleration = biological − chronological (years)")
            fig2.update_yaxes(title="subjects")
            _plotly_theme(fig2, height=380, legend=False)
            st.plotly_chart(fig2, use_container_width=True)
            caption(
                "The spread of biological-age gaps. Mass to the <b>right</b> of zero = aging "
                "faster than the calendar; mass to the <b>left</b> = aging slower. The width of "
                "this distribution is the population's variation in pace of aging."
            )

    if risk.size:
        st.markdown(
            "<div class='gq-card'><h4>Same model, different question: mortality risk</h4>",
            unsafe_allow_html=True,
        )
        order = np.argsort(chrono) if chrono is not None else np.argsort(pred)
        xvals = (chrono if chrono is not None else pred)[order]
        fig3 = go.Figure(
            go.Scatter(
                x=xvals,
                y=100 * risk[order],
                mode="markers",
                marker=dict(size=6, color=UP, opacity=0.6, line=dict(width=0)),
                hovertemplate="age %{x:.0f} → %{y:.1f}% 10-yr risk<extra></extra>",
            )
        )
        fig3.update_xaxes(title="chronological age (years)")
        fig3.update_yaxes(title="predicted 10-year mortality risk (%)")
        _plotly_theme(fig3, height=320, legend=False)
        st.plotly_chart(fig3, use_container_width=True)
        st.markdown("</div>", unsafe_allow_html=True)
        caption(
            "PhenoAge is built on a mortality model, so it also yields a 10-year mortality "
            "risk for each subject. This is why a clock's <b>age accuracy</b> and its "
            "<b>mortality prediction</b> are related but not the same thing."
        )

    out = pd.DataFrame({"biological_age_phenoage": pred})
    if chrono is not None:
        out.insert(0, "chronological_age", chrono)
        out["age_acceleration"] = res["age_acceleration"]
    if risk.size:
        out["mortality_risk_10yr"] = risk
    st.download_button(
        "⬇ Download results (CSV)",
        out.to_csv(index=False),
        file_name="phenoage_results.csv",
        mime="text/csv",
    )

    st.markdown(
        "<div class='gq-note' style='margin-top:8px'><b>Keep in mind.</b> PhenoAge is a research "
        "instrument for populations, not a medical diagnosis. On the simulated cohort the numbers "
        "illustrate the method; on your own data, results depend on assay units matching the "
        "expected columns.</div>",
        unsafe_allow_html=True,
    )


with tab_clock:
    render_clock()


# ======================================================================= #
# RESILIENCE
# ======================================================================= #


def render_resilience():
    st.markdown("### Resilience · early-warning signals of aging")
    st.markdown(
        "<div class='gq-card'><div class='gq-analysis'>"
        "Complex systems that are close to breaking down get <b>slower to recover</b> from small "
        "knocks and their signals become <b>more erratic and more synchronised</b>. Ecologists and "
        "physicists call these <b>critical-slowing-down</b> early-warning signals; the same idea "
        "applies to an aging body. GeroQuery looks for two of them across age: rising "
        "<b>variance</b> (the state wobbles more) and rising <b>cross-correlation</b> between "
        "markers (they move together more). When both climb with age, resilience is eroding."
        "</div></div>",
        unsafe_allow_html=True,
    )

    src = st.radio(
        "Data source",
        ["Simulated example cohort", "Upload my own CSV"],
        horizontal=True,
        key="res_src",
    )
    df = None
    if src == "Real NHANES cohort":
        df = svc.store.get_dataset("clinical_nhanes_slice")
        st.markdown(
            "<span class='gq-real'>REAL DATA</span> "
            f"<span class='gq-note'>{len(df):,} US adults, NHANES 2017-2018, six "
            "routine blood markers. On this cohort the variance signal replicates "
            "and the cross-correlation signal does not.</span>",
            unsafe_allow_html=True,
        )
    elif src == "Method-validation fixture":
        df = svc.store.get_dataset("clinical_synthetic_csd")
        st.markdown(
            "<span class='gq-sim'>SYNTHETIC ON PURPOSE</span> "
            "<span class='gq-note'>720 subjects with critical slowing down planted "
            "by construction. Not evidence about people — it exists so the "
            "estimator can be checked against a known answer. Compare its verdict "
            "with the real cohort's.</span>",
            unsafe_allow_html=True,
        )
    else:
        st.caption(
            "Upload a cross-sectional CSV with an `age` column and numeric biomarker columns."
        )
        up = st.file_uploader("Upload cohort CSV", type=["csv"], key="res_upload")
        if up is not None:
            df = pd.read_csv(up)

    n_strata = st.slider(
        "Age strata",
        3,
        12,
        6,
        help="How many age bands to split the cohort into before measuring the "
        "early-warning signals.",
    )

    if df is None:
        st.info("Pick the simulated cohort or upload a CSV, then compute resilience.")
        return
    if not st.button("Compute resilience", type="primary"):
        return

    biomarkers = [c for c in df.columns if c not in ("subject_id", "age", "sex")]
    try:
        res = svc.resilience_csd(data=df, biomarker_cols=biomarkers, n_strata=n_strata)
    except Exception as exc:  # noqa: BLE001
        st.error(str(exc))
        return

    declines = res["resilience_declines"]
    verdict_color = UP if declines else CONF["robust"]
    verdict = "Resilience declines with age" if declines else "No clear resilience decline detected"
    glyph = "⚠" if declines else "✓"
    if declines:
        lead = (
            "Both variance and cross-correlation rise across age bands — the system is slower "
            "to recover and its markers move together more, the signature of eroding resilience."
        )
    else:
        lead = "The early-warning signals did not both rise consistently with age in this dataset."
    st.markdown(
        f"<div class='gq-card'><div class='gq-verdict'>"
        f"<div class='gq-glyph' style='color:{verdict_color}'>{glyph}</div>"
        f"<div class='gq-verdict-main'><div>{pill(verdict, verdict_color)}</div>"
        f"<div class='lead'>{lead}</div></div></div></div>",
        unsafe_allow_html=True,
    )

    frame = pd.DataFrame(
        {
            "age": res["strata_midpoints"],
            "variance": res["variance"],
            "cross_correlation": res["cross_correlation"],
        }
    )
    c1, c2 = st.columns(2)
    with c1:
        f1 = go.Figure(
            go.Scatter(
                x=frame["age"],
                y=frame["variance"],
                mode="lines+markers",
                line=dict(color=UP, width=2.5),
                marker=dict(size=8, color=UP),
            )
        )
        f1.update_xaxes(title="age (band midpoint, years)")
        f1.update_yaxes(title="variance of health state")
        _plotly_theme(f1, height=320, legend=False)
        st.plotly_chart(f1, use_container_width=True)
        caption(
            "<b>Variance vs age.</b> How much the overall biomarker state wobbles within each "
            "age band. A rising line means older bands are less stable — a core early-warning "
            "signal."
        )
    with c2:
        f2 = go.Figure(
            go.Scatter(
                x=frame["age"],
                y=frame["cross_correlation"],
                mode="lines+markers",
                line=dict(color=ACCENT, width=2.5),
                marker=dict(size=8, color=ACCENT),
            )
        )
        f2.update_xaxes(title="age (band midpoint, years)")
        f2.update_yaxes(title="mean cross-correlation between markers")
        _plotly_theme(f2, height=320, legend=False)
        st.plotly_chart(f2, use_container_width=True)
        caption(
            "<b>Cross-correlation vs age.</b> How synchronised the markers are. When markers "
            "increasingly move together with age, the system has lost independent buffering "
            "capacity — the second early-warning signal."
        )

    with st.expander("Method, assumptions & limitations (read this)"):
        st.markdown(f"**Method:** `{res['method']}`")
        for a in res["assumptions"]:
            st.markdown(f"- {a}")
        st.markdown(
            "This is a **cross-sectional proxy**: it compares different people at different ages, "
            "not the same person recovering over time. True critical-slowing-down is a "
            "*longitudinal* phenomenon; with repeated within-person measurements the recovery-rate "
            "metric applies directly. The signal is real dynamical-systems theory (Scheffer et al., "
            "2009; Gijzel et al., 2017; Pyrkov et al., 2021), applied here as an age-stratified "
            "approximation."
        )


with tab_res:
    render_resilience()


# ======================================================================= #
# ABOUT & DATA
# ======================================================================= #

with tab_about:
    st.markdown("### About GeroQuery")
    st.markdown(
        "<div class='gq-card'><div class='gq-analysis'>"
        "GeroQuery unifies the scattered evidence on how genes relate to aging into one "
        "gene-first search experience, and adds two biomarker tools — a real biological-age clock "
        "and a resilience/early-warning analysis. It exists to save aging researchers from "
        "manually cross-referencing a dozen separate databases to answer one simple question."
        "</div></div>",
        unsafe_allow_html=True,
    )

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            "<div class='gq-card'><h4>What is real</h4>"
            "<div class='gq-note'>"
            "• <b>Gene aging effects</b> — measured, not asserted. 485,905 "
            "Hedges' <i>g</i> estimates over 46,091 genes from 31 checksum-pinned "
            "GEO DataSets, pooled by random effects, with intervals wide enough "
            "to say <i>we cannot tell</i>.<br><br>"
            "• <b>Curated flags</b> — the ingested HAGR releases: GenAge, CellAge, "
            "LongevityMap (including its null results), GenDR.<br><br>"
            "• <b>Interventions</b> — DrugAge and GenDR, one record per compound "
            "and organism, median over the significant experiments only.<br><br>"
            "• <b>Clinical cohorts</b> — NHANES 2017-2018 from CDC, checksum "
            "verified.<br><br>"
            "• <b>The clocks</b> — the published Levine PhenoAge coefficients plus "
            "236 wrapped biolearn and pyaging clocks, validated on real GEO "
            "methylation."
            "</div></div>",
            unsafe_allow_html=True,
        )
    with c2:
        st.markdown(
            "<div class='gq-card'><h4>What is simulated (and labelled)</h4>"
            "<div class='gq-note'>"
            "• <b>The method-validation fixture</b> is "
            "<span class='gq-sim'>SYNTHETIC ON PURPOSE</span> — 720 subjects with "
            "critical slowing down planted by construction. Not a stand-in for "
            "missing data: it is the only way to ask whether the estimator "
            "recovers an effect <i>known</i> to be there, which no real dataset "
            "can answer. Its own dataset id, and a test asserts it never merges "
            "with the NHANES table.<br><br>"
            "Earlier versions of this project shipped fabricated GEO study accessions, effect "
            "sizes, and p-values to fake a meta-analysis. <b>Those have been removed entirely.</b> "
            "GeroQuery no longer invents evidence."
            "</div></div>",
            unsafe_allow_html=True,
        )

    st.markdown("<div class='gq-card'><h4>Data sources &amp; licences</h4>", unsafe_allow_html=True)
    src_df = pd.DataFrame(svc.sources())
    show = src_df[["name", "omics", "federated", "cacheable"]].copy()
    show["license"] = [s["name"] for s in src_df["license"]]
    st.dataframe(show, use_container_width=True, hide_index=True)
    st.markdown("</div>", unsafe_allow_html=True)
    caption(
        "Open, redistributable curated data is bundled; large or controlled sources "
        "(UK Biobank, protected GTEx) are federate/link-only and never re-hosted — enforced "
        "in code by each adapter's licence gate."
    )

    st.markdown("<div class='gq-card'><h4>All references</h4>", unsafe_allow_html=True)
    for r in svc.references():
        st.markdown(
            f"<div class='gq-note' style='margin-bottom:4px'>{r['citation']} "
            f"<a href='{r['url']}' target='_blank'>↗</a></div>",
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)
