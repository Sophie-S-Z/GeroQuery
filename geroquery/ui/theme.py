"""Visual system for the GeroQuery dashboard.

Two decisions drive everything here.

**The typeface pairing is the argument.** EB Garamond carries the prose, because
this surface asks you to *read* an interpretation before you act on it — a
confidence interval spanning zero is a sentence, not a number. Courier Prime
carries every measured quantity: effect sizes, intervals, gene identifiers,
sample counts. The split is not decorative. A reader should be able to tell, from
the shape of the glyphs alone, which parts of the page were measured and which
were written.

**Dark is the default.** The use scene is a researcher reading interval estimates
for long stretches, usually indoors, often at night. Both themes are supplied and
switch at runtime; the palettes are tuned separately rather than one being an
inversion of the other, because a serif at 16px needs different contrast headroom
on dark than on light.

Every colour pair below meets WCAG AA for its role (body ≥4.5:1, large ≥3:1).
"""

from __future__ import annotations

from dataclasses import dataclass

FONT_IMPORT = (
    "@import url('https://fonts.googleapis.com/css2?"
    "family=EB+Garamond:ital,wght@0,400;0,500;0,600;1,400&"
    "family=Courier+Prime:ital,wght@0,400;0,700;1,400&display=swap');"
)

SERIF = "'EB Garamond', 'Iowan Old Style', Georgia, serif"
MONO = "'Courier Prime', 'SFMono-Regular', Consolas, monospace"


@dataclass(frozen=True)
class Palette:
    name: str
    bg: str  # page
    surface: str  # raised panel
    surface_alt: str  # inset / table stripe
    line: str  # hairline
    ink: str  # body text
    muted: str  # secondary text
    accent: str  # interactive / brand
    up: str  # increases with age
    down: str  # decreases with age
    null: str  # no detectable effect
    real: str  # measured-data badge
    sim: str  # constructed-data badge
    plot_grid: str


DARK = Palette(
    name="dark",
    bg="#12110f",
    surface="#1a1917",
    surface_alt="#211f1c",
    line="#33302b",
    ink="#f2ede3",
    muted="#a49c8e",
    accent="#d4a24c",
    up="#e2705f",
    down="#6aa9e0",
    null="#8f8779",
    real="#7bc49a",
    sim="#d9a441",
    plot_grid="#2b2823",
)

LIGHT = Palette(
    name="light",
    bg="#faf7f0",
    surface="#ffffff",
    surface_alt="#f4efe4",
    line="#ddd5c4",
    ink="#241f18",
    muted="#6b6355",
    accent="#8a5a12",
    up="#b23a29",
    down="#1f5f96",
    null="#7c7466",
    real="#1f6b45",
    sim="#8a5a12",
    plot_grid="#e6dfd0",
)

PALETTES = {"dark": DARK, "light": LIGHT}


def css(p: Palette) -> str:
    """The whole stylesheet, parameterised by palette."""
    return f"""
<style>
{FONT_IMPORT}

:root {{
  --bg:{p.bg}; --surface:{p.surface}; --surface-alt:{p.surface_alt};
  --line:{p.line}; --ink:{p.ink}; --muted:{p.muted}; --accent:{p.accent};
  --up:{p.up}; --down:{p.down}; --null:{p.null};
  --real:{p.real}; --sim:{p.sim};
  --serif:{SERIF}; --mono:{MONO};
}}

/* ---- base ---------------------------------------------------------- */
.stApp {{ background:var(--bg); }}
html, body, [class*="css"], .stMarkdown, p, li, span, div {{
  font-family:var(--serif); color:var(--ink);
}}
.block-container {{ padding-top:2.2rem; padding-bottom:5rem; max-width:1120px; }}

h1,h2,h3,h4 {{
  font-family:var(--serif); color:var(--ink); font-weight:600;
  letter-spacing:-0.01em;
}}
h1 {{ font-size:2.6rem; line-height:1.1; margin:0 0 .3rem 0; }}
h2 {{ font-size:1.6rem; margin:2.6rem 0 .9rem 0; }}
h3 {{ font-size:1.22rem; margin:2.2rem 0 .8rem 0; }}
p {{ font-size:1.06rem; line-height:1.62; }}
a {{ color:var(--accent); text-decoration:none; border-bottom:1px solid transparent; }}
a:hover {{ border-bottom-color:var(--accent); }}
code, .gq-num {{ font-family:var(--mono); font-size:.92em; }}

/* Everything measured is monospace. That is the tell. */
.gq-measure {{ font-family:var(--mono); font-variant-numeric:tabular-nums; }}

/* ---- masthead ------------------------------------------------------ */
.gq-mast {{ border-bottom:1px solid var(--line); padding-bottom:1.6rem; margin-bottom:.4rem; }}
.gq-mast h1 {{ font-weight:500; }}
.gq-mast .lede {{
  font-size:1.12rem; color:var(--muted); max-width:62ch; margin:.55rem 0 0 0; line-height:1.6;
}}
.gq-mast .meta {{
  font-family:var(--mono); font-size:.76rem; color:var(--muted);
  letter-spacing:.02em; margin-top:1rem;
}}

/* ---- tabs: the brief asked for room ------------------------------- */
.stTabs [data-baseweb="tab-list"] {{
  gap:2.6rem; border-bottom:1px solid var(--line);
  margin:2rem 0 2.2rem 0; padding-bottom:0;
}}
.stTabs [data-baseweb="tab"] {{
  font-family:var(--serif); font-size:1.06rem; font-weight:500;
  color:var(--muted); background:transparent;
  padding:.7rem .2rem 1rem .2rem; height:auto;
}}
.stTabs [data-baseweb="tab"]:hover {{ color:var(--ink); }}
.stTabs [aria-selected="true"] {{ color:var(--ink); }}
.stTabs [data-baseweb="tab-highlight"] {{ background:var(--accent); height:2px; }}
.stTabs [data-baseweb="tab-border"] {{ display:none; }}

/* ---- panels -------------------------------------------------------- */
.gq-panel {{
  background:var(--surface); border:1px solid var(--line);
  border-radius:14px; padding:1.5rem 1.7rem; margin-bottom:1.1rem;
}}
.gq-panel h4 {{
  font-family:var(--mono); font-size:.72rem; font-weight:700;
  letter-spacing:.11em; text-transform:uppercase; color:var(--muted);
  margin:0 0 1rem 0;
}}

/* ---- the verdict --------------------------------------------------- */
.gq-verdict {{ margin:.2rem 0 .1rem 0; }}
.gq-verdict .claim {{ font-size:1.5rem; line-height:1.32; font-weight:500; }}
.gq-verdict .because {{
  color:var(--muted); font-size:1rem; line-height:1.6; margin-top:.55rem; max-width:66ch;
}}

/* ---- the interval, as a figure ------------------------------------ */
.gq-interval {{ font-family:var(--mono); font-size:1.42rem; letter-spacing:-.01em; }}
.gq-interval .ci {{ font-size:.92rem; color:var(--muted); margin-left:.5rem; }}

/* ---- badges -------------------------------------------------------- */
.gq-badge {{
  display:inline-block; font-family:var(--mono); font-size:.66rem; font-weight:700;
  letter-spacing:.1em; text-transform:uppercase;
  padding:.3rem .6rem; border-radius:5px; border:1px solid transparent;
}}
.gq-real {{ color:var(--real); border-color:var(--real); }}
.gq-sim  {{ color:var(--sim);  border-color:var(--sim);  }}
.gq-tag  {{ color:var(--muted); border-color:var(--line); }}

/* ---- data tables --------------------------------------------------- */
.gq-table {{ width:100%; border-collapse:collapse; font-family:var(--mono); font-size:.83rem; }}
.gq-table th {{
  text-align:left; font-weight:700; color:var(--muted); text-transform:uppercase;
  letter-spacing:.07em; font-size:.68rem; padding:.55rem .7rem;
  border-bottom:1px solid var(--line);
}}
.gq-table td {{ padding:.5rem .7rem; border-bottom:1px solid var(--line); color:var(--ink); }}
.gq-table tr:last-child td {{ border-bottom:none; }}
.gq-table td.num {{ text-align:right; font-variant-numeric:tabular-nums; }}

/* ---- notes --------------------------------------------------------- */
.gq-note {{ color:var(--muted); font-size:.93rem; line-height:1.6; }}
.gq-caveat {{
  font-family:var(--serif); font-style:italic; color:var(--muted);
  font-size:.98rem; line-height:1.6; max-width:66ch;
}}

/* ---- controls ------------------------------------------------------ */
.stButton > button {{
  font-family:var(--mono); font-size:.78rem; font-weight:700; letter-spacing:.07em;
  text-transform:uppercase; border-radius:7px;
  border:1px solid var(--accent); background:transparent; color:var(--accent);
  padding:.6rem 1.2rem; transition:background .16s ease, color .16s ease;
}}
.stButton > button:hover {{ background:var(--accent); color:var(--bg); border-color:var(--accent);
  }}
.stButton > button:focus-visible {{ outline:2px solid var(--accent); outline-offset:2px; }}

.stTextInput input, .stSelectbox div[data-baseweb="select"] > div {{
  font-family:var(--mono); font-size:.9rem;
  background:var(--surface); border-color:var(--line); color:var(--ink);
  border-radius:8px;
}}
.stTextInput input::placeholder {{ color:var(--muted); }}
.stRadio label, .stCheckbox label, .stSlider label {{ font-family:var(--serif); color:var(--ink); }}

/* ---- sidebar ------------------------------------------------------- */
section[data-testid="stSidebar"] {{ background:var(--surface); border-right:1px solid var(--line);
  }}
section[data-testid="stSidebar"] * {{ color:var(--ink); }}
section[data-testid="stSidebar"] h2 {{ font-size:1.05rem; margin:1.2rem 0 .5rem 0; }}
section[data-testid="stSidebar"] .gq-note {{ color:var(--muted); }}

/* ---- misc ---------------------------------------------------------- */
hr {{ border:none; border-top:1px solid var(--line); margin:2.4rem 0; }}
[data-testid="stMetricValue"] {{ font-family:var(--mono); color:var(--ink); }}
[data-testid="stMetricLabel"] {{ font-family:var(--serif); color:var(--muted); }}
#MainMenu, footer {{ visibility:hidden; }}
</style>
"""


def plotly_layout(p: Palette, height: int = 340) -> dict:
    """Chart chrome that matches the page rather than fighting it."""
    return dict(
        template="none",
        height=height,
        margin=dict(l=8, r=14, t=26, b=8),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Courier Prime, monospace", size=12, color=p.muted),
        hoverlabel=dict(
            font=dict(family="Courier Prime, monospace", size=12),
            bgcolor=p.surface,
            bordercolor=p.line,
        ),
        xaxis=dict(gridcolor=p.plot_grid, zerolinecolor=p.line, linecolor=p.line),
        yaxis=dict(gridcolor=p.plot_grid, zerolinecolor=p.line, linecolor=p.line),
        showlegend=False,
    )
