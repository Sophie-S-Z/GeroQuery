"""The Streamlit dashboard, exercised headlessly.

These exist because three shipped bugs got past a green test suite: the gene
explorer showed "isn't in the curated aging knowledge base yet" for *every* gene,
the resilience tab could never load a cohort, and the references linked to the
wrong papers. None of it was caught, because nothing ran the UI.

`AppTest` runs the real script. It is the only test in this repo that would have
caught any of the three.
"""

from __future__ import annotations

import pytest

pytest.importorskip("streamlit", reason="dashboard extra not installed")

from streamlit.testing.v1 import AppTest  # noqa: E402

APP = "geroquery/ui/streamlit_app.py"
TIMEOUT = 300


@pytest.fixture(scope="module")
def app():
    return AppTest.from_file(APP, default_timeout=TIMEOUT).run()


def _body(at: AppTest) -> str:
    return " ".join(m.value for m in at.markdown)


def test_the_app_renders_without_raising(app):
    assert not app.exception, [e.value for e in app.exception]


def test_all_four_tabs_are_present(app):
    assert [t.label for t in app.tabs] == [
        "Gene explorer",
        "Aging clock",
        "Resilience",
        "About & data",
    ]


def test_a_gene_with_data_gets_a_verdict_not_a_dead_end(app):
    """The regression that shipped: `report["knowledge"]` was hardcoded None, so
    every gene fell through to a 'not in the knowledge base' warning and the rest
    of the explorer was unreachable."""
    body = _body(app)
    assert "curated aging knowledge base" not in body
    # CDKN1A is the default query: human CI excludes zero, mouse spans it.
    assert "rises with age" in body
    assert "no detectable change" in body


def test_the_verdict_is_read_from_the_interval_not_the_point_estimate(app):
    """A pooled effect can be positive while its interval spans zero. The page
    must say 'no detectable change' in that case, or it is asserting an effect
    the data does not support."""
    body = _body(app)
    assert "the interval spans zero" in body or "interval lies entirely" in body


def test_masthead_counts_are_computed_not_typed(app):
    """A hardcoded corpus size on a page whose argument is provenance goes stale
    the first time the panel changes."""
    from geroquery.api.service import GeroService

    signatures = GeroService().store.query_signatures()
    assert f"{len(signatures):,} effect sizes" in _body(app)


@pytest.mark.parametrize(
    ("cohort", "expected"),
    [
        ("Real NHANES cohort", "No resilience decline detected"),
        ("Method-validation fixture", "Resilience declines with age"),
    ],
)
def test_resilience_runs_and_separates_the_two_cohorts(cohort, expected):
    """The shipped bug: the radio's options were rewritten but its branches were
    not, so `src` never matched and no cohort could ever load.

    The pair also pins the science — the estimator must find the planted effect
    and decline to claim it on real data.
    """
    at = AppTest.from_file(APP, default_timeout=TIMEOUT).run()
    at.radio(key="res_src").set_value(cohort).run()
    compute = next(b for b in at.button if "early-warning" in b.label.lower())
    compute.click().run()

    assert not at.exception, [e.value for e in at.exception]
    body = _body(at)
    assert expected in body
    assert "Kendall" in body  # the trend table rendered


def test_the_clock_runs_on_the_real_nhanes_cohort():
    at = AppTest.from_file(APP, default_timeout=TIMEOUT).run()
    run = next(b for b in at.button if "phenoage" in b.label.lower())
    run.click().run()
    assert not at.exception, [e.value for e in at.exception]
    labels = [m.label for m in at.metric]
    assert "r with chronological age" in labels


def test_both_themes_render_and_carry_the_brief_typography():
    for dark, present, absent in ((True, "#12110f", "#faf7f0"), (False, "#faf7f0", "#12110f")):
        at = AppTest.from_file(APP, default_timeout=TIMEOUT).run()
        at.toggle[0].set_value(dark).run()
        assert not at.exception, [e.value for e in at.exception]
        style = " ".join(m.value for m in at.markdown if "<style>" in m.value)
        assert present in style and absent not in style
        assert "EB+Garamond" in style and "Courier+Prime" in style


def test_constructed_data_is_badged_wherever_it_appears():
    at = AppTest.from_file(APP, default_timeout=TIMEOUT).run()
    at.radio(key="res_src").set_value("Method-validation fixture").run()
    assert "constructed" in _body(at)
