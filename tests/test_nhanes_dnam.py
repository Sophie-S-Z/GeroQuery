"""Tests for the NHANES 1999-2002 cross-layer adapter.

Everything here runs offline against the committed 600-row real sample except
the ``live`` tests, which are the only thing that would notice NCHS re-issuing
a file.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from geroquery.sources import nhanes_dnam as nd
from geroquery.sources.manifest import (
    NHANES_1999_2002,
    NHANES_DNAM,
    NHANES_MORTALITY,
)


@pytest.fixture(scope="module")
def sample() -> pd.DataFrame:
    return nd.load_sample()


# --- the two coding traps ----------------------------------------------------


def test_dnam_sex_coding_is_the_inverse_of_nhanes_sex_coding():
    """XY_Estimation is 1=female, 2=male; RIAGENDR is 1=male, 2=female.

    Both are plausible-looking 1/2 integers, so reusing one map for the other
    fails silently and mislabels every subject. This pins the inversion.
    """
    assert nd.SEX_LABELS == {1: "male", 2: "female"}
    assert nd.DNAM_SEX_LABELS == {1: "female", 2: "male"}
    assert nd.DNAM_SEX_LABELS[1] != nd.SEX_LABELS[1]


def test_dnam_and_reported_sex_agree_on_the_committed_sample(sample):
    """~97.6% agreement is the check working; ~2.4% would mean it is inverted."""
    qc = nd.sex_discordance(sample)
    assert qc["n_compared"] == len(sample)
    assert qc["rate"] < 0.10, (
        "DNAm-predicted sex disagrees with reported sex on more than 10% of "
        "samples. Either the coding was flipped again or the upstream file "
        "changed."
    )


def test_the_two_cycles_use_different_variable_names_for_the_same_marker():
    """Creatinine and ALP are renamed between 1999-2000 and 2001-2002.

    Using one cycle's names for both drops the other cycle's subjects entirely,
    and does it without raising, so the failure looks like a smaller cohort
    rather than a bug.
    """
    by_name = {cycle.name: cycle for cycle in nd.CYCLES}
    assert by_name["1999-2000"].markers["creatinine"] == ("bio", "LBXSCR")
    assert by_name["2001-2002"].markers["creatinine"] == ("bio", "LBDSCR")
    assert by_name["1999-2000"].markers["alp"] == ("bio", "LBXSAPSI")
    assert by_name["2001-2002"].markers["alp"] == ("bio", "LBDSAPSI")


def test_both_cycles_are_represented_in_the_built_cohort(sample):
    """The direct consequence of the map above being right."""
    cycles = set(sample["cycle"].unique())
    assert cycles == {"1999-2000", "2001-2002"}
    counts = sample["cycle"].value_counts()
    assert counts.min() > 0.25 * len(sample), (
        "One cycle contributes almost nothing, which is what a wrong variable "
        "name looks like from the outside."
    )


# --- the mortality file layout ----------------------------------------------


def test_mortality_layout_parses_a_literal_record():
    """Two records, one death and one survivor, sliced from a literal.

    Written as text rather than a fixture file so the column positions are
    visible in the test that pins them. The survivor's line is deliberately
    short: the published files do not pad trailing blanks, and slicing an
    unpadded 46-character line truncates permth_exm to one digit.
    """
    text = (
        "00002         1100600                     177177\n"
        "00005         10   ..                     244244"
    )
    frame = nd.parse_mortality(text)
    assert list(frame["SEQN"]) == [2, 5]

    died, alive = frame.iloc[0], frame.iloc[1]
    assert died["eligstat"] == nd.ELIGIBLE
    assert died["mortstat"] == nd.DECEASED
    assert died["ucod_leading"] == "006"  # Alzheimer's disease
    assert died["permth_exm"] == 177

    assert alive["mortstat"] == 0
    assert alive["ucod_leading"] is None
    assert alive["permth_exm"] == 244


def test_mortality_layout_tolerates_unpadded_short_lines():
    short = "00005         10   ..                     244244"[:46]
    frame = nd.parse_mortality(short)
    # Truncated input must not silently yield a plausible small number.
    assert pd.isna(frame.loc[0, "permth_exm"]) or frame.loc[0, "permth_exm"] != 244


def test_every_leading_cause_code_has_a_label():
    assert set(nd.UCOD_LEADING_LABELS) == {f"{i:03d}" for i in range(1, 11)}


# --- the frame contract ------------------------------------------------------


def test_resilience_frame_excludes_every_clock_and_cell_fraction(sample):
    """The health state must not silently acquire twelve clock columns.

    ``ResilienceService`` infers its biomarker list by excluding known
    non-biomarker columns. This frame carries clocks, cell fractions, and
    mortality variables, so handing it over unnarrowed would put an epigenetic
    clock into the health state whose variance the CSD estimator measures.
    """
    narrowed = nd.resilience_frame(sample)
    assert list(narrowed.columns) == [
        "subject_id",
        "age",
        "sex",
        *nd.MARKERS,
        "survey_weight",
    ]
    for column in (*nd.CLOCK_COLUMNS, *nd.CELL_FRACTIONS, "died", "followup_years"):
        assert column not in narrowed.columns


def test_six_marker_health_state_matches_the_2017_2018_cohort():
    """The same six markers, so the resilience estimator sees the state it was
    validated on. A different set would make the two cohorts incomparable."""
    from geroquery.sources import nhanes

    assert nd.MARKERS == nhanes.MARKERS
    assert nd.PHENOAGE_MARKERS == nhanes.PHENOAGE_MARKERS


def test_age_acceleration_is_orthogonal_to_chronological_age(sample):
    """The residual definition, not the raw difference.

    A raw ``clock - age`` difference correlates with age by construction for
    every published clock, so any association it shows with an age-related
    outcome is partly that artefact.
    """
    acceleration = nd.age_acceleration(sample, "HorvathAge")
    ok = acceleration.notna()
    assert acceleration[ok].mean() == pytest.approx(0.0, abs=1e-8)
    assert np.corrcoef(sample.loc[ok, "age"], acceleration[ok])[0, 1] == pytest.approx(
        0.0, abs=1e-8
    )


def test_age_acceleration_refuses_clocks_not_measured_in_years():
    """DunedinPoAm is a rate and HorvathTelo is a length; neither has a
    meaningful difference from a chronological age."""
    frame = nd.load_sample()
    for clock in ("DunedinPoAm", "HorvathTelo"):
        assert clock not in nd.AGE_LIKE_CLOCKS
        with pytest.raises(Exception, match="not measured in years"):
            nd.age_acceleration(frame, clock)


def test_telomere_length_shortens_with_age(sample):
    """The sign is the check. DNAmTL rising with age would mean the join is
    scrambled or the column was misread."""
    ok = sample[["age", "HorvathTelo"]].notna().all(axis=1)
    assert np.corrcoef(sample.loc[ok, "age"], sample.loc[ok, "HorvathTelo"])[0, 1] < -0.2


def test_committed_sample_covers_the_full_column_contract(sample):
    for column in (
        "subject_id",
        "cycle",
        "age",
        "sex",
        "dnam_sex",
        "died",
        "followup_years",
        "survey_weight",
        *nd.MARKERS,
        *nd.CLOCK_COLUMNS,
    ):
        assert column in sample.columns, column
    assert sample["age"].min() >= nd.MIN_AGE
    assert sample["age"].max() <= nd.AGE_TOPCODE
    assert set(sample["died"].unique()) <= {0, 1}
    assert (sample["followup_years"] >= 0).all()


def test_caveats_name_the_age_floor_and_the_crp_assay_change():
    """The two things most likely to be forgotten when reading a result off
    this cohort."""
    joined = " ".join(nd.CAVEATS)
    assert "50" in joined and "dynamic range" in joined
    assert "high-sensitivity" in joined
    assert "perturb" in joined


def test_adapter_reports_which_mode_produced_the_frame(tmp_path):
    """A hazard ratio from the 600-row sample must never be reported as the
    cohort result, so the mode is returned rather than logged."""
    source = nd.NhanesDnamSource(data_dir=tmp_path)
    nd.write_sample(nd.load_sample(), tmp_path)
    frame, mode = source.crosslayer_frame(prefer_full=False)
    assert mode == "sample"
    assert len(frame) > 0


# --- provenance --------------------------------------------------------------


def test_every_new_artifact_is_pinned_with_a_digest_and_a_size():
    artifacts = [*NHANES_1999_2002.values(), *NHANES_MORTALITY.values(), NHANES_DNAM]
    assert len(artifacts) == 11
    for artifact in artifacts:
        assert len(artifact.sha256) == 64
        assert artifact.n_bytes > 0
        assert artifact.url.startswith("https://")
        assert artifact.attribution


@pytest.mark.live
def test_live_upstream_bytes_still_match_the_pinned_digests():
    """The only thing that notices NCHS re-issuing a file."""
    from geroquery.sources.fetch import fetch_artifact

    for artifact in (NHANES_DNAM, NHANES_MORTALITY["MORT_1999"]):
        path = fetch_artifact(artifact, allow_network=True)
        assert path.stat().st_size == artifact.n_bytes


@pytest.mark.live
def test_live_full_cohort_has_the_expected_shape():
    """Sanity numbers from the 2026-08-05 build. A drift here means either the
    upstream changed or the join did."""
    frame = nd.load_full(allow_network=True)
    assert len(frame) == 2517
    assert int(frame["died"].sum()) == 1350
    assert frame["age"].min() == 50.0
    assert nd.sex_discordance(frame)["n_discordant"] == 60
