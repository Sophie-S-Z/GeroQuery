"""GEO DataSets adapter: SOFT parsing, age bands, and contrast construction.

All offline. The rules encoded here decide which real samples become a signature
row, so they are tested against the actual level vocabulary the panel contains —
``"20 - 39 y"``, ``"E16.5"``, ``"P7"``, ``"aged (24 months)"``, ``"ad libitum"``
— rather than against tidy invented strings.
"""

from __future__ import annotations

import pytest

from geroquery.exceptions import SourceError
from geroquery.sources import geo
from geroquery.sources.manifest import GEO_AGING_PANEL

# A miniature but structurally faithful GDS: two age levels, a controlled
# variable with a control arm, a sex variable, sample sources, and a table with
# an annotated, an ambiguous, and an unannotated probe.
SOFT = """^DATASET = GDS9999
!dataset_title = Age effect on the widget
!dataset_platform = GPL999
!dataset_platform_organism = Homo sapiens
!dataset_sample_organism = Homo sapiens
!dataset_value_type = count
!dataset_sample_count = 8
!dataset_reference_series = GSE9999
!dataset_pubmed_id = 12345678
!dataset_update_date = Jan 01 2015
^SUBSET = GDS9999_1
!subset_dataset_id = GDS9999
!subset_description = 20 - 39 y
!subset_sample_id = GSM1,GSM2,GSM3,GSM4
!subset_type = age
^SUBSET = GDS9999_2
!subset_dataset_id = GDS9999
!subset_description = 60 - 79 y
!subset_sample_id = GSM5,GSM6,GSM7,GSM8
!subset_type = age
^SUBSET = GDS9999_3
!subset_dataset_id = GDS9999
!subset_description = control
!subset_sample_id = GSM1,GSM2,GSM3,GSM5,GSM6,GSM7
!subset_type = disease state
^SUBSET = GDS9999_4
!subset_dataset_id = GDS9999
!subset_description = widgetosis
!subset_sample_id = GSM4,GSM8
!subset_type = disease state
^SUBSET = GDS9999_5
!subset_dataset_id = GDS9999
!subset_description = male
!subset_sample_id = GSM1,GSM2,GSM3,GSM4,GSM5,GSM6,GSM7,GSM8
!subset_type = gender
^DATASET = GDS9999
#ID_REF = Platform reference identifier
#IDENTIFIER = identifier
#GSM1 = Value for GSM1: 26 year old; src: human frontal cortex
#GSM5 = Value for GSM5: 70 year old; src: human frontal cortex
!dataset_table_begin
ID_REF\tIDENTIFIER\tGSM1\tGSM2\tGSM3\tGSM4\tGSM5\tGSM6\tGSM7\tGSM8\tGene ID
p1\tAAA\t100\t110\t105\t108\t200\t210\t205\t208\t111
p2\tAAA\t10\t11\t10\t11\t20\t21\t20\t21\t111
p3\tBBB///CCC\t50\t51\t52\t53\t54\t55\t56\t57\t222///333
p4\tDDD\t70\t71\t72\t73\t74\t75\t76\t77\tnull
!dataset_table_end
"""


@pytest.fixture(scope="module")
def header() -> geo.GdsHeader:
    return geo.parse_header(SOFT)


# ---- URLs -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("accession", "shard"),
    [("GDS707", "GDSnnn"), ("GDS156", "GDSnnn"), ("GDS1803", "GDS1nnn"), ("GDS5226", "GDS5nnn")],
)
def test_gds_directory_shards_by_leading_digits(accession, shard):
    assert geo.gds_directory(accession) == shard
    assert f"/{shard}/{accession}/" in geo.soft_url(accession)


def test_soft_url_matches_the_pinned_manifest_urls():
    """The adapter and the manifest must agree, or a verified download would be
    a different file from the one the adapter goes looking for."""
    for accession, artifact in GEO_AGING_PANEL.items():
        assert geo.soft_url(accession) == artifact.url


def test_gds_directory_rejects_a_non_gds_accession():
    with pytest.raises(ValueError, match="Not a GDS accession"):
        geo.gds_directory("GSE1572")


# ---- age vocabulary -------------------------------------------------------


@pytest.mark.parametrize(
    ("description", "years"),
    [
        ("106 y", 106.0),
        ("20 - 39 y", 29.5),
        ("4 to 13 years", 8.5),
        ("6 mo", 0.5),
        ("24 m", 2.0),
        ("12 weeks", 12 / 52.18),
        ("day 14", 14 / 365.25),
        ("postnatal day 5", 5 / 365.25),
        ("P7", 7 / 365.25),
        ("aged (24 months)", 2.0),
    ],
)
def test_age_in_years_parses_the_real_level_vocabulary(description, years):
    assert geo.age_in_years(description) == pytest.approx(years, rel=1e-6)


@pytest.mark.parametrize("description", ["E16.5", "Embryonic day 14", "embryonic", "young", "old"])
def test_prenatal_and_qualitative_levels_have_no_numeric_age(description):
    assert geo.age_in_years(description) is None


@pytest.mark.parametrize(
    ("description", "species", "group"),
    [
        ("26 years", "human", "young"),
        ("70 years", "human", "old"),
        ("50 years", "human", None),  # deliberately in the gap between the bands
        ("10 years", "human", None),  # paediatric: excluded, not "young"
        ("2 mo", "mouse", "young"),
        ("24 mo", "mouse", "old"),
        ("12 mo", "mouse", None),
        ("day 14", "mouse", None),  # developmental, not young adult
        ("young", "human", "young"),
        ("aged", "mouse", "old"),
        ("adult", "mouse", None),
        ("middle-aged", "human", None),
    ],
)
def test_age_group_bands(description, species, group):
    assert geo.age_group(description, species) == group


@pytest.mark.parametrize(
    "description",
    [
        "control",
        "Control (healthy, age-matched)",
        "wild type",
        "ad libitum",
        "normal diet",
        "baseline",
        "non-smoker",
        "wild type control",
    ],
)
def test_control_levels_are_recognised(description):
    assert geo.is_control_level(description)


@pytest.mark.parametrize(
    "description",
    [
        "schizophrenia",
        "calorie restricted",
        "IQGAP2 knockout",
        "severe presbycusis",
        "4 h after 1st session",
        "smoker",
    ],
)
def test_non_control_levels_are_not_mistaken_for_controls(description):
    assert not geo.is_control_level(description)


# ---- parsing --------------------------------------------------------------


def test_parse_header_reads_metadata_subsets_and_sample_sources(header):
    assert header.accession == "GDS9999"
    assert header.species == "human"
    assert header.series_id == "GSE9999"
    assert header.pubmed_id == "12345678"
    assert header.sample_count == 8
    assert header.variables == ["age", "disease state", "gender"]
    assert len(header.by_variable("age")) == 2
    assert header.sample_sources["GSM5"] == "human frontal cortex"


def test_parse_table_keeps_annotation_and_coerces_sample_values(header):
    table = geo.parse_table(SOFT)
    assert list(table.index) == ["p1", "p2", "p3", "p4"]
    assert table.loc["p1", "GSM1"] == 100.0
    assert set(table.columns) >= {"IDENTIFIER", "Gene ID", "GSM1", "GSM8"}


def test_parse_table_refuses_text_with_no_table():
    with pytest.raises(SourceError, match="no dataset table"):
        geo.parse_table("^DATASET = GDS1\n!dataset_title = nothing here\n")


def test_probe_annotation_prefers_entrez_and_drops_ambiguous_probes():
    table = geo.parse_table(SOFT)
    keys, identifiers = geo.probe_annotation(table)
    # p3 measures two genes; p4 has no Entrez id so it falls back to its symbol.
    assert keys.to_dict() == {"p1": "111", "p2": "111", "p4": "SYM:DDD"}
    assert identifiers.loc["111", "symbol"] == "AAA"


# ---- contrast construction ------------------------------------------------


def test_build_contrasts_restricts_to_the_control_arm(header):
    contrasts, skipped = geo.build_contrasts(header)
    assert not skipped
    assert len(contrasts) == 1
    contrast = contrasts[0]
    assert contrast.young_samples == ("GSM1", "GSM2", "GSM3")  # GSM4 is a case
    assert contrast.old_samples == ("GSM5", "GSM6", "GSM7")  # GSM8 is a case
    assert contrast.restrictions == ("disease state = control",)
    assert contrast.tissue == "brain"
    assert contrast.sex == "male"
    assert contrast.age_range == "30y vs 70y"  # midpoints of the two level ranges
    assert contrast.study_id == "GEO:GDS9999"


def test_a_dataset_with_no_control_arm_is_skipped_with_a_reason():
    soft = SOFT.replace("!subset_description = control", "!subset_description = mild widgetosis")
    contrasts, skipped = geo.build_contrasts(geo.parse_header(soft))
    assert contrasts == []
    assert "no recognisable control level" in skipped[0]


def test_groups_below_the_minimum_are_skipped_with_counts():
    soft = SOFT.replace(
        "!subset_sample_id = GSM1,GSM2,GSM3,GSM5,GSM6,GSM7",
        "!subset_sample_id = GSM1,GSM5,GSM6,GSM7",
    )
    contrasts, skipped = geo.build_contrasts(geo.parse_header(soft))
    assert contrasts == []
    assert "young=1" in skipped[0]


def test_stratified_variables_produce_one_contrast_per_level_with_distinct_ids():
    soft = SOFT + ""
    soft = soft.replace(
        "^DATASET = GDS9999\n#ID_REF",
        "^SUBSET = GDS9999_6\n"
        "!subset_dataset_id = GDS9999\n"
        "!subset_description = left widget\n"
        "!subset_sample_id = GSM1,GSM2,GSM3,GSM5,GSM6,GSM7\n"
        "!subset_type = tissue\n"
        "^SUBSET = GDS9999_7\n"
        "!subset_dataset_id = GDS9999\n"
        "!subset_description = right widget\n"
        "!subset_sample_id = GSM4,GSM8\n"
        "!subset_type = tissue\n"
        "^DATASET = GDS9999\n#ID_REF",
    )
    contrasts, skipped = geo.build_contrasts(geo.parse_header(soft))
    assert [c.study_id for c in contrasts] == ["GEO:GDS9999:left_widget"]
    assert any("right widget" in reason for reason in skipped)


def test_sex_is_unspecified_when_the_dataset_declares_no_sex_variable():
    soft = SOFT.replace("!subset_type = gender", "!subset_type = other")
    contrasts, _ = geo.build_contrasts(geo.parse_header(soft))
    assert contrasts[0].sex == geo.SEX_UNSPECIFIED


def test_unmodelled_organism_is_skipped_rather_than_coerced():
    soft = SOFT.replace(
        "!dataset_sample_organism = Homo sapiens", "!dataset_sample_organism = Danio rerio"
    )
    contrasts, skipped = geo.build_contrasts(geo.parse_header(soft))
    assert contrasts == []
    assert "not modelled" in skipped[0]


# ---- tissue vocabulary ----------------------------------------------------


@pytest.mark.parametrize(
    ("text", "tissue"),
    [
        ("human frontal cortex", "brain"),
        ("vastus lateralis biopsy, 69 yr old male", "skeletal muscle"),
        ("bone marrow adipocyte", "adipose tissue"),
        ("hematopoietic stem cell (HSC)", "bone marrow"),
        ("monocyte-derived dendritic cells", "blood"),
        ("hepatocellular carcinoma at 6 months", "liver"),
        ("some entirely unfamiliar structure", None),
    ],
)
def test_normalize_tissue(text, tissue):
    assert geo.normalize_tissue(text) == tissue
