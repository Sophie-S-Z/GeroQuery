"""M1 idmap — resolution correctness against known ground truth."""

import pytest

from geroquery.exceptions import GeneNotFoundError, TissueNotFoundError, UnknownSpeciesError
from geroquery.idmap import get_resolver

CDKN2A_HUMAN = "ENSG00000147889"
CDKN2A_MOUSE = "ENSMUSG00000044303"


@pytest.fixture(scope="module")
def r():
    return get_resolver()


@pytest.mark.parametrize(
    "query",
    ["CDKN2A", "cdkn2a", "ENSG00000147889", "1029", "P42771", "p16", "INK4A"],
)
def test_symbol_alias_ensembl_entrez_uniprot_round_trip(r, query):
    g = r.resolve_gene(query, species="human")
    assert g.canonical_id == CDKN2A_HUMAN
    assert g.symbol == "CDKN2A"
    assert g.matched_from == query


def test_deprecated_id_redirects(r):
    assert r.resolve_gene("ENSG00000147888").canonical_id == CDKN2A_HUMAN


def test_cross_species_ambiguity_is_surfaced_not_dropped(r):
    g = r.resolve_gene("CDKN2A")  # no species
    assert g.canonical_id == CDKN2A_HUMAN  # human preferred, deterministic
    assert CDKN2A_MOUSE in g.ambiguous_candidates


def test_species_filter_scopes_result(r):
    assert r.resolve_gene("CDKN2A", species="mouse").canonical_id == CDKN2A_MOUSE


def test_batch_resolution_marks_unknowns_none(r):
    out = r.resolve_batch(["TP53", "___nope___", "SIRT1"], species="human")
    assert out["TP53"] is not None and out["SIRT1"] is not None
    assert out["___nope___"] is None


def test_orthologs_span_species(r):
    ids = {g.canonical_id for g in r.orthologs(CDKN2A_HUMAN)}
    assert {CDKN2A_HUMAN, CDKN2A_MOUSE} <= ids


def test_unknown_gene_raises(r):
    with pytest.raises(GeneNotFoundError):
        r.resolve_gene("definitely_not_a_gene")


@pytest.mark.parametrize(
    "age,species,expected",
    [(61.25, "human", 0.5), (2.0, "mouse", 0.5), (1.9, "rat", 0.5)],
)
def test_fractional_age_math(r, age, species, expected):
    assert r.fractional_age(age, species) == pytest.approx(expected, abs=0.02)


def test_fractional_age_unknown_species_raises(r):
    with pytest.raises(UnknownSpeciesError):
        r.fractional_age(50, "griffin")


@pytest.mark.parametrize(
    "label,expected",
    [("whole blood", "UBERON:0000178"), ("liver", "UBERON:0002107"), ("plasma", "UBERON:0001969")],
)
def test_tissue_mapping(r, label, expected):
    assert r.map_tissue(label).term_id == expected


def test_unknown_tissue_raises(r):
    with pytest.raises(TissueNotFoundError):
        r.map_tissue("zzzznotarealtissue")
