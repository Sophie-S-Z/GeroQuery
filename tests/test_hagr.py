"""HAGR parsers: GenAge, CellAge, LongevityMap, DrugAge, GenDR.

Offline, against literal excerpts of the real release files — including the one
malformed GenDR line that is actually in the current release.
"""

from __future__ import annotations

import pytest

from geroquery.sources import hagr

GENAGE_HUMAN = """GenAge ID,symbol,name,entrez gene id,uniprot,why
1,GHR,growth hormone receptor,2690,GHR_HUMAN,mammal
2,CDKN2A,cyclin dependent kinase inhibitor 2A,1029,P42771,"functional,putative"
3,,unnamed,999,,mammal
"""

GENAGE_MODELS = """GenAge ID,symbol,name,organism,entrez gene id,avg lifespan change (max obsv),\
lifespan effect,longevity influence
1,aak-2,AMP-Activated Kinase,Caenorhabditis elegans,181727,,Increase and Decrease,Pro-Longevity
2,Trp53,transformation related protein 53,Mus musculus,22059,20,Increase,Pro-Longevity
"""

CELLAGE = """Entrez ID\tGene symbol\tGene name\tCancer Cell\tType of senescence\t\
Senescence Effect\tReference
1029\tCDKN2A\tcyclin dependent kinase inhibitor 2A\tNo\tReplicative\tInduces\t26583757
"""

LONGEVITYMAP = """id,Association,Population,Variant(s),Gene(s),PubMed,
1,significant,Japanese,rs2802292,"FOXO3, FOXO3B",18765803,
2,non-significant,Dutch,HLA-B40,HLA-B,1859103,
"""

DRUGAGE = """compound_name,species,strain,dosage,age_at_initiation,treatment_duration,\
avg_lifespan_change_percent,avg_lifespan_significance,max_lifespan_change_percent,\
max_lifespan_significance,gender,weight_change_percent,weight_change_significance,ITP,pubmed_id
Rapamycin,Mus musculus,UM-HET3,14 ppm,,,10,S,,NA,Male,,NA,Yes,19587680
Rapamycin,Mus musculus,UM-HET3,14 ppm,,,20,S,,NA,Female,,NA,Yes,19587680
Rapamycin,Mus musculus,UM-HET3,4 ppm,,,3,NS,,NA,Male,,NA,No,23594965
Placebo,Mus musculus,UM-HET3,none,,,1,NS,,NA,Male,,NA,No,19587680
"""

# The current GenDR release genuinely contains a line like this: an unquoted
# field runs many "symbol;organism" pairs together on one row.
GENDR = """GenDR ID,gene symbol,species,entrez gene id,gene name
1,SIR2,Saccharomyces cerevisiae,851520,Silent Information Regulator 2
2,Rrm1,Mus musculus,20133,ribonucleotide reductase M1
3,His4:CG33871;Drosophila melanogaster,His4:CG33873;Drosophila melanogaster,,broken row
"""


# ---- gene tables ----------------------------------------------------------


def test_genage_human_carries_the_evidence_code_into_the_assertion():
    out = hagr.parse_genage_human(GENAGE_HUMAN)
    assert [a.symbol for a in out] == ["GHR", "CDKN2A"]  # the empty-symbol row is dropped
    assert out[0].entrez == "2690"
    assert out[0].species == "human"
    assert "mammal" in out[0].assertion
    assert "functional,putative" in out[1].assertion


def test_genage_models_leaves_non_mammalian_species_unresolved():
    out = hagr.parse_genage_models(GENAGE_MODELS)
    worm, mouse = out
    assert worm.organism == "Caenorhabditis elegans"
    assert worm.species is None  # parsed and counted, never labelled mammalian
    assert mouse.species == "mouse"
    assert "pro-longevity" in mouse.assertion
    assert "max observed change 20%" in mouse.assertion


def test_cellage_records_the_senescence_effect_and_links_to_pubmed():
    (entry,) = hagr.parse_cellage(CELLAGE)
    assert entry.database == "CellAge"
    assert entry.assertion == "induces senescence (replicative)"
    assert entry.url == "https://pubmed.ncbi.nlm.nih.gov/26583757/"


def test_longevitymap_splits_multi_gene_rows_and_keeps_null_results():
    out = hagr.parse_longevitymap(LONGEVITYMAP)
    assert [a.symbol for a in out] == ["FOXO3", "FOXO3B", "HLA-B"]
    assert "significant longevity association in Japanese" in out[0].assertion
    # The null result is retained: a list built only from hits cannot say that a
    # gene was tested and did not replicate.
    assert out[2].assertion.startswith("non-significant")


def test_gendr_drops_the_malformed_upstream_row():
    out = hagr.parse_gendr(GENDR)
    assert [a.symbol for a in out] == ["SIR2", "Rrm1"]
    assert all(";" not in a.organism for a in out)


@pytest.mark.parametrize(
    ("fields", "ok"),
    [
        (("CDKN2A",), True),
        (("Rrm1", "Mus musculus"), True),
        (("His4:CG33873;Drosophila melanogaster",), False),
        (("",), False),
        (("x" * 65,), False),
    ],
)
def test_is_well_formed(fields, ok):
    assert hagr.is_well_formed(*fields) is ok


# ---- interventions --------------------------------------------------------


def test_drugage_medians_only_the_significant_experiments():
    out = {i.intervention_id: i for i in hagr.parse_drugage(DRUGAGE)}
    rapamycin = out["DRUGAGE:rapamycin:mus_musculus"]
    # 10 and 20 are significant; the 3 is not and must not pull the median down.
    assert rapamycin.lifespan_effect_pct == pytest.approx(15.0)
    assert rapamycin.n_experiments == 3
    assert rapamycin.source == "ITP"  # at least one row is flagged as an ITP study


def test_drugage_reports_none_when_nothing_reached_significance():
    """None means "tested, no significant effect" — deliberately not 0.0, which
    would read as a measured null effect of exactly zero."""
    out = {i.intervention_id: i for i in hagr.parse_drugage(DRUGAGE)}
    placebo = out["DRUGAGE:placebo:mus_musculus"]
    assert placebo.lifespan_effect_pct is None
    assert placebo.n_experiments == 1
    assert placebo.source == "DrugAge"


def test_drug_interventions_never_claim_gene_targets():
    """DrugAge records no gene targets; asserting them by hand is exactly the
    fabricated edge this ingestion exists to remove."""
    assert all(i.gene_symbols == () for i in hagr.parse_drugage(DRUGAGE))


def test_gendr_interventions_link_the_genes_of_each_organism():
    out = {i.intervention_id: i for i in hagr.gendr_interventions(hagr.parse_gendr(GENDR))}
    assert set(out) == {
        "GENDR:dietary_restriction:saccharomyces_cerevisiae",
        "GENDR:dietary_restriction:mus_musculus",
    }
    mouse = out["GENDR:dietary_restriction:mus_musculus"]
    assert mouse.gene_symbols == ("Rrm1",)
    assert mouse.itype == "dietary"
    assert mouse.lifespan_effect_pct is None
