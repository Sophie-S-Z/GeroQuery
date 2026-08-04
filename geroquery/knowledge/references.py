"""Real, verifiable literature references used across the curated knowledge base.

Every reference is a real publication. Where we are confident of the PubMed ID
it is recorded and the link points straight at the article; otherwise the link
is a PubMed title search that resolves to the same paper. We never invent a
numeric identifier — an unverifiable PMID would itself be a form of false
evidence, which this project explicitly avoids.
"""

from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import quote_plus


@dataclass(frozen=True)
class Reference:
    key: str
    authors: str
    year: int
    journal: str
    title: str
    pmid: str | None = None

    @property
    def short(self) -> str:
        """A compact in-text citation, e.g. 'Krishnamurthy et al., 2004'."""
        return f"{self.authors}, {self.year}"

    @property
    def url(self) -> str:
        if self.pmid:
            return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"
        return f"https://pubmed.ncbi.nlm.nih.gov/?term={quote_plus(self.title)}"

    def to_dict(self) -> dict:
        return {
            "key": self.key,
            "authors": self.authors,
            "year": self.year,
            "journal": self.journal,
            "title": self.title,
            "pmid": self.pmid,
            "citation": f"{self.authors} ({self.year}). {self.title}. {self.journal}.",
            "url": self.url,
        }


_REFS: list[Reference] = [
    Reference(
        "lopez_otin_2013",
        "López-Otín et al.",
        2013,
        "Cell",
        "The hallmarks of aging",
        "23746838",
    ),
    Reference(
        "lopez_otin_2023",
        "López-Otín et al.",
        2023,
        "Cell",
        "Hallmarks of aging: An expanding universe",
        "36599349",
    ),
    Reference(
        "krishnamurthy_2004",
        "Krishnamurthy et al.",
        2004,
        "J Clin Invest",
        "Ink4a/Arf expression is a biomarker of aging",
        "15520862",
    ),
    Reference(
        "baker_2011",
        "Baker et al.",
        2011,
        "Nature",
        "Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders",
        "22048312",
    ),
    Reference(
        "freund_2012",
        "Freund et al.",
        2012,
        "Mol Biol Cell",
        "Lamin B1 loss is a senescence-associated biomarker",
        "22496421",
    ),
    Reference(
        "shimi_2011",
        "Shimi et al.",
        2011,
        "Genes Dev",
        "The role of nuclear lamin B1 in cell proliferation and senescence",
        "22156207",
    ),
    Reference(
        "tanaka_2018",
        "Tanaka et al.",
        2018,
        "Aging Cell",
        "Plasma proteomic signature of age in healthy humans",
        "29785772",
    ),
    Reference(
        "lehallier_2019",
        "Lehallier et al.",
        2019,
        "Nat Med",
        "Undulating changes in human plasma proteome profiles across the lifespan",
        "31806903",
    ),
    Reference(
        "willcox_2008",
        "Willcox et al.",
        2008,
        "PNAS",
        "FOXO3A genotype is strongly associated with human longevity",
        "18765803",
    ),
    Reference(
        "flachsbart_2009",
        "Flachsbart et al.",
        2009,
        "PNAS",
        "Association of FOXO3A variation with human longevity confirmed in German centenarians",
        "19196970",
    ),
    Reference(
        "suh_2008",
        "Suh et al.",
        2008,
        "PNAS",
        "Functionally significant insulin-like growth factor I receptor mutations in centenarians",
        "18316725",
    ),
    Reference(
        "harrison_2009",
        "Harrison et al.",
        2009,
        "Nature",
        "Rapamycin fed late in life extends lifespan in genetically heterogeneous mice",
        "19587680",
    ),
    Reference(
        "kuroo_1997",
        "Kuro-o et al.",
        1997,
        "Nature",
        "Mutation of the mouse klotho gene leads to a syndrome resembling ageing",
        "9363890",
    ),
    Reference(
        "kurosu_2005",
        "Kurosu et al.",
        2005,
        "Science",
        "Suppression of aging in mice by the hormone Klotho",
        "16123266",
    ),
    Reference(
        "junnila_2013",
        "Junnila et al.",
        2013,
        "Nat Rev Endocrinol",
        "The GH/IGF-1 axis in ageing and longevity",
        "23726003",
    ),
    Reference(
        "martin_montalvo_2013",
        "Martín-Montalvo et al.",
        2013,
        "Nat Commun",
        "Metformin improves healthspan and lifespan in mice",
        "23900241",
    ),
    Reference(
        "harrison_2014",
        "Harrison et al.",
        2014,
        "Aging Cell",
        "Acarbose, 17-alpha-estradiol, and nordihydroguaiaretic acid extend "
        "mouse lifespan preferentially in males",
        "24245565",
    ),
    Reference(
        "xu_2018",
        "Xu et al.",
        2018,
        "Nat Med",
        "Senolytics improve physical function and increase lifespan in old age",
        "29988130",
    ),
    Reference(
        "mills_2016",
        "Mills et al.",
        2016,
        "Cell Metab",
        "Long-term administration of nicotinamide mononucleotide mitigates "
        "age-associated physiological decline in mice",
        "28068222",
    ),
    Reference(
        "levine_2018",
        "Levine et al.",
        2018,
        "Aging (Albany NY)",
        "An epigenetic biomarker of aging for lifespan and healthspan",
        "29676998",
    ),
    Reference(
        "levine_2013",
        "Levine",
        2013,
        "J Gerontol A Biol Sci Med Sci",
        "Modeling the rate of senescence: can estimated biological age predict "
        "mortality more accurately than chronological age?",
        "23213031",
    ),
    Reference(
        "liu_2018",
        "Liu et al.",
        2018,
        "PLoS Med",
        "A new aging measure captures morbidity and mortality risk across "
        "diverse subpopulations from NHANES IV",
        "30130351",
    ),
    Reference(
        "scheffer_2009",
        "Scheffer et al.",
        2009,
        "Nature",
        "Early-warning signals for critical transitions",
        "19727193",
    ),
    Reference(
        "gijzel_2017",
        "Gijzel et al.",
        2017,
        "J Gerontol A Biol Sci Med Sci",
        "Dynamical resilience indicators in time series of self-rated health "
        "correspond to frailty and health in older adults",
        "28329147",
    ),
    Reference(
        "pyrkov_2021",
        "Pyrkov et al.",
        2021,
        "Nat Commun",
        "Longitudinal analysis of blood markers reveals progressive loss of "
        "resilience and predicts human lifespan limit",
        "34039909",
    ),
    Reference(
        "peters_2015",
        "Peters et al.",
        2015,
        "Nat Commun",
        "The transcriptional landscape of age in human peripheral blood",
        "26490707",
    ),
    Reference(
        "de_magalhaes_2009",
        "de Magalhães et al.",
        2009,
        "Aging Cell",
        "Meta-analysis of age-related gene expression profiles identifies "
        "common signatures of aging",
        "19419974",
    ),
    Reference(
        "tacutu_2018",
        "Tacutu et al.",
        2018,
        "Nucleic Acids Res",
        "Human Ageing Genomic Resources: new and updated databases",
        "29121237",
    ),
]

REFERENCES: dict[str, Reference] = {r.key: r for r in _REFS}


def ref(key: str) -> Reference:
    return REFERENCES[key]
