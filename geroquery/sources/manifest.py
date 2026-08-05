"""Pinned manifest of external data artifacts.

Every byte of real upstream data GeroQuery ingests is declared here first: exact
URL, release, expected SHA-256, expected size, and licence. Nothing is fetched
that is not in this table, and nothing that fails its checksum is ever used.

Why pin rather than just download: NHANES and friends silently re-publish files
under the same path. Without a checksum, an upstream edit becomes an unexplained
change in your results months later, with no way to tell which. With one, the
fetch fails loudly and you go look.

Checksums below were computed from the actual bytes served on the date in
``VERIFIED_ON``. If a fetch starts failing, the upstream file changed — confirm
against the CDC release notes and bump both ``sha256`` and ``VERIFIED_ON``
rather than deleting the check.
"""

from __future__ import annotations

from dataclasses import dataclass

# Bumped whenever any artifact below changes. Recorded alongside derived data so
# a result can be traced to the exact inputs that produced it.
MANIFEST_VERSION = "2026.3"

# Date the checksums were last confirmed against live upstream bytes.
VERIFIED_ON = "2026-08-05"


@dataclass(frozen=True)
class RemoteArtifact:
    """One pinned remote file."""

    key: str
    url: str
    sha256: str
    n_bytes: int
    release: str
    license: str
    attribution: str
    description: str

    @property
    def filename(self) -> str:
        """Basename used for the on-disk cache entry."""
        return self.url.rsplit("/", 1)[-1]


# --- NHANES 2017-2018 (release "J") ----------------------------------------
#
# URL pattern note, learned the hard way: the human-facing
# ``/Nchs/Nhanes/2017-2018/DEMO_J.XPT`` path returns an HTML page, not XPORT
# bytes, so ``pandas.read_sas`` fails with "Header record is not an XPORT file".
# The path below is the one that actually serves the data file.
_NHANES_BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"

_NHANES_LICENSE = "US public domain (NCHS)"
_NHANES_ATTRIBUTION = (
    "Centers for Disease Control and Prevention (CDC), National Center for Health "
    "Statistics (NCHS). National Health and Nutrition Examination Survey Data, "
    "2017-2018. Hyattsville, MD: U.S. Department of Health and Human Services."
)

NHANES_2017_2018: dict[str, RemoteArtifact] = {
    "DEMO_J": RemoteArtifact(
        key="DEMO_J",
        url=f"{_NHANES_BASE}/DEMO_J.xpt",
        sha256="c0b46e0345ea19404928656277c8b0d10b0cca348a9b2fe4fc3c67e8b7ee73ec",
        n_bytes=3412720,
        release="NHANES 2017-2018",
        license=_NHANES_LICENSE,
        attribution=_NHANES_ATTRIBUTION,
        description="Demographics: age (RIDAGEYR, topcoded at 80), sex (RIAGENDR).",
    ),
    "BIOPRO_J": RemoteArtifact(
        key="BIOPRO_J",
        url=f"{_NHANES_BASE}/BIOPRO_J.xpt",
        sha256="5bcd5722c1892883b96a9d7fed0befadabacae313f800fadc0133cd5dd00c4c6",
        n_bytes=2106080,
        release="NHANES 2017-2018",
        license=_NHANES_LICENSE,
        attribution=_NHANES_ATTRIBUTION,
        description="Standard biochemistry profile: albumin, creatinine, glucose.",
    ),
    "HSCRP_J": RemoteArtifact(
        key="HSCRP_J",
        url=f"{_NHANES_BASE}/HSCRP_J.xpt",
        sha256="19109fab0661ba23352bdf17dd9828bbd5dd8be34759888a2a64d9e33397bfe6",
        n_bytes=202000,
        release="NHANES 2017-2018",
        license=_NHANES_LICENSE,
        attribution=_NHANES_ATTRIBUTION,
        description="High-sensitivity C-reactive protein (LBXHSCRP, mg/L).",
    ),
    "CBC_J": RemoteArtifact(
        key="CBC_J",
        url=f"{_NHANES_BASE}/CBC_J.xpt",
        sha256="00f964098dc91272e415344554cbf1b627ecaeb2eda201518527446a7d81e742",
        n_bytes=1476320,
        release="NHANES 2017-2018",
        license=_NHANES_LICENSE,
        attribution=_NHANES_ATTRIBUTION,
        description="Complete blood count: lymphocyte percent, red cell distribution width.",
    ),
}

# --- NHANES 1999-2002: the cross-layer cohort -------------------------------
#
# 2017-2018 has the best clinical panel but no methylation. 1999-2002 is the
# only NHANES with both: NCHS assayed DNA methylation on a subsample of adults
# aged 50+ and published the derived clocks in July 2024. That makes these two
# cycles the only place where a clock, a health state, and a death are all
# observed on the same person, which is the whole reason this block exists.
#
# Two cycles, so two of everything. Variable names are NOT stable across them —
# creatinine is LBXSCR in 1999-2000 and LBDSCR in 2001-2002, alkaline
# phosphatase LBXSAPSI then LBDSAPSI. The per-cycle map lives in
# :mod:`geroquery.sources.nhanes_dnam`; getting it wrong yields a silently
# smaller cohort rather than an error, so it is pinned by a test.
_NHANES_1999_BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles"
_NHANES_2001_BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles"

_NHANES_1999_2002_ATTRIBUTION = (
    "Centers for Disease Control and Prevention (CDC), National Center for Health "
    "Statistics (NCHS). National Health and Nutrition Examination Survey Data, "
    "1999-2002. Hyattsville, MD: U.S. Department of Health and Human Services."
)

# key, url base, filename, sha256, bytes, description
_NHANES_1999_2002_RAW: tuple[tuple[str, str, str, str, int, str], ...] = (
    (
        "DEMO_1999",
        _NHANES_1999_BASE,
        "DEMO.xpt",
        "a17bd1ebfd6007b8dc93f16259cb090678eb6c467a39363b7e932bba15febf89",
        11500560,
        "Demographics 1999-2000: age (RIDAGEYR, topcoded at 85), sex, 4-year MEC weight.",
    ),
    (
        "LAB18_1999",
        _NHANES_1999_BASE,
        "LAB18.xpt",
        "36dc6f02045ceb61e146885bb04400fd580c9c35cc3bfc4269e314f0e06fe9c0",
        2223120,
        "Biochemistry profile 1999-2000: albumin, creatinine (LBXSCR), glucose, ALP.",
    ),
    (
        "LAB25_1999",
        _NHANES_1999_BASE,
        "LAB25.xpt",
        "f69ef53e50fb3eea6a9e03b0e1167f201ecdb6ad485793002f4f7bee100fbd45",
        1487520,
        "Complete blood count 1999-2000: lymphocyte percent, RDW, MCV, WBC.",
    ),
    (
        "LAB11_1999",
        _NHANES_1999_BASE,
        "LAB11.xpt",
        "9ed9904dfc804b50c2c08b763b43d7df702913fc010c0f877e10f9f195ee183a",
        469040,
        "C-reactive protein 1999-2000 (LBXCRP, mg/dL) - standard assay, not hs-CRP.",
    ),
    (
        "DEMO_2001",
        _NHANES_2001_BASE,
        "DEMO_B.xpt",
        "6458dc2307cee244ab6ba26f70ec78c7b842ffc68e668bc43f0ab6b9737937d9",
        3273520,
        "Demographics 2001-2002: age, sex, 4-year MEC weight.",
    ),
    (
        "L40_2001",
        _NHANES_2001_BASE,
        "L40_B.xpt",
        "16f907f3a945cc24d053a9492a3f90edd3e0335f70ba047bbd9952ba6a00e182",
        2448480,
        "Biochemistry profile 2001-2002: albumin, creatinine (LBDSCR), glucose, ALP.",
    ),
    (
        "L25_2001",
        _NHANES_2001_BASE,
        "L25_B.xpt",
        "4bce7e5622c4b79b9965d7869875e05dae547c4b222fd89db47b8efb430b85f3",
        1671760,
        "Complete blood count 2001-2002: lymphocyte percent, RDW, MCV, WBC.",
    ),
    (
        "L11_2001",
        _NHANES_2001_BASE,
        "L11_B.xpt",
        "1a78e83a59aaa4a04aeb8ed7abeaeb4b87e9950601ac428c52fae9ce53e9bbbc",
        446240,
        "C-reactive protein 2001-2002 (LBXCRP, mg/dL) - standard assay, not hs-CRP.",
    ),
)

NHANES_1999_2002: dict[str, RemoteArtifact] = {
    key: RemoteArtifact(
        key=key,
        url=f"{base}/{filename}",
        sha256=sha256,
        n_bytes=n_bytes,
        release="NHANES 1999-2002",
        license=_NHANES_LICENSE,
        attribution=_NHANES_1999_2002_ATTRIBUTION,
        description=description,
    )
    for key, base, filename, sha256, n_bytes, description in _NHANES_1999_2002_RAW
}

# --- NHANES DNA methylation epigenetic biomarkers ---------------------------
#
# Published 2024-07-31. 4,449 rows, of which 2,532 carry measurements: adults
# aged 50+ from the 1999-2000 and 2001-2002 cycles, whole blood, EPIC array.
#
# What is public is the *derived* layer — NCHS ran the clocks themselves and
# released the per-participant outputs (12 clocks, 6 cell fractions, the GrimAge
# component predictors, and a 4-year DNAm survey weight). The CpG-level betas
# are Research Data Center only.
#
# That split is a feature here, not a limitation. It means this path needs no
# biolearn, no pyaging, no torch, and no second interpreter — and it sidesteps
# the normalization mismatch documented in docs/RESULTS_METHYLATION_CLOCKS.md,
# because these are the survey's own batch-corrected values rather than ours.
NHANES_DNAM = RemoteArtifact(
    key="nhanes_dnam",
    url="https://wwwn.cdc.gov/nchs/data/nhanes/dnam/dnmepi.sas7bdat",
    sha256="583e0660eafac5b4b2600a18bc4ca6e8a342fd1e87a78c1b03063b14223d956f",
    n_bytes=1245184,
    release="NHANES 1999-2002 DNAm epigenetic biomarkers (released 2024-07-31)",
    license=_NHANES_LICENSE,
    attribution=_NHANES_1999_2002_ATTRIBUTION,
    description=(
        "Per-participant DNAm clocks (Horvath, Hannum, SkinBlood, PhenoAge, GrimAge, "
        "GrimAge2, DunedinPoAm, HorvathTelo, Zhang, Lin, Weidner, VidalBralo), "
        "cell-type fractions, and the DNAm 4-year survey weight."
    ),
)

# --- NCHS public-use linked mortality files ---------------------------------
#
# Follow-up through 2019-12-31, linked to the survey by SEQN. This is the hard
# outcome the repo has never had: every result until now has been a measurement
# validated against another measurement.
#
# Public-use caveat that has to be carried into any result computed from these:
# NCHS perturbs a subset of records to prevent re-identification, substituting
# synthetic follow-up time or cause of death. Aggregate estimates are designed
# to survive that; individual records are not trustworthy on their own.
_MORTALITY_BASE = "https://ftp.cdc.gov/pub/HEALTH_STATISTICS/NCHS/datalinkage/linked_mortality"
_MORTALITY_ATTRIBUTION = (
    "National Center for Health Statistics. NCHS Data Linked to NDI Mortality Files "
    "(public-use, follow-up through December 31, 2019). Hyattsville, MD."
)

NHANES_MORTALITY: dict[str, RemoteArtifact] = {
    key: RemoteArtifact(
        key=key,
        url=f"{_MORTALITY_BASE}/{filename}",
        sha256=sha256,
        n_bytes=n_bytes,
        release="NCHS 2019 public-use linked mortality file",
        license=_NHANES_LICENSE,
        attribution=_MORTALITY_ATTRIBUTION,
        description=description,
    )
    for key, filename, sha256, n_bytes, description in (
        (
            "MORT_1999",
            "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
            "562bd367107add0b2b40fed36a062d8b194fca31e11aab897a0efad4bebc8b26",
            487666,
            "Mortality follow-up for NHANES 1999-2000 (9,965 records).",
        ),
        (
            "MORT_2001",
            "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
            "2788c8c8a02995ea5686882b5fe65b4b9cd415016031525f8a3c1a3646c47f1b",
            540362,
            "Mortality follow-up for NHANES 2001-2002 (11,039 records).",
        ),
    )
}

# --- NCBI GEO DataSets: the aging-signature panel ---------------------------
#
# Each entry is a GDS "full SOFT" file: the curated value matrix plus the joined
# platform annotation, so probes arrive carrying a gene symbol and an Entrez id.
#
# The panel is not hand-picked. It is every GDS that declares an ``age`` subset
# variable (a GEO DataSets query, 189 records) and survives the contrast rules in
# :mod:`geroquery.sources.geo` — an adult young-vs-old split in a modelled
# species, restricted to the control arm of every other subset variable, with at
# least three samples per group. Whatever passes, ships; the panel changes only
# by changing a rule, not by preferring a result.
#
# 27 distinct GEO Series are represented by 31 datasets: GEO sometimes splits one
# experiment across two array halves (GDS287/GDS288, GDS355/GDS356,
# GDS2961/GDS2962, GDS472/GDS473). Those pairs share subjects, so ``series_id``
# is carried onto every Study row and the meta-analysis caveat is documented
# rather than left for a reader to discover.
_GEO_PANEL: tuple[tuple[str, str, int, str, str, str, str, str | None, int, str], ...] = (
    # accession, sha256, bytes, GEO update date, organism, platform, series, PMID, n, title
    (
        "GDS156",
        "e067441c5ba7724c63e0135e628ec7b9e11d319f5267b77e1ec6baa7d48f550e",
        3469683,
        "Apr 17 2003",
        "Homo sapiens",
        "GPL91",
        "GSE80",
        "12204100",
        12,
        "Muscle function and aging (HG-U95A)",
    ),
    (
        "GDS1803",
        "dcf9c413699f890a53f8b3c65fd7a1d47ce14c2fbf961f356c006109cd2fe16d",
        9470248,
        "Oct 13 2006",
        "Mus musculus",
        "GPL1261",
        "GSE4332",
        "15967997",
        8,
        "Age effect on hematopoietic stem cells",
    ),
    (
        "GDS2019",
        "6fc4223d7ef28784de0f8b62766f05d907f1e57f7b4230bc120acf6ff45b5733",
        3592679,
        "Oct 10 2012",
        "Mus musculus",
        "GPL81",
        "GSE3129",
        "19943135",
        14,
        "Age effect on livers of long-lived Snell dwarf mutants",
    ),
    (
        "GDS2612",
        "147e17d9339c8b051a4e7887d61548dd9ca854c279add316fbf7ea155571ae13",
        6872128,
        "Apr 10 2007",
        "Mus musculus",
        "GPL339",
        "GSE6323",
        "17381838",
        15,
        "Caloric restriction effect on aged skeletal muscle",
    ),
    (
        "GDS287",
        "0485d058ca305fad7a98475a0e070b1ac3b64a291b53e6c08c3e33ffac89a11a",
        5627183,
        "Nov 02 2004",
        "Homo sapiens",
        "GPL96",
        "GSE362",
        "12783983",
        15,
        "Muscle function and aging - male (HG-U133A)",
    ),
    (
        "GDS288",
        "2c48f25653894b89b86b0c227759ad8c8e5a35f786205442c0d025d51992e330",
        4044281,
        "Nov 02 2004",
        "Homo sapiens",
        "GPL97",
        "GSE362",
        "12783983",
        15,
        "Muscle function and aging - male (HG-U133B)",
    ),
    (
        "GDS2929",
        "4c650c797d0b9c6aa0a73cb97be771fe36e7a5539b230c5e03b6498cc945c72b",
        10570095,
        "Jun 25 2008",
        "Mus musculus",
        "GPL1261",
        "GSE6591",
        "17726092",
        15,
        "Aging lungs and genetic background",
    ),
    (
        "GDS2961",
        "6ec78687fca2ff14fda70db99366a3d41e0cf9a238385545bb23695898506d2d",
        2481060,
        "Jan 04 2008",
        "Mus musculus",
        "GPL738",
        "GSE7829",
        "17499630",
        67,
        "Male and female thymi response to aging and caloric restriction (A)",
    ),
    (
        "GDS2962",
        "bec2419eb45cdd206150a7677da7e4a8d6d231604348328b62ea304f09445a41",
        2512465,
        "May 14 2015",
        "Mus musculus",
        "GPL782",
        "GSE7829",
        "17499630",
        67,
        "Male and female thymi response to aging and caloric restriction (B)",
    ),
    (
        "GDS2972",
        "01843e36b0e71d89730ef22db20fc72c30d2a122bd3d287982db948fba157ceb",
        3888578,
        "Oct 10 2012",
        "Mus musculus",
        "GPL81",
        "GSE8146",
        "17316780",
        20,
        "Vitamin E supplementation effect on aged heart",
    ),
    (
        "GDS2973",
        "5256bbd9e8ccb1c677ec459ac10e2301c48be6137aab2885169aa325c28e0a75",
        11354648,
        "Oct 03 2008",
        "Mus musculus",
        "GPL1261",
        "GSE8150",
        "17316780",
        20,
        "Vitamin E supplementation effect on aged brain",
    ),
    (
        "GDS3182",
        "84b83a338f656476e8e3b73e092a3381e16f3d5bd96ce052f01394eaacc44254",
        15837748,
        "Aug 22 2008",
        "Homo sapiens",
        "GPL570",
        "GSE9103",
        None,
        37,
        "Young and aged skeletal muscles response to long-term vigorous endurance exercise",
    ),
    (
        "GDS355",
        "4ae0e8031fae67bc06f110107314ec9ad250c3ca0d40146facf5f7da04bf1828",
        1757299,
        "Jun 23 2003",
        "Mus musculus",
        "GPL75",
        "GSE459",
        None,
        15,
        "Calorie restriction and aging (Mu11K-A)",
    ),
    (
        "GDS356",
        "5eec399a40fc8250022ee44155a17eea4f0a4e7d8f9ee5e371e35aee29b8c126",
        889138,
        "Jun 23 2003",
        "Mus musculus",
        "GPL76",
        "GSE459",
        None,
        15,
        "Calorie restriction and aging (Mu11K-B)",
    ),
    (
        "GDS3942",
        "52084f53058f07a433244c060aa6a9d13064aa0ded83fea9b916f3cc590d4b16",
        14763680,
        "Dec 20 2011",
        "Homo sapiens",
        "GPL570",
        "GSE32719",
        "22123971",
        27,
        "Aging effect on bone marrow hematopoietic stem cells",
    ),
    (
        "GDS3976",
        "8cd7cb310826330741a12c0b8e833d6aaa7ff32641f885829dde583fb777b05c",
        11885562,
        "Jun 01 2015",
        "Mus musculus",
        "GPL1261",
        "GSE27686",
        "21549326",
        16,
        "Premature Hematopoietic Aging model: bone marrow hematopoietic stem cells",
    ),
    (
        "GDS4522",
        "7c81a971c5c709c841fea1303c9d3e13588ad917f26d43abe69abd351cf870f0",
        17037828,
        "Sep 27 2013",
        "Homo sapiens",
        "GPL570",
        "GSE21935",
        "21538462",
        42,
        "Schizophrenia: postmortem superior temporal cortex",
    ),
    (
        "GDS4523",
        "3ee7cd36c540e2e7ddbe0405a068f31746fe11b14396cc02011c198508ea06f7",
        17291515,
        "Nov 05 2013",
        "Homo sapiens",
        "GPL570",
        "GSE17612",
        "19255580",
        51,
        "Schizophrenia: postmortem anterior prefrontal cortex",
    ),
    (
        "GDS472",
        "f5986ad188a97bbab0542aac79b89ec6038b90410c77e7c9de2cd333d4545d51",
        5349269,
        "Oct 10 2012",
        "Homo sapiens",
        "GPL96",
        "GSE674",
        "15036396",
        15,
        "Muscle function and aging - female (HG-U133A)",
    ),
    (
        "GDS473",
        "4a4883b2100cc2e4e5c99a0fa567d13f6f10587907bdead946acbbd437f092dc",
        3774287,
        "Oct 10 2012",
        "Homo sapiens",
        "GPL97",
        "GSE674",
        "15036396",
        15,
        "Muscle function and aging - female (HG-U133B)",
    ),
    (
        "GDS4858",
        "0f6676f9a6496ec3c3b786c68c711e6a0475e0f10ac46902a1b640333a58ce6d",
        13603848,
        "Apr 22 2014",
        "Homo sapiens",
        "GPL570",
        "GSE38718",
        "23418191",
        22,
        "Skeletal muscles from men and women of various ages",
    ),
    (
        "GDS4874",
        "c51223f10e286788071dbbf71019c25cd7d09ee3008c2257377d6a0b18a3a987",
        10902128,
        "Feb 26 2016",
        "Mus musculus",
        "GPL1261",
        "GSE46646",
        "23951254",
        12,
        "IQGAP2 knockout model of hepatocellular carcinoma: time course",
    ),
    (
        "GDS4892",
        "1ae95f27aded7c19e15b069e5d3b2fda70bb253b6be22889f43b22922b8f9330",
        11475844,
        "Jul 08 2014",
        "Mus musculus",
        "GPL1261",
        "GSE50821",
        "24797481",
        14,
        "Age effect on skeletal muscle precursor cells",
    ),
    (
        "GDS4904",
        "c2736a75a19da1e6fa51dc4a297fce471f51ccf068fa32405651c1db89b469f6",
        10890228,
        "Jul 02 2014",
        "Mus musculus",
        "GPL1261",
        "GSE52550",
        "24280126",
        12,
        "PGC-1alpha deficiency effect on aged gastrocnemius muscle",
    ),
    (
        "GDS5204",
        "c1db88fc831b1466ec1d3774110d710b32356aa3503ebb27ded3c3d0206b9d57",
        17518950,
        "Oct 01 2014",
        "Homo sapiens",
        "GPL570",
        "GSE53890",
        "24670762",
        41,
        "Age effect on normal adult brain: frontal cortical region",
    ),
    (
        "GDS5216",
        "0d6255849fb4de2151b8e25b75d29fa9f34d1b77caa66895855259f26179d4eb",
        17412444,
        "Sep 17 2014",
        "Homo sapiens",
        "GPL570",
        "GSE25941",
        "22302958",
        36,
        "Age effect on the skeletal muscle",
    ),
    (
        "GDS5217",
        "e795e1fd37104f821ce4123d4c4f237741df32af3d89b245e6f6db679e904b96",
        24886302,
        "Sep 17 2014",
        "Homo sapiens",
        "GPL570",
        "GSE28392",
        "22302958",
        70,
        "Resistance exercise effect on MHC I and MHC IIa muscle fibers of young and old women",
    ),
    (
        "GDS5218",
        "9a81dc8f321fd003f0f7e1e4b008cb109c0453682b469080f80e3789953b3bb8",
        33472006,
        "Jun 08 2015",
        "Homo sapiens",
        "GPL570",
        "GSE28422",
        "22302958",
        110,
        "Resistance exercise effect on skeletal muscles of young and old adults",
    ),
    (
        "GDS5226",
        "8ed1a68a75ee5e292b5cbad79fcf5b87421eab777426076f01f857151da038b2",
        10075204,
        "Jun 08 2015",
        "Mus musculus",
        "GPL6246",
        "GSE25905",
        "21545734",
        18,
        "Age effect on adipocytes of the bone marrow and epididymis",
    ),
    (
        "GDS5286",
        "7130ebe3217d252ad5a586fe0ead49c0e56bcde6c75d8ed64c6f571afc12292d",
        5580710,
        "Oct 22 2014",
        "Homo sapiens",
        "GPL571",
        "GSE58015",
        "25191744",
        9,
        "Monocyte-derived dendritic cells from young and aged donors",
    ),
    (
        "GDS707",
        "a04e662be8d6c94edbe79792dc6a66a4d7d1d3f297c0b25a6d8d90861d4986cc",
        4756323,
        "Oct 10 2012",
        "Homo sapiens",
        "GPL8300",
        "GSE1572",
        "15190254",
        30,
        "Aging brain: frontal cortex expression profiles at various ages",
    ),
)

GEO_LICENSE = "NCBI GEO - US public domain, attribute"
GEO_ATTRIBUTION = (
    "National Center for Biotechnology Information (NCBI), Gene Expression Omnibus (GEO). "
    "Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and "
    "hybridization array data repository. Nucleic Acids Res. 2002;30(1):207-10."
)


def _gds_url(accession: str) -> str:
    """FTP location of a GDS full SOFT file. ``GDS707`` shards under ``GDSnnn``."""
    digits = accession[len("GDS") :]
    shard = f"GDS{digits[:-3]}nnn" if len(digits) > 3 else "GDSnnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/datasets/{shard}/{accession}/soft/{accession}_full.soft.gz"


GEO_AGING_PANEL: dict[str, RemoteArtifact] = {
    accession: RemoteArtifact(
        key=accession,
        url=_gds_url(accession),
        sha256=sha256,
        n_bytes=n_bytes,
        release=f"GEO DataSets, curated {updated}",
        license=GEO_LICENSE,
        attribution=GEO_ATTRIBUTION,
        description=f"{title} [{organism}, {platform}, {series}, n={samples}]"
        + (f" PMID:{pubmed}" if pubmed else ""),
    )
    for accession, sha256, n_bytes, updated, organism, platform, series, pubmed, samples, title in (
        _GEO_PANEL
    )
}

# Series id per dataset, for the non-independence caveat above.
GEO_SERIES_BY_ACCESSION: dict[str, str] = {row[0]: row[6] for row in _GEO_PANEL}
GEO_PUBMED_BY_ACCESSION: dict[str, str | None] = {row[0]: row[7] for row in _GEO_PANEL}


# --- GEO DNA-methylation series: the clock validation panel -----------------
#
# Two hand-verified 450K blood series. These exist so the 236 wired aging clocks
# have real data to run on: before them, "236 clocks" meant 236 objects that had
# been shown not to crash.
#
# Unlike the expression panel these are GEO *Series*, not DataSets — no GDS
# exists for 450K methylation. Free-text age parsing is therefore unavoidable,
# which is acceptable for two hand-checked series in a way it would not be for
# thirty; ``sources/methylation.py`` validates each series' characteristic keys
# at parse time so an upstream format change fails loudly.
_METHYLATION_RAW: tuple[tuple[str, str, int, str, str], ...] = (
    (
        "GSE64495",
        "30114290d102c8b8f78110f9b5341580614c03b3c933389cbb71d456decb69f9",
        275855955,
        "Illumina 450K, whole blood, n=113, ages 0-94",
        "Primary clock validation set. Ships the authors' own Horvath-clock output "
        "per sample, so a wrapper can be checked against a published number rather "
        "than only against chronological age.",
    ),
    (
        "GSE30870",
        "86f1f7e6144bf398429ad5c528a16ff7bee56a7291893ec7861d243822644e30",
        82331350,
        "Illumina 450K, whole blood and cord blood, n=40",
        "Newborns vs nonagenarians. An extreme-contrast check: a clock that cannot "
        "separate cord blood from 90-year-old blood is broken in a way no "
        "correlation coefficient would hide.",
    ),
)


def _gse_matrix_url(accession: str) -> str:
    """FTP location of a GEO Series matrix. ``GSE30870`` shards under ``GSE30nnn``."""
    digits = accession[len("GSE") :]
    shard = f"GSE{digits[:-3]}nnn" if len(digits) > 3 else "GSEnnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{shard}/{accession}"
        f"/matrix/{accession}_series_matrix.txt.gz"
    )


METHYLATION_PANEL: dict[str, RemoteArtifact] = {
    accession: RemoteArtifact(
        key=accession,
        url=_gse_matrix_url(accession),
        sha256=sha256,
        n_bytes=n_bytes,
        release=f"GEO Series matrix ({platform})",
        license=GEO_LICENSE,
        attribution=GEO_ATTRIBUTION,
        description=description,
    )
    for accession, sha256, n_bytes, platform, description in _METHYLATION_RAW
}

# --- HAGR curated aging knowledge ------------------------------------------
#
# The Human Ageing Genomic Resources set. These replace what used to be a
# twelve-gene hand-written table: real curated assertions, with the same
# checksum discipline as everything else. Each artifact is a small zip; the
# member named below is the one that is parsed.
_HAGR_LICENSE = (
    "HAGR - free for non-commercial use, attribution required "
    "(https://genomics.senescence.info/legal.html)"
)
_HAGR_ATTRIBUTION = (
    "Human Ageing Genomic Resources (HAGR), de Magalhaes JP et al. "
    "https://genomics.senescence.info/"
)

# key, url, member file inside the zip, sha256, bytes, description
_HAGR_RAW: tuple[tuple[str, str, str, str, int, str], ...] = (
    (
        "genage_human",
        "https://genomics.senescence.info/genes/human_genes.zip",
        "genage_human.csv",
        "cd4b1d63c9a3fb9574da8453487ad528622b2ae70b6afcde0c8b372b08c52c7b",
        9465,
        "GenAge human: genes with evidence of a role in human ageing.",
    ),
    (
        "genage_models",
        "https://genomics.senescence.info/genes/models_genes.zip",
        "genage_models.csv",
        "12c08f92de4cb464894bb739a1926720213c1b6850533f773d8e2baa0a437d1c",
        48796,
        "GenAge model organisms: genes whose manipulation changes lifespan.",
    ),
    (
        "cellage",
        "https://genomics.senescence.info/cells/cellAge.zip",
        "cellage3.tsv",
        "943506acb3715a8d507459d69bc9a64c3e7863c809c62457171de03dc0062845",
        21540,
        "CellAge: genes experimentally shown to induce, inhibit, or regulate senescence.",
    ),
    (
        "longevitymap",
        "https://genomics.senescence.info/longevity/longevity_genes.zip",
        "longevity.csv",
        "5226be2961140bf1b7776a23b7eb7aede4c7a5e896f3bc21f33e1ddd30625763",
        24769,
        "LongevityMap: human genetic variants tested for association with longevity.",
    ),
    (
        "drugage",
        "https://genomics.senescence.info/drugs/dataset.zip",
        "drugage.csv",
        "d12f59717de60a207748f53bdbcdb484ed07e30d8608eb8e4a101e2a08d7fa80",
        45365,
        "DrugAge: compounds tested for an effect on model-organism lifespan.",
    ),
    (
        "gendr",
        "https://genomics.senescence.info/diet/dataset.zip",
        "gendr_manipulations.csv",
        "b6480eddc4c1e033d21d92e057b9277da373b817588f897844e47ef331433d66",
        8209,
        "GenDR: genes required for the lifespan extension of dietary restriction.",
    ),
)

HAGR: dict[str, RemoteArtifact] = {
    key: RemoteArtifact(
        key=key,
        url=url,
        sha256=sha256,
        n_bytes=n_bytes,
        release="HAGR current release",
        license=_HAGR_LICENSE,
        attribution=_HAGR_ATTRIBUTION,
        description=description,
    )
    for key, url, _member, sha256, n_bytes, description in _HAGR_RAW
}

# Which file inside each zip carries the data.
HAGR_MEMBERS: dict[str, str] = {key: member for key, _u, member, _s, _n, _d in _HAGR_RAW}


# AnAge: maximum-lifespan estimates, used to normalize chronological age to
# fractional lifespan so an "old mouse" and an "old human" are comparable.
# Ingested rather than transcribed: the nine numbers that used to live in
# idmap/data/anage.json were correct, but nothing recorded which AnAge release
# they came from or when they would go stale.
ANAGE = RemoteArtifact(
    key="anage",
    url="https://genomics.senescence.info/species/dataset.zip",
    sha256="e3ddb66e32e973a79932859ba53013e8f60d957c6ec01c6eb573e3ea3018d630",
    n_bytes=169150,
    release="HAGR AnAge current release",
    license=_HAGR_LICENSE,
    attribution=_HAGR_ATTRIBUTION,
    description="AnAge: maximum longevity and life-history traits across species.",
)

# Flat lookup across every collection, for the generic fetch CLI.
MANIFEST: dict[str, RemoteArtifact] = {
    **NHANES_2017_2018,
    **NHANES_1999_2002,
    **NHANES_MORTALITY,
    **GEO_AGING_PANEL,
    **METHYLATION_PANEL,
    **HAGR,
    "anage": ANAGE,
    "nhanes_dnam": NHANES_DNAM,
}


def get_artifact(key: str) -> RemoteArtifact:
    """Look up a pinned artifact, with a listing of valid keys on failure."""
    try:
        return MANIFEST[key]
    except KeyError:
        raise KeyError(f"Unknown artifact {key!r}. Known artifacts: {sorted(MANIFEST)}") from None
