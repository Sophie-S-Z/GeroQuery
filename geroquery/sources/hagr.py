"""Human Ageing Genomic Resources (HAGR) adapters — real curated knowledge.

Five HAGR databases replace what used to be a hand-written table of twelve genes
and seven interventions in ``etl/build_fixtures.py``. The facts in that table
were real; the table was not sourced, so nothing could tell you when it went out
of date or where a claim came from. These are the actual releases, checksum
pinned like every other artifact, with per-row provenance back to the database
and, where HAGR supplies one, the PubMed id.

* **GenAge (human)** — genes with evidence of a role in human ageing.
* **GenAge (model organisms)** — genes whose manipulation changes lifespan.
* **CellAge** — genes shown to induce, inhibit, or regulate cellular senescence.
* **LongevityMap** — human variants tested for association with longevity,
  *including the ones that found nothing*. Keeping the non-significant entries is
  the point: a gene list built only from positive associations cannot tell you
  that a gene was looked at and did not replicate.
* **DrugAge** — compounds tested for an effect on model-organism lifespan.
* **GenDR** — genes required for dietary restriction to extend lifespan.

Two boundaries are enforced here rather than left to the reader:

**Species.** GeroQuery models human and mouse. HAGR's model-organism tables are
mostly yeast, worm, and fly. Those rows are parsed and counted but their species
is left unresolved, so nothing silently arrives in the store labelled as if it
were mammalian.

**Gene-intervention links.** DrugAge records no gene targets. The old table
asserted them by hand (rapamycin -> MTOR). Inventing that edge is exactly the
kind of plausible-looking claim this project should not make, so drug
interventions ship with no linked genes and only GenDR — which does record
gene-level dependence — populates that field.
"""

from __future__ import annotations

import csv
import io
import zipfile
from dataclasses import dataclass
from pathlib import Path

from ..exceptions import SourceError
from .base import Capabilities, License, SourceAdapter
from .manifest import HAGR, HAGR_MEMBERS

HAGR_BASE = "https://genomics.senescence.info"

DATABASE_URLS = {
    "GenAge": f"{HAGR_BASE}/genes/",
    "GenAge (models)": f"{HAGR_BASE}/genes/models.html",
    "CellAge": f"{HAGR_BASE}/cells/",
    "LongevityMap": f"{HAGR_BASE}/longevity/",
    "DrugAge": f"{HAGR_BASE}/drugs/",
    "GenDR": f"{HAGR_BASE}/diet/",
}

ORGANISM_TO_SPECIES = {"Homo sapiens": "human", "Mus musculus": "mouse"}

# DrugAge marks a lifespan result significant ("S"), not significant ("NS"), or
# untested ("NA"). Only "S" rows contribute to a reported effect size.
SIGNIFICANT = "S"


@dataclass(frozen=True)
class GeneAssertion:
    """One curated claim about one gene, before identifier resolution."""

    database: str
    symbol: str
    entrez: str | None
    organism: str
    assertion: str
    url: str

    @property
    def species(self) -> str | None:
        return ORGANISM_TO_SPECIES.get(self.organism)


@dataclass(frozen=True)
class LifespanIntervention:
    """One compound or regimen tested against lifespan in one organism."""

    intervention_id: str
    name: str
    itype: str
    source: str
    organism: str
    # None means "tested, and no significant lifespan effect was reported" —
    # distinct from a small positive number, and worth keeping.
    lifespan_effect_pct: float | None
    n_experiments: int
    gene_symbols: tuple[str, ...]
    url: str


def _read_member(path: Path, member: str) -> str:
    try:
        with zipfile.ZipFile(path) as archive:
            return archive.read(member).decode("utf-8", "replace")
    except (KeyError, zipfile.BadZipFile) as exc:
        raise SourceError(
            f"Could not read {member!r} from {path.name}.",
            detail={"path": str(path), "member": member, "error": str(exc)},
        ) from exc


def _rows(text: str, delimiter: str = ",") -> list[dict[str, str]]:
    return [
        {(k or "").strip(): (v or "").strip() for k, v in row.items()}
        for row in csv.DictReader(io.StringIO(text), delimiter=delimiter)
    ]


def _pubmed(reference: str) -> str | None:
    ref = reference.strip().split(",")[0].split(";")[0].strip()
    return f"https://pubmed.ncbi.nlm.nih.gov/{ref}/" if ref.isdigit() else None


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in text.lower()).strip("_")


# Longest real gene symbol is well under this; longest organism binomial too.
_MAX_FIELD = 64


def is_well_formed(*fields: str) -> bool:
    """Reject fields that carry the marks of a broken CSV row.

    GenDR's current release contains one physical line where an unquoted field
    runs many ``symbol;organism`` pairs together, so a naive parse yields a gene
    called ``His4:CG33873;Drosophila melanogaster`` and an organism to match.
    Rows like that are dropped and counted rather than admitted: a fabricated
    gene symbol in the curated table would join to nothing and look like a
    coverage gap instead of a parse failure.
    """
    return all(
        field and ";" not in field and len(field) <= _MAX_FIELD and "\t" not in field
        for field in fields
    )


def _float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# --- per-database parsers ---------------------------------------------------


def parse_genage_human(text: str) -> list[GeneAssertion]:
    """``GenAge ID, symbol, name, entrez gene id, uniprot, why``.

    ``why`` is HAGR's evidence code — ``mammal``, ``model``, ``cell``,
    ``functional``, ``putative``, ``upstream``, ``downstream``, possibly several
    comma-separated. It is carried into the assertion verbatim: "putative" and
    "mammal" are very different strengths of claim and flattening them to
    "ageing-associated" would lose that.
    """
    out = []
    for row in _rows(text):
        symbol = row.get("symbol", "")
        if not is_well_formed(symbol):
            continue
        evidence = row.get("why", "") or "unspecified"
        out.append(
            GeneAssertion(
                database="GenAge",
                symbol=symbol,
                entrez=row.get("entrez gene id") or None,
                organism="Homo sapiens",
                assertion=f"human ageing-associated (evidence: {evidence})",
                url=DATABASE_URLS["GenAge"],
            )
        )
    return out


def parse_genage_models(text: str) -> list[GeneAssertion]:
    """``... organism, entrez gene id, avg lifespan change (max obsv), lifespan effect,
    longevity influence``."""
    out = []
    for row in _rows(text):
        symbol = row.get("symbol", "")
        organism = row.get("organism", "")
        if not is_well_formed(symbol, organism):
            continue
        influence = row.get("longevity influence", "") or "Unannotated"
        effect = row.get("lifespan effect", "") or "unspecified"
        change = row.get("avg lifespan change (max obsv)", "")
        detail = f"{influence.lower()}; lifespan effect: {effect.lower()}"
        if change:
            detail += f"; max observed change {change}%"
        out.append(
            GeneAssertion(
                database="GenAge (models)",
                symbol=symbol,
                entrez=row.get("entrez gene id") or None,
                organism=organism,
                assertion=detail,
                url=DATABASE_URLS["GenAge (models)"],
            )
        )
    return out


def parse_cellage(text: str) -> list[GeneAssertion]:
    """``Entrez ID, Gene symbol, Gene name, Cancer Cell, Type of senescence,
    Senescence Effect, Reference`` (tab separated)."""
    out = []
    for row in _rows(text, delimiter="\t"):
        symbol = row.get("Gene symbol", "")
        if not is_well_formed(symbol):
            continue
        effect = (row.get("Senescence Effect", "") or "regulates").lower()
        kind = (row.get("Type of senescence", "") or "unclear").lower()
        out.append(
            GeneAssertion(
                database="CellAge",
                symbol=symbol,
                entrez=row.get("Entrez ID") or None,
                organism="Homo sapiens",
                assertion=f"{effect} senescence ({kind})",
                url=_pubmed(row.get("Reference", "")) or DATABASE_URLS["CellAge"],
            )
        )
    return out


def parse_longevitymap(text: str) -> list[GeneAssertion]:
    """``id, Association, Population, Variant(s), Gene(s), PubMed``.

    One row can name several genes; each becomes its own assertion, because the
    store joins curated knowledge one gene at a time.
    """
    out = []
    for row in _rows(text):
        genes = [g.strip() for g in row.get("Gene(s)", "").split(",") if g.strip()]
        association = (row.get("Association", "") or "unspecified").lower()
        population = row.get("Population", "") or "unspecified population"
        variant = row.get("Variant(s)", "") or "unspecified variant"
        for symbol in (g for g in genes if is_well_formed(g)):
            out.append(
                GeneAssertion(
                    database="LongevityMap",
                    symbol=symbol,
                    entrez=None,  # LongevityMap indexes by symbol only.
                    organism="Homo sapiens",
                    assertion=f"{association} longevity association in {population} ({variant})",
                    url=_pubmed(row.get("PubMed", "")) or DATABASE_URLS["LongevityMap"],
                )
            )
    return out


def parse_gendr(text: str) -> list[GeneAssertion]:
    """``GenDR ID, gene symbol, species, entrez gene id, gene name``."""
    out = []
    for row in _rows(text):
        symbol = row.get("gene symbol", "")
        organism = row.get("species", "")
        if not is_well_formed(symbol, organism):
            continue
        out.append(
            GeneAssertion(
                database="GenDR",
                symbol=symbol,
                entrez=row.get("entrez gene id") or None,
                organism=organism,
                assertion="required for the lifespan extension of dietary restriction",
                url=DATABASE_URLS["GenDR"],
            )
        )
    return out


def parse_drugage(text: str) -> list[LifespanIntervention]:
    """Aggregate DrugAge experiments into one record per compound and organism.

    DrugAge is a table of *experiments*: the same compound appears many times at
    different doses, strains, and sexes. The median of the significant average
    lifespan changes is reported, and ``n_experiments`` says how many rows stand
    behind it, so a result from one experiment is distinguishable from a result
    that replicated.
    """
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in _rows(text):
        name, organism = row.get("compound_name", ""), row.get("species", "")
        if is_well_formed(organism) and name:
            grouped.setdefault((name, organism), []).append(row)

    out = []
    for (name, organism), rows in sorted(grouped.items()):
        changes = sorted(
            value
            for row in rows
            if row.get("avg_lifespan_significance") == SIGNIFICANT
            and (value := _float(row.get("avg_lifespan_change_percent", ""))) is not None
        )
        median = None
        if changes:
            mid = len(changes) // 2
            middle = changes[mid] if len(changes) % 2 else (changes[mid - 1] + changes[mid]) / 2
            median = round(middle, 4)
        # NIA Interventions Testing Program provenance, where DrugAge records it.
        source = "ITP" if any(r.get("ITP", "").lower() == "yes" for r in rows) else "DrugAge"
        out.append(
            LifespanIntervention(
                intervention_id=f"DRUGAGE:{_slug(name)}:{_slug(organism)}",
                name=name,
                itype="drug",
                source=source,
                organism=organism,
                lifespan_effect_pct=median,
                n_experiments=len(rows),
                gene_symbols=(),  # DrugAge records no gene targets. See module docstring.
                url=DATABASE_URLS["DrugAge"],
            )
        )
    return out


def gendr_interventions(assertions: list[GeneAssertion]) -> list[LifespanIntervention]:
    """One dietary-restriction record per organism, linked to its GenDR genes."""
    by_organism: dict[str, list[str]] = {}
    for item in assertions:
        if item.database == "GenDR" and item.organism:
            by_organism.setdefault(item.organism, []).append(item.symbol)
    return [
        LifespanIntervention(
            intervention_id=f"GENDR:dietary_restriction:{_slug(organism)}",
            name="dietary restriction",
            itype="dietary",
            source="GenDR",
            organism=organism,
            # GenDR asserts gene-level dependence, not a lifespan percentage.
            lifespan_effect_pct=None,
            n_experiments=len(symbols),
            gene_symbols=tuple(sorted(set(symbols))),
            url=DATABASE_URLS["GenDR"],
        )
        for organism, symbols in sorted(by_organism.items())
    ]


_PARSERS = {
    "genage_human": parse_genage_human,
    "genage_models": parse_genage_models,
    "cellage": parse_cellage,
    "longevitymap": parse_longevitymap,
    "gendr": parse_gendr,
}


def load_assertions(paths: dict[str, Path]) -> list[GeneAssertion]:
    """Parse every gene-level HAGR database whose artifact is present."""
    out: list[GeneAssertion] = []
    for key, parser in _PARSERS.items():
        path = paths.get(key)
        if path is not None:
            out.extend(parser(_read_member(path, HAGR_MEMBERS[key])))
    return out


def load_interventions(paths: dict[str, Path]) -> list[LifespanIntervention]:
    """Parse DrugAge and derive the GenDR dietary-restriction records."""
    out: list[LifespanIntervention] = []
    if (path := paths.get("drugage")) is not None:
        out.extend(parse_drugage(_read_member(path, HAGR_MEMBERS["drugage"])))
    if (path := paths.get("gendr")) is not None:
        out.extend(gendr_interventions(parse_gendr(_read_member(path, HAGR_MEMBERS["gendr"]))))
    return out


class HagrSource(SourceAdapter):
    """Adapter over the HAGR curated aging databases."""

    name = "hagr"

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("annotation", "intervention"),
            species=("human", "mouse"),
            federated=False,
            cacheable=True,
            notes=(
                "GenAge, CellAge, LongevityMap, DrugAge, GenDR — checksum-pinned "
                f"releases ({', '.join(sorted(HAGR))})."
            ),
        )

    def license(self) -> License:
        artifact = HAGR["genage_human"]
        return License(artifact.license, redistributable=True, attribution=artifact.attribution)
