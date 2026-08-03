"""GTEx Portal API v2 adapter — open summary tier, queried live.

Serves three things, all from the open API:

* ``lookup_gene`` — symbol or Ensembl id -> versioned GENCODE id (GTEx keys
  everything on ``ENSG…​.17``, not the bare Ensembl id).
* ``median_expression`` — median TPM per tissue, with UBERON ontology ids.
* ``sample_expression`` — the per-sample TPM distribution for one tissue.

**What this adapter deliberately does not do: age stratification.**
``/expression/geneExpression`` accepts ``attributeSubset=ageBracket``, and it is
tempting to read that as "GTEx will give me expression by age". It does not — the
open endpoint returns a single record with ``subsetGroup: null`` and the full
undivided sample list. Donor age is part of GTEx's protected tier and reaching it
requires dbGaP authorization, which is what the separate ``gtex-protected``
federate-only adapter represents. GTEx therefore contributes *tissue context*
here, not an aging signature. Verified against the live API on 2026-08-03.

Requests must pass ``datasetId``. Without it every endpoint returns HTTP 200 with
an empty ``data`` list — a silent empty result, not an error, which is the kind of
thing that turns into "why does this gene have no expression anywhere".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .. import config
from ..exceptions import SourceError
from .base import Capabilities, License, SourceAdapter

API_ROOT = "https://gtexportal.org/api/v2"

# The release the open API actually serves. `gtex_v10` is accepted by the query
# validator but returns nothing, so it is not a safe default.
DEFAULT_DATASET_ID = "gtex_v8"

_HEADERS = {
    "User-Agent": "GeroQuery/1.0 (+https://github.com/praneel/geroquery)",
    "Accept": "application/json",
}


@dataclass(frozen=True)
class GtexGene:
    gencode_id: str
    symbol: str
    ensembl_id: str
    chromosome: str | None = None
    description: str | None = None


@dataclass(frozen=True)
class TissueExpression:
    tissue_id: str
    uberon_id: str | None
    median_tpm: float
    unit: str
    gencode_id: str
    symbol: str
    dataset_id: str

    def to_dict(self) -> dict:
        return {
            "tissue_id": self.tissue_id,
            "uberon_id": self.uberon_id,
            "median_tpm": self.median_tpm,
            "unit": self.unit,
            "gencode_id": self.gencode_id,
            "symbol": self.symbol,
            "dataset_id": self.dataset_id,
        }


class GtexOpenSource(SourceAdapter):
    """Live adapter over the GTEx Portal open summary API."""

    name = "gtex-open"

    def __init__(
        self,
        dataset_id: str = DEFAULT_DATASET_ID,
        allow_network: bool | None = None,
        timeout: float | None = None,
    ):
        self.dataset_id = dataset_id
        self._allow_network = allow_network
        self._timeout = timeout

    # ---- adapter surface -------------------------------------------------

    def capabilities(self) -> Capabilities:
        return Capabilities(
            source_name=self.name,
            omics=("transcriptome",),
            species=("human",),
            federated=True,
            cacheable=True,
            notes=(
                f"GTEx Portal API v2 ({self.dataset_id}) open summary tier: median TPM per "
                "tissue with UBERON ids, and per-sample TPM distributions. Donor age is "
                "NOT available here — age-linked GTEx is dbGaP-controlled, see "
                "'gtex-protected'."
            ),
        )

    def license(self) -> License:
        return License(
            "open-summary",
            redistributable=True,
            attribution="GTEx Portal API v2 (GTEx Consortium). Open access summary statistics.",
        )

    @property
    def allow_network(self) -> bool:
        return config.ALLOW_NETWORK if self._allow_network is None else self._allow_network

    @property
    def timeout(self) -> float:
        return config.HTTP_TIMEOUT if self._timeout is None else self._timeout

    # ---- transport -------------------------------------------------------

    def _get(self, path: str, params: dict[str, Any]) -> list[dict]:
        if not self.allow_network:
            raise SourceError(
                f"GTEx is a live source and network access is disabled. "
                f"Set GEROQUERY_ALLOW_NETWORK=1 to query {path}.",
                detail={"source": self.name, "path": path},
            )
        import httpx

        try:
            resp = httpx.get(
                f"{API_ROOT}{path}",
                params=params,
                headers=_HEADERS,
                timeout=self.timeout,
                follow_redirects=True,
            )
            resp.raise_for_status()
            payload = resp.json()
        except Exception as exc:
            raise SourceError(
                f"GTEx request to {path} failed: {exc}",
                detail={"source": self.name, "path": path, "params": params},
            ) from exc
        data = payload.get("data")
        return data if isinstance(data, list) else []

    # ---- queries ---------------------------------------------------------

    def lookup_gene(self, gene: str) -> GtexGene:
        """Resolve a symbol or Ensembl id to GTEx's versioned GENCODE id."""
        rows = self._get("/reference/gene", {"geneId": gene})
        if not rows:
            raise SourceError(f"GTEx has no gene matching {gene!r}.", detail={"query": gene})
        # Prefer an exact symbol match; GTEx will also return alias hits.
        wanted = gene.strip().upper()
        row = next(
            (r for r in rows if str(r.get("geneSymbolUpper", "")).upper() == wanted), rows[0]
        )
        gencode = str(row["gencodeId"])
        return GtexGene(
            gencode_id=gencode,
            symbol=row.get("geneSymbol", gene),
            ensembl_id=gencode.split(".")[0],
            chromosome=row.get("chromosome"),
            description=row.get("description"),
        )

    def median_expression(
        self, gene: str, tissue_ids: list[str] | None = None
    ) -> list[TissueExpression]:
        """Median TPM across tissues for one gene."""
        resolved = self.lookup_gene(gene)
        params: dict[str, Any] = {
            "gencodeId": resolved.gencode_id,
            "datasetId": self.dataset_id,
        }
        if tissue_ids:
            params["tissueSiteDetailId"] = tissue_ids
        rows = self._get("/expression/medianGeneExpression", params)
        out = [
            TissueExpression(
                tissue_id=r["tissueSiteDetailId"],
                uberon_id=r.get("ontologyId"),
                median_tpm=float(r["median"]),
                unit=r.get("unit", "TPM"),
                gencode_id=r["gencodeId"],
                symbol=r.get("geneSymbol", resolved.symbol),
                dataset_id=r.get("datasetId", self.dataset_id),
            )
            for r in rows
            if r.get("median") is not None
        ]
        out.sort(key=lambda t: -t.median_tpm)
        return out

    def sample_expression(self, gene: str, tissue_id: str) -> dict:
        """Per-sample TPM values for one gene in one tissue.

        No donor ages accompany these values — see the module docstring. Returned
        with ``age_stratified: False`` so a caller cannot mistake the distribution
        for an age-resolved one.
        """
        resolved = self.lookup_gene(gene)
        rows = self._get(
            "/expression/geneExpression",
            {
                "gencodeId": resolved.gencode_id,
                "tissueSiteDetailId": tissue_id,
                "datasetId": self.dataset_id,
            },
        )
        values: list[float] = []
        uberon = None
        unit = "TPM"
        for r in rows:
            values.extend(float(v) for v in r.get("data", []))
            uberon = r.get("ontologyId", uberon)
            unit = r.get("unit", unit)
        return {
            "gene": resolved.symbol,
            "gencode_id": resolved.gencode_id,
            "tissue_id": tissue_id,
            "uberon_id": uberon,
            "unit": unit,
            "n_samples": len(values),
            "values": values,
            "dataset_id": self.dataset_id,
            "age_stratified": False,
            "note": (
                "GTEx open API returns donor ages only in the protected tier; these "
                "samples are unstratified. Age-resolved GTEx requires dbGaP access."
            ),
        }

    def fetch_signature(self, gene: str, tissue_id: str | None = None) -> dict:
        """Uniform entry point used by the federated query path."""
        if tissue_id:
            return self.sample_expression(gene, tissue_id)
        return {
            "gene": gene,
            "dataset_id": self.dataset_id,
            "median_expression": [t.to_dict() for t in self.median_expression(gene)],
            "age_stratified": False,
        }
