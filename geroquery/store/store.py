"""M4 store — analytical storage & query.

Storage split:
  * phenotype/clinical datasets (the simulated example cohort, or any uploaded
    matrix persisted for reuse) -> Parquet, read directly by pandas/DuckDB.
  * small relational metadata (curated flags, interventions, dataset registry)
    -> SQLite.

Callers see only the small interface below. Swapping SQLite for Postgres, or
local Parquet for S3/HF, would not change this surface.

Note: earlier versions kept a partitioned Parquet table of fabricated aging
"signatures". That has been removed — the gene-level aging evidence GeroQuery
serves is real, curated, and cited, and lives in :mod:`geroquery.knowledge`,
not in a generated data table.
"""

from __future__ import annotations

import hashlib
import sqlite3
from pathlib import Path

import pandas as pd

from .. import DATA_VERSION, config
from ..exceptions import GeroQueryError
from ..models import CuratedFlag, Intervention
from ..sources import CuratedKnowledgeSource, InterventionSource, LocalEvidenceSource


class DatasetNotFoundError(GeroQueryError):
    code = "dataset_not_found"
    http_status = 404


# The clinical example cohort shipped for trying the clock/resilience tools.
_EXAMPLE_COHORT_CSV = "example_cohort_simulated.csv"
_EXAMPLE_COHORT_ID = "example_cohort_simulated"


class GeroStore:
    def __init__(self, data_home: Path | None = None):
        self.data_home = Path(data_home or config.DATA_HOME)
        self.ds_dir = self.data_home / "datasets"
        self.meta_db = self.data_home / "metadata.sqlite"

    # ---- build -----------------------------------------------------------

    def is_built(self) -> bool:
        return self.meta_db.exists() and self.ds_dir.exists()

    def ensure_built(self) -> GeroStore:
        if not self.is_built():
            self.build()
        return self

    def build(self, cohort_csv: Path | None = None) -> GeroStore:
        """(Re)materialize the store from the bundled source adapters."""
        self.data_home.mkdir(parents=True, exist_ok=True)
        self.ds_dir.mkdir(parents=True, exist_ok=True)

        # Licence gate: only persist what may be redistributed.
        LocalEvidenceSource().assert_cacheable()

        # Clinical / phenotype datasets -> one parquet each, registered below.
        cohort_csv = cohort_csv or (config.SOURCES_DATA / _EXAMPLE_COHORT_CSV)
        dataset_registry = []
        if cohort_csv.exists():
            df = pd.read_csv(cohort_csv)
            out = self.ds_dir / f"{_EXAMPLE_COHORT_ID}.parquet"
            df.to_parquet(out, index=False)
            dataset_registry.append(
                {
                    "dataset_id": _EXAMPLE_COHORT_ID,
                    "kind": "clinical",
                    "n_rows": len(df),
                    "columns": ",".join(df.columns),
                    "path": out.name,
                    "description": "SIMULATED example cohort — nine PhenoAge clinical markers, "
                    "age and sex. A transparent worked example for the clock and resilience "
                    "tools; not real patient data.",
                }
            )

        self._build_metadata(dataset_registry)
        return self

    def _build_metadata(self, dataset_registry: list[dict]):
        curated = CuratedKnowledgeSource().flags()
        interventions = InterventionSource().interventions()

        if self.meta_db.exists():
            self.meta_db.unlink()
        con = sqlite3.connect(self.meta_db)
        try:
            con.executescript(
                "CREATE TABLE curated_knowledge (gene_id TEXT, database TEXT, assertion TEXT,"
                " url TEXT);"
                " CREATE TABLE interventions (intervention_id TEXT PRIMARY KEY, name TEXT,"
                " itype TEXT, source TEXT, organism TEXT, lifespan_effect_pct REAL,"
                " linked_gene_ids TEXT, url TEXT);"
                " CREATE TABLE datasets (dataset_id TEXT PRIMARY KEY, kind TEXT, n_rows INTEGER,"
                " columns TEXT, path TEXT, description TEXT);"
                " CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);"
                " CREATE INDEX idx_curated_gene ON curated_knowledge(gene_id);"
            )
            con.executemany(
                "INSERT INTO curated_knowledge VALUES (:gene_id,:database,:assertion,:url)",
                [
                    {
                        "gene_id": c.gene_id,
                        "database": c.database,
                        "assertion": c.assertion,
                        "url": c.url,
                    }
                    for c in curated
                ],
            )
            con.executemany(
                "INSERT INTO interventions VALUES (:intervention_id,:name,:itype,:source,:organism,"
                ":lifespan_effect_pct,:linked_gene_ids,:url)",
                [
                    {
                        "intervention_id": iv.intervention_id,
                        "name": iv.name,
                        "itype": iv.itype,
                        "source": iv.source,
                        "organism": iv.organism,
                        "lifespan_effect_pct": iv.lifespan_effect_pct,
                        "linked_gene_ids": "|".join(iv.linked_gene_ids),
                        "url": iv.url,
                    }
                    for iv in interventions
                ],
            )
            con.executemany(
                "INSERT INTO datasets VALUES "
                "(:dataset_id,:kind,:n_rows,:columns,:path,:description)",
                dataset_registry,
            )
            con.execute("INSERT INTO meta VALUES ('data_version', ?)", (self._compute_version(),))
            con.commit()
        finally:
            con.close()

    def _compute_version(self) -> str:
        """DATA_VERSION plus a short content hash of the bundled source files."""
        h = hashlib.sha256()
        for name in ("curated_knowledge.csv", "interventions.csv", _EXAMPLE_COHORT_CSV):
            p = config.SOURCES_DATA / name
            if p.exists():
                h.update(name.encode())
                h.update(p.read_bytes())
        return f"{DATA_VERSION}+{h.hexdigest()[:12]}"

    # ---- relational metadata --------------------------------------------

    def _meta_con(self) -> sqlite3.Connection:
        self.ensure_built()
        con = sqlite3.connect(self.meta_db)
        con.row_factory = sqlite3.Row
        return con

    def curated_flags(self, gene_id: str) -> list[CuratedFlag]:
        con = self._meta_con()
        try:
            rows = con.execute(
                "SELECT gene_id, database, assertion, url FROM curated_knowledge WHERE gene_id = ?",
                (gene_id,),
            ).fetchall()
        finally:
            con.close()
        return [CuratedFlag(**dict(r)) for r in rows]

    def interventions(
        self, name: str | None = None, gene_id: str | None = None
    ) -> list[Intervention]:
        con = self._meta_con()
        try:
            if name is not None:
                rows = con.execute(
                    "SELECT * FROM interventions WHERE lower(name) = ?", (name.lower(),)
                ).fetchall()
            else:
                rows = con.execute("SELECT * FROM interventions ORDER BY name").fetchall()
        finally:
            con.close()
        out = []
        for r in rows:
            linked = tuple(x for x in (r["linked_gene_ids"] or "").split("|") if x)
            if gene_id is not None and gene_id not in linked:
                continue
            out.append(
                Intervention(
                    intervention_id=r["intervention_id"],
                    name=r["name"],
                    itype=r["itype"],
                    source=r["source"],
                    organism=r["organism"],
                    lifespan_effect_pct=r["lifespan_effect_pct"],
                    linked_gene_ids=linked,
                    url=r["url"],
                )
            )
        return out

    def list_datasets(self) -> list[dict]:
        con = self._meta_con()
        try:
            return [
                dict(r)
                for r in con.execute("SELECT * FROM datasets ORDER BY dataset_id").fetchall()
            ]
        finally:
            con.close()

    def get_dataset(self, dataset_id: str) -> pd.DataFrame:
        con = self._meta_con()
        try:
            row = con.execute(
                "SELECT path FROM datasets WHERE dataset_id = ?", (dataset_id,)
            ).fetchone()
        finally:
            con.close()
        if row is None:
            raise DatasetNotFoundError(f"Unknown dataset {dataset_id!r}.", detail=dataset_id)
        return pd.read_parquet(self.ds_dir / row["path"])

    def version(self) -> str:
        con = self._meta_con()
        try:
            row = con.execute("SELECT value FROM meta WHERE key = 'data_version'").fetchone()
        finally:
            con.close()
        return row["value"] if row else DATA_VERSION
