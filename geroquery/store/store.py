"""M4 store — analytical storage & query.

Storage split, per the PRD:
  * large harmonized matrices/signatures  -> Parquet (partitioned) read by DuckDB
  * small relational metadata (studies, curated flags, interventions, datasets)
    -> SQLite

Callers see only the small interface below; whether a query hits Parquet via
DuckDB or a SQLite row is hidden. Swapping SQLite for Postgres, or local Parquet
for S3/HF, would not change this surface.
"""

from __future__ import annotations

import hashlib
import sqlite3
from pathlib import Path

import duckdb
import pandas as pd

from .. import DATA_VERSION, config
from ..exceptions import GeroQueryError
from ..models import AgingSignature, CuratedFlag, Intervention, Study
from ..sources import CuratedKnowledgeSource, InterventionSource, LocalSignatureSource


class DatasetNotFoundError(GeroQueryError):
    code = "dataset_not_found"
    http_status = 404


class GeroStore:
    def __init__(self, data_home: Path | None = None):
        self.data_home = Path(data_home or config.DATA_HOME)
        self.sig_dir = self.data_home / "signatures"  # partitioned parquet
        self.ds_dir = self.data_home / "datasets"
        self.meta_db = self.data_home / "metadata.sqlite"

    # ---- build -----------------------------------------------------------

    def is_built(self) -> bool:
        return self.meta_db.exists() and self.sig_dir.exists()

    def ensure_built(self) -> GeroStore:
        if not self.is_built():
            self.build()
        return self

    def build(self, clinical_csv: Path | None = None) -> GeroStore:
        """(Re)materialize the store from the bundled source adapters."""
        self.data_home.mkdir(parents=True, exist_ok=True)
        self.sig_dir.mkdir(parents=True, exist_ok=True)
        self.ds_dir.mkdir(parents=True, exist_ok=True)

        sig_source = LocalSignatureSource()
        # Licence gate: only persist what may be redistributed.
        sig_source.assert_cacheable()

        signatures = sig_source.signatures()
        sig_df = pd.DataFrame([s.to_dict() for s in signatures])

        con = duckdb.connect()
        try:
            con.register("sig", sig_df)
            # Partitioned Parquet by species + omic_layer (predicate pushdown).
            con.execute(
                f"COPY (SELECT * FROM sig) TO '{self.sig_dir.as_posix()}' "
                "(FORMAT PARQUET, PARTITION_BY (species, omic_layer), OVERWRITE_OR_IGNORE 1)"
            )
        finally:
            con.close()

        # Clinical / phenotype datasets -> one parquet each, registered below.
        clinical_csv = clinical_csv or (config.SOURCES_DATA / "clinical_nhanes_slice.csv")
        dataset_registry = []
        if clinical_csv.exists():
            df = pd.read_csv(clinical_csv)
            out = self.ds_dir / "clinical_nhanes_slice.parquet"
            df.to_parquet(out, index=False)
            dataset_registry.append(
                {
                    "dataset_id": "clinical_nhanes_slice",
                    "kind": "clinical",
                    "n_rows": len(df),
                    "columns": ",".join(df.columns),
                    "path": out.name,
                    "description": "NHANES-style clinical marker slice (demo) with age strata.",
                }
            )

        self._build_metadata(sig_source, dataset_registry)
        return self

    def _build_metadata(self, sig_source: LocalSignatureSource, dataset_registry: list[dict]):
        studies = sig_source.studies()
        curated = CuratedKnowledgeSource().flags()
        interventions = InterventionSource().interventions()

        if self.meta_db.exists():
            self.meta_db.unlink()
        con = sqlite3.connect(self.meta_db)
        try:
            con.executescript(
                "CREATE TABLE studies (study_id TEXT PRIMARY KEY, source TEXT, omic_layer TEXT,"
                " species TEXT, tissue TEXT, sample_size INTEGER, processing_method TEXT,"
                " license TEXT, url TEXT, version TEXT);"
                " CREATE TABLE curated_knowledge (gene_id TEXT, database TEXT, assertion TEXT,"
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
                "INSERT INTO studies VALUES (:study_id,:source,:omic_layer,:species,:tissue,"
                ":sample_size,:processing_method,:license,:url,:version)",
                [
                    s.to_dict()
                    | {
                        "tissue": s.tissue,
                        "sample_size": s.sample_size,
                        "processing_method": s.processing_method,
                        "license": s.license,
                        "url": s.url,
                        "version": s.version,
                    }
                    for s in studies
                ],
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
        """DATA_VERSION plus a short content hash of the signature parquet."""
        h = hashlib.sha256()
        for p in sorted(self.sig_dir.rglob("*.parquet")):
            h.update(p.name.encode())
            h.update(str(p.stat().st_size).encode())
        return f"{DATA_VERSION}+{h.hexdigest()[:12]}"

    # ---- query -----------------------------------------------------------

    def _sig_glob(self) -> str:
        return (self.sig_dir / "**" / "*.parquet").as_posix()

    def query_signatures(
        self,
        gene_id: str | None = None,
        species: str | None = None,
        tissue: str | None = None,
        omic_layer: str | None = None,
        sex: str | None = None,
    ) -> list[AgingSignature]:
        self.ensure_built()
        clauses, params = [], []
        for col, val in [
            ("gene_id", gene_id),
            ("species", species),
            ("omic_layer", omic_layer),
            ("tissue", tissue),
            ("sex", sex),
        ]:
            if val is not None:
                clauses.append(f"{col} = ?")
                params.append(val)
        where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
        sql = (
            f"SELECT gene_id, study_id, omic_layer, species, tissue, sex, age_range, "
            f"effect_size, direction, p_value, q_value, standard_error, source "
            f"FROM read_parquet('{self._sig_glob()}', hive_partitioning=true){where} "
            f"ORDER BY gene_id, omic_layer, study_id"
        )
        con = duckdb.connect()
        try:
            rows = con.execute(sql, params).fetchall()
        finally:
            con.close()
        cols = [
            "gene_id",
            "study_id",
            "omic_layer",
            "species",
            "tissue",
            "sex",
            "age_range",
            "effect_size",
            "direction",
            "p_value",
            "q_value",
            "standard_error",
            "source",
        ]
        return [AgingSignature(**dict(zip(cols, r, strict=True))) for r in rows]

    # ---- relational metadata --------------------------------------------

    def _meta_con(self) -> sqlite3.Connection:
        self.ensure_built()
        con = sqlite3.connect(self.meta_db)
        con.row_factory = sqlite3.Row
        return con

    def list_studies(self) -> list[Study]:
        con = self._meta_con()
        try:
            rows = con.execute("SELECT * FROM studies ORDER BY study_id").fetchall()
        finally:
            con.close()
        return [Study(**dict(r)) for r in rows]

    def get_study(self, study_id: str) -> Study | None:
        con = self._meta_con()
        try:
            row = con.execute("SELECT * FROM studies WHERE study_id = ?", (study_id,)).fetchone()
        finally:
            con.close()
        return Study(**dict(row)) if row else None

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
