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

import glob as glob_module
import hashlib
import re
import sqlite3
from pathlib import Path

import duckdb
import pandas as pd

from .. import DATA_VERSION, config
from ..exceptions import GeroQueryError
from ..models import AgingSignature, CuratedFlag, Intervention, Study
from ..sources import (
    CuratedKnowledgeSource,
    InterventionSource,
    LocalSignatureSource,
    NhanesClinicalSource,
    nhanes,
)


class DatasetNotFoundError(GeroQueryError):
    code = "dataset_not_found"
    http_status = 404


class GeroStore:
    def __init__(self, data_home: Path | None = None):
        self.data_home = Path(data_home or config.DATA_HOME)
        self.sig_dir = self.data_home / "signatures"  # partitioned parquet
        self.ds_dir = self.data_home / "datasets"
        self.meta_db = self.data_home / "metadata.sqlite"
        # Once built, stays built for this process. Without the latch, every
        # query stat()'d the filesystem twice before doing any work.
        self._built = False
        # One long-lived DuckDB connection per store. Reconnecting per query
        # re-opened and re-scanned the Parquet footers each time.
        self._duck: duckdb.DuckDBPyConnection | None = None

    # ---- build -----------------------------------------------------------

    def is_built(self) -> bool:
        if self._built:
            return True
        self._built = self.meta_db.exists() and self.sig_dir.exists()
        return self._built

    def ensure_built(self) -> GeroStore:
        if not self.is_built():
            self.build()
        return self

    def _duck_con(self) -> duckdb.DuckDBPyConnection:
        if self._duck is None:
            self._duck = duckdb.connect()
        return self._duck

    def close(self) -> None:
        """Release the DuckDB connection. Safe to call more than once."""
        if self._duck is not None:
            self._duck.close()
            self._duck = None

    def build(
        self,
        clinical_csv: Path | None = None,
        prefer_full_clinical: bool = True,
        prefer_full_signatures: bool = True,
    ) -> GeroStore:
        """(Re)materialize the store from the bundled source adapters.

        Both ``prefer_full_*`` flags select between the full table built by
        ``make data`` and the committed offline slice of the same real data.
        Tests pin them to the slice so a result does not depend on whether the
        developer happens to have run the download.
        """
        self.data_home.mkdir(parents=True, exist_ok=True)
        self.sig_dir.mkdir(parents=True, exist_ok=True)
        self.ds_dir.mkdir(parents=True, exist_ok=True)

        sig_source = LocalSignatureSource(prefer_full=prefer_full_signatures)
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
        dataset_registry = self._register_clinical(clinical_csv, prefer_full_clinical)

        self._build_metadata(sig_source, dataset_registry)
        # Drop any connection held over from a previous build: it may have cached
        # metadata for Parquet files this build just replaced.
        self.close()
        self._built = True
        return self

    def _write_dataset(
        self, dataset_id: str, df: pd.DataFrame, kind: str, description: str
    ) -> dict:
        out = self.ds_dir / f"{dataset_id}.parquet"
        df.to_parquet(out, index=False)
        return {
            "dataset_id": dataset_id,
            "kind": kind,
            "n_rows": len(df),
            "columns": ",".join(df.columns),
            "path": out.name,
            "description": description,
        }

    def _register_clinical(
        self, clinical_csv: Path | None = None, prefer_full: bool = True
    ) -> list[dict]:
        """Materialize the clinical datasets, real and synthetic, kept separate.

        ``clinical_nhanes_slice`` is real NHANES 2017-2018. ``clinical_synthetic_csd``
        is the generated method-validation fixture. They are deliberately two
        datasets rather than one: the synthetic one has a critical-slowing-down
        signal planted in it, and a caller that cannot tell them apart will read
        that planted signal as a finding about people.
        """
        registry: list[dict] = []

        if clinical_csv is not None:
            df = pd.read_csv(clinical_csv)
            registry.append(
                self._write_dataset(
                    "clinical_nhanes_slice",
                    df,
                    "clinical",
                    f"Clinical marker table loaded from {clinical_csv.name}.",
                )
            )
        else:
            source = NhanesClinicalSource()
            source.assert_cacheable()  # NHANES is public domain; this must pass.
            df, mode = source.clinical_frame(prefer_full=prefer_full)
            caveat = (
                "Full verified cohort."
                if mode == "full"
                else "OFFLINE SAMPLE of real rows — estimates from it are not the "
                "reported cohort result. Run `make data` for the full cohort."
            )
            registry.append(
                self._write_dataset(
                    "clinical_nhanes_slice",
                    df,
                    "clinical",
                    f"REAL {nhanes.RELEASE} clinical markers, joined on SEQN "
                    f"({mode} mode, n={len(df)}). {caveat} Age topcoded at 80; "
                    f"survey weights carried but not applied.",
                )
            )

        synthetic_csv = config.SOURCES_DATA / "clinical_synthetic_csd.csv"
        if synthetic_csv.exists():
            df = pd.read_csv(synthetic_csv)
            registry.append(
                self._write_dataset(
                    "clinical_synthetic_csd",
                    df,
                    "synthetic",
                    "SYNTHETIC method-validation fixture. Critical slowing down is planted "
                    "in it by construction (latent-factor variance grows with age). Use to "
                    "check that estimators recover a known effect — never as evidence "
                    "about human biology.",
                )
            )
        return registry

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
                " license TEXT, url TEXT, version TEXT, series_id TEXT);"
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
                ":sample_size,:processing_method,:license,:url,:version,:series_id)",
                # to_dict() drops None values, so every optional column is filled
                # back in explicitly — a missing key is a bind error, not a NULL.
                [
                    s.to_dict()
                    | {
                        "tissue": s.tissue,
                        "sample_size": s.sample_size,
                        "processing_method": s.processing_method,
                        "license": s.license,
                        "url": s.url,
                        "version": s.version,
                        "series_id": s.series_id,
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

    # Partition values come from a closed vocabulary; anything else is refused
    # rather than interpolated into a glob path.
    _PARTITION_VALUE = re.compile(r"^[A-Za-z0-9_.-]+$")

    def _sig_glob(self, species: str | None = None, omic_layer: str | None = None) -> str:
        """Glob narrowed to the relevant hive partitions.

        The data is partitioned ``species=<s>/omic_layer=<o>``. Naming the
        partition directories directly means DuckDB opens only those Parquet
        files. Leaving it as ``**`` and relying on the WHERE clause makes pruning
        depend on the filter value being visible at plan time, which it is not
        when the value arrives as a bound ``?`` parameter — so every file gets
        opened and its footer read on every query.
        """
        parts = []
        for name, value in (("species", species), ("omic_layer", omic_layer)):
            if value is None:
                parts.append("*")
            elif self._PARTITION_VALUE.match(value):
                parts.append(f"{name}={value}")
            else:
                raise GeroQueryError(f"Invalid {name} filter {value!r}.", detail={name: value})
        return (self.sig_dir.joinpath(*parts) / "*.parquet").as_posix()

    def query_signatures(
        self,
        gene_id: str | None = None,
        species: str | None = None,
        tissue: str | None = None,
        omic_layer: str | None = None,
        sex: str | None = None,
    ) -> list[AgingSignature]:
        self.ensure_built()
        glob = self._sig_glob(species=species, omic_layer=omic_layer)
        if not glob_module.glob(glob):
            # No partition matches the requested species/omic_layer combination.
            # read_parquet raises on an empty glob, so answer directly.
            return []

        # species/omic_layer are already enforced by the partition path above;
        # repeating them in WHERE would only re-filter rows we know match.
        clauses, params = [], []
        for col, val in [
            ("gene_id", gene_id),
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
            f"FROM read_parquet('{glob}', hive_partitioning=true){where} "
            f"ORDER BY gene_id, omic_layer, study_id"
        )
        rows = self._duck_con().execute(sql, params).fetchall()
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
