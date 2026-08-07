"""M4 store — round-trip, filtered query, versioning, datasets."""

import pytest

from geroquery.store import DatasetNotFoundError, GeroStore


def test_query_signatures_filters(store):
    sigs = store.query_signatures(
        gene_id="ENSG00000147889", species="human", omic_layer="transcriptome"
    )
    # 14 of the panel's 16 human contrasts carry a CDKN2A probe. Pinned rather
    # than bounded: the panel is fixed by the manifest, so a change here means
    # the manifest, the contrast rules, or the probe annotation moved.
    assert len(sigs) == 14
    assert all(s.species == "human" and s.omic_layer == "transcriptome" for s in sigs)
    # Deliberately not "all up". On real data CDKN2A splits 6 up / 8 down; see
    # test_api.test_p16_does_not_replicate_across_the_real_panel.
    assert {s.direction for s in sigs} == {"up", "down"}


def test_partition_pruning_by_tissue(store):
    brain = store.query_signatures(gene_id="ENSG00000147889", tissue="brain")
    assert brain and all(s.tissue == "brain" for s in brain)


def test_query_glob_narrows_to_the_requested_partitions(store):
    """species/omic_layer must be resolved in the glob path, not left to a bound
    `?` parameter — otherwise DuckDB opens every Parquet file on every query."""
    everything = store._sig_glob()
    assert everything.endswith("/*/*/*.parquet")

    narrowed = store._sig_glob(species="human", omic_layer="transcriptome")
    assert narrowed.endswith("/species=human/omic_layer=transcriptome/*.parquet")

    partial = store._sig_glob(species="mouse")
    assert partial.endswith("/species=mouse/*/*.parquet")


def test_query_rejects_partition_values_that_could_escape_the_glob(store):
    from geroquery.exceptions import GeroQueryError

    with pytest.raises(GeroQueryError, match="Invalid species"):
        store.query_signatures(species="../../etc")


def test_query_returns_empty_for_a_partition_that_does_not_exist(store):
    """mouse has no proteome partition. An empty glob must answer [], not raise
    out of read_parquet."""
    assert store.query_signatures(species="mouse", omic_layer="proteome") == []


def test_filtered_and_unfiltered_queries_agree(store):
    """Narrowing the glob must not change results, only the files opened."""
    narrowed = store.query_signatures(gene_id="ENSMUSG00000044303", species="mouse")
    broad = [
        s
        for s in store.query_signatures(gene_id="ENSMUSG00000044303")
        if s.species == "mouse" and s.omic_layer == "transcriptome"
    ]
    assert narrowed and [s.to_dict() for s in narrowed] == [s.to_dict() for s in broad]


def test_duckdb_connection_is_reused_across_queries(store):
    store.query_signatures(gene_id="ENSG00000147889")
    first = store._duck
    store.query_signatures(gene_id="ENSG00000113368")
    assert first is not None and store._duck is first


def test_is_built_latches_after_the_first_check(store):
    assert store.is_built() is True
    assert store._built is True


def test_studies_and_provenance(store):
    studies = store.list_studies()
    assert len(studies) == 32  # 32 contrasts over 31 GEO DataSets
    one = store.get_study(studies[0].study_id)
    assert one is not None
    assert one.source == "GEO"
    assert one.url.startswith("https://www.ncbi.nlm.nih.gov/")
    assert one.version.startswith("GEO DataSets, curated ")


def test_studies_carry_the_series_id_that_exposes_shared_subjects(store):
    """GEO splits some experiments across two array halves, so two study_ids can
    share subjects. A random-effects pool assumes independence; the series id is
    what lets a caller check that instead of trusting the study count."""
    studies = store.list_studies()
    assert all(s.series_id for s in studies)
    series = [s.series_id for s in studies]
    assert len(set(series)) == 27 < len(series)
    shared = {s for s in series if series.count(s) > 1}
    assert "GSE362" in shared  # GDS287 / GDS288: HG-U133A and HG-U133B, same men


def test_curated_and_interventions(store):
    flags = store.curated_flags("ENSG00000147889")  # CDKN2A
    assert {f.database for f in flags} >= {"GenAge", "CellAge"}
    assert all(f.url.startswith("https://") for f in flags)

    rapa = {i.organism: i for i in store.interventions(name="rapamycin")}
    assert rapa["Mus musculus"].lifespan_effect_pct == 13.0
    assert rapa["Mus musculus"].source == "ITP"  # DrugAge flags the ITP studies
    # DrugAge records no gene targets, so none are asserted. The old fixture
    # linked rapamycin to MTOR by hand; that edge was never in the data.
    assert all(i.linked_gene_ids == () for i in rapa.values())


def test_dietary_restriction_is_the_one_intervention_with_real_gene_links(store):
    dr = [i for i in store.interventions(name="dietary restriction")]
    assert dr, "GenDR should contribute a dietary-restriction record per organism"
    mouse = next(i for i in dr if i.organism == "Mus musculus")
    assert mouse.source == "GenDR" and mouse.itype == "dietary"
    assert len(mouse.linked_gene_ids) >= 5


def test_dataset_load_and_missing(store):
    df = store.get_dataset("clinical_nhanes_slice")
    assert df.shape[0] == 600  # committed NHANES sample; see the `store` fixture
    with pytest.raises(DatasetNotFoundError):
        store.get_dataset("no_such_dataset")


def test_clinical_dataset_is_real_nhanes(store):
    """The headline clinical dataset must be real measurements, not generated.

    Guards the property that actually matters: subject ids are NHANES SEQNs, and
    the values look like NHANES rather than like a normal draw around a constant.
    """
    df = store.get_dataset("clinical_nhanes_slice")
    assert df["subject_id"].str.startswith("NHANES:").all()
    assert df["age"].min() >= 20 and df["age"].max() == 80  # topcoded
    # Real hs-CRP is strongly right-skewed; the synthetic generator's is not.
    assert (df["crp"] > 0).all()
    assert df["crp"].skew() > 2.0


def test_real_and_synthetic_clinical_are_separate_datasets(store):
    """The planted-effect fixture must never be reachable under the real id.

    If these ever collapse into one dataset, a caller reads a critical-slowing-down
    signal that was engineered into the generator as a finding about people.
    """
    ids = {d["dataset_id"] for d in store.list_datasets()}
    assert {"clinical_nhanes_slice", "clinical_synthetic_csd"} <= ids

    by_id = {d["dataset_id"]: d for d in store.list_datasets()}
    assert "REAL" in by_id["clinical_nhanes_slice"]["description"]
    assert "SYNTHETIC" in by_id["clinical_synthetic_csd"]["description"]

    synthetic = store.get_dataset("clinical_synthetic_csd")
    assert synthetic.shape[0] == 720
    assert not synthetic["subject_id"].str.startswith("NHANES:").any()


def test_version_is_reproducible(tmp_path):
    """Two independent builds from the same fixtures produce the same version
    and identical query results."""
    s1 = GeroStore(data_home=tmp_path / "a").build()
    s2 = GeroStore(data_home=tmp_path / "b").build()
    assert s1.version() == s2.version()
    q1 = s1.query_signatures(gene_id="ENSG00000113368")
    q2 = s2.query_signatures(gene_id="ENSG00000113368")
    assert [x.to_dict() for x in q1] == [x.to_dict() for x in q2]


# ---- the store path is data, not syntax (bug #25) --------------------------


@pytest.mark.parametrize(
    "dirname",
    [
        "project [old]",  # glob character class
        "O'Brien lab",  # terminates a SQL string literal
    ],
)
def test_the_store_works_from_a_path_containing_glob_or_sql_metacharacters(tmp_path, dirname):
    """A checkout directory is not a pattern and not SQL.

    `_sig_glob` interpolated the raw path into both a `glob.glob()` pattern and
    a `read_parquet('...')` literal. A `[` or `]` anywhere above the store made
    the glob match nothing, so `query_signatures` took its empty-glob early exit
    and every gene answered HTTP 200 with `n_signatures: 0` — the dashboard
    reporting the whole corpus as never measured, with no error anywhere. An
    apostrophe instead closed the SQL literal and every query died in the
    DuckDB parser.

    Both are ordinary things to have in a Windows or macOS home directory.
    """
    home = tmp_path / dirname
    home.mkdir()
    store = GeroStore(data_home=home).build(
        prefer_full_clinical=False, prefer_full_signatures=False
    )

    rows = store.query_signatures(species="human")
    assert rows, "the store built but answered nothing — the glob silently missed"
    assert all(r.species == "human" for r in rows)

    # The narrowed path (both partition columns pinned) has to work too.
    narrowed = store.query_signatures(species="human", omic_layer="transcriptome")
    assert narrowed


def test_concurrent_signature_queries_each_get_their_own_results(store):
    """FastAPI runs sync endpoints on a threadpool, so this is the real shape.

    One `DuckDBPyConnection` was shared by every request. Two simultaneous
    `GET /v1/gene/{id}/signature` calls land on different threadpool workers and
    interleave `execute`/`fetchall` on the same handle, so a request can be
    handed the other request's result set. Each query here asks for a different
    species and checks it got its own answer back.
    """
    import concurrent.futures

    expected = {sp: store.query_signatures(species=sp) for sp in ("human", "mouse")}
    assert all(expected.values()), "fixture needs rows for both species"

    def ask(species):
        rows = store.query_signatures(species=species)
        return species, [r.species for r in rows]

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(ask, sp) for sp in ("human", "mouse") * 12]
        for future in concurrent.futures.as_completed(futures):
            species, got = future.result()
            assert got, f"{species} query came back empty under concurrency"
            assert set(got) == {species}, f"{species} query received another query's rows"
            assert len(got) == len(expected[species])
