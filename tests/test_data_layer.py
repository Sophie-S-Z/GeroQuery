"""Manifest, checksum-verifying fetch, and the real NHANES adapter.

Everything here runs offline. The one test that needs live upstream bytes is
marked ``live`` and is excluded from the default run.
"""

from __future__ import annotations

import hashlib
import shutil

import pytest

from geroquery.exceptions import ChecksumMismatchError, NetworkDisabledError, SourceError
from geroquery.sources import fetch as fetch_mod
from geroquery.sources import nhanes
from geroquery.sources.manifest import MANIFEST, NHANES_2017_2018, RemoteArtifact, get_artifact

# ---- manifest -------------------------------------------------------------


def test_every_artifact_is_fully_pinned():
    """A manifest entry without a real digest is worse than no manifest — it
    looks like verification while verifying nothing."""
    assert MANIFEST
    for key, art in MANIFEST.items():
        assert art.key == key
        assert art.url.startswith("https://")
        assert len(art.sha256) == 64 and int(art.sha256, 16) >= 0  # valid hex digest
        assert art.n_bytes > 0
        assert art.license and art.attribution and art.description


def test_nhanes_urls_use_the_datafiles_path():
    """The human-facing /Nchs/Nhanes/<years>/ path serves HTML, not XPORT, so a
    URL drifting back to it produces 'Header record is not an XPORT file'."""
    for art in NHANES_2017_2018.values():
        assert "/Nchs/Data/Nhanes/Public/2017/DataFiles/" in art.url
        assert art.url.endswith(".xpt")


def test_get_artifact_rejects_unknown_key():
    with pytest.raises(KeyError, match="Unknown artifact"):
        get_artifact("NOT_A_REAL_KEY")


# ---- fetch ----------------------------------------------------------------


@pytest.fixture
def fake_artifact(tmp_path):
    """An artifact whose bytes we control, so verification can be tested exactly."""
    payload = b"geroquery-fetch-layer-test-payload\n" * 16
    art = RemoteArtifact(
        key="FAKE",
        url="https://example.invalid/fake.bin",
        sha256=hashlib.sha256(payload).hexdigest(),
        n_bytes=len(payload),
        release="test",
        license="test",
        attribution="test",
        description="test artifact",
    )
    return art, payload, tmp_path / "cache"


def test_fetch_returns_cached_bytes_without_network(fake_artifact):
    art, payload, cache = fake_artifact
    cache.mkdir()
    fetch_mod.cache_path(art, cache).write_bytes(payload)

    path = fetch_mod.fetch_artifact(art, directory=cache, allow_network=False)
    assert path.read_bytes() == payload


def test_fetch_refuses_to_reach_out_when_network_disabled(fake_artifact):
    art, _payload, cache = fake_artifact
    with pytest.raises(NetworkDisabledError) as exc:
        fetch_mod.fetch_artifact(art, directory=cache, allow_network=False)
    assert "make data" in str(exc.value)


def test_corrupt_cache_entry_is_not_trusted(fake_artifact):
    """A cache hit is re-verified, not assumed good. A truncated file must not be
    silently handed back as valid data."""
    art, payload, cache = fake_artifact
    cache.mkdir()
    fetch_mod.cache_path(art, cache).write_bytes(payload[:10])

    assert fetch_mod.is_cached(art, cache) is False
    with pytest.raises(NetworkDisabledError):
        fetch_mod.fetch_artifact(art, directory=cache, allow_network=False)


def test_verify_reports_size_and_digest_mismatch_distinctly(fake_artifact):
    art, payload, cache = fake_artifact
    cache.mkdir()
    path = fetch_mod.cache_path(art, cache)

    path.write_bytes(payload[:-1])
    with pytest.raises(ChecksumMismatchError, match="bytes"):
        fetch_mod.verify(path, art)

    # Right length, wrong content: only the digest can catch this.
    path.write_bytes(b"x" * len(payload))
    with pytest.raises(ChecksumMismatchError, match="SHA-256"):
        fetch_mod.verify(path, art)


def test_failed_verification_leaves_no_cache_entry(fake_artifact, monkeypatch):
    """An interrupted or tampered download must not leave a file that a later
    offline run would pick up as a valid cache hit."""
    art, _payload, cache = fake_artifact

    def _bad_download(artifact, dest, timeout):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"wrong bytes entirely")

    monkeypatch.setattr(fetch_mod, "_download", _bad_download)
    with pytest.raises(ChecksumMismatchError):
        fetch_mod.fetch_artifact(art, directory=cache, allow_network=True)

    assert not fetch_mod.cache_path(art, cache).exists()
    assert not list(cache.glob("*.part"))


def test_successful_download_is_verified_then_placed(fake_artifact, monkeypatch):
    art, payload, cache = fake_artifact

    def _good_download(artifact, dest, timeout):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(payload)

    monkeypatch.setattr(fetch_mod, "_download", _good_download)
    path = fetch_mod.fetch_artifact(art, directory=cache, allow_network=True)
    assert path.read_bytes() == payload
    assert fetch_mod.is_cached(art, cache)


# ---- NHANES adapter -------------------------------------------------------


def test_marker_map_matches_the_clock_feature_contract():
    """The six markers are exactly what the clinical clocks consume. If these
    drift apart, clock application silently loses features."""
    from geroquery.clocks.registry import CLINICAL_FEATURES

    assert set(nhanes.MARKERS) == set(CLINICAL_FEATURES)


def test_offline_sample_is_real_nhanes_rows():
    df = nhanes.load_sample()
    assert len(df) == nhanes.SAMPLE_SIZE
    assert df["subject_id"].str.fullmatch(r"NHANES:\d+").all()
    assert set(nhanes.MARKERS) <= set(df.columns)
    assert "survey_weight" in df.columns
    assert df["age"].min() >= nhanes.MIN_AGE
    assert df["age"].max() <= nhanes.AGE_TOPCODE
    assert df[list(nhanes.MARKERS)].notna().all().all()  # complete cases only


def test_adapter_falls_back_to_sample_and_says_so(tmp_path, monkeypatch):
    """With no cache and no network, the adapter must return real sample rows and
    report mode='sample' — a caller has to be able to tell which cohort a number
    came from."""
    monkeypatch.setattr(fetch_mod.config, "CACHE_HOME", tmp_path / "empty")
    monkeypatch.setattr(fetch_mod.config, "ALLOW_NETWORK", False)

    # Point at a data dir holding only the committed sample. Without this the
    # test passes or fails on whether the developer has run `make data`: the
    # adapter prefers a prebuilt clinical_nhanes_full.csv, which is git-ignored
    # and therefore present locally but absent in CI.
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    shutil.copy(nhanes.sample_path(), nhanes.sample_path(data_dir))

    frame, mode = nhanes.NhanesClinicalSource(data_dir).clinical_frame()
    assert mode == "sample"
    assert len(frame) == nhanes.SAMPLE_SIZE


def test_adapter_prefers_a_prebuilt_full_table(tmp_path):
    """`make data` writes the full table; the adapter must use it rather than
    re-parsing XPORT on every store build."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    sample = nhanes.load_sample()
    shutil.copy(nhanes.sample_path(), nhanes.sample_path(data_dir))
    sample.head(11).to_csv(nhanes.full_path(data_dir), index=False)

    frame, mode = nhanes.NhanesClinicalSource(data_dir=data_dir).clinical_frame()
    assert mode == "full"
    assert len(frame) == 11


def test_missing_sample_raises_with_a_repair_instruction(tmp_path):
    with pytest.raises(SourceError, match="--write-sample"):
        nhanes.load_sample(tmp_path)


def test_nhanes_licence_permits_caching():
    """NHANES is US public domain, so the store's licence gate must pass. This is
    the counterpart to the federate-only sources that it must refuse."""
    nhanes.NhanesClinicalSource().assert_cacheable()


def test_caveats_are_carried_with_the_data():
    joined = " ".join(nhanes.CAVEATS).lower()
    assert "topcod" in joined
    assert "cross-sectional" in joined
    assert "weight" in joined


@pytest.mark.live
def test_full_download_matches_pinned_checksums():
    """Live: the pinned digests still describe what CDC serves today."""
    frame = nhanes.load_full(allow_network=True)
    assert len(frame) == 4895
    assert frame["age"].between(20, 80).all()
