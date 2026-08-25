"""The QC annotation bundle: one file, one fingerprint, atomic write.

The five derived tables (CDS exons, transcript metadata, UTR exons, APPRIS gene bodies,
all-GTF gene bodies) used to be five pickles reused whenever they existed. Existence is not
validity: a changed GTF served the previous run's coordinates, and the numbers that came
out looked entirely plausible. The bundle is keyed by the CONTENT of both inputs, the size
filters, the builder's own source and a schema version.

Synthetic inputs throughout -- a four-line GTF and a two-line APPRIS table -- so the whole
file runs in a second and asserts on identity rather than on biology.
"""
from __future__ import annotations

import os
import pickle
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "common" / "ribo_seq_qc"))

#: A minimal GENCODE-shaped GTF. CDS/UTR rows must carry `gene_type "protein_coding"` and
#: a `transcript_id` present in the APPRIS table, which is exactly what the parser filters
#: on -- so this fixture also documents the two attributes it requires.
_ATTR = ('gene_id "G1"; transcript_id "T1.1"; gene_type "protein_coding"; gene_name "A";')
GTF = "\n".join([
    'chr1\tsrc\tgene\t100\t900\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding";',
    'chr1\tsrc\tCDS\t200\t400\t.\t+\t0\t' + _ATTR,
    'chr1\tsrc\tUTR\t100\t199\t.\t+\t.\t' + _ATTR,
    'chr1\tsrc\tUTR\t401\t900\t.\t+\t.\t' + _ATTR,
]) + "\n"
APPRIS = "T1.1|G1|-|-|A-201|A|801|UTR5:1-100|CDS:101-301|UTR3:302-801|\n"


@pytest.fixture
def env(tmp_path, monkeypatch):
    """A private output root and a private copy of both inputs."""
    gtf = tmp_path / "in.gtf"
    gtf.write_text(GTF)
    appris = tmp_path / "in.tsv"
    appris.write_text(APPRIS)
    monkeypatch.setenv("RIBOFLOW_PAPER_OUT", str(tmp_path / "out"))
    monkeypatch.setenv("RIBOFLOW_PAPER_GTF", str(gtf))
    monkeypatch.setenv("RIBOFLOW_PAPER_APPRIS", str(appris))
    import config
    return {"config": config, "gtf": gtf, "appris": appris, "root": tmp_path}


def rebuilt(config):
    """Call `load_bundle` and report whether it had to rebuild."""
    path = Path(config.bundle_path())
    before = path.stat().st_mtime_ns if path.exists() else None
    config.load_bundle()
    return before != path.stat().st_mtime_ns


def test_the_bundle_is_one_file_with_every_payload(env):
    config = env["config"]
    payloads = config.load_bundle()
    assert set(payloads) == set(config.BUNDLE_PAYLOADS)
    assert len(list(Path(config.cache_dir()).glob("*.pkl"))) == 1, \
        "five separate pickles are what went stale independently"


def test_unchanged_inputs_are_reused(env):
    config = env["config"]
    assert rebuilt(config), "the first call must build"
    assert not rebuilt(config), "a matching fingerprint must be reused"


def test_changed_gtf_content_rebuilds(env):
    config = env["config"]
    config.load_bundle()
    env["gtf"].write_text(GTF + GTF.splitlines()[1] + "\n")
    assert rebuilt(config)


def test_changed_appris_content_rebuilds(env):
    config = env["config"]
    config.load_bundle()
    env["appris"].write_text(APPRIS + APPRIS)
    assert rebuilt(config)


def test_identical_content_at_a_moved_path_is_reused(env, monkeypatch):
    """Keyed by CONTENT, so relocating an input is not a reason to re-parse a GTF."""
    config = env["config"]
    config.load_bundle()
    moved = env["root"] / "moved.tsv"
    shutil.copy(env["appris"], moved)
    monkeypatch.setenv("RIBOFLOW_PAPER_APPRIS", str(moved))
    assert not rebuilt(config)


def test_a_changed_builder_rebuilds(env):
    """The parse rule is part of the identity: same inputs, different code, different
    tables. Simulated by perturbing the fingerprint's builder term directly, so the test
    never edits a file inside the repository."""
    config = env["config"]
    config.load_bundle()
    original = config.bundle_fingerprint
    config.bundle_fingerprint = lambda *a, **k: "a-different-builder-digest"
    try:
        assert rebuilt(config)
    finally:
        config.bundle_fingerprint = original


def test_a_corrupt_bundle_rebuilds_rather_than_half_loading(env):
    config = env["config"]
    config.load_bundle()
    Path(config.bundle_path()).write_bytes(b"not a pickle at all")
    assert rebuilt(config)
    assert set(config.load_bundle()) == set(config.BUNDLE_PAYLOADS)


def test_a_schema_bump_rebuilds(env):
    config = env["config"]
    config.load_bundle()
    path = Path(config.bundle_path())
    with path.open("rb") as handle:
        document = pickle.load(handle)
    document["schema_version"] = config.CACHE_SCHEMA_VERSION + 100
    with path.open("wb") as handle:
        pickle.dump(document, handle)
    assert rebuilt(config)


def test_the_write_is_atomic(env):
    """A failed write must leave the previous bundle intact and no temp file behind.

    `os.replace` is what guarantees it: a reader sees the old bundle or the new one, never
    a truncated pickle. Simulated by making the pickle step raise mid-write."""
    config = env["config"]
    config.load_bundle()
    path = Path(config.bundle_path())
    good = path.read_bytes()

    import pickle as pickle_module
    original = pickle_module.dump

    def explode(*args, **kwargs):
        raise RuntimeError("disk full")

    pickle_module.dump = explode
    try:
        with pytest.raises(RuntimeError):
            config._write_bundle({k: None for k in config.BUNDLE_PAYLOADS},
                                 str(env["gtf"]), str(env["appris"]))
    finally:
        pickle_module.dump = original

    assert path.read_bytes() == good, "the previous bundle must survive a failed write"
    assert not list(Path(config.cache_dir()).glob("*.tmp")), "no temp file left behind"


def test_every_loader_reads_the_same_bundle(env):
    """Five public names, one validated source."""
    config = env["config"]
    for loader, payload in ((config.load_annotation, "appris_cds"),
                            (config.load_appris_meta, "appris_meta"),
                            (config.load_appris_utr, "appris_utr"),
                            (config.load_gene_bodies, "appris_gene_body"),
                            (config.load_all_gene_bodies, "all_gene_bodies")):
        assert loader() is config.load_bundle()[payload] or \
            loader().equals(config.load_bundle()[payload]), payload
    assert len(list(Path(config.cache_dir()).glob("*.pkl"))) == 1
