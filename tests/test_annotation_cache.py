"""The annotation cache: reused only when every input that shapes it is unchanged.

One bundle is shared across 24 coverage builds, so the only thing standing between that and
a wrong answer is the fingerprint. These tests assert it moves when it must and holds when
it must not, and that a reused bundle equals a freshly built one.
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "coverage"))

import annotation_cache as ac  # noqa: E402


@pytest.fixture
def annotation(inputs):
    return dict(gtf=inputs.gtf, appris=inputs.appris, regions=None,
                left_span=35, right_span=10)


def digest(annotation, **overrides):
    merged = dict(annotation, **overrides)
    return ac.fingerprint(merged["gtf"], merged["appris"], merged["regions"],
                          merged["left_span"], merged["right_span"])[0]


# ── the fingerprint ──────────────────────────────────────────────────────────

def test_the_same_inputs_give_the_same_fingerprint(annotation):
    assert digest(annotation) == digest(annotation)


def test_changed_gtf_CONTENT_changes_the_fingerprint(annotation, tmp_path, inputs):
    """By content, not by path or mtime: a rebuilt GTF at the same path with different
    exons must not be served from a stale cache."""
    import conftest
    other = conftest.build_synthetic_gtf(tmp_path / "other.gtf",
                                         exons=conftest.VARIANT_EXONS,
                                         geometry=conftest.VARIANT_GEOMETRY)
    assert digest(annotation, gtf=other) != digest(annotation)


def test_a_byte_identical_copy_at_a_different_path_keeps_the_fingerprint(annotation,
                                                                         tmp_path):
    copy = tmp_path / "copied.gtf"
    copy.write_bytes(Path(annotation["gtf"]).read_bytes())
    assert digest(annotation, gtf=copy) == digest(annotation), \
        "identity is by content, so a move must not invalidate the cache"


def test_each_span_parameter_changes_the_fingerprint(annotation):
    assert digest(annotation, left_span=30) != digest(annotation)
    assert digest(annotation, right_span=15) != digest(annotation)


def test_supplying_a_regions_bed_changes_the_fingerprint(annotation, inputs):
    assert digest(annotation, regions=inputs.regions) != digest(annotation)


def test_the_fingerprint_covers_the_code_that_builds_the_bundle(annotation, monkeypatch,
                                                                tmp_path):
    """A changed transcript_coords.py can change the bundle without any input changing,
    so the module digests are part of the fingerprint."""
    assert "transcript_coords.py" in ac.SOURCE_MODULES
    assert "transcript_regions.py" in ac.SOURCE_MODULES
    before = digest(annotation)
    monkeypatch.setattr(ac, "CACHE_VERSION", ac.CACHE_VERSION + 1)
    assert digest(annotation) != before, "the schema version must be in the fingerprint"


def test_the_recorded_inputs_name_content_not_paths(annotation):
    _digest, inputs_record = ac.fingerprint(
        annotation["gtf"], annotation["appris"], None, 35, 10)
    for name, record in inputs_record.items():
        assert "sha256" in record and "bytes" in record, name
        assert "/" not in record["name"], name


# ── reuse ────────────────────────────────────────────────────────────────────

def test_a_matching_cache_is_reused(annotation, tmp_path):
    cache = tmp_path / "annotation.pkl"
    first, reused = ac.load_or_build(cache, **annotation)
    assert reused is False and cache.exists()
    second, reused = ac.load_or_build(cache, **annotation)
    assert reused is True
    assert list(second["transcripts"]["transcript_id"]) == \
        list(first["transcripts"]["transcript_id"])
    assert second["n_positions"] == first["n_positions"]


def test_a_changed_input_rebuilds_rather_than_reusing(annotation, tmp_path):
    import conftest
    cache = tmp_path / "annotation.pkl"
    ac.load_or_build(cache, **annotation)
    other = conftest.build_synthetic_gtf(tmp_path / "other.gtf", exons=conftest.VARIANT_EXONS,
                                         geometry=conftest.VARIANT_GEOMETRY)
    _bundle, reused = ac.load_or_build(cache, **dict(annotation, gtf=other))
    assert reused is False, "a stale cache must never be served"


def test_a_corrupt_cache_is_rebuilt_not_raised(annotation, tmp_path):
    cache = tmp_path / "annotation.pkl"
    cache.write_bytes(b"not a pickle")
    _bundle, reused = ac.load_or_build(cache, **annotation)
    assert reused is False


def test_the_cache_is_written_atomically(annotation, tmp_path):
    """A half-written bundle at the final path would be read as complete by the next run."""
    cache = tmp_path / "annotation.pkl"
    ac.load_or_build(cache, **annotation)
    assert not list(tmp_path.glob("*.tmp-*")), "a temporary file was left behind"
    with open(cache, "rb") as handle:
        stored = pickle.load(handle)
    assert set(stored) == {"fingerprint", "inputs", "version", "bundle"}
    assert stored["version"] == ac.CACHE_VERSION


def test_no_cache_path_still_builds(annotation):
    bundle, reused = ac.load_or_build(None, **annotation)
    assert reused is False and bundle["n_positions"] > 0


def test_the_bundle_carries_everything_the_builder_needs(annotation):
    bundle, _ = ac.load_or_build(None, **annotation)
    for key in ("headers", "coords", "cds_table", "transcripts", "exons", "n_positions",
                "regions", "ribo_bins", "region_summary", "stop_ids",
                "index_of_id", "index_of_base"):
        assert key in bundle, key


# ── the cache changes nothing about the output ───────────────────────────────

def test_a_cached_build_equals_an_uncached_one(inputs, tmp_path):
    """The whole justification for sharing one bundle across 24 samples: the coverage file
    must be identical whether the annotation was parsed for this sample or reused."""
    import h5py
    import build_shared_coverage
    import conftest

    cache = tmp_path / "annotation.pkl"
    uncached, _ = build_shared_coverage.build(
        conftest.build_config(inputs, output=tmp_path / "uncached"))
    built, first_report = build_shared_coverage.build(
        conftest.build_config(inputs, output=tmp_path / "cached",
                              annotation_cache=cache))
    reused, second_report = build_shared_coverage.build(
        conftest.build_config(inputs, output=tmp_path / "reused",
                              annotation_cache=cache))

    assert first_report["annotation_cache_reused"] is False
    assert second_report["annotation_cache_reused"] is True, "the second build rebuilt it"

    def contents(path):
        datasets, attributes = {}, {}
        with h5py.File(path, "r") as handle:
            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    values = obj[()]
                    datasets[name] = (str(obj.dtype),
                                      values.tolist() if obj.dtype.kind != "O"
                                      else [str(v) for v in values])
            handle.visititems(visit)
            attributes = {k: str(handle.attrs[k]) for k in handle.attrs
                          if k != "created_utc"}
        return datasets, attributes

    baseline = contents(uncached)
    assert len(baseline[0]) == 11, "the comparison must cover the whole file"
    assert contents(built) == baseline
    assert contents(reused) == baseline
