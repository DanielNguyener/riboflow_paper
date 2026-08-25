"""The coverage HDF5 container: schema, atomicity, overflow, derived regions and slicing.

Synthetic and BAM-free. Everything the writer guards against is exercised as a FAILURE
path, because a guard that is never tripped in a test is a guard nobody knows works.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
COVERAGE_DIR = REPO / "code" / "coverage"


def _load(name):
    if str(COVERAGE_DIR) not in sys.path:
        sys.path.insert(0, str(COVERAGE_DIR))
    spec = importlib.util.spec_from_file_location(name, COVERAGE_DIR / ("%s.py" % name))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def cs():
    return _load("coverage_schema")


@pytest.fixture(scope="module")
def bsc():
    return _load("build_shared_coverage")


# ── a two-transcript cohort ──────────────────────────────────────────────────
#   TA  20 nt, UTR5 [0,4) CDS [4,16) UTR3 [16,20)
#   TB  10 nt, CDS [0,10)

def _transcripts():
    return pd.DataFrame([
        dict(transcript_id="ENSTA.1", gene_id="ENSGA.1", gene_name="AGENE",
             transcript_len=20, cds_start=4, cds_end=16, coverage_offset=0),
        dict(transcript_id="ENSTB.1", gene_id="ENSGB.1", gene_name="BGENE",
             transcript_len=10, cds_start=0, cds_end=10, coverage_offset=20),
    ])


N_POSITIONS = 30
PROVENANCE = {"parameters": {"paper_cds_trim": 2, "genome_uniqueness": "NH==1",
                             "txome_uniqueness": "MAPQ>=42",
                             "psite_placement": "cigar_aware",
                             "appris_principal_ranks_consumed": False}}


def _writer(cs, path, **overrides):
    kwargs = dict(sample="SYN", transcripts=_transcripts(), provenance=PROVENANCE,
                  paper_cds_trim=2)
    kwargs.update(overrides)
    return cs.CoverageWriter(path, **kwargs)


def _fill(writer, arrays=None):
    arrays = arrays or {}
    for signal in ["genome_psite", "txome_psite", "genome_footprint", "txome_footprint"]:
        writer.write_signal(signal, arrays.get(
            signal, np.zeros(N_POSITIONS, dtype=np.int64)))


@pytest.fixture
def built(cs, tmp_path):
    path = tmp_path / "SYN.shared_coverage.h5"
    psite = np.zeros(N_POSITIONS, dtype=np.int64)
    psite[5] = 3                 # TA, inside its CDS
    psite[25] = 7                # TB, position 5 of its CDS
    with _writer(cs, path) as writer:
        _fill(writer, {"genome_psite": psite})
    return path


# ── round trip ───────────────────────────────────────────────────────────────

def test_a_written_file_validates(cs, built):
    assert cs.validate_file(built) == []


def test_the_file_holds_only_the_arrays_and_the_cds_bounds(cs, built):
    import h5py
    with h5py.File(built, "r") as handle:
        assert set(handle) == {"transcripts", "coverage"}
        assert set(handle["transcripts"]) == set(cs.TRANSCRIPT_COLUMNS)
        assert set(handle["coverage"]) == set(cs.SIGNALS)
        assert set(handle.attrs) == set(cs.ROOT_ATTRS)
        assert handle["transcripts"]["transcript_id"].dtype.kind == "S", \
            "identifiers are fixed-width bytes, not a variable-length string heap"


def test_tracks_are_sliced_by_offset_and_length(cs, built):
    with cs.open_coverage(built) as coverage:
        track = coverage.get_track(0, "genome_psite")
        assert len(track) == 20
        assert int(track.sum()) == 3 and int(track.argmax()) == 5
        assert len(coverage.get_track(1, "genome_psite")) == 10


def test_regions_are_derived_from_the_cds_bounds(cs, built):
    with cs.open_coverage(built) as coverage:
        assert coverage.regions_of(0) == {"UTR5": (0, 4), "CDS": (4, 16), "UTR3": (16, 20)}
        assert coverage.regions_of(1) == {"CDS": (0, 10)}


def test_a_transcript_without_a_cds_has_no_regions(cs):
    assert cs.regions_from_cds(20, cs.NO_CDS, cs.NO_CDS) == {}
    assert cs.regions_from_cds(20, 0, 20) == {"CDS": (0, 20)}


def test_slice_region_uses_the_files_own_trim(cs, built):
    with cs.open_coverage(built) as coverage:
        assert coverage.trim == 2
        assert coverage.slice_region(0, "CDS") == (6, 14)       # trimmed by the file's 2
        assert coverage.slice_region(0, "CDS", trim=0) == (4, 16)


def test_slice_region_collapses_rather_than_inverting(cs, built):
    """A CDS shorter than twice the trim must not produce end < start."""
    with cs.open_coverage(built) as coverage:
        start, end = coverage.slice_region(1, "CDS", trim=99)
        assert end == start


def test_a_missing_region_is_reported_not_invented(cs, built):
    with cs.open_coverage(built) as coverage:
        with pytest.raises(cs.SchemaError) as excinfo:
            coverage.slice_region(1, "UTR5")
        assert "UTR5" in str(excinfo.value)


def test_the_dtype_is_int32(cs, built):
    with cs.open_coverage(built) as coverage:
        assert coverage.get_track(0, "genome_psite").dtype == np.int32


# ── derived counts and the coverage states ───────────────────────────────────

def test_the_coverage_states_are_distinguishable(cs):
    assert cs.describe_coverage_state(0, 0) == cs.STATE_NO_READS
    assert cs.describe_coverage_state(7, 0) == cs.STATE_OUTSIDE_SLICE
    assert cs.describe_coverage_state(3, 3) == cs.STATE_COVERED


def test_event_counts_are_computed_from_the_track(cs, built):
    with cs.open_coverage(built) as coverage:
        assert coverage.event_counts(0)["genome_psite"] == 3
        assert coverage.event_counts(1)["genome_psite"] == 7
        assert coverage.event_counts(1)["txome_psite"] == 0


def test_reads_outside_the_slice_are_not_confused_with_no_reads(cs, built):
    """Transcript B has 7 P-site events but none inside a heavily trimmed CDS."""
    with cs.open_coverage(built) as coverage:
        start, end = coverage.slice_region(1, "CDS", trim=99)
        assert cs.describe_coverage_state(
            coverage.event_counts(1)["genome_psite"],
            int(coverage.get_track(1, "genome_psite")[start:end].sum())
        ) == cs.STATE_OUTSIDE_SLICE


def test_cds_window_sums_match_the_writer_side_definition(cs, built):
    """The concordance keys derive from these: P-sites on the untrimmed CDS, footprints
    on the trimmed interior. Both are exact integer sums."""
    with cs.open_coverage(built) as coverage:
        values = coverage.signal("genome_psite")
        assert coverage.cds_window_sums(values, 0).tolist() == [3, 7]
        assert coverage.cds_window_sums(values, 2).tolist() == [0, 7]   # TA's event is at 5
        assert coverage.cds_window_sums(values, 6).tolist() == [0, 0]
        assert coverage.cds_window_sums(values, 99).tolist() == [0, 0]  # past the end


def test_window_sums_treat_a_missing_cds_as_zero(cs):
    values = np.ones(30, dtype=np.int32)
    sums = cs.window_sums(values, np.array([0, 20]), np.array([4, cs.NO_CDS]),
                          np.array([16, cs.NO_CDS]))
    assert sums.tolist() == [12, 0]


# ── identifier resolution ────────────────────────────────────────────────────

def test_gene_resolves_versioned_and_unversioned(cs, built):
    with cs.open_coverage(built) as coverage:
        assert coverage.resolve_gene("ENSGA.1") == 0
        assert coverage.resolve_gene("ENSGA") == 0


def test_transcript_resolves_versioned_and_unversioned(cs, built):
    with cs.open_coverage(built) as coverage:
        assert coverage.index_of_transcript("ENSTB.1") == 1
        assert coverage.index_of_transcript("ENSTB") == 1


def test_an_unknown_gene_is_an_error(cs, built):
    with cs.open_coverage(built) as coverage:
        with pytest.raises(cs.SchemaError):
            coverage.resolve_gene("ENSGZZZ")


def test_a_transcript_id_from_the_wrong_gene_is_rejected(cs, built):
    with cs.open_coverage(built) as coverage:
        with pytest.raises(cs.SchemaError) as excinfo:
            coverage.resolve_gene("ENSGA", transcript_id="ENSTB.1")
        assert "belongs to gene" in str(excinfo.value)


def test_transcript_info_is_identity_and_geometry_only(cs, built):
    with cs.open_coverage(built) as coverage:
        info = coverage.transcript_info(0)
        assert info == {"transcript_id": "ENSTA.1", "gene_id": "ENSGA.1",
                        "gene_name": "AGENE", "transcript_len": 20,
                        "cds_start": 4, "cds_end": 16, "cds_len": 12}


# ── atomicity and refusal ────────────────────────────────────────────────────

def test_an_interrupted_write_leaves_nothing_at_the_final_path(cs, tmp_path):
    path = tmp_path / "X.h5"
    with pytest.raises(RuntimeError):
        with _writer(cs, path) as writer:
            writer.write_signal("genome_psite", np.zeros(N_POSITIONS, dtype=np.int64))
            raise RuntimeError("simulated crash")
    assert not path.exists()
    assert list(tmp_path.glob("X.h5.tmp-*")) == []


def test_finalize_refuses_when_a_signal_was_never_written(cs, tmp_path):
    path = tmp_path / "Y.h5"
    writer = _writer(cs, path)
    writer.write_signal("genome_psite", np.zeros(N_POSITIONS, dtype=np.int64))
    with pytest.raises(cs.SchemaError) as excinfo:
        writer.finalize()
    assert "never written" in str(excinfo.value)
    writer.abort()
    assert not path.exists()


def test_a_stale_temporary_file_is_refused_not_reused(cs, tmp_path):
    stale = tmp_path / "W.h5.tmp-999-20260101T000000"
    stale.write_bytes(b"partial")
    with pytest.raises(cs.SchemaError) as excinfo:
        _writer(cs, tmp_path / "W.h5")
    assert "stale" in str(excinfo.value)


def test_a_transcript_table_missing_a_column_is_refused(cs, tmp_path):
    with pytest.raises(cs.SchemaError) as excinfo:
        _writer(cs, tmp_path / "C.h5", transcripts=_transcripts().drop(columns="cds_end"))
    assert "cds_end" in str(excinfo.value)


# ── overflow and shape guards ────────────────────────────────────────────────

def test_int32_overflow_is_a_hard_failure(cs, tmp_path):
    writer = _writer(cs, tmp_path / "O.h5")
    values = np.zeros(N_POSITIONS, dtype=np.int64)
    values[3] = 2 ** 31
    with pytest.raises(cs.SchemaError) as excinfo:
        writer.write_signal("genome_psite", values)
    assert "overflow" in str(excinfo.value)
    writer.abort()


def test_negative_coverage_is_rejected(cs, tmp_path):
    writer = _writer(cs, tmp_path / "N.h5")
    values = np.zeros(N_POSITIONS, dtype=np.int64)
    values[1] = -1
    with pytest.raises(cs.SchemaError):
        writer.write_signal("genome_psite", values)
    writer.abort()


def test_a_wrongly_sized_signal_is_rejected(cs, tmp_path):
    writer = _writer(cs, tmp_path / "S.h5")
    with pytest.raises(cs.SchemaError) as excinfo:
        writer.write_signal("genome_psite", np.zeros(N_POSITIONS + 1, dtype=np.int64))
    assert "shape" in str(excinfo.value)
    writer.abort()


def test_an_unknown_signal_name_is_rejected(cs, tmp_path):
    writer = _writer(cs, tmp_path / "U.h5")
    with pytest.raises(cs.SchemaError):
        writer.write_signal("not_a_signal", np.zeros(N_POSITIONS, dtype=np.int64))
    writer.abort()


# ── validator ────────────────────────────────────────────────────────────────

def test_validate_rejects_a_missing_file(cs, tmp_path):
    assert cs.validate_file(tmp_path / "nope.h5")


def test_validate_rejects_a_non_hdf5_file(cs, tmp_path):
    path = tmp_path / "junk.h5"
    path.write_bytes(b"not hdf5")
    assert cs.validate_file(path)


def test_validate_rejects_an_earlier_schema(cs, built):
    import h5py
    with h5py.File(built, "r+") as handle:
        handle.attrs["schema"] = "riboflow_paper/shared-coverage/2"
    problems = cs.validate_file(built)
    assert any("schema" in p for p in problems)


def test_validate_catches_cds_bounds_outside_the_transcript(cs, built):
    import h5py
    with h5py.File(built, "r+") as handle:
        handle["transcripts"]["cds_end"][0] = 21          # past a 20 nt transcript
    problems = cs.validate_file(built)
    assert any("CDS bounds" in p for p in problems)


def test_validate_catches_a_half_missing_cds(cs, built):
    import h5py
    with h5py.File(built, "r+") as handle:
        handle["transcripts"]["cds_start"][1] = cs.NO_CDS
    problems = cs.validate_file(built)
    assert any("cds_start == cds_end == -1" in p for p in problems)


def test_validate_catches_unsorted_transcript_order(cs, built):
    import h5py
    with h5py.File(built, "r+") as handle:
        handle["transcripts"]["transcript_id"][0] = b"ZZZZ.9"
    problems = cs.validate_file(built)
    assert any("sorted" in p for p in problems)


def test_validate_catches_a_broken_offset(cs, built):
    import h5py
    with h5py.File(built, "r+") as handle:
        handle["transcripts"]["coverage_offset"][1] = 19
    problems = cs.validate_file(built)
    assert any("running sum" in p for p in problems)


def test_opening_an_invalid_file_raises_rather_than_failing_later(cs, built):
    import h5py
    with h5py.File(built, "r+") as handle:
        handle["transcripts"]["cds_end"][0] = 21
    with pytest.raises(cs.SchemaError):
        cs.open_coverage(built)


# ── accumulators ─────────────────────────────────────────────────────────────

def test_points_accumulate_into_int32(bsc):
    counts = bsc.accumulate_points(np.array([0, 0, 3], dtype=np.int64), 5)
    assert counts.dtype == np.int32
    assert counts.tolist() == [2, 0, 0, 1, 0]


def test_intervals_accumulate_as_depth(bsc):
    depth = bsc.accumulate_intervals(np.array([0, 2]), np.array([3, 5]), 6)
    assert depth.tolist() == [1, 1, 2, 1, 1, 0]


def test_an_interval_escaping_its_transcript_is_caught(bsc):
    """The difference array must balance; an end past the coordinate would leak depth."""
    with pytest.raises(bsc.BuildError) as excinfo:
        bsc.accumulate_intervals(np.array([0]), np.array([6]), 5)
    assert "outside the coordinate" in str(excinfo.value)


def test_a_missing_cds_is_stored_as_minus_one(bsc, cs):
    regions = pd.DataFrame([dict(transcript_index=0, label="CDS", start=4, end=16)])
    starts, ends = bsc._cds_windows(regions, 2)
    assert starts.tolist() == [4, cs.NO_CDS] and ends.tolist() == [16, cs.NO_CDS]


# ── the plotting-side contract ───────────────────────────────────────────────

def _plotter():
    panels = REPO / "code" / "panels"
    for directory in (str(panels), str(COVERAGE_DIR)):
        if directory not in sys.path:
            sys.path.insert(0, directory)
    spec = importlib.util.spec_from_file_location(
        "plot_transcript_coverage", panels / "plot_transcript_coverage.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_panel_generator_renders_from_a_coverage_file(built, tmp_path):
    module = _plotter()
    tracks = module.load_tracks(built, gene_id="ENSGA", region="whole")
    assert tracks["transcript_id"] == "ENSTA.1"
    assert tracks["gene_name"] == "AGENE"
    assert len(tracks["x"]) == 20
    for key in ("g_ps", "t_ps", "g_fp", "t_fp"):
        assert len(tracks["raw"][key]) == 20
    assert tracks["states"]["g_ps"] == "covered"
    assert tracks["states"]["t_ps"] == "no_reads_assigned"

    import matplotlib.pyplot as plt
    figure, axes = module.plot_coverage(tracks, signal="both")
    try:
        assert len(axes) == 2
        figure.savefig(tmp_path / "panel.png", dpi=40)
    finally:
        plt.close(figure)
    assert (tmp_path / "panel.png").exists()


def test_the_panel_generator_uses_the_files_own_trim_for_the_cds_view(built):
    module = _plotter()
    tracks = module.load_tracks(built, gene_id="ENSGA", region="cds")
    # CDS is [4, 16) and the file's trim is 2, so the CDS-RELATIVE axis runs [2, 10)
    assert (tracks["x_start"], tracks["x_end"]) == (2, 10)
    assert tracks["slice"] == [6, 14]
    assert len(tracks["raw"]["g_ps"]) == 8


def test_the_panel_generator_refuses_an_ambiguous_gene(cs, tmp_path):
    """Two transcripts of one gene, no --transcript-id: it must stop and list them."""
    transcripts = _transcripts()
    transcripts.loc[1, "gene_id"] = transcripts.loc[0, "gene_id"]      # same gene now
    path = tmp_path / "AMBIG.shared_coverage.h5"
    with _writer(cs, path, sample="AMBIG", transcripts=transcripts) as writer:
        _fill(writer)
    module = _plotter()
    schema_module, _metrics = module._import_coverage_modules()
    with pytest.raises(schema_module.SchemaError) as excinfo:
        module.load_tracks(path, gene_id="ENSGA")
    message = str(excinfo.value)
    assert "ENSTA.1" in message and "ENSTB.1" in message
    assert "not guessed" in message
