"""The generic gene-ID coverage plotter.

Give it one coverage HDF5 and a gene ID and it draws genome-versus-transcriptome coverage
in the shared transcript coordinate, with the transcript's regions marked. It is the tool
a reader uses on their own data; the paper's Figure 3 A and B are two invocations of it.

Two properties are tested here rather than assumed:

  * it reads ONE file and builds nothing -- no BAM, no GTF, no pipeline invocation;
  * an unusable input fails with a message that says what to do, not a traceback.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

import plot_transcript_coverage as plotter
from conftest import GENE_IDS, TX_MINUS, TX_PLUS, TX_SHORT

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "code" / "panels" / "plot_transcript_coverage.py"


# ── it is a reader, not a pipeline ───────────────────────────────────────────

def test_the_plotter_never_imports_a_bam_library():
    """A plotting program that can open a BAM is a plotting program that will one day
    quietly recompute something."""
    tree = ast.parse(SOURCE.read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    for forbidden in ("pysam", "pyranges", "build_shared_coverage",
                      "build_cohort_coverage"):
        assert forbidden not in imported, "the plotter imports %s" % forbidden


def test_the_plotter_takes_no_bam_or_annotation_argument():
    text = SOURCE.read_text()
    for flag in ("--genome-bam", "--transcriptome-bam", "--bam", "--gtf", "--appris"):
        assert '"%s"' % flag not in text, "the plotter accepts %s" % flag


def test_the_plotter_never_runs_the_builder():
    """It may PRINT the build command as guidance. It must not execute it."""
    tree = ast.parse(SOURCE.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            assert not (node.value.id == "subprocess"), "the plotter shells out"
        if isinstance(node, ast.Name):
            assert node.id not in ("system", "popen"), node.id


def test_the_coverage_file_flag_is_required_and_named_coverage_h5():
    parser_text = SOURCE.read_text()
    assert '"--coverage-h5", required=True' in parser_text


# ── pre-flight validation ────────────────────────────────────────────────────

def test_a_missing_file_is_refused_with_the_build_command(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        plotter.check_coverage_file(tmp_path / "absent.h5")
    message = str(excinfo.value)
    assert "does not exist" in message
    assert "build_shared_coverage.py" in message, "the message must say how to make one"


def test_a_file_that_is_not_hdf5_is_refused(tmp_path):
    fake = tmp_path / "fake.h5"
    fake.write_text("definitely not HDF5")
    with pytest.raises(SystemExit) as excinfo:
        plotter.check_coverage_file(fake)
    assert "not a usable coverage file" in str(excinfo.value)


def test_a_valid_file_reports_its_identity(coverage_path):
    identity = plotter.check_coverage_file(coverage_path)
    assert identity["sample"] == "SYN"
    assert identity["assay"] == "ribo"
    assert identity["routes"] == ["genome", "transcriptome"]
    assert identity["psite_placement"] == "cigar_aware"
    assert identity["coordinate_system"] == "transcript_5p_to_3p"


def test_the_expected_sample_is_enforced_when_given(coverage_path):
    """Plotting the wrong cell line is a mistake that produces a perfectly nice figure."""
    plotter.check_coverage_file(coverage_path, expect_sample="SYN")
    with pytest.raises(SystemExit) as excinfo:
        plotter.check_coverage_file(coverage_path, expect_sample="HeLa")
    assert "but 'HeLa' was expected" in str(excinfo.value)


# ── selection by identifier ──────────────────────────────────────────────────

def test_a_gene_id_resolves_to_its_transcript(coverage_path):
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS])
    assert tracks["transcript_id"] == TX_PLUS


def test_an_unversioned_gene_id_resolves(coverage_path):
    tracks = plotter.load_tracks(coverage_path,
                                 gene_id=GENE_IDS[TX_PLUS].split(".", 1)[0])
    assert tracks["transcript_id"] == TX_PLUS


def test_a_transcript_id_resolves(coverage_path):
    tracks = plotter.load_tracks(coverage_path, transcript_id=TX_MINUS)
    assert tracks["transcript_id"] == TX_MINUS


def test_neither_identifier_is_an_error(coverage_path):
    with pytest.raises(SystemExit) as excinfo:
        plotter.load_tracks(coverage_path)
    assert "--gene-id" in str(excinfo.value)


def test_an_unknown_gene_is_an_error_naming_it(coverage_path):
    import coverage_schema
    with pytest.raises(coverage_schema.SchemaError) as excinfo:
        plotter.load_tracks(coverage_path, gene_id="ENSG_NOT_HERE")
    assert "ENSG_NOT_HERE" in str(excinfo.value)


def test_a_gene_name_is_not_a_selection_key(coverage_path):
    """Gene NAME is display metadata: neither unique nor stable. Accepting it would make
    the selection silently ambiguous."""
    import coverage_schema
    with pytest.raises(coverage_schema.SchemaError):
        plotter.load_tracks(coverage_path, gene_id="PLUSGENE")


# ── the region overlay ───────────────────────────────────────────────────────

def test_the_default_overlay_is_the_canonical_regions(coverage_path):
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS])
    assert tracks["overlay"] == "canonical"
    assert sorted(row[0] for row in plotter.overlay_intervals(tracks)) == \
        ["CDS", "UTR3", "UTR5"]


def test_the_overlay_can_be_switched_off(coverage_path):
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS], overlay="none")
    assert tracks["overlay"] == "none"
    assert plotter.overlay_intervals(tracks) == []


def test_the_regions_tile_the_transcript_in_order(coverage_path):
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS])
    intervals = plotter.overlay_intervals(tracks)
    assert intervals[0][1] == 0
    assert intervals[-1][2] == tracks["transcript_len"]
    for (_a, _s1, end), (_b, start, _e2) in zip(intervals, intervals[1:]):
        assert end == start, "the regions must be contiguous"


def test_the_cds_region_is_the_stored_bounds(coverage_path):
    import coverage_schema
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS])
    with coverage_schema.open_coverage(coverage_path) as coverage:
        index = coverage.index_of_transcript(TX_PLUS)
        assert tracks["regions"]["CDS"] == (int(coverage.cds_start[index]),
                                            int(coverage.cds_end[index]))


def test_asking_for_an_overlay_the_file_cannot_supply_is_an_error():
    """A silent downgrade would produce a plot labelled with regions that were never
    stored."""
    with pytest.raises(SystemExit) as excinfo:
        plotter.resolve_overlay("canonical", {})
    assert "no CDS bounds" in str(excinfo.value)


def test_auto_falls_back_to_none_without_regions():
    assert plotter.resolve_overlay("auto", {"CDS": (0, 10)}) == "canonical"
    assert plotter.resolve_overlay("auto", {}) == "none"


# ── the plotted window ───────────────────────────────────────────────────────

def test_the_whole_transcript_is_the_default_window(coverage_path):
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS])
    assert tracks["slice"] == [0, 200]
    assert tracks["x"][0] == 0 and tracks["x"][-1] == 199


def test_the_cds_window_uses_the_files_own_trim(coverage_path):
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS], region="cds")
    assert tracks["trim"] == 15
    # CDS-relative axis: the window runs [trim, cds_len - trim)
    assert tracks["x_start"] == 15
    assert tracks["x_end"] == 96 - 15


def test_a_transcript_whose_cds_does_not_survive_the_trim_is_refused(coverage_path):
    with pytest.raises(SystemExit) as excinfo:
        plotter.load_tracks(coverage_path, transcript_id=TX_SHORT, region="cds")
    assert "does not survive" in str(excinfo.value)


def test_normalization_is_explicit_and_does_not_touch_the_raw_arrays(coverage_path):
    raw = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS])
    scaled = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS],
                                 normalize="max")
    assert (raw["values"]["g_ps"] == raw["raw"]["g_ps"]).all()
    assert (scaled["raw"]["g_ps"] == raw["raw"]["g_ps"]).all(), \
        "normalization must not alter the recorded raw counts"
    assert scaled["values"]["g_ps"].max() == pytest.approx(1.0)


def test_correlations_are_computed_on_the_raw_counts(coverage_path):
    """Correlating the NORMALIZED values would silently change what the number means."""
    import compute_coverage_concordance as ccc
    tracks = plotter.load_tracks(coverage_path, gene_id=GENE_IDS[TX_PLUS],
                                 normalize="max")
    correlations = plotter.annotate_correlations(tracks)
    assert correlations["psite"]["spearman"] == ccc._spear(tracks["raw"]["g_ps"],
                                                           tracks["raw"]["t_ps"])


# ── end to end ───────────────────────────────────────────────────────────────

def test_the_cli_renders_and_records_what_it_did(coverage_path, tmp_path):
    """`render()` returns the record describing the figure it just drew.

    That record used to be written beside the panel as a `.inputs.json` sidecar. Nothing
    downstream read those files, so they are no longer produced -- but the facts in them
    still matter, so they are asserted here against the returned value instead.
    """
    output = tmp_path / "gapdh"
    record = plotter.render(["--coverage-h5", str(coverage_path),
                             "--gene-id", GENE_IDS[TX_PLUS],
                             "--output", str(output), "--format", "png",
                             "--annotate-correlation"])
    assert output.with_suffix(".png").exists()
    assert not output.with_suffix(".inputs.json").exists(), \
        "the panel must not write a JSON sidecar any more"

    assert record["resolved"]["transcript_id"] == TX_PLUS
    assert record["region_overlay"] == "canonical"
    assert record["coverage_identity"]["psite_placement"] == "cigar_aware"
    assert record["correlations"]["psite"]["spearman"] is not None
    assert "/" not in record["coverage_file"], "the record should not embed a path"
    assert record["outputs"] == [str(output.with_suffix(".png"))]


def test_the_cli_refuses_to_overwrite_without_force(coverage_path, tmp_path):
    output = tmp_path / "gapdh"
    argv = ["--coverage-h5", str(coverage_path), "--gene-id", GENE_IDS[TX_PLUS],
            "--output", str(output), "--format", "png"]
    assert plotter.main(argv) == 0
    with pytest.raises(SystemExit):
        plotter.main(argv)
    assert plotter.main(argv + ["--force"]) == 0


def test_the_cli_reports_a_missing_coverage_file_before_anything_else(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        plotter.main(["--coverage-h5", str(tmp_path / "nope.h5"),
                      "--gene-id", "ENSGPLUS", "--output", str(tmp_path / "x")])
    assert "does not exist" in str(excinfo.value)
    assert not list(tmp_path.glob("x.*")), "nothing should have been drawn"
