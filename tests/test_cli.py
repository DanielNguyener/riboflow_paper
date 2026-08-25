"""Every entry point: `--help` works, and a missing input fails early and says what.

A pipeline that runs for twenty minutes and then dies on a typo in a path has wasted
twenty minutes. Each program here is checked to refuse before it computes, and to name
the flag rather than raise a traceback from inside a library.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

#: Every program with a command-line interface, and nothing else.
#:
#: `transcript_coords`, `transcript_regions`, `annotation_cache` and `psite_placement` are
#: NOT here: no workflow runs them as programs, so they are libraries and their argparse
#: layers are gone. A module that only advertises `--help` is an interface to maintain
#: with no caller to justify it.
ENTRY_POINTS = [
    "code/make_tables.py",
    "code/make_panels.py",
    "code/assemble_figures.py",
    "code/make_figures.py",
    "code/coverage/build_shared_coverage.py",
    "code/coverage/build_cohort_coverage.py",
    "code/coverage/compute_coverage_concordance.py",
    "code/coverage/coverage_schema.py",
    "code/coverage/validate_assignment_policy.py",
    "code/panels/plot_transcript_coverage.py",
    "code/alignment_fate/build_gene_read_partition.py",
    "code/alignment_fate/build_gene_partition_data.py",
    "code/alignment_fate/build_locus_data.py",
    "code/panels/plot_gene_read_partition.py",
    "code/ribo_rna/count_transcript_reads.py",
    "code/ribo_rna/build_count_matrices.py",
    "code/te_route/plot_te_route_panels.py",
    "code/common/build_orf_catalog.py",
]

#: CLIs that exist only so an orchestrator can spawn them as subprocesses. They are given a
#: complete argument list by `make_tables.py` or `make_panels.py` and are not documented for
#: direct use, so they are not held to the bare-`--help` contract above.
ORCHESTRATED = [
    "code/panels/plot_cds_periodicity_difference.py",
    "code/panels/plot_fig05_plos_panels.py",
    "code/panels/plot_gene_partition.py",
    "code/panels/plot_locus_coverage.py",
    "code/panels/plot_multimap_biotype.py",
    "code/panels/plot_nonselected_isoform_reach.py",
    "code/panels/plot_per_transcript_concordance.py",
    "code/panels/plot_pooled_concordance.py",
    "code/panels/plot_read_id_union.py",
    "code/panels/plot_readlen_psite_selection.py",
    "code/panels/plot_route_read_counts.py",
    "code/read_taxonomy/compute_concordance.py",
    "code/read_taxonomy/compute_reach.py",
    "code/read_taxonomy/compute_taxonomy.py",
    "code/read_taxonomy/compute_tie_biotype.py",
    "code/read_taxonomy/run_read_taxonomy.py",
    "code/ribo_seq_qc/determine_offset_method.py",
    "code/ribo_seq_qc/run_pipeline.py",
    "code/ribo_seq_qc/run_transcriptome_qc.py",
]


def bare_environment():
    """The caller's environment with every `RIBOFLOW_PAPER_*` variable removed.

    These programs are configurable entirely from the environment, which is the documented
    way to run them. So a test that asks "does this refuse with no arguments?" while the
    real GTF and APPRIS are exported is not asking about the code -- it is asking about the
    operator's shell, and it gets a different answer depending on who runs it. Worse, the
    answer it got with those variables set was `build_orf_catalog.py` running to completion
    and writing into the repository's own `results/`.

    Stripping the prefix here makes every verdict below a property of the program.
    """
    environment = dict(os.environ)
    for name in [k for k in environment if k.startswith("RIBOFLOW_PAPER_")]:
        del environment[name]
    return environment


def run(relative, *args):
    return subprocess.run([sys.executable, str(REPO / relative), *args],
                          capture_output=True, text=True, cwd=str(REPO),
                          env=bare_environment())


@pytest.mark.parametrize("program", ENTRY_POINTS)
def test_every_entry_point_exists(program):
    assert (REPO / program).exists(), program


@pytest.mark.parametrize("program", ENTRY_POINTS)
def test_every_entry_point_has_help(program):
    result = run(program, "--help")
    assert result.returncode == 0, result.stderr[-2000:]
    assert "usage:" in result.stdout.lower()


@pytest.mark.parametrize("program", ENTRY_POINTS)
def test_no_entry_point_runs_with_no_arguments(program):
    """A pipeline stage with a safe zero-argument default is a stage that can quietly do
    the wrong thing. Every one of these must be told what to work on."""
    result = run(program)
    assert result.returncode != 0, "%s did something with no arguments" % program


# ── build_shared_coverage ────────────────────────────────────────────────────

REQUIRED_BUILD_FLAGS = ["--sample", "--genome-bam", "--transcriptome-bam",
                        "--gtf", "--appris", "--qc-genome", "--qc-txome"]


@pytest.mark.parametrize("flag", REQUIRED_BUILD_FLAGS)
def test_build_shared_coverage_demands_each_required_flag(flag):
    argv = []
    for name in REQUIRED_BUILD_FLAGS:
        if name != flag:
            argv += [name, "x"]
    result = run("code/coverage/build_shared_coverage.py", *argv)
    assert result.returncode != 0
    assert flag in result.stderr, result.stderr[-800:]


def test_build_shared_coverage_names_every_missing_input_at_once(tmp_path):
    result = run("code/coverage/build_shared_coverage.py",
                 "--sample", "X",
                 "--genome-bam", str(tmp_path / "a.bam"),
                 "--transcriptome-bam", str(tmp_path / "b.bam"),
                 "--gtf", str(tmp_path / "c.gtf"),
                 "--appris", str(tmp_path / "d.tsv"),
                 "--qc-genome", str(tmp_path / "e.csv"),
                 "--qc-txome", str(tmp_path / "f.csv"))
    assert result.returncode != 0
    message = result.stderr + result.stdout
    assert "do not exist" in message
    for flag in ("--genome-bam", "--transcriptome-bam", "--gtf", "--appris",
                 "--qc-genome", "--qc-txome"):
        assert flag in message, "%s was not listed" % flag


def test_build_shared_coverage_has_no_placement_choice():
    result = run("code/coverage/build_shared_coverage.py", "--help")
    assert "--psite-placement" not in result.stdout
    assert "reference_offset" not in result.stdout


def test_build_shared_coverage_does_not_record_paths_by_default():
    result = run("code/coverage/build_shared_coverage.py", "--help")
    assert "--record-input-paths" in result.stdout


# ── make_tables ──────────────────────────────────────────────────────────────

def test_make_tables_requires_bams():
    result = run("code/make_tables.py", "--all")
    assert result.returncode != 0
    assert "--bams is required" in result.stderr


def test_make_tables_rejects_a_bams_path_that_is_not_a_directory(tmp_path):
    not_a_directory = tmp_path / "file.txt"
    not_a_directory.write_text("")
    result = run("code/make_tables.py", "--bams", str(not_a_directory), "--all")
    assert result.returncode != 0
    assert "not a directory" in result.stderr


def test_make_tables_rejects_an_unknown_stage(tmp_path):
    result = run("code/make_tables.py", "--bams", str(tmp_path), "--stages", "nope")
    assert result.returncode != 0
    assert "unknown stage" in result.stderr


def test_make_tables_requires_a_stage_selection(tmp_path):
    result = run("code/make_tables.py", "--bams", str(tmp_path))
    assert result.returncode != 0
    assert "--all or --stages" in result.stderr


def test_make_tables_lists_its_stages_in_help():
    result = run("code/make_tables.py", "--help")
    for stage in ("coverage", "concordance", "taxonomy", "reach", "multimap_biotype",
                  "te_counts", "te_normalize", "te_stats", "gene_partition", "locus"):
        assert stage in result.stdout


def _make_tables():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "make_tables", REPO / "code" / "make_tables.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_stage_runs_after_everything_it_declares_a_dependency_on():
    """The stage table is the only schedule. If a stage names a dependency, the order it
    is listed in has to honour it -- `reach` reads both the taxonomy and the
    alignment-concordance masters and cannot be scheduled before either."""
    make_tables = _make_tables()
    order = make_tables.STAGE_ORDER
    for name, _run, needs, _anno, _outputs in make_tables.STAGES:
        for dependency in needs:
            assert dependency in order, "%s depends on unknown stage %s" % (name, dependency)
            assert order.index(dependency) < order.index(name), \
                "%s is scheduled before its dependency %s" % (name, dependency)
    assert order[0] == "annotation"
    assert set(order) == set(make_tables.STAGE_RUN)


def test_asking_for_a_stage_pulls_in_its_dependencies():
    make_tables = _make_tables()
    assert make_tables.required_stages({"reach"}) == [
        "annotation", "taxonomy", "alignment_concordance", "reach"]
    assert make_tables.required_stages({"annotation"}) == ["annotation"]


def test_every_declared_output_has_a_shipped_counterpart_and_vice_versa():
    """`data/` mirrors the output root exactly, which is what makes `--compare` a prefix
    swap. A stage that declares an output nothing ships, or a shipped file neither a stage
    nor the exemption list accounts for, means one of the two drifted.

    Every shipped table has a stage behind it, except the third-party inputs
    `EXTERNAL_INPUTS` names with their source: a shipped file no stage builds and no
    source is recorded for cannot be regenerated, and a table nobody can regenerate is not
    a result."""
    make_tables = _make_tables()
    accounted = set(make_tables.OUTPUTS) | set(make_tables.EXTERNAL_INPUTS)
    shipped = {str(p.relative_to(REPO / "data"))
               for p in (REPO / "data").rglob("*") if p.is_file() and p.suffix != ".md"}
    assert accounted == shipped, (
        "declared but not shipped: %s; shipped but not declared: %s"
        % (sorted(accounted - shipped), sorted(shipped - accounted)))


def test_every_shipped_table_now_has_a_stage_behind_it():
    """The named-exemption list is gone: nothing under `data/` is hand-produced any more."""
    make_tables = _make_tables()
    assert not hasattr(make_tables, "SHIPPED_NOT_GENERATED_HERE")


def test_the_shipped_path_is_derived_not_tabulated():
    make_tables = _make_tables()
    for relative in make_tables.OUTPUTS:
        assert make_tables.shipped_for(relative) == REPO / "data" / relative


# ── the annotation inputs have no guessed default ────────────────────────────

def test_the_gtf_has_no_default_and_says_so(monkeypatch):
    sys.path.insert(0, str(REPO / "code" / "common" / "ribo_seq_qc"))
    import config
    monkeypatch.delenv("RIBOFLOW_PAPER_GTF", raising=False)
    with pytest.raises(config.AnnotationError) as excinfo:
        config.gtf_path()
    message = str(excinfo.value)
    assert "RIBOFLOW_PAPER_GTF" in message
    assert "--gtf" in message


def test_the_appris_table_has_no_default_and_says_so(monkeypatch):
    sys.path.insert(0, str(REPO / "code" / "common" / "ribo_seq_qc"))
    import config
    monkeypatch.delenv("RIBOFLOW_PAPER_APPRIS", raising=False)
    with pytest.raises(config.AnnotationError) as excinfo:
        config.appris_path()
    assert "RIBOFLOW_PAPER_APPRIS" in str(excinfo.value)


def test_the_bam_root_has_no_default_and_says_so(monkeypatch):
    sys.path.insert(0, str(REPO / "code" / "common"))
    import bam_inputs
    monkeypatch.delenv("RIBOFLOW_PAPER_BAMS", raising=False)
    with pytest.raises(bam_inputs.InputError) as excinfo:
        bam_inputs.bams_root()
    message = str(excinfo.value)
    assert "RIBOFLOW_PAPER_BAMS" in message
    assert "make_tables.py" in message, "the message should say how to set it"


def test_the_output_root_does_have_a_default(monkeypatch):
    """Output is different from input: writing to `results/` is a safe default, and
    demanding it every time would be noise."""
    sys.path.insert(0, str(REPO / "code" / "common"))
    import bam_inputs
    monkeypatch.delenv("RIBOFLOW_PAPER_OUT", raising=False)
    assert bam_inputs.output_root() == REPO / "results"


def test_the_output_root_is_overridable(monkeypatch, tmp_path):
    sys.path.insert(0, str(REPO / "code" / "common"))
    import bam_inputs
    monkeypatch.setenv("RIBOFLOW_PAPER_OUT", str(tmp_path))
    assert bam_inputs.output_root() == tmp_path


# ── the other BAM-reading programs ───────────────────────────────────────────

def test_count_transcript_reads_names_its_missing_inputs(tmp_path):
    """The annotation and the QC tables are inputs too, and a missing one must name itself.

    The CDS this program counts into comes from the GTF and the read-length window from the
    QC master, so neither is optional and neither has a guessed default."""
    result = run("code/ribo_rna/count_transcript_reads.py", "--sample", "X",
                 "--ribo-genome-bam", str(tmp_path / "a"),
                 "--ribo-txome-bam", str(tmp_path / "b"),
                 "--rna-genome-bam", str(tmp_path / "c"),
                 "--rna-txome-bam", str(tmp_path / "d"),
                 "--gtf", str(tmp_path / "g"), "--appris", str(tmp_path / "p"),
                 "--qc-genome", str(tmp_path / "qg"), "--qc-txome", str(tmp_path / "qt"))
    assert result.returncode != 0
    message = result.stderr + result.stdout
    assert "do not exist" in message
    for flag in ("--rna-txome-bam", "--gtf", "--qc-txome"):
        assert flag in message


def test_count_transcript_reads_has_no_region_option():
    """Every count it produces is a CDS count; there is no whole-transcript mode to pick.

    The option existed, defaulted to `whole`, and the published table was built with it --
    which is how a coding-density figure came to count UTR reads."""
    result = run("code/ribo_rna/count_transcript_reads.py", "--help")
    assert result.returncode == 0
    # `--regions` (the actual-regions BED) is a different, still-valid flag, so match the
    # whole option token rather than a prefix of it. Nor is there a `whole` mode left to
    # choose: argparse accepts unambiguous abbreviations, so `--region whole` would quietly
    # bind to `--regions` if the choice string survived anywhere.
    assert not re.search(r"--region\b", result.stdout)
    assert "whole" not in result.stdout


def test_compute_coverage_concordance_requires_a_coverage_directory():
    result = run("code/coverage/compute_coverage_concordance.py")
    assert result.returncode != 0
    assert "--coverage" in result.stderr


def test_coverage_schema_validate_reports_a_missing_file(tmp_path):
    result = run("code/coverage/coverage_schema.py", "--validate",
                 str(tmp_path / "nope.h5"))
    assert result.returncode != 0
    assert "does not exist" in (result.stdout + result.stderr)


def test_coverage_schema_validate_rejects_a_file_that_is_not_hdf5(tmp_path):
    fake = tmp_path / "fake.h5"
    fake.write_text("this is not HDF5")
    result = run("code/coverage/coverage_schema.py", "--validate", str(fake))
    assert result.returncode != 0
    assert "cannot open as HDF5" in (result.stdout + result.stderr)


# ── --stages accepts the names it advertises ─────────────────────────────────
# Regression: STAGE_RUN is keyed by stage name while STAGES is a list of tuples, so
# validating a name against STAGES rejected every valid stage. These go through the real
# command line, because the bug lived in argument parsing and unit-testing the table
# would not have caught it.

def stages_validate(stages, cwd):
    """`--validate` is read-only, so it is safe to point --bams at any directory."""
    return subprocess.run(
        [sys.executable, str(REPO / "code" / "make_tables.py"),
         "--bams", str(cwd), "--stages", stages, "--validate",
         "--output", str(Path(cwd) / "out")],
        capture_output=True, text=True, cwd=str(REPO))


def test_a_single_valid_stage_is_accepted(tmp_path):
    result = stages_validate("qc", tmp_path)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "unknown stage" not in (result.stdout + result.stderr)
    assert "qc" in result.stdout


def test_several_comma_separated_stages_are_accepted(tmp_path):
    result = stages_validate("qc,coverage,concordance", tmp_path)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "unknown stage" not in (result.stdout + result.stderr)
    for stage in ("qc", "coverage", "concordance"):
        assert stage in result.stdout


def test_every_declared_stage_name_is_accepted(tmp_path):
    """Whatever `--help` advertises must parse. Derived from the table, never re-listed."""
    make_tables = _make_tables()
    for stage in make_tables.STAGE_ORDER:
        result = stages_validate(stage, tmp_path)
        assert result.returncode == 0, "%s was rejected: %s" % (stage, result.stderr[-800:])


def test_an_unknown_stage_is_rejected_and_named(tmp_path):
    result = stages_validate("not_a_stage", tmp_path)
    assert result.returncode != 0
    assert "not_a_stage" in result.stderr, "the message must name the offending stage"


def test_one_bad_stage_rejects_the_whole_request(tmp_path):
    result = stages_validate("qc,not_a_stage", tmp_path)
    assert result.returncode != 0
    assert "not_a_stage" in result.stderr


def test_requested_stages_pull_in_their_dependencies_in_order(tmp_path):
    """`reach` reads the taxonomy and alignment-concordance masters."""
    make_tables = _make_tables()
    assert make_tables.required_stages({"reach"}) == [
        "annotation", "taxonomy", "alignment_concordance", "reach"]
    result = stages_validate("reach", tmp_path)
    assert result.returncode == 0, result.stderr[-2000:]
    line = [l for l in result.stdout.splitlines() if "Stages that would run" in l]
    assert line, result.stdout
    order = [s.strip() for s in line[0].split(":", 1)[1].split(",")]
    for earlier, later in (("annotation", "reach"), ("taxonomy", "reach"),
                           ("alignment_concordance", "reach")):
        assert order.index(earlier) < order.index(later), order


def test_validate_is_read_only_and_needs_no_generated_output(tmp_path):
    """A --validate that leaves directories behind has written to the tree it was only
    asked to inspect."""
    output = tmp_path / "out"
    result = stages_validate("qc,coverage", tmp_path)
    assert result.returncode == 0, result.stderr[-2000:]
    assert not output.exists(), "--validate created %s" % output
    assert sorted(p.name for p in tmp_path.iterdir()) == []


def test_every_main_block_is_accounted_for():
    """Every `__main__` under `code/` is either a public entry point or an orchestrated one.

    Two lists, because they carry different obligations. `ENTRY_POINTS` is what a user
    runs, and the tests above hold it to the help and no-argument contracts.
    `ORCHESTRATED` is invoked only by `make_tables.py` / `make_panels.py` with a full
    argument list, so a bare `--help` is not part of its contract.

    The check is that the union is exact. A new CLI must be classified explicitly, and a
    file left in either list after its CLI is deleted fails here rather than lingering."""
    found = set()
    for path in sorted((REPO / "code").rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        if "__main__" in path.read_text():
            found.add(str(path.relative_to(REPO)))
    classified = set(ENTRY_POINTS) | set(ORCHESTRATED)
    assert not (set(ENTRY_POINTS) & set(ORCHESTRATED)), "a file cannot be in both lists"
    assert found == classified, (
        "unclassified CLI(s): %s; classified but no longer a CLI: %s"
        % (sorted(found - classified), sorted(classified - found)))
