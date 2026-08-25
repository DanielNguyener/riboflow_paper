"""What a coverage file says about how it was made.

A coverage HDF5 is meant to be shareable, and a shared file is only useful if it answers,
without reference to the machine that produced it:

    which sample, which assay, which two routes
    which annotation, by content and not by path
    which P-site rule
    which code, by content
    what command

All of it is one JSON string in the root attribute `provenance`; everything here is
asserted on a file the real builder actually wrote.
"""
from __future__ import annotations

import json

import pytest

import coverage_schema
from conftest import SAMPLE, TRIM, build_config


@pytest.fixture
def provenance(coverage):
    return coverage.provenance


# ── identity ─────────────────────────────────────────────────────────────────

def test_the_schema_and_its_version_are_recorded(coverage):
    assert coverage.handle.attrs["schema"] == coverage_schema.SCHEMA
    assert coverage.schema_version == coverage_schema.SCHEMA_VERSION


def test_the_sample_assay_and_routes_are_recorded(coverage):
    assert coverage.sample == SAMPLE
    assert coverage.assay == "ribo"
    assert tuple(coverage.routes) == ("genome", "transcriptome")


def test_each_signal_name_maps_to_one_route_and_measure():
    """The four dataset names are the contract; the mapping is a constant, not an attr."""
    for name in coverage_schema.SIGNALS:
        assert coverage_schema.SIGNAL_ROUTE[name] in ("genome", "transcriptome")
        assert coverage_schema.SIGNAL_MEASURE[name] in ("psite", "footprint")


def test_the_coordinate_convention_is_stated_on_the_file(coverage, provenance):
    assert coverage.coordinate_system == "transcript_5p_to_3p"
    assert provenance["parameters"]["stop_codon_assignment"] == "utr3"
    assert provenance["parameters"]["exon_source"] == "gencode_exon_features"


def test_identity_is_available_as_one_call(coverage):
    identity = coverage.identity()
    assert identity["sample"] == SAMPLE
    assert identity["psite_placement"] == "cigar_aware"
    assert identity["paper_cds_trim"] == TRIM
    assert identity["n_transcripts"] == 3


# ── the P-site policy ────────────────────────────────────────────────────────

def test_the_psite_rule_is_recorded_in_two_places(coverage, provenance):
    import psite_placement
    assert coverage.psite_placement == psite_placement.PSITE_PLACEMENT
    assert provenance["parameters"]["psite_placement"] == "cigar_aware"


def test_the_two_assignment_policies_are_recorded_separately(provenance):
    """P-site and footprint assignment are different rules. Recording one
    `assignment_rule` would misdescribe half the file."""
    policies = provenance["assignment_policies"]
    assert policies["psite"]["rule"] == "first_exon_overlap"
    assert policies["footprint"]["rule"] == "max_exon_overlap"
    assert policies["psite"] != policies["footprint"]


def test_the_parameters_that_change_the_numbers_are_all_recorded(provenance):
    for key in ("paper_cds_trim", "genome_uniqueness", "txome_uniqueness",
                "psite_placement", "reference_name"):
        assert key in provenance["parameters"], key


def test_no_appris_rank_is_invented(provenance):
    """No formal APPRIS principal rank reaches this pipeline, and a placeholder would be
    worse than an absence."""
    assert provenance["parameters"]["appris_principal_ranks_consumed"] is False
    assert "appris_category" not in json.dumps(provenance)


# ── inputs, by content ───────────────────────────────────────────────────────

def test_every_input_is_identified_by_digest(provenance):
    for name, record in provenance["inputs"].items():
        assert "name" in record and "bytes" in record, name
        assert "sha256" in record or "index_sha256" in record, name


def test_input_paths_are_not_recorded_by_default(provenance):
    """An absolute path names the machine that built the file, the layout of a private
    project and often the operator. The digest is what identifies the data."""
    for name, record in provenance["inputs"].items():
        assert "path" not in record, "%s leaked a filesystem path" % name
        assert "/" not in record["name"], name


def test_the_recorded_command_carries_no_paths_by_default(provenance):
    assert provenance["generation"]["paths_redacted"] is True
    assert provenance["generation"]["command"].startswith("code/coverage/")


@pytest.mark.parametrize("argv,expected", [
    (["build_shared_coverage.py", "--sample", "HeLa"],
     "code/coverage/build_shared_coverage.py --sample HeLa"),
    (["b.py", "--genome-bam", "/home/me/secret/HeLa.bam"],
     "code/coverage/b.py --genome-bam HeLa.bam"),
    # the `--flag=value` form carries the path in the value half
    (["b.py", "--report=/Users/me/project/out.json"],
     "code/coverage/b.py --report=out.json"),
    (["b.py", "--output", "/scratch/user/results/coverage"],
     "code/coverage/b.py --output coverage"),
])
def test_a_recorded_command_never_names_a_directory(argv, expected):
    """Asserted on a controlled argv rather than the ambient one: the recorded command is
    whatever `sys.argv` held, which under a test runner is the runner's own command line."""
    record = coverage_schema.invocation(argv)
    assert record["command"] == expected
    assert record["paths_redacted"] is True


def test_paths_survive_when_they_are_asked_for():
    record = coverage_schema.invocation(
        ["b.py", "--report=/Users/me/out.json", "--gtf", "/ref/g.gtf"], record_paths=True)
    assert "/Users/me/out.json" in record["command"]
    assert "/ref/g.gtf" in record["command"]
    assert record["paths_redacted"] is False


def test_paths_can_be_recorded_deliberately(inputs, tmp_path):
    """Opt-in, for someone who wants them and knows what they are sharing."""
    import build_shared_coverage
    path, _ = build_shared_coverage.build(build_config(
        inputs, record_input_paths=True, output=tmp_path / "with-paths"))
    with coverage_schema.open_coverage(path) as handle:
        record = handle.provenance
    assert record["generation"]["paths_redacted"] is False
    assert "path" in record["inputs"]["gtf"]
    assert str(inputs.gtf) == record["inputs"]["gtf"]["path"]


def test_the_whole_provenance_blob_is_free_of_absolute_paths(provenance):
    """The check that matters: nothing anywhere in the record names a home directory."""
    blob = json.dumps(provenance)
    for marker in ("/Users/", "/home/", "/private/var", "/tmp/"):
        assert marker not in blob, "provenance contains %r" % marker


# ── code identity ────────────────────────────────────────────────────────────

def test_the_code_version_covers_every_module_that_defines_a_number(provenance):
    modules = provenance["code_version"]["modules"]
    assert set(modules) == set(coverage_schema.CODE_VERSION_MODULES)
    assert all(digest for digest in modules.values()), "a module was not found"


def test_the_provenance_is_one_root_attribute(coverage):
    """Everything about how the file was made is one JSON string on the root, so a plain
    h5py reader sees it without walking groups."""
    blob = coverage.handle.attrs["provenance"]
    assert json.loads(blob if isinstance(blob, str) else blob.decode())["code_version"]


def test_the_code_version_changes_when_the_code_changes(tmp_path):
    """Distinguishes 'built by different code' from 'built from different inputs', which
    input digests alone cannot do."""
    import shutil
    fake = tmp_path / "code"
    fake.mkdir()
    source = coverage_schema.Path(coverage_schema.__file__).parent
    for name in coverage_schema.CODE_VERSION_MODULES:
        shutil.copy(source / name, fake / name)
    before = coverage_schema.code_version(fake)["combined_sha256"]
    (fake / "psite_placement.py").write_text("# changed\n")
    after = coverage_schema.code_version(fake)["combined_sha256"]
    assert before != after


def test_a_missing_module_is_recorded_as_absent_not_skipped(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    version = coverage_schema.code_version(empty)
    assert set(version["modules"]) == set(coverage_schema.CODE_VERSION_MODULES)
    assert all(v is None for v in version["modules"].values())


# ── software and counts ──────────────────────────────────────────────────────

def test_the_library_versions_are_recorded(provenance):
    for library in ("python", "numpy", "pandas", "pysam", "h5py", "scipy"):
        assert library in provenance["software"], library


def test_the_counts_describe_the_file_that_was_written(coverage, provenance):
    assert provenance["counts"]["n_transcripts"] == coverage.n_transcripts
    assert provenance["counts"]["n_positions"] == coverage.n_positions


def test_the_region_summary_records_how_the_stop_codon_was_handled(provenance):
    summary = provenance["regions"]
    assert summary["stop_codon_assignment"] == "utr3"
    # These synthetic transcripts have no annotated stop codon, so nothing is relocated
    # -- and that is recorded rather than left to be inferred from a zero.
    assert summary["n_stop_relocated"] == 0
    assert summary["n_no_annotated_stop"] == 3
