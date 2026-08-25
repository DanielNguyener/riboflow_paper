"""The cohort driver: manifest validation, opt-in scope, reuse, and the checksum manifest.

The cohort product is roughly a gigabyte and hours of I/O. Three properties matter:

  * a bad manifest is reported in full, before any compute starts;
  * a whole-cohort run is asked for, never inferred;
  * a built file can be recognised later without transferring it.
"""
from __future__ import annotations

import csv

import pytest

import build_cohort_coverage as cohort
from conftest import SAMPLE

SCHEMA = cohort.SCHEMA_VERSION


def write_manifest(path, rows):
    columns = ["schema_version", "sample_id"] + list(cohort.BAM_COLUMNS)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    return path


def good_row(inputs, sample=SAMPLE):
    """A manifest row pointing at the synthetic BAMs. The RNA columns reuse the ribo
    files: this driver only validates that they exist."""
    return {
        "schema_version": SCHEMA, "sample_id": sample,
        "ribo_genome_bam": str(inputs.genome_bam),
        "ribo_genome_bai": str(inputs.genome_bam) + ".bai",
        "ribo_txome_bam": str(inputs.txome_bam),
        "ribo_txome_bai": str(inputs.txome_bam) + ".bai",
        "rna_genome_bam": str(inputs.genome_bam),
        "rna_genome_bai": str(inputs.genome_bam) + ".bai",
        "rna_txome_bam": str(inputs.txome_bam),
        "rna_txome_bai": str(inputs.txome_bam) + ".bai",
    }


# ── manifest validation ──────────────────────────────────────────────────────

def test_a_valid_manifest_is_read(inputs, tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    rows = cohort.read_manifest(path)
    assert [r["sample_id"] for r in rows] == [SAMPLE]


def test_an_empty_manifest_is_rejected(tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [])
    with pytest.raises(cohort.CohortError) as excinfo:
        cohort.read_manifest(path)
    assert "no rows" in str(excinfo.value)


def test_a_wrong_schema_version_is_rejected(inputs, tmp_path):
    row = good_row(inputs)
    row["schema_version"] = "something/else/9"
    path = write_manifest(tmp_path / "m.tsv", [row])
    with pytest.raises(cohort.CohortError) as excinfo:
        cohort.read_manifest(path)
    assert "schema_version" in str(excinfo.value)


def test_a_duplicate_sample_id_is_rejected(inputs, tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs), good_row(inputs)])
    with pytest.raises(cohort.CohortError) as excinfo:
        cohort.read_manifest(path)
    assert "duplicate" in str(excinfo.value)


def test_an_empty_sample_id_is_rejected(inputs, tmp_path):
    row = good_row(inputs)
    row["sample_id"] = ""
    path = write_manifest(tmp_path / "m.tsv", [row])
    with pytest.raises(cohort.CohortError) as excinfo:
        cohort.read_manifest(path)
    assert "empty sample_id" in str(excinfo.value)


def test_every_problem_is_reported_at_once(inputs, tmp_path):
    """Reporting one problem per run turns a five-minute fix into five runs."""
    bad = good_row(inputs, "A")
    bad["schema_version"] = "wrong"
    worse = good_row(inputs, "")
    worse["schema_version"] = "wrong"
    path = write_manifest(tmp_path / "m.tsv", [bad, worse])
    with pytest.raises(cohort.CohortError) as excinfo:
        cohort.read_manifest(path)
    message = str(excinfo.value)
    assert message.count("schema_version") == 2
    assert "empty sample_id" in message


# ── file validation ──────────────────────────────────────────────────────────

def test_validation_passes_when_every_file_is_present(inputs, tmp_path):
    rows = cohort.read_manifest(write_manifest(tmp_path / "m.tsv", [good_row(inputs)]))
    assert cohort.validate(rows, None) == []


def test_a_missing_bam_is_named(inputs, tmp_path):
    row = good_row(inputs)
    row["ribo_txome_bam"] = str(tmp_path / "absent.bam")
    rows = cohort.read_manifest(write_manifest(tmp_path / "m.tsv", [row]))
    problems = cohort.validate(rows, None, cohort.COVERAGE_COLUMNS)
    assert len(problems) == 1
    assert "ribo_txome_bam" in problems[0] and "absent.bam" in problems[0]


def test_an_empty_bam_is_caught_as_well_as_a_missing_one(inputs, tmp_path):
    """A zero-byte file exists. It is still not a BAM."""
    hollow = tmp_path / "hollow.bam"
    hollow.write_bytes(b"")
    row = good_row(inputs)
    row["ribo_genome_bam"] = str(hollow)
    rows = cohort.read_manifest(write_manifest(tmp_path / "m.tsv", [row]))
    problems = cohort.validate(rows, None, cohort.COVERAGE_COLUMNS)
    assert any("is empty" in p for p in problems)


def test_a_relative_manifest_path_resolves_against_the_bam_root(inputs, tmp_path):
    row = good_row(inputs)
    row["ribo_genome_bam"] = inputs.genome_bam.name
    rows = cohort.read_manifest(write_manifest(tmp_path / "m.tsv", [row]))
    resolved = cohort.resolve(rows[0], "ribo_genome_bam", inputs.root)
    assert resolved == inputs.root / inputs.genome_bam.name
    assert resolved.exists()


def test_an_absolute_manifest_path_is_left_alone(inputs, tmp_path):
    rows = cohort.read_manifest(write_manifest(tmp_path / "m.tsv", [good_row(inputs)]))
    assert cohort.resolve(rows[0], "ribo_genome_bam", "/elsewhere") == inputs.genome_bam


# ── scope is opt-in ──────────────────────────────────────────────────────────

def test_selecting_nothing_is_an_error_not_a_whole_cohort_run(inputs, tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    with pytest.raises(SystemExit) as excinfo:
        cohort.main(["--manifest", str(path), "--gtf", str(inputs.gtf),
                     "--appris", str(inputs.appris),
                     "--qc-genome", str(inputs.qc_genome),
                     "--qc-txome", str(inputs.qc_txome),
                     "--output", str(tmp_path / "out")])
    assert "no implicit whole-cohort run" in str(excinfo.value)


def test_samples_and_all_together_are_refused(inputs, tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    with pytest.raises(SystemExit) as excinfo:
        cohort.main(["--manifest", str(path), "--all", "--samples", SAMPLE,
                     "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                     "--qc-genome", str(inputs.qc_genome),
                     "--qc-txome", str(inputs.qc_txome),
                     "--output", str(tmp_path / "out")])
    assert "not both" in str(excinfo.value)


def test_an_unknown_sample_is_refused(inputs, tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    with pytest.raises(SystemExit) as excinfo:
        cohort.main(["--manifest", str(path), "--samples", "NOPE",
                     "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                     "--qc-genome", str(inputs.qc_genome),
                     "--qc-txome", str(inputs.qc_txome),
                     "--output", str(tmp_path / "out")])
    assert "not in the manifest" in str(excinfo.value)


@pytest.mark.parametrize("omit", ["gtf", "appris", "qc-genome", "qc-txome"])
def test_each_required_build_input_is_demanded_by_name(inputs, tmp_path, omit):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    argv = ["--manifest", str(path), "--all", "--output", str(tmp_path / "out")]
    for flag, value in (("gtf", inputs.gtf), ("appris", inputs.appris),
                        ("qc-genome", inputs.qc_genome), ("qc-txome", inputs.qc_txome)):
        if flag != omit:
            argv += ["--" + flag, str(value)]
    with pytest.raises(SystemExit) as excinfo:
        cohort.main(argv)
    assert "--%s is required" % omit in str(excinfo.value)


def test_validate_reports_and_exits_without_building(inputs, tmp_path, capsys):
    path = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    out = tmp_path / "out"
    assert cohort.main(["--manifest", str(path), "--bams", str(inputs.root),
                        "--validate", "--output", str(out)]) == 0
    assert "VALID" in capsys.readouterr().out
    assert not out.exists(), "--validate must not build anything"


# ── a real cohort run ────────────────────────────────────────────────────────

@pytest.fixture
def built_cohort(inputs, tmp_path):
    manifest = write_manifest(tmp_path / "m.tsv", [good_row(inputs)])
    out = tmp_path / "cohort"
    code = cohort.main(["--manifest", str(manifest), "--samples", SAMPLE,
                        "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                        "--regions", str(inputs.regions),
                        "--qc-genome", str(inputs.qc_genome),
                        "--qc-txome", str(inputs.qc_txome),
                        "--output", str(out), "--gzip-level", "1"])
    return code, out, manifest, inputs


def test_the_cohort_driver_builds_and_reports_success(built_cohort):
    code, out, _manifest, _inputs = built_cohort
    assert code == 0
    assert (out / ("%s.shared_coverage.h5" % SAMPLE)).exists()


def test_a_checksum_manifest_is_written(built_cohort):
    _code, out, _manifest, _inputs = built_cohort
    checksums = out / "coverage_checksums.tsv"
    assert checksums.exists()
    rows = list(csv.DictReader(open(checksums), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["sample_id"] == SAMPLE
    assert len(row["sha256"]) == 64
    assert len(row["provenance_sha256"]) == 64
    assert row["schema_version"].endswith("/3")
    assert int(row["n_transcripts"]) == 3


def test_the_recorded_checksum_matches_the_file_on_disk(built_cohort):
    """A checksum manifest that does not describe the files beside it is worse than none."""
    _code, out, _manifest, _inputs = built_cohort
    row = next(csv.DictReader(open(out / "coverage_checksums.tsv"), delimiter="\t"))
    path = out / row["filename"]
    assert cohort.sha256_file(path) == row["sha256"]
    assert path.stat().st_size == int(row["bytes"])


def test_skip_existing_does_not_rebuild(built_cohort):
    """Reuse is by explicit request, and it must not silently rewrite the file it claims
    to be reusing."""
    _code, out, manifest, inputs = built_cohort
    path = out / ("%s.shared_coverage.h5" % SAMPLE)
    before = path.read_bytes()
    code = cohort.main(["--manifest", str(manifest), "--samples", SAMPLE,
                        "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                        "--regions", str(inputs.regions),
                        "--qc-genome", str(inputs.qc_genome),
                        "--qc-txome", str(inputs.qc_txome),
                        "--output", str(out), "--skip-existing"])
    assert code == 0
    assert path.read_bytes() == before


def test_without_skip_existing_the_file_is_rebuilt(built_cohort):
    _code, out, manifest, inputs = built_cohort
    path = out / ("%s.shared_coverage.h5" % SAMPLE)
    before = path.stat().st_mtime_ns
    code = cohort.main(["--manifest", str(manifest), "--samples", SAMPLE,
                        "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                        "--regions", str(inputs.regions),
                        "--qc-genome", str(inputs.qc_genome),
                        "--qc-txome", str(inputs.qc_txome),
                        "--output", str(out), "--gzip-level", "1"])
    assert code == 0
    assert path.stat().st_mtime_ns != before


def test_a_rebuild_from_the_same_inputs_has_the_same_provenance_digest(built_cohort):
    """`created_utc` differs between runs, so the files are not byte-identical -- but the
    inputs, parameters and code are, and that is what the provenance digest records."""
    _code, out, manifest, inputs = built_cohort
    first = next(csv.DictReader(open(out / "coverage_checksums.tsv"), delimiter="\t"))
    cohort.main(["--manifest", str(manifest), "--samples", SAMPLE,
                 "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                 "--regions", str(inputs.regions),
                 "--qc-genome", str(inputs.qc_genome),
                 "--qc-txome", str(inputs.qc_txome),
                 "--output", str(out), "--gzip-level", "1"])
    second = next(csv.DictReader(open(out / "coverage_checksums.tsv"), delimiter="\t"))
    assert first["provenance_sha256"] == second["provenance_sha256"]


def test_a_changed_annotation_changes_the_provenance_digest(built_cohort, tmp_path):
    """The same path with different content is exactly what a path comparison cannot see
    and a digest must."""
    from conftest import VARIANT_EXONS, VARIANT_GEOMETRY, build_synthetic_gtf

    _code, out, manifest, inputs = built_cohort
    first = next(csv.DictReader(open(out / "coverage_checksums.tsv"), delimiter="\t"))
    build_synthetic_gtf(inputs.gtf, exons=VARIANT_EXONS, geometry=VARIANT_GEOMETRY)
    cohort.main(["--manifest", str(manifest), "--samples", SAMPLE,
                 "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                 "--regions", str(inputs.regions),
                 "--qc-genome", str(inputs.qc_genome),
                 "--qc-txome", str(inputs.qc_txome),
                 "--output", str(out), "--gzip-level", "1"])
    second = next(csv.DictReader(open(out / "coverage_checksums.tsv"), delimiter="\t"))
    assert first["provenance_sha256"] != second["provenance_sha256"]


def test_a_failed_sample_is_reported_and_the_exit_status_is_nonzero(inputs, tmp_path):
    """One bad BAM must not discard twenty good hours, but it must not pass silently."""
    broken = tmp_path / "broken.bam"
    broken.write_bytes(b"not a bam at all")
    row = good_row(inputs)
    row["ribo_genome_bam"] = str(broken)
    row["ribo_genome_bai"] = str(broken)
    manifest = write_manifest(tmp_path / "m.tsv", [row])
    code = cohort.main(["--manifest", str(manifest), "--samples", SAMPLE,
                        "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
                        "--qc-genome", str(inputs.qc_genome),
                        "--qc-txome", str(inputs.qc_txome),
                        "--output", str(tmp_path / "out")])
    assert code == 1
