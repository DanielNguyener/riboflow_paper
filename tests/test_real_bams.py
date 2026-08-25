"""Acceptance: the pipeline on REAL alignments, from an empty output directory.

Everything else in this suite runs on a synthetic cohort, which proves the logic but not
that the logic survives contact with GENCODE, a 300 MB BAM and a real APPRIS table. This
module runs the actual chain

    BAM + GTF + APPRIS  ->  transcript coordinate  ->  coverage HDF5
                        ->  concordance tables     ->  example vectors  ->  a panel

on one real sample, in a temporary directory, and asserts the invariants that only real
data can violate.

SKIPPED unless all three inputs are supplied:

    RIBOFLOW_PAPER_BAMS=/path/to/riboflow/output \\
    RIBOFLOW_PAPER_GTF=/path/to/gencode.v34.annotation.gtf.gz \\
    RIBOFLOW_PAPER_APPRIS=/path/to/appris_human_v2_transcript_lengths.tsv \\
        python -m pytest tests/test_real_bams.py -q

Optional: RIBOFLOW_PAPER_REGIONS (the actual-regions BED, used as a cross-check) and
RIBOFLOW_PAPER_SAMPLE (default HeLa).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"

BAMS = os.environ.get("RIBOFLOW_PAPER_BAMS")
GTF = os.environ.get("RIBOFLOW_PAPER_GTF")
APPRIS = os.environ.get("RIBOFLOW_PAPER_APPRIS")
REGIONS = os.environ.get("RIBOFLOW_PAPER_REGIONS")
SAMPLE = os.environ.get("RIBOFLOW_PAPER_SAMPLE", "HeLa")

GENOME_TEMPLATE = "{s}/genome/alignment_ribo/merged/{s}.post_dedup.bam"
TXOME_TEMPLATE = "{s}/transcriptome/alignment_ribo/merged/{s}.transcriptome.post_dedup.bam"

SKIP_REASON = (
    "the real-BAM acceptance run needs RIBOFLOW_PAPER_BAMS, RIBOFLOW_PAPER_GTF and "
    "RIBOFLOW_PAPER_APPRIS. Run it with:\n"
    "    RIBOFLOW_PAPER_BAMS=/path/to/riboflow/output \\\n"
    "    RIBOFLOW_PAPER_GTF=/path/to/gencode.v34.annotation.gtf.gz \\\n"
    "    RIBOFLOW_PAPER_APPRIS=/path/to/appris_human_v2_transcript_lengths.tsv \\\n"
    "        python -m pytest tests/test_real_bams.py -q")


def _missing():
    if not (BAMS and GTF and APPRIS):
        return SKIP_REASON
    for label, path in (("BAM tree", BAMS), ("GTF", GTF), ("APPRIS", APPRIS)):
        if not Path(path).exists():
            return "%s does not exist: %s" % (label, path)
    for label, template in (("genome BAM", GENOME_TEMPLATE),
                            ("txome BAM", TXOME_TEMPLATE)):
        if not (Path(BAMS) / template.format(s=SAMPLE)).exists():
            return "the %s for %s is not under %s" % (label, SAMPLE, BAMS)
    return None


pytestmark = [pytest.mark.bams,
              pytest.mark.skipif(_missing() is not None, reason=_missing() or "")]


def genome_bam():
    return Path(BAMS) / GENOME_TEMPLATE.format(s=SAMPLE)


def txome_bam():
    return Path(BAMS) / TXOME_TEMPLATE.format(s=SAMPLE)


def qc_master(kind):
    """The shipped QC masters carry the published read-length window and offsets."""
    return REPO / "data" / "ribo_seq_qc" / kind / "tables" / "readlen_window_qc.csv"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """One real coverage build, from an EMPTY output directory. No prebuilt cache and no
    prebuilt CDS table: the annotation is parsed here, from the GTF, exactly as a clean
    checkout would have to."""
    sys.path.insert(0, str(CODE / "coverage"))
    import build_shared_coverage
    import coverage_schema

    output = tmp_path_factory.mktemp("real")
    argv = ["--sample", SAMPLE,
            "--genome-bam", str(genome_bam()),
            "--transcriptome-bam", str(txome_bam()),
            "--gtf", GTF, "--appris", APPRIS,
            "--qc-genome", str(qc_master("genome")),
            "--qc-txome", str(qc_master("transcriptome")),
            "--output", str(output),
            "--report", str(output / "report.json")]
    if REGIONS:
        argv += ["--regions", REGIONS]
    assert build_shared_coverage.main(argv) == 0
    path = output / ("%s.shared_coverage.h5" % SAMPLE)
    report = json.loads((output / "report.json").read_text())
    return path, report, output, coverage_schema


# ── the build ────────────────────────────────────────────────────────────────

def test_a_real_coverage_file_is_produced_and_validates(built):
    path, _report, _output, coverage_schema = built
    assert path.exists()
    assert coverage_schema.validate_file(path) == []


def test_the_annotation_was_parsed_here_not_reused(built):
    """The point of the acceptance run: nothing came from a previous build."""
    _path, report, output, _schema = built
    assert report["n_transcripts"] > 10_000
    assert not list(output.glob("**/*.pkl")), "a cache was written into the output tree"


def test_a_prebuilt_cds_table_is_not_accepted_here():
    """Supplying one would leave the step that most needs proving unproven, so it is
    rejected rather than honoured."""
    assert os.environ.get("RIBOFLOW_PAPER_CDS_TABLE") is None, (
        "unset RIBOFLOW_PAPER_CDS_TABLE: this run builds the CDS table from the GTF on "
        "purpose")


def test_the_file_records_the_real_inputs_by_digest(built):
    path, _report, _output, coverage_schema = built
    with coverage_schema.open_coverage(path) as coverage:
        inputs = coverage.provenance["inputs"]
    assert inputs["gtf"]["name"] == Path(GTF).name
    assert len(inputs["gtf"]["sha256"]) == 64
    assert inputs["genome_bam"]["name"] == genome_bam().name
    for name, record in inputs.items():
        assert "path" not in record, "%s leaked a filesystem path" % name


def test_the_cds_bounds_lie_inside_every_transcript(built):
    """On real GENCODE this is the check that catches a wrong exon set: the builder
    asserts spliced exon length == reference length, and the stored CDS must fit."""
    path, _report, _output, coverage_schema = built
    with coverage_schema.open_coverage(path) as coverage:
        has_cds = coverage.cds_start >= 0
        assert has_cds.sum() > 10_000
        assert (coverage.cds_end[has_cds] <= coverage.transcript_len[has_cds]).all()
        assert (coverage.cds_start[has_cds] < coverage.cds_end[has_cds]).all()


def test_the_stop_codon_is_relocated_into_the_utr3_on_real_annotation(built):
    """Real transcripts have annotated stop codons; the synthetic ones do not, so this is
    where the relocation rule meets the case it was written for."""
    path, _report, _output, coverage_schema = built
    with coverage_schema.open_coverage(path) as coverage:
        summary = coverage.provenance["regions"]
    assert summary["stop_codon_assignment"] == "utr3"
    assert summary["n_stop_relocated"] > 10_000


def test_every_route_has_real_coverage(built):
    path, _report, _output, coverage_schema = built
    with coverage_schema.open_coverage(path) as coverage:
        for signal in coverage_schema.SIGNALS:
            total = int(coverage.signal(signal).sum(dtype=np.int64))
            assert total > 100_000, "%s is implausibly small: %d" % (signal, total)


def test_the_psite_rule_is_cigar_aware_and_recorded(built):
    path, _report, _output, coverage_schema = built
    with coverage_schema.open_coverage(path) as coverage:
        assert coverage.psite_placement == "cigar_aware"


def test_no_psite_lands_on_a_base_its_read_does_not_cover(built):
    """The property the CIGAR-aware rule exists for, on real spliced alignments."""
    import pysam
    sys.path.insert(0, str(CODE / "coverage"))
    sys.path.insert(0, str(CODE / "common"))
    import bam_inputs
    import psite_placement

    offsets = psite_placement.load_offsets(qc_master("genome"), SAMPLE)
    checked = spliced = 0
    handle = pysam.AlignmentFile(str(genome_bam()), "rb")
    try:
        for read in handle.fetch(until_eof=True):
            if not bam_inputs.is_unique_genome_read(read):
                continue
            offset = offsets.get(read.query_length)
            if offset is None:
                continue
            position = psite_placement.place(read, offset)
            if position is None:
                continue
            covered = {ref for _q, ref in read.get_aligned_pairs(matches_only=True)}
            assert position in covered, read.query_name
            checked += 1
            if "N" in psite_placement.cigar_signature(read):
                spliced += 1
            if checked >= 200_000:
                break
    finally:
        handle.close()
    assert checked > 1000
    assert spliced > 0, "no spliced read was seen, so the check would be vacuous"


# ── the downstream consumers ─────────────────────────────────────────────────

def run(program, *args):
    return subprocess.run([sys.executable, str(CODE / program), *args],
                          capture_output=True, text=True, cwd=str(REPO))


def test_concordance_runs_from_the_hdf5_without_touching_a_bam(built):
    path, _report, output, _schema = built
    result = run("coverage/compute_coverage_concordance.py",
                 "--coverage", str(path.parent),
                 "--output", str(output / "concordance"))
    assert result.returncode == 0, result.stderr[-3000:]
    for name in ("region_concordance_per_sample.tsv", "region_coverage_per_sample.tsv"):
        assert (output / "concordance" / name).exists()


def test_the_generic_plotter_renders_from_the_real_file(built, tmp_path):
    """GAPDH by gene ID, the way a reader would plot their own gene."""
    path, _report, _output, _schema = built
    result = run("panels/plot_transcript_coverage.py",
                 "--coverage-h5", str(path), "--gene-id", "ENSG00000111640",
                 "--expect-sample", SAMPLE,
                 "--output", str(tmp_path / "gapdh"), "--format", "png")
    assert result.returncode == 0, result.stderr[-3000:]
    assert (tmp_path / "gapdh.png").exists()
    assert not (tmp_path / "gapdh.inputs.json").exists(), \
        "the panel must not write a JSON sidecar any more"
    # Asserted on stdout because this runs the CLI as a subprocess. The in-process test
    # `test_plot_coverage.py::test_the_cli_renders_and_records_what_it_did` checks the
    # same facts against the record `render()` returns.
    assert "-> ENST" in result.stdout and "GAPDH" in result.stdout, result.stdout[-2000:]
    assert "region overlay (canonical)" in result.stdout, result.stdout[-2000:]


def test_the_plotter_refuses_the_wrong_sample(built, tmp_path):
    path, _report, _output, _schema = built
    result = run("panels/plot_transcript_coverage.py",
                 "--coverage-h5", str(path), "--gene-id", "ENSG00000111640",
                 "--expect-sample", "NOT_THIS_ONE",
                 "--output", str(tmp_path / "wrong"), "--format", "png")
    assert result.returncode != 0
    assert "was expected" in (result.stderr + result.stdout)
