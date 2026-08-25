"""`qc_core`: the calculations the genome and transcriptome QC steps share.

Both routes build a different `(length, rel_pos)` table and then run IDENTICAL code on it.
That shared half lives in `qc_core`, and if it drifted the genome-versus-transcriptome
comparison would become a comparison of QC tables. The unit tests below pin the shared
calculations; the last test runs the real transcriptome step end to end on the synthetic BAM,
which is what proves the wiring as well as the arithmetic.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "ribo_seq_qc"))

import qc_core  # noqa: E402


# ── the dependency guard ─────────────────────────────────────────────────────

def test_a_present_package_passes():
    qc_core.require("os", "sys")


def test_a_missing_package_is_an_actionable_error_not_an_install():
    with pytest.raises(SystemExit) as excinfo:
        qc_core.require("a_package_that_does_not_exist")
    message = str(excinfo.value)
    assert "a_package_that_does_not_exist" in message
    assert "requirements.txt" in message, "the error must say how to fix it"


def test_no_qc_script_installs_anything_at_runtime():
    """A script that pip-installs into whatever interpreter is running can silently change
    the versions a published number was produced with."""
    for script in sorted((REPO / "code" / "ribo_seq_qc").glob("*.py")):
        text = script.read_text()
        assert '"-m", "pip"' not in text, script.name
        assert "pip install" not in text or script.name == "qc_core.py", script.name


# ── read-length selection ────────────────────────────────────────────────────

def test_the_window_expands_from_the_peak_until_85_percent_is_captured():
    counts = {26: 10, 27: 100, 28: 60, 29: 20, 30: 10}
    lengths, lo, hi, captured = qc_core.select_read_lengths(counts)
    assert lo == 27 and hi == 29
    assert lengths == [27, 28, 29]
    assert captured / sum(counts.values()) > 0.85


def test_expansion_always_takes_the_richer_neighbour():
    counts = {25: 5, 26: 40, 27: 100, 28: 5}
    lengths, lo, hi, _ = qc_core.select_read_lengths(counts)
    assert lo == 26 and hi == 27, "26 has 40 reads and 28 has 5, so 26 is taken first"
    assert lengths == [26, 27]


def test_only_lengths_21_to_40_are_considered():
    """TE_model slices `data_tmp.iloc[6:26]` off a distribution starting at 15 nt, i.e.
    21..40. Counts outside that window are invisible: they are neither reachable by
    expansion nor part of the 85 % denominator."""
    assert (qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN) == (21, 40)

    counts = {15: 10_000, 30: 100, 45: 10_000}
    lengths, lo, hi, captured = qc_core.select_read_lengths(counts)
    assert (lo, hi) == (30, 30), "the 15 nt and 45 nt mass is outside the search range"
    assert captured == 100, "and outside the denominator too"

    edges = {21: 50, 40: 50}
    lengths, lo, hi, captured = qc_core.select_read_lengths(edges)
    assert (lo, hi) == (21, 40) and captured == 100
    assert lengths == list(range(21, 41)), "the interval is inclusive of both bounds"


def test_a_tie_for_the_mode_resolves_to_the_shortest_tied_length():
    """The original reads `.values[0]` off a frame ordered by ascending `read_length`, so
    the FIRST maximum wins. A dict-order `max()` would be at the mercy of insertion order."""
    counts = {26: 100, 30: 100, 34: 100}
    _lengths, lo, hi, _captured = qc_core.select_read_lengths(counts)
    assert lo <= 26 <= hi
    # built in a different insertion order, the answer must not move
    shuffled = {34: 100, 26: 100, 30: 100}
    assert qc_core.select_read_lengths(shuffled) == qc_core.select_read_lengths(counts)


def test_a_tie_between_neighbours_takes_the_longer_length():
    """`counts[mmax + 1] >= counts[mmin - 1]` -- the `>=` sends ties upward.

    Sized so exactly ONE expansion happens: 80 of 100 is under the threshold, 90 is over,
    so the single step taken is the one the tie decides.
    """
    counts = {29: 10, 30: 80, 31: 10}
    _lengths, lo, hi, captured = qc_core.select_read_lengths(counts)
    assert (lo, hi) == (30, 31), "29 and 31 are tied at 10, so the longer one is taken"
    assert captured == 90


def test_capturing_exactly_85_percent_expands_once_more():
    """The loop is `while value <= pct_85`, not `<`. A window sitting exactly on the
    threshold is not finished."""
    counts = {29: 85, 30: 15}                       # peak is exactly 85 % of 100
    _lengths, lo, hi, captured = qc_core.select_read_lengths(counts)
    assert (lo, hi) == (29, 30), "85 <= 85.0, so it expands"
    assert captured == 100

    just_over = {29: 86, 30: 14}
    _lengths, lo, hi, captured = qc_core.select_read_lengths(just_over)
    assert (lo, hi) == (29, 29), "86 > 85.0, so it stops"
    assert captured == 86


def test_the_returned_interval_is_contiguous_and_inclusive():
    """Every length between lo and hi is present, including ones with zero reads: the
    interval is a range, not the set of lengths that happened to be observed."""
    counts = {25: 100, 26: 0, 27: 0, 28: 60, 29: 40}
    lengths, lo, hi, _captured = qc_core.select_read_lengths(counts)
    assert lengths == list(range(lo, hi + 1))
    assert 26 in lengths and 27 in lengths, "empty bins inside the interval are kept"


def test_the_peak_alone_can_already_satisfy_the_target():
    _lengths, lo, hi, _ = qc_core.select_read_lengths({30: 1000, 31: 1})
    assert (lo, hi) == (30, 30)


def test_an_empty_cds_histogram_is_an_error_not_a_silent_window():
    with pytest.raises(ValueError, match="no CDS-assigned reads"):
        qc_core.select_read_lengths({})
    with pytest.raises(ValueError):
        qc_core.select_read_lengths({15: 10_000})   # all outside 21-40


# ── the CDS-assigned histogram the selection runs on ─────────────────────────

def test_non_cds_reads_are_excluded_from_the_transcriptome_histogram(tmp_path):
    """Selection is on `get_length_dist("CDS")`, so a read has to land in RiboPy's
    extended-boundary CDS -- narrower than the annotated CDS by `right_span + 1` at the
    5' end and `left_span - 3` at the 3' end."""
    import conftest
    import pysam
    import region_lib as rl

    ref = "ENSTX|UTR5:1-100|CDS:101-400|UTR3:401-600|"
    start_site, stop_site = 100, 400                  # 0-based, half-open
    core_lo = start_site + rl.DEFAULT_RIGHT_SPAN + 1  # 111
    core_hi = stop_site - rl.DEFAULT_LEFT_SPAN        # 365

    placements = [
        (0, "utr5"), (core_lo - 1, "just before the core"), (core_lo, "first core base"),
        (250, "mid CDS"), (core_hi - 1, "last core base"), (core_hi, "just after"),
        (500, "utr3"),
    ]
    records = [dict(ref=ref, pos=pos, cigar="28M", mapq=42, name=label.replace(" ", "_"))
               for pos, label in placements]
    path = conftest._write_bam(tmp_path / "cds.bam", [(ref, 600)], records)

    with pysam.AlignmentFile(str(path), "rb") as bam:
        hist = qc_core.cds_length_hist_transcriptome(bam, 21, 40)

    assert hist == {28: 3}, (
        "only the three 5' ends inside [%d, %d) count; got %s" % (core_lo, core_hi, hist))


def test_the_transcriptome_histogram_respects_the_length_gate(tmp_path):
    import conftest
    import pysam

    ref = "ENSTX|UTR5:1-100|CDS:101-400|UTR3:401-600|"
    records = [dict(ref=ref, pos=200, cigar="%dM" % n, mapq=42, name="len%d" % n)
               for n in (20, 21, 30, 40, 41)]
    path = conftest._write_bam(tmp_path / "lens.bam", [(ref, 600)], records)
    with pysam.AlignmentFile(str(path), "rb") as bam:
        hist = qc_core.cds_length_hist_transcriptome(bam, 21, 40)
    assert sorted(hist) == [21, 30, 40], "20 and 41 nt are outside the gate"


def test_a_genome_read_is_counted_once_under_overlapping_annotation():
    """Two APPRIS isoforms can put the same genomic base in both their CDS cores. RiboPy
    counts one alignment once; a per-transcript join would count it twice and inflate
    exactly the bins under dense annotation."""
    pr = pytest.importorskip("pyranges")
    reads = pd.DataFrame({"Chromosome": ["chr1"], "pos5": [1000], "Strand": ["+"],
                          "length": [29]})
    overlapping = pr.PyRanges(pd.DataFrame({
        "Chromosome": ["chr1", "chr1", "chr1"],
        "Start": [900, 950, 990], "End": [1100, 1050, 1010],
        "Strand": ["+", "+", "+"]}))
    assert qc_core.cds_length_hist_genome(reads, overlapping, 21, 40) == {29: 1}


def test_genome_reads_outside_the_cds_core_and_off_strand_are_excluded():
    pr = pytest.importorskip("pyranges")
    reads = pd.DataFrame({
        "Chromosome": ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2"],
        "pos5":       [900,    1000,   1099,   899,    1100,   1000],
        "Strand":     ["+",    "+",    "+",    "+",    "+",    "+"],
        "length":     [29,     29,     29,     29,     29,     29]})
    core = pr.PyRanges(pd.DataFrame({"Chromosome": ["chr1"], "Start": [900],
                                     "End": [1100], "Strand": ["+"]}))
    assert qc_core.cds_length_hist_genome(reads, core, 21, 40) == {29: 3}, \
        "900, 1000 and 1099 are inside [900, 1100); 899 is before it, 1100 is past the " \
        "half-open end, and chr2 is elsewhere"

    wrong_strand = reads.assign(Strand="-")
    assert qc_core.cds_length_hist_genome(wrong_strand, core, 21, 40) == {}, \
        "the CDS core is on +, so nothing on - can be inside it"


def test_both_routes_select_identically_from_equivalent_cds_histograms(tmp_path):
    """The routes differ only in how a 5' end is tested for CDS membership. Given
    histograms that agree, the interval must agree -- that is what makes the
    genome-versus-transcriptome comparison a comparison of alignments, not of QC code."""
    import conftest
    import pysam
    import region_lib as rl

    ref = "ENSTX|UTR5:1-100|CDS:101-400|UTR3:401-600|"
    core_lo = 100 + rl.DEFAULT_RIGHT_SPAN + 1
    plan = {26: 20, 27: 40, 28: 100, 29: 50, 30: 10}

    records, genome_rows = [], []
    for length, n in plan.items():
        for i in range(n):
            records.append(dict(ref=ref, pos=core_lo + (i % 40), cigar="%dM" % length,
                                mapq=42, name="t%d_%d" % (length, i)))
            genome_rows.append({"Chromosome": "chr1", "pos5": 5000 + (i % 40),
                                "Strand": "+", "length": length})
    path = conftest._write_bam(tmp_path / "eq.bam", [(ref, 600)], records)
    with pysam.AlignmentFile(str(path), "rb") as bam:
        txome_hist = qc_core.cds_length_hist_transcriptome(bam, 21, 40)

    pr = pytest.importorskip("pyranges")
    core = pr.PyRanges(pd.DataFrame({"Chromosome": ["chr1"], "Start": [5000],
                                     "End": [5100], "Strand": ["+"]}))
    genome_hist = qc_core.cds_length_hist_genome(pd.DataFrame(genome_rows), core, 21, 40)

    assert txome_hist == genome_hist == plan
    assert qc_core.select_read_lengths(txome_hist) == \
        qc_core.select_read_lengths(genome_hist)


# ── the metagene ─────────────────────────────────────────────────────────────

def frame(pairs):
    return pd.DataFrame([{"length": l, "rel_pos": float(p)} for l, p in pairs])


def test_the_metagene_counts_only_inside_the_window():
    reads = frame([(30, -12), (30, -12), (30, 0), (30, 99), (30, -999)])
    counts = qc_core.metagene_counts(reads, 50, 30)
    assert counts[30][-12] == 2
    assert counts[30][0] == 1
    assert 99 not in counts[30] and -999 not in counts[30]


def test_reads_without_a_coding_reference_are_excluded():
    reads = pd.DataFrame([{"length": 30, "rel_pos": -12.0},
                          {"length": 30, "rel_pos": float("nan")}])
    assert qc_core.metagene_counts(reads, 50, 30)[30][-12] == 1


# ── the output table ─────────────────────────────────────────────────────────

def test_every_observed_length_gets_a_row_selected_or_not(tmp_path):
    length_counts = {28: 5, 30: 95}
    phase2 = {30: {"psite_offset": 12, "frame0_pct": 80.0, "frame1_pct": 10.0,
                   "frame2_pct": 10.0, "periodic": True, "n_first10": 10,
                   "n_first10_frame0": 8, "n_first10_frame1": 1, "n_first10_frame2": 1}}
    qc = qc_core.window_qc_table(length_counts, 100, phase2, str(tmp_path), "SYN", 50.0)

    assert list(qc["read_length"]) == [28, 30]
    assert list(qc["in_phase1"]) == [False, True]
    unselected = qc[qc["read_length"] == 28].iloc[0]
    assert pd.isna(unselected["psite_offset"]), "an unselected length has no offset"
    assert unselected["periodic"] is False or unselected["periodic"] == False  # noqa: E712


def test_the_published_column_order_is_preserved(tmp_path):
    """Downstream loaders and the shipped tables index these by name, but a reordering
    would still show up as a diff against every published CSV."""
    qc = qc_core.window_qc_table({30: 10}, 10, {}, str(tmp_path), "SYN", 50.0)
    assert list(qc.columns) == [
        "read_length", "n_reads", "pct_reads", "in_phase1",
        "psite_offset", "frame0_pct", "frame1_pct", "frame2_pct", "periodic",
        "n_first10", "n_first10_frame0", "n_first10_frame1", "n_first10_frame2"]


# ── end to end ───────────────────────────────────────────────────────────────

def test_the_transcriptome_step_runs_end_to_end(tmp_path):
    """The real script, on the synthetic BAM: proves the wiring into `qc_core`, not just
    the arithmetic."""
    import conftest

    bam = tmp_path / "SYN.transcriptome.bam"
    conftest.build_txome_bam(bam)
    out = tmp_path / "out"
    env = dict(os.environ,
               PYTHONPATH=os.pathsep.join([str(REPO / "code" / "common"),
                                           str(REPO / "code" / "common" / "ribo_seq_qc"),
                                           str(REPO / "code" / "ribo_seq_qc")]),
               RIBOFLOW_PAPER_QC_TX_OUT=str(out), MPLBACKEND="Agg")
    script = str(REPO / "code" / "ribo_seq_qc" / "01t_readlen_psite_qc_transcriptome.py")
    result = subprocess.run(
        [sys.executable, script, "--sample", "SYN", "--bam", str(bam)],
        capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stderr[-3000:]

    staging = out / "tables" / "_staging"
    qc = pd.read_csv(staging / "SYN_readlen_window_qc.csv")
    assert list(qc.columns)[:4] == ["read_length", "n_reads", "pct_reads", "in_phase1"]
    assert qc["in_phase1"].any(), "no read length was selected"
    # One QC table per sample, and only one: the selected lengths and their offsets are
    # columns of this file, so a second `psite_shifts` table held nothing new.
    assert sorted(p.name for p in staging.glob("SYN_*.csv")) == \
        ["SYN_readlen_window_qc.csv"]

    # The metagenes are diagnostics, so the default run writes none of them.
    assert not (out / "plots").exists(), \
        "the default run should produce tables only, no diagnostic plots"

    # `--plots` still draws them, into a second output root so the tables above are not
    # overwritten (the step refuses nothing, but re-running would just redo the work).
    plotted = tmp_path / "out_plots"
    result = subprocess.run(
        [sys.executable, script, "--sample", "SYN", "--bam", str(bam), "--plots"],
        capture_output=True, text=True,
        env=dict(env, RIBOFLOW_PAPER_QC_TX_OUT=str(plotted)))
    assert result.returncode == 0, result.stderr[-3000:]
    for plot in ("SYN_preshift.pdf", "SYN_postshift.pdf"):
        assert (plotted / "plots" / "metagene" / plot).exists(), plot
