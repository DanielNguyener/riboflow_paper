"""The CDS-assigned counting rule behind Figure 4.

Everything here runs on the synthetic three-transcript cohort in `conftest.py`, through the
real annotation bundle and the real BAM readers -- no mock of either. The geometry is known
exactly, so every expected count below is written down by hand rather than snapshotted.

WHY THE SYNTHETIC READS MAKE THIS TEST SHARP
--------------------------------------------
`genome_read(tid, rel)` and `txome_read(tid, rel)` place a read so that its P-SITE lands on
cds-relative `rel`. Its 5' nucleotide therefore sits at `rel - OFFSET`. The fixture includes
a read at `rel = 5`, whose 5' end is at cds-relative -7 -- inside the 5'UTR. Counting by
P-site would include it; counting by 5' end must not. One read decides the question.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from conftest import (CDS_START, GEOMETRY, OFFSET, OTHER_LEN, READ_LEN, SAMPLE,
                      TX_MINUS, TX_PLUS, TX_SHORT, cds_rel_to_genomic)

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"
for _entry in (CODE / "ribo_rna", CODE / "coverage", CODE / "common",
               CODE / "common" / "ribo_seq_qc"):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))


@pytest.fixture
def rrl():
    import ribo_rna_lib

    return ribo_rna_lib


@pytest.fixture
def bundle(inputs, rrl):
    built, _reused = rrl.load_bundle(inputs.gtf, inputs.appris, inputs.regions,
                                     inputs.root / "cache.pkl")
    return built


@pytest.fixture
def universe(bundle, inputs, rrl):
    ids, _report = rrl.build_universe(bundle, inputs.txome_bam)
    return ids


# ── the universe ─────────────────────────────────────────────────────────────

def test_the_universe_is_the_annotation_intersected_with_the_reference(bundle, inputs, rrl):
    """Not `length_filtered`, not "has a UTR" -- what both sides can actually name."""
    ids, report = rrl.build_universe(bundle, inputs.txome_bam)
    assert ids == sorted([TX_PLUS, TX_MINUS, TX_SHORT])
    assert report["n_universe"] == 3
    assert report["n_annotated_not_in_reference"] == 0
    assert report["n_reference_not_annotated"] == 0


def test_the_reference_name_is_read_for_its_transcript_id_and_nothing_else(inputs, rrl):
    """The `|CDS:..|` field in a reference name is never parsed here.

    Taking the CDS from the header on one route and from the annotation on the other would
    compare two definitions of a CDS, not two alignment routes."""
    ids = rrl.transcript_ids_in_reference(inputs.txome_bam)
    assert ids == {TX_PLUS, TX_MINUS, TX_SHORT}
    source = (CODE / "ribo_rna" / "ribo_rna_lib.py").read_text()
    assert "CDS:" not in source.split('"""', 2)[2], \
        "ribo_rna_lib must not parse a CDS span out of a reference name"


def test_a_transcript_with_no_5_prime_utr_keeps_its_cds(bundle, universe, rrl):
    """A CDS starting at transcript coordinate 0 is a transcript with no annotated 5'UTR.

    The superseded helper required all three regions, which silently dropped 691 real
    transcripts from the published universe -- every one of them for a missing UTR, not a
    missing CDS."""
    spans = rrl.transcript_cds_spans(bundle, universe)
    faked = bundle["regions"].copy()
    mask = (faked["label"] == "CDS") & (faked["transcript_id"] == TX_PLUS)
    faked.loc[mask, "start"] = 0
    faked.loc[mask, "end"] = spans[TX_PLUS][1] - spans[TX_PLUS][0]
    no_utr5 = dict(bundle, regions=faked)

    rebuilt = rrl.transcript_cds_spans(no_utr5, universe)
    assert rebuilt[TX_PLUS][0] == 0
    assert set(rebuilt) == set(spans), "no transcript is dropped for lacking a UTR"


def test_the_two_routes_describe_the_same_cds(bundle, universe, rrl):
    """Genomic CDS length == transcript-coordinate CDS length, per transcript.

    If they ever differed, a genome count and a transcriptome count of the same transcript
    would be measurements of different intervals."""
    spans = rrl.transcript_cds_spans(bundle, universe)
    genomic = rrl.genome_cds_intervals(bundle, universe).df
    per_transcript = (genomic["End"] - genomic["Start"]).groupby(
        genomic["transcript_id"]).sum()
    for tid in universe:
        assert int(per_transcript[tid]) == spans[tid][1] - spans[tid][0], tid
        assert int(per_transcript[tid]) == GEOMETRY[tid][3]


def test_a_disagreeing_genomic_cds_is_refused_not_absorbed(bundle, universe, rrl):
    broken = bundle["cds_table"].copy()
    row = broken.index[broken["transcript_id"] == TX_PLUS][0]
    broken.loc[row, "exon_len"] = int(broken.loc[row, "exon_len"]) + 3
    with pytest.raises(rrl.RiboRnaError, match="differs from their"):
        rrl.genome_cds_intervals(dict(bundle, cds_table=broken), universe)


# ── what a count is ──────────────────────────────────────────────────────────

def test_the_genome_route_counts_the_5_prime_nucleotide_and_never_a_p_site(
        bundle, universe, inputs, rrl):
    """The rel-5 read decides this: its P-site is in the CDS, its 5' end is in the UTR5.

    Expected, by hand, over the synthetic genome BAM:

      TX_PLUS   rel 20 (x2, two selected lengths) and rel 40 -> 5' ends at cds_rel 8 and 28
                plus the junction-spanning read, whose 5' base is genomic 1050 = cds_rel 50
      TX_MINUS  rel 30 and rel 50 -> cds_rel 18 and 38
      TX_SHORT  rel 12 -> cds_rel 0
      excluded  rel 5 (5' end at cds_rel -7, in the UTR5), the intergenic read, the
                NH:i:2 multimapper, the unselected read length, secondary, supplementary
    """
    cds_pr = rrl.genome_cds_intervals(bundle, universe)
    counts, assigned, ambiguous, retained = rrl.count_genome_cds(
        inputs.genome_bam, cds_pr, stranded=True, read_lengths={READ_LEN, 28})

    assert counts == {TX_PLUS: 4, TX_MINUS: 2, TX_SHORT: 1}
    assert assigned == 7
    assert ambiguous == 0
    assert retained == 9, "6 in-window plan reads + junction + intergenic + short_tx"


def test_the_transcriptome_route_counts_the_same_way(bundle, universe, inputs, rrl):
    """Same rule, same read population, on the other reference.

    The fixture's `clip5` read straddles the CDS start with its 5' end in the UTR5 and is
    excluded; `clip3` straddles the CDS end with its 5' end still in the CDS and is kept.
    A rule that asked "does the read OVERLAP the CDS" would keep both."""
    spans = rrl.transcript_cds_spans(bundle, universe)
    counts, assigned, retained = rrl.count_txome_cds(
        inputs.txome_bam, spans, read_lengths={READ_LEN, 28})

    assert counts == {TX_PLUS: 4, TX_MINUS: 2, TX_SHORT: 1}
    assert assigned == 7
    assert retained == 11, "the MAPQ-41 read and the unselected length are both out"


def test_the_5_prime_end_of_a_reverse_read_is_its_high_coordinate(rrl):
    """`reference_end - 1`, not `reference_start`.

    TX_MINUS runs 3'->5' along the genome, so a minus-strand read placed by
    `reference_start` would be counted `length - 1` bases away from where it starts."""
    genomic = cds_rel_to_genomic(TX_MINUS, 18)
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr2", "LN": 10000}]})
    read = pysam.AlignedSegment(header)
    read.reference_id = 0
    read.reference_start = genomic - READ_LEN + 1
    read.cigarstring = "%dM" % READ_LEN
    read.flag = 16
    assert rrl._five_prime(read) == genomic


def test_no_p_site_offset_reaches_this_module():
    """Figure 4 must be unable to move when the P-site rule changes.

    Figure 3 applies the offset through `psite_placement.place`; Figure 4 shares that
    module's read-length loader and nothing else from it."""
    source = (CODE / "ribo_rna" / "ribo_rna_lib.py").read_text()
    assert "psite" not in source.lower().replace("p-site", "")
    counting = (CODE / "ribo_rna" / "count_transcript_reads.py").read_text()
    assert "load_selected_lengths" in counting
    assert "load_offsets" not in counting, "the offsets must not even be loaded here"
    assert "psite_placement.place" not in counting


def test_figure_3_and_figure_4_read_the_same_window_from_the_same_loader(inputs):
    """One loader, one `in_phase1` filter -- so the two figures describe one read set."""
    import psite_placement

    offsets = psite_placement.load_offsets(inputs.qc_genome, SAMPLE)
    lengths = psite_placement.load_selected_lengths(inputs.qc_genome, SAMPLE)
    assert lengths == sorted(offsets)
    assert OTHER_LEN not in lengths, "an unselected length must not reach either figure"


# ── ambiguity ────────────────────────────────────────────────────────────────

def _join_frame(rows):
    return pd.DataFrame(rows, columns=["read_idx", "transcript_id"])


def test_a_read_in_two_genes_is_excluded_and_tallied(rrl):
    """APPRIS names one transcript per gene, so two transcripts means two genes and the
    annotation does not say which one the read came from."""
    frame = _join_frame([(0, TX_PLUS), (1, TX_PLUS), (1, TX_MINUS), (2, TX_MINUS)])
    counts, assigned, ambiguous = rrl._resolve_overlaps(frame)
    assert counts == {TX_PLUS: 1, TX_MINUS: 1}
    assert (assigned, ambiguous) == (2, 1)


def test_exclusion_does_not_depend_on_row_order(rrl):
    """`keep="first"` made the answer a property of the join's row order.

    Reversing the rows changed which transcript a straddling read was credited to, without
    changing a single input read. Both steps here are set operations, so they cannot."""
    rows = [(0, TX_PLUS), (1, TX_PLUS), (1, TX_MINUS), (2, TX_MINUS), (3, TX_SHORT)]
    forward = rrl._resolve_overlaps(_join_frame(rows))
    backward = rrl._resolve_overlaps(_join_frame(list(reversed(rows))))
    shuffled = rrl._resolve_overlaps(
        _join_frame([rows[i] for i in (2, 0, 4, 1, 3)]))
    assert forward == backward == shuffled


def test_two_cds_exons_of_one_transcript_are_one_assignment(rrl):
    """A duplicate row for the SAME transcript is not ambiguity, it is one read."""
    frame = _join_frame([(0, TX_PLUS), (0, TX_PLUS)])
    counts, assigned, ambiguous = rrl._resolve_overlaps(frame)
    assert (counts, assigned, ambiguous) == ({TX_PLUS: 1}, 1, 0)


def test_no_transcript_is_chosen_lexicographically():
    source = (CODE / "ribo_rna" / "ribo_rna_lib.py").read_text()
    assert 'keep="first"' not in source
    assert "keep='first'" not in source
    assert "idxmax" not in source


# ── the read filters ─────────────────────────────────────────────────────────

def test_a_genome_read_without_an_nh_tag_is_refused(inputs, rrl, tmp_path):
    """No MAPQ fallback. A BAM that cannot report multiplicity cannot be filtered on it,
    and quietly substituting a threshold would put multimappers in a unique-read count."""
    import bam_inputs

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 10000}]}
    path = tmp_path / "no_nh.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        segment = pysam.AlignedSegment(out.header)
        segment.query_name = "untagged"
        segment.reference_id = 0
        segment.reference_start = 1000
        segment.cigarstring = "30M"
        segment.query_sequence = "A" * 30
        segment.query_qualities = pysam.qualitystring_to_array("I" * 30)
        segment.mapping_quality = 255       # would pass any MAPQ rule
        segment.flag = 0
        out.write(segment)
    pysam.index(str(path))

    with pytest.raises(bam_inputs.InputError, match="no NH tag"):
        rrl.count_genome_cds(path, None, stranded=True, read_lengths={30})


def test_a_paired_alignment_is_a_validation_failure(rrl, tmp_path):
    """One count per alignment double-counts a paired fragment. This cohort's RNA-seq is
    single-end; a paired library needs an explicit fragment rule, not a silent guess."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 10000}]}
    path = tmp_path / "paired.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        segment = pysam.AlignedSegment(out.header)
        segment.query_name = "a_fragment"
        segment.reference_id = 0
        segment.reference_start = 1000
        segment.cigarstring = "30M"
        segment.query_sequence = "A" * 30
        segment.query_qualities = pysam.qualitystring_to_array("I" * 30)
        segment.mapping_quality = 255
        segment.set_tag("NH", 1)
        segment.flag = 1                     # paired
        out.write(segment)
    pysam.index(str(path))

    with pytest.raises(rrl.RiboRnaError, match="paired-end"):
        rrl.assert_single_end(path)


def test_the_single_end_bams_pass_validation(inputs, rrl):
    for path in (inputs.genome_bam, inputs.txome_bam):
        rrl.assert_single_end(path)


def test_each_route_applies_its_own_uniqueness_rule(inputs, rrl):
    """Genome `NH == 1`, transcriptome `MAPQ >= 42` -- never one shared cut-off.

    The fixture's genome multimapper carries MAPQ 3 and the transcriptome's low-confidence
    read carries MAPQ 41; each is excluded by the rule its own route uses."""
    import bam_inputs

    def names(path, predicate):
        handle = pysam.AlignmentFile(str(path), "rb")
        try:
            return {r.query_name for r in handle.fetch(until_eof=True) if predicate(r)}
        finally:
            handle.close()

    assert "excl_multimapper" not in names(inputs.genome_bam,
                                           bam_inputs.is_unique_genome_read)
    assert "excl_low_mapq" not in names(inputs.txome_bam, bam_inputs.is_unique_txome_read)
    assert bam_inputs.txome_min_mapq() == 42


# ── the statistics ───────────────────────────────────────────────────────────

def test_the_correlation_has_no_cpm_and_no_normalisation(rrl):
    """Both statistics are computed on RAW counts.

    Dividing a route by its own library size rescales away the very quantity the route
    comparison is about -- how much of the library each reference places in a CDS."""
    tids = ["t%d" % i for i in range(200)]
    rng = np.random.default_rng(0)
    ribo = {t: int(v) for t, v in zip(tids, rng.integers(0, 500, len(tids)))}
    rna = {t: int(v) for t, v in zip(tids, rng.integers(0, 500, len(tids)))}

    result = rrl.correlate(ribo, rna, tids)
    assert set(result) == {"n_transcripts", "spearman_rho", "pearson_log2_raw"}
    assert result["n_transcripts"] == 200

    # Spearman is rank-based, so a per-assay rescale cannot move it -- which is exactly why
    # it is the headline. Pearson on log2(raw+1) is NOT scale-free, and that is the point:
    # the genome route's lower CDS depth enters it directly.
    scaled = {t: v * 10 for t, v in ribo.items()}
    assert rrl.correlate(scaled, rna, tids)["spearman_rho"] == pytest.approx(
        result["spearman_rho"])
    assert rrl.correlate(scaled, rna, tids)["pearson_log2_raw"] != pytest.approx(
        result["pearson_log2_raw"])


def test_an_empty_route_is_nan_not_a_number(rrl):
    tids = ["a", "b", "c"]
    result = rrl.correlate({}, {"a": 1, "b": 2, "c": 3}, tids)
    assert np.isnan(result["spearman_rho"])
    assert np.isnan(result["pearson_log2_raw"])


# ── the interface that was removed ───────────────────────────────────────────

REMOVED_SYMBOLS = ("EXPECTED_UNIVERSE", "log2cpm", "ccc(", "bland_altman",
                   "pearson_raw", "ribo_denom", "rna_denom",
                   "count_genome_bam", "count_txome_bam",
                   "assert_uniqueness_not_shared")


@pytest.mark.parametrize("symbol", REMOVED_SYMBOLS)
def test_no_superseded_symbol_survives_anywhere(symbol):
    """There must never be an old implementation and a new one side by side.

    `region_lib.build_genome_exon_table` is not in this list: it survives
    because it builds RiboPy's extended-boundary CDS core for the read-length selection,
    a different region over a different transcript set from the canonical CDS Figure 4
    counts into. It stopped being Figure 4's universe helper; it did not stop existing."""
    hits = []
    for path in sorted((CODE).rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        if symbol in path.read_text():
            hits.append(str(path.relative_to(REPO)))
    assert not hits, "%r survives in %s" % (symbol, hits)


def test_the_superseded_route_program_is_gone():
    """Its universe build and its counting both moved into `count_transcript_reads.py`.
    A second copy could only disagree with the first."""
    assert not (CODE / "ribo_rna" / "compute_ribo_rna_route.py").exists()
    assert sorted(p.name for p in (CODE / "ribo_rna").glob("*.py")) == [
        "build_count_matrices.py", "count_transcript_reads.py", "ribo_rna_lib.py"]


def test_no_whole_transcript_mode_survives():
    for path in sorted((CODE / "ribo_rna").glob("*.py")):
        text = path.read_text()
        assert not re.search(r'"whole"', text), path
        assert "_staging_whole" not in text, path
        assert "SUPPLEMENTARY" not in text, path


# ── end to end, through the program ──────────────────────────────────────────

def test_the_program_writes_both_tables(inputs, tmp_path):
    """The route rows are always produced; the count table only when asked for."""
    counts_out = tmp_path / "counts.tsv"
    route_out = tmp_path / "route.tsv"
    result = subprocess.run(
        [sys.executable, str(CODE / "ribo_rna" / "count_transcript_reads.py"),
         "--sample", SAMPLE,
         "--ribo-genome-bam", str(inputs.genome_bam),
         "--ribo-txome-bam", str(inputs.txome_bam),
         "--rna-genome-bam", str(inputs.genome_bam),
         "--rna-txome-bam", str(inputs.txome_bam),
         "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
         "--regions", str(inputs.regions),
         "--annotation-cache", str(tmp_path / "cache.pkl"),
         "--qc-genome", str(inputs.qc_genome), "--qc-txome", str(inputs.qc_txome),
         "--route-output", str(route_out), "--counts-output", str(counts_out)],
        capture_output=True, text=True, cwd=str(REPO))
    assert result.returncode == 0, result.stdout + result.stderr

    counts = pd.read_csv(counts_out, sep="\t")
    assert list(counts.columns) == [
        "transcript_id", "genome_ribo_reads", "genome_rna_reads",
        "txome_ribo_reads", "txome_rna_reads"]
    assert sorted(counts["transcript_id"]) == sorted([TX_PLUS, TX_MINUS, TX_SHORT])
    assert int(counts["genome_ribo_reads"].sum()) == 7
    assert int(counts["txome_ribo_reads"].sum()) == 7

    route = pd.read_csv(route_out, sep="\t")
    assert list(route["route"]) == ["genome", "transcriptome"]
    assert set(route["region"]) == {"cds"}
    for column in ("spearman_rho", "pearson_log2_raw", "ribo_cds_frac",
                   "ribo_ambiguous_excluded", "rna_ambiguous_excluded"):
        assert column in route.columns
    assert "pearson_r" not in route.columns, "the CPM metric is gone, not renamed"


def test_the_count_table_is_optional(inputs, tmp_path):
    """Twenty-four samples contribute route rows; one contributes a worked example."""
    route_out = tmp_path / "route.tsv"
    result = subprocess.run(
        [sys.executable, str(CODE / "ribo_rna" / "count_transcript_reads.py"),
         "--sample", SAMPLE,
         "--ribo-genome-bam", str(inputs.genome_bam),
         "--ribo-txome-bam", str(inputs.txome_bam),
         "--rna-genome-bam", str(inputs.genome_bam),
         "--rna-txome-bam", str(inputs.txome_bam),
         "--gtf", str(inputs.gtf), "--appris", str(inputs.appris),
         "--regions", str(inputs.regions),
         "--annotation-cache", str(tmp_path / "cache.pkl"),
         "--qc-genome", str(inputs.qc_genome), "--qc-txome", str(inputs.qc_txome),
         "--route-output", str(route_out)],
        capture_output=True, text=True, cwd=str(REPO))
    assert result.returncode == 0, result.stdout + result.stderr
    assert route_out.exists()
    assert not list(tmp_path.glob("*counts*"))
