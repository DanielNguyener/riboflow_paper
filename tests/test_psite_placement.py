"""P-site placement: the offset is walked along the READ, not along the reference.

Every expectation here is computed by hand from the CIGAR and written in the test, so a
change in behaviour shows up as a disagreement with arithmetic rather than with a golden
value nobody can re-derive.

The rule under test (`psite_placement.place`) is: return the reference position of the
aligned read base `offset` steps from the read's 5' end, or None if the read has fewer
aligned bases than that. `get_aligned_pairs(matches_only=True)` is ascending in QUERY
order, so on the minus strand -- where the read's 5' end is the base at the HIGHEST
reference coordinate -- the walk runs from the end of the list.

Two consequences are deliberate and tested:

  * insertions and soft clips consume query but have no reference position, so they are
    not aligned bases and do not advance the walk;
  * deletions and splice gaps consume reference but have no query position, so the walk
    steps over them -- which is what keeps a junction-spanning read's P-site out of the
    intron.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
COVERAGE = REPO / "code" / "coverage"
if str(COVERAGE) not in sys.path:
    sys.path.insert(0, str(COVERAGE))

import psite_placement as pp                                       # noqa: E402

pysam = pytest.importorskip("pysam")

CONTIG = "chr1"
CONTIG_LENGTH = 100_000
START = 1000


def header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]})


def read(cigar, *, pos=START, reverse=False, mapq=255, nh=1):
    """An AlignedSegment with the given CIGAR, built in memory (no BAM needed)."""
    segment = pysam.AlignedSegment(header())
    segment.query_name = "r"
    segment.reference_id = 0
    segment.reference_start = pos
    segment.cigarstring = cigar
    n_query = sum(length for op, length in segment.cigartuples
                  if op in (0, 1, 4, 7, 8))
    segment.query_sequence = "A" * n_query
    segment.query_qualities = pysam.qualitystring_to_array("I" * n_query)
    segment.mapping_quality = mapq
    if nh is not None:
        segment.set_tag("NH", int(nh))
    segment.flag = 16 if reverse else 0
    return segment


def aligned_reference_positions(segment):
    return [ref for _query, ref in segment.get_aligned_pairs(matches_only=True)]


# ── the six alignment shapes, forward strand ─────────────────────────────────

def test_ordinary_alignment_places_offset_bases_from_the_start():
    # 30M at 1000: aligned bases are ref 1000..1029, one per read base.
    assert pp.place(read("30M"), 12) == 1012


def test_offset_zero_is_the_first_aligned_base():
    assert pp.place(read("30M"), 0) == 1000
    assert pp.place(read("5S25M"), 0) == 1000, "the walk starts at the first ALIGNED base"


def test_soft_clipping_does_not_advance_the_walk():
    # 5S25M at 1000: the 5 clipped bases have no reference position, so aligned base 12
    # is ref 1012 -- the same answer as 30M, because clipping removed read bases that
    # were never placed anywhere.
    assert pp.place(read("5S25M"), 12) == 1012


def test_insertions_do_not_advance_the_walk():
    # 10M2I18M at 1000: aligned bases are ref 1000..1009 then ref 1010..1027. The two
    # inserted bases have no reference position, so aligned base 12 is ref 1012.
    assert pp.place(read("10M2I18M"), 12) == 1012


def test_deletions_are_stepped_over():
    # 10M2D20M at 1000: aligned bases are ref 1000..1009 then ref 1012..1031 -- ref 1010
    # and 1011 are deleted and carry no read base. Aligned base 12 is therefore ref 1014,
    # NOT ref 1012.
    assert pp.place(read("10M2D20M"), 12) == 1014


def test_a_junction_is_stepped_over_instead_of_landing_in_the_intron():
    # 10M940N20M at 1000: exon 1 is ref 1000..1009, exon 2 is ref 1950..1969.
    # Aligned base 12 is the third base of exon 2 -> ref 1952.
    # Walking 12 positions along the REFERENCE would give ref 1012, which is intronic and
    # which the read does not cover at all. That is the failure this rule exists to avoid.
    segment = read("10M940N20M")
    assert pp.place(segment, 12) == 1952
    assert 1012 not in aligned_reference_positions(segment)


def test_a_read_entirely_inside_one_exon_is_unaffected_by_the_junction():
    # The same spliced read, offset 5: still inside exon 1.
    assert pp.place(read("10M940N20M"), 5) == 1005


# ── the same shapes, reverse strand ──────────────────────────────────────────

def test_reverse_ordinary_alignment_walks_from_the_high_coordinate_end():
    # 30M at 1000 reverse: the read's 5' end is ref 1029, and the walk runs downward.
    assert pp.place(read("30M", reverse=True), 12) == 1017
    assert pp.place(read("30M", reverse=True), 0) == 1029


def test_reverse_soft_clipping():
    # 5S25M at 1000 reverse: aligned bases are ref 1000..1024, 5' end is ref 1024.
    assert pp.place(read("5S25M", reverse=True), 12) == 1012
    assert pp.place(read("5S25M", reverse=True), 0) == 1024


def test_reverse_insertion():
    # 10M2I18M at 1000 reverse: 28 aligned bases, ref 1000..1009 and 1010..1027.
    # 5' end is ref 1027; 12 aligned bases down is index 27 - 12 = 15 -> ref 1015.
    assert pp.place(read("10M2I18M", reverse=True), 12) == 1015


def test_reverse_deletion():
    # 10M2D20M at 1000 reverse: 30 aligned bases, ref 1000..1009 and 1012..1031.
    # 5' end is ref 1031; index 29 - 12 = 17 -> the 8th base of the second block -> 1019.
    assert pp.place(read("10M2D20M", reverse=True), 12) == 1019


def test_reverse_junction_is_stepped_over():
    # 20M940N10M at 1000 reverse: exon 1 ref 1000..1019, exon 2 ref 1960..1969.
    # The read's 5' end is ref 1969. Twelve aligned bases down crosses the junction:
    # index 29 - 12 = 17 -> ref 1017.
    # Walking the reference downward from reference_end - 1 would give 1969 - 12 = 1957,
    # which is intronic.
    segment = read("20M940N10M", reverse=True)
    assert pp.place(segment, 12) == 1017
    assert 1957 not in aligned_reference_positions(segment)


# ── the placement is never off the read ──────────────────────────────────────

SHAPES = ["30M", "5S25M", "25M5S", "10M2I18M", "10M2D20M", "10M940N20M",
          "20M940N10M", "3S10M2D15M2S", "8M100N8M50N14M"]


@pytest.mark.parametrize("cigar", SHAPES)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("offset", [0, 1, 5, 11, 12, 19])
def test_the_psite_is_always_a_base_the_read_covers(cigar, reverse, offset):
    segment = read(cigar, reverse=reverse)
    covered = set(aligned_reference_positions(segment))
    position = pp.place(segment, offset)
    if position is None:
        assert len(covered) <= offset
    else:
        assert position in covered


@pytest.mark.parametrize("cigar", SHAPES)
@pytest.mark.parametrize("reverse", [False, True])
def test_consecutive_offsets_step_one_aligned_base_at_a_time(cigar, reverse):
    """Offsets 0..n-1 must enumerate the read's aligned bases in 5'->3' order."""
    segment = read(cigar, reverse=reverse)
    positions = aligned_reference_positions(segment)
    expected = list(reversed(positions)) if reverse else positions
    walked = [pp.place(segment, k) for k in range(len(expected))]
    assert walked == expected


# ── the undefined case ───────────────────────────────────────────────────────

def test_an_offset_past_the_last_aligned_base_is_undefined():
    assert pp.place(read("30M"), 29) == 1029
    assert pp.place(read("30M"), 30) is None
    assert pp.place(read("30M"), 999) is None


def test_a_heavily_clipped_read_can_be_undefined():
    # 25S5M has only 5 aligned bases, so an offset of 11 has nowhere to land.
    assert pp.place(read("25S5M"), 4) == 1004
    assert pp.place(read("25S5M"), 5) is None
    assert pp.place(read("25S5M"), 11) is None


def test_an_undefined_placement_is_none_on_both_strands():
    assert pp.place(read("25S5M", reverse=True), 11) is None


def test_a_negative_offset_is_rejected_rather_than_wrapping():
    """Python would happily index pairs[-1]; that would silently place the P-site at the
    3' end of the read."""
    with pytest.raises(pp.PlacementError):
        pp.place(read("30M"), -1)


# ── CIGAR helpers ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("cigar,expected", [
    ("30M", "M"), ("5S25M", "MS"), ("10M2I18M", "IM"), ("10M2D20M", "DM"),
    ("10M940N20M", "MN"), ("3S10M2D15M2S", "DMS"),
])
def test_cigar_signature(cigar, expected):
    assert pp.cigar_signature(read(cigar)) == expected


@pytest.mark.parametrize("cigar,pure", [
    ("30M", True), ("30=", True), ("15=15X", True),
    ("5S25M", False), ("10M2I18M", False), ("10M940N20M", False),
])
def test_is_pure_match(cigar, pure):
    assert pp.is_pure_match(read(cigar)) is pure


# ── the recorded policy ──────────────────────────────────────────────────────

def test_there_is_exactly_one_placement_rule():
    """The retired reference-offset rule must not come back by accident."""
    assert pp.PSITE_PLACEMENT == "cigar_aware"
    assert not hasattr(pp, "place_reference_offset")
    assert not hasattr(pp, "PLACEMENT_RULES")
    assert not hasattr(pp, "DEFAULT_PLACEMENT")


def test_place_takes_no_rule_argument():
    import inspect
    assert list(inspect.signature(pp.place).parameters) == ["read", "offset"]


def test_no_builder_offers_a_placement_choice():
    """A `--psite-placement` flag anywhere would mean the rule is still selectable."""
    for name in ("build_shared_coverage.py", "build_cohort_coverage.py"):
        text = (COVERAGE / name).read_text()
        assert "--psite-placement" not in text, name
        assert "reference_offset" not in text, name


# ── summarize_placements over a real BAM ─────────────────────────────────────

def _write_bam(path, records):
    head = {"HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(path), "wb", header=head) as out:
        for i, segment in enumerate(sorted(records, key=lambda s: s.reference_start)):
            segment.query_name = "read%d" % i
            out.write(segment)
    pysam.index(str(path))
    return path


def test_summarize_placements_counts_each_shape(tmp_path):
    records = [read("30M", pos=1000), read("30M", pos=2000),
               read("10M940N20M", pos=3000), read("10M2D20M", pos=5000),
               read("10M2I18M", pos=6000), read("5S25M", pos=7000),
               read("25S5M", pos=8000),                       # undefined at offset 12
               read("30M", pos=9000, nh=2, mapq=1)]           # a multimapper, excluded
    bam = _write_bam(tmp_path / "s.bam", records)

    result = pp.summarize_placements(bam, {30: 12, 35: 12})
    counts = result["counts"]
    assert counts["considered"] == 7, "the NH:i:2 read must be excluded"
    assert counts["pure_match"] == 2
    assert counts["spliced"] == 1
    assert counts["deletion"] == 1
    assert counts["insertion"] == 1
    assert counts["soft_clipped"] == 2                          # 5S25M and 25S5M
    assert counts["undefined"] == 1                             # 25S5M has 5 aligned bases
    assert counts["placed"] == 6
    assert result["undefined_by_cigar_signature"] == {"MS": 1}
    assert result["psite_placement"] == "cigar_aware"
    assert result["bam"] == "s.bam", "the report must not embed the full path"


def test_summarize_placements_skips_lengths_without_an_offset(tmp_path):
    bam = _write_bam(tmp_path / "s.bam", [read("30M", pos=1000), read("40M", pos=2000)])
    result = pp.summarize_placements(bam, {30: 12})
    assert result["counts"]["considered"] == 1


def test_summarize_placements_skips_secondary_and_supplementary(tmp_path):
    secondary = read("30M", pos=2000)
    secondary.flag |= 256
    supplementary = read("30M", pos=3000)
    supplementary.flag |= 2048
    bam = _write_bam(tmp_path / "s.bam",
                     [read("30M", pos=1000), secondary, supplementary])
    assert pp.summarize_placements(bam, {30: 12})["counts"]["considered"] == 1


# ── load_offsets ─────────────────────────────────────────────────────────────

QC_HEADER = "sample,read_length,in_phase1,psite_offset\n"


def test_load_offsets_selects_only_phase1_rows(tmp_path):
    qc = tmp_path / "qc.csv"
    qc.write_text(QC_HEADER + "HeLa,26,True,11\nHeLa,27,True,11\nHeLa,31,False,12\n"
                              "A549,26,True,12\n")
    assert pp.load_offsets(qc, "HeLa") == {26: 11, 27: 11}


def test_load_offsets_names_the_missing_sample(tmp_path):
    qc = tmp_path / "qc.csv"
    qc.write_text(QC_HEADER + "HeLa,26,True,11\n")
    with pytest.raises(pp.PlacementError) as excinfo:
        pp.load_offsets(qc, "Nope")
    assert "Nope" in str(excinfo.value) and "HeLa" in str(excinfo.value)


def test_load_offsets_names_the_missing_column(tmp_path):
    qc = tmp_path / "qc.csv"
    qc.write_text("sample,read_length\nHeLa,26\n")
    with pytest.raises(pp.PlacementError) as excinfo:
        pp.load_offsets(qc, "HeLa")
    assert "in_phase1" in str(excinfo.value)


def test_load_offsets_rejects_a_sample_with_no_phase1_lengths(tmp_path):
    qc = tmp_path / "qc.csv"
    qc.write_text(QC_HEADER + "HeLa,26,False,11\n")
    with pytest.raises(pp.PlacementError) as excinfo:
        pp.load_offsets(qc, "HeLa")
    assert "in_phase1" in str(excinfo.value)
