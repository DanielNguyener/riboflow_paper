"""The shared transcript coordinate: exon map, ordering, reconciliation, round-trip.

Synthetic throughout -- no GTF download, no BAM. The geometry exercises every
rule the coordinate depends on: both strands, a splice junction, a single-exon transcript,
and a transcript whose exons appear in the GTF in the wrong order.

The map is built from complete `exon` features rather than the region-split CDS + UTR
caches so that the spliced length equals the reference length for every transcript.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
COVERAGE_DIR = REPO / "code" / "coverage"


@pytest.fixture(scope="module")
def tc():
    return __import__("transcript_coords")


# ── synthetic geometry ───────────────────────────────────────────────────────
# TP  chr1 '+'  two exons, 12 + 8 = 20 nt   -- splice junction at tx 12
# TM  chr2 '-'  two exons, 10 + 5 = 15 nt   -- minus strand, so 5'->3' is DESCENDING
# TS  chr3 '+'  one exon, 6 nt              -- single exon
GEOMETRY = {
    "ENSTP.1": ("chr1", "+", [(100, 112), (200, 208)], 20),
    "ENSTM.2": ("chr2", "-", [(500, 510), (300, 305)], 15),
    "ENSTS.3": ("chr3", "+", [(50, 56)], 6),
}


def _features(shuffle=False):
    """The shape `parse_gtf_features` returns: {tid: {"exon": [...], "CDS": [...]}}."""
    out = {}
    for tid, (chrom, strand, exons, _length) in GEOMETRY.items():
        rows = [(chrom, start, end, strand) for start, end in exons]
        if shuffle:
            rows = list(reversed(rows))
        out[tid] = {"exon": rows, "CDS": rows[:1]}
    return out


def _headers():
    return {
        tid: {
            "transcript_id": tid, "gene_id": "ENSG%s" % tid[4],
            "transcript_name": "%s-201" % tid, "gene_name": "GENE%s" % tid[4],
            "transcript_len": length,
        }
        for tid, (_c, _s, _e, length) in GEOMETRY.items()
    }


@pytest.fixture
def coords(tc):
    return tc.build_transcript_coords(_features(), _headers())


# ── ordering and geometry ────────────────────────────────────────────────────

def test_transcripts_are_in_sorted_id_order(coords):
    """Storage order matters: the pooled-Pearson reconstruction relies on the
    covered subset in storage order equalling sorted(covered)."""
    ids = coords["transcripts"]["transcript_id"].tolist()
    assert ids == sorted(ids)


def test_exons_are_ordered_five_to_three_on_the_plus_strand(coords):
    exons = coords["exons"]
    plus = exons[exons["transcript_index"] == 1]          # ENSTP.1 sorts second
    assert plus["g_start"].tolist() == [100, 200]
    assert plus["tx_start"].tolist() == [0, 12]
    assert plus["tx_end"].tolist() == [12, 20]


def test_exons_are_ordered_five_to_three_on_the_minus_strand(coords):
    """On '-' the 5' exon has the HIGHER genomic coordinate."""
    exons = coords["exons"]
    minus = exons[exons["transcript_index"] == 0]         # ENSTM.2 sorts first
    assert minus["g_start"].tolist() == [500, 300]
    assert minus["tx_start"].tolist() == [0, 10]


def test_gtf_row_order_does_not_matter(tc):
    """Exons are ordered by coordinate and strand, not by their order in the file."""
    a = tc.build_transcript_coords(_features(shuffle=False), _headers())
    b = tc.build_transcript_coords(_features(shuffle=True), _headers())
    assert a["exons"].equals(b["exons"])


def test_coverage_offsets_are_the_running_sum(coords):
    transcripts = coords["transcripts"]
    expected = np.concatenate([[0], np.cumsum(transcripts["transcript_len"])[:-1]])
    assert transcripts["coverage_offset"].tolist() == expected.tolist()
    assert coords["n_positions"] == int(transcripts["transcript_len"].sum())


# ── the reconciliation gate ──────────────────────────────────────────────────

def test_reconciliation_failure_is_fatal_and_names_the_transcript(tc):
    headers = _headers()
    headers["ENSTP.1"]["transcript_len"] = 21             # one nt too long
    with pytest.raises(tc.CoordinateError) as excinfo:
        tc.build_transcript_coords(_features(), headers)
    message = str(excinfo.value)
    assert "ENSTP.1" in message
    assert "spliced" in message


def test_a_transcript_with_no_exon_feature_is_fatal(tc):
    features = _features()
    del features["ENSTS.3"]
    with pytest.raises(tc.CoordinateError) as excinfo:
        tc.build_transcript_coords(features, _headers())
    assert "ENSTS.3" in str(excinfo.value)


def test_exons_on_two_chromosomes_are_rejected(tc):
    features = _features()
    features["ENSTP.1"]["exon"][1] = ("chrX", 200, 208, "+")
    with pytest.raises(tc.CoordinateError) as excinfo:
        tc.build_transcript_coords(features, _headers())
    assert "chromosomes" in str(excinfo.value)


def test_overlapping_exons_are_rejected(tc):
    features = _features()
    features["ENSTP.1"]["exon"] = [("chr1", 100, 112, "+"), ("chr1", 105, 113, "+")]
    headers = _headers()
    headers["ENSTP.1"]["transcript_len"] = 20
    with pytest.raises(tc.CoordinateError) as excinfo:
        tc.build_transcript_coords(features, headers)
    assert "overlapping" in str(excinfo.value)


# ── validation of the assembled tables ───────────────────────────────────────

def test_validate_rejects_unsorted_storage_order(tc, coords):
    transcripts = coords["transcripts"].iloc[::-1].reset_index(drop=True)
    with pytest.raises(tc.CoordinateError) as excinfo:
        tc._validate(transcripts, coords["exons"])
    assert "sorted" in str(excinfo.value)


def test_validate_rejects_a_broken_offset(tc, coords):
    transcripts = coords["transcripts"].copy()
    transcripts.loc[1, "coverage_offset"] = 999
    with pytest.raises(tc.CoordinateError) as excinfo:
        tc._validate(transcripts, coords["exons"])
    assert "coverage_offset" in str(excinfo.value)


def test_validate_rejects_non_contiguous_exons(tc, coords):
    exons = coords["exons"].copy()
    exons.loc[exons.index[-1], "tx_start"] += 1
    with pytest.raises(tc.CoordinateError):
        tc._validate(coords["transcripts"], exons)


# ──────────────────── region overlays ────────────────────
# Region overlays: header vs BED conventions, the stop-codon relocation, ribopy bins.

@pytest.fixture(scope="module")
def tr():
    return __import__("transcript_regions")


# Modelled on the real GAPDH record: L=1348, UTR5:1-151, CDS:152-1159, UTR3:1160-1348.
# The BED for it is UTR5 [0,151) CDS [151,1156) UTR3 [1156,1348).
GAPDH_NAME = ("ENST00000396861.5|ENSG00000111640.15|OTTHUMG00000137379.3|"
              "OTTHUMT00000268060.1|GAPDH-205|GAPDH|1348|UTR5:1-151|CDS:152-1159|"
              "UTR3:1160-1348|")
# A stopless transcript, modelled on ENST00000625258.1: CDS runs to the transcript end.
STOPLESS_NAME = ("ENSTSTOPLESS.1|ENSGX.1|-|-|X-201|XGENE|192|UTR5:1-61|CDS:62-192|")
# No 5'UTR and no 3'UTR in the header, CDS spans the whole transcript.
WHOLE_NAME = "ENSTWHOLE.1|ENSGY.1|-|-|Y-201|YGENE|492|CDS:1-492|"


@pytest.fixture
def appris(tmp_path):
    path = tmp_path / "lengths.tsv"
    path.write_text("%s\t1348\n%s\t192\n%s\t492\n"
                    % (GAPDH_NAME, STOPLESS_NAME, WHOLE_NAME))
    return path


@pytest.fixture
def headers(tr, appris):
    return tr.parse_reference_headers(appris)


# ── header parsing ───────────────────────────────────────────────────────────

def test_header_fields_are_parsed(headers):
    gapdh = headers["ENST00000396861.5"]
    assert gapdh["gene_id"] == "ENSG00000111640.15"
    assert gapdh["transcript_name"] == "GAPDH-205"
    assert gapdh["gene_name"] == "GAPDH"
    assert gapdh["transcript_len"] == 1348
    assert gapdh["hdr_utr5"] == (1, 151)
    assert gapdh["hdr_cds"] == (152, 1159)      # 1-based inclusive, stop INCLUDED
    assert gapdh["hdr_utr3"] == (1160, 1348)


def test_absent_header_regions_are_none_not_zero_length(headers):
    assert headers["ENSTWHOLE.1"]["hdr_utr5"] is None
    assert headers["ENSTWHOLE.1"]["hdr_utr3"] is None


def test_a_header_without_a_cds_is_rejected(tr, tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("ENSTNOCDS.1|ENSGZ.1|-|-|Z-201|ZGENE|100|UTR5:1-100|\t100\n")
    with pytest.raises(tr.RegionError) as excinfo:
        tr.parse_reference_headers(path)
    assert "CDS" in str(excinfo.value)


# ── the stop-codon relocation ────────────────────────────────────────────────

def test_stop_codon_is_moved_from_cds_into_utr3(tr, headers):
    """The measured rule: header CDS 152-1159 (stop in) -> normalized [151, 1156) (stop out)."""
    regions = tr.derive_normalized_regions(headers["ENST00000396861.5"], True)
    assert regions == {"UTR5": (0, 151), "CDS": (151, 1156), "UTR3": (1156, 1348)}
    assert regions["CDS"][1] - regions["CDS"][0] == 1005      # the published GAPDH cds_len


def test_without_an_annotated_stop_nothing_is_relocated(tr, headers):
    regions = tr.derive_normalized_regions(headers["ENSTSTOPLESS.1"], False)
    assert regions == {"UTR5": (0, 61), "CDS": (61, 192)}
    assert "UTR3" not in regions                              # nothing to relocate into


def test_a_stop_only_utr3_is_created_when_the_header_has_none(tr, headers):
    """495 real transcripts do exactly this: no header UTR3, but the relocated stop makes one."""
    regions = tr.derive_normalized_regions(headers["ENSTWHOLE.1"], True)
    assert regions == {"CDS": (0, 489), "UTR3": (489, 492)}


# ── tiling and normalization ─────────────────────────────────────────────────

def test_regions_tile_the_transcript_exactly(tr, headers):
    for tid, has_stop in (("ENST00000396861.5", True), ("ENSTSTOPLESS.1", False),
                          ("ENSTWHOLE.1", True)):
        regions = tr.derive_normalized_regions(headers[tid], has_stop)
        tr.check_tiling(tid, regions, headers[tid]["transcript_len"])


def test_a_gap_in_the_tiling_is_fatal(tr):
    with pytest.raises(tr.RegionError) as excinfo:
        tr.check_tiling("T", {"UTR5": (0, 10), "CDS": (12, 20)}, 20)
    assert "contiguous" in str(excinfo.value)


def test_regions_not_reaching_the_end_are_fatal(tr):
    with pytest.raises(tr.RegionError):
        tr.check_tiling("T", {"CDS": (0, 15)}, 20)


# ── the BED cross-check ──────────────────────────────────────────────────────

def _bed_line(name, start, end, label):
    return "%s\t%d\t%d\t%s\t0\t+\n" % (name, start, end, label)


def test_a_matching_bed_is_accepted_and_marks_the_source(tr, headers, tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text(
        _bed_line(GAPDH_NAME, 0, 151, "UTR5")
        + _bed_line(GAPDH_NAME, 151, 1156, "CDS")
        + _bed_line(GAPDH_NAME, 1156, 1348, "UTR3"))
    only = {"ENST00000396861.5": headers["ENST00000396861.5"]}
    rows, summary = tr.build_regions(only, {"ENST00000396861.5"},
                                     tr.parse_actual_regions_bed(bed))
    assert summary["n_checked_against_bed"] == 1
    assert {r["label"]: (r["start"], r["end"]) for r in rows} == {
        "UTR5": (0, 151), "CDS": (151, 1156), "UTR3": (1156, 1348)}
    assert all(r["source"] == "bed" for r in rows)
    # both raw conventions are preserved verbatim, and they differ at the CDS end
    cds = next(r for r in rows if r["label"] == "CDS")
    assert cds["raw_header_end_1based"] == 1159
    assert cds["raw_bed_end"] == 1156


def test_a_disagreeing_bed_is_fatal(tr, headers, tmp_path):
    """The BED is a cross-check, not a fallback: disagreement means the rule is wrong here."""
    bed = tmp_path / "regions.bed"
    bed.write_text(
        _bed_line(GAPDH_NAME, 0, 151, "UTR5")
        + _bed_line(GAPDH_NAME, 151, 1150, "CDS")          # wrong end
        + _bed_line(GAPDH_NAME, 1150, 1348, "UTR3"))
    only = {"ENST00000396861.5": headers["ENST00000396861.5"]}
    with pytest.raises(tr.RegionError) as excinfo:
        tr.build_regions(only, {"ENST00000396861.5"}, tr.parse_actual_regions_bed(bed))
    assert "disagree" in str(excinfo.value)


def test_an_unknown_bed_label_is_rejected(tr, tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text(_bed_line(GAPDH_NAME, 0, 10, "PROMOTER"))
    with pytest.raises(tr.RegionError) as excinfo:
        tr.parse_actual_regions_bed(bed)
    assert "PROMOTER" in str(excinfo.value)


# ── ribopy bins ──────────────────────────────────────────────────────────────

def test_ribo_bins_match_the_region_lib_formula(tr, headers):
    """Boundaries are a verbatim port of region_lib.classify, using the STOP-INCLUSIVE end."""
    left, right = 35, 10
    rows = tr.build_ribo_region_bins({"ENST00000396861.5": headers["ENST00000396861.5"]},
                                     left, right)
    bins = {r["label"]: (r["start"], r["end"]) for r in rows}
    start_site, stop_site = 151, 1159                      # header CDS end, stop-inclusive
    assert bins["UTR5_OUTER"] == (0, start_site - left)
    assert bins["START_WINDOW"] == (start_site - left, start_site + right + 1)
    assert bins["CDS_CORE"] == (start_site + right + 1, stop_site - left)
    assert bins["STOP_WINDOW"] == (stop_site - left, stop_site + right + 1)
    assert bins["UTR3_OUTER"] == (stop_site + right + 1, 1348)


def test_ribo_bins_keep_the_historical_aliases(tr, headers):
    rows = tr.build_ribo_region_bins({"ENST00000396861.5": headers["ENST00000396861.5"]},
                                     35, 10)
    assert {r["label"]: r["ribopy_alias"] for r in rows} == {
        "UTR5_OUTER": "UTR5", "START_WINDOW": "UTR5J", "CDS_CORE": "CDS",
        "STOP_WINDOW": "UTR3J", "UTR3_OUTER": "UTR3"}


def test_ribo_bins_are_clipped_not_inverted_on_a_short_transcript(tr, headers):
    """Spans wider than the transcript must not produce end < start."""
    rows = tr.build_ribo_region_bins({"ENSTSTOPLESS.1": headers["ENSTSTOPLESS.1"]}, 300, 300)
    assert rows, "expected at least one clipped bin"
    assert all(r["end"] > r["start"] for r in rows)
    assert all(0 <= r["start"] and r["end"] <= 192 for r in rows)


def test_negative_spans_are_rejected(tr, headers):
    with pytest.raises(tr.RegionError):
        tr.build_ribo_region_bins(headers, -1, 10)
