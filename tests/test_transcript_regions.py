"""Region overlays: header vs BED conventions, the stop-codon relocation, ribopy bins.

The behaviour under test:
the reference header's CDS INCLUDES the stop codon, `actual_regions.bed`'s CDS EXCLUDES it
and assigns those 3 nt to UTR3, and the two agree on the CDS start. The normalized form
follows the BED, and is derivable from the header alone given whether a stop is annotated.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

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
def tr():
    return _load("transcript_regions")


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


def test_the_heuristic_agrees_with_the_gtf_rule_on_these_cases(tr, headers):
    """The GTF stop_codon feature is authoritative; the heuristic is the no-GTF fallback."""
    heuristic = tr.heuristic_stop_codon_ids(headers)
    assert "ENST00000396861.5" in heuristic      # has a header UTR3
    assert "ENSTWHOLE.1" in heuristic            # no UTR3 but CDS length 492 % 3 == 0
    assert "ENSTSTOPLESS.1" not in heuristic     # no UTR3 and CDS length 131 % 3 != 0


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


def test_bin_provenance_records_how_the_spans_were_chosen(tr):
    provenance = tr.ribo_bin_provenance(35, 10, "cli_default")
    assert provenance["algorithm"] == "ribopy_get_extended_boundaries"
    assert provenance["parameter_source"] == "cli_default"
    assert provenance["stop_site_source"] == "header_cds_end_stop_inclusive"


def test_negative_spans_are_rejected(tr, headers):
    with pytest.raises(tr.RegionError):
        tr.build_ribo_region_bins(headers, -1, 10)
