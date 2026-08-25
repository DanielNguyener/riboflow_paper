"""The shared transcript coordinate: exon map, ordering, reconciliation, round-trip.

Synthetic throughout -- no GTF download, no BAM. The geometry exercises every
rule the coordinate depends on: both strands, a splice junction, a single-exon transcript,
and a transcript whose exons appear in the GTF in the wrong order.

The map is built from complete `exon` features rather than the region-split CDS + UTR
caches so that the spliced length equals the reference length for every transcript.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
COVERAGE_DIR = REPO / "code" / "coverage"


def _load(name):
    """Load a coverage module by path. Ordinary deterministic import -- no sys.modules
    pre-registration, which code in this repository does not use."""
    if str(COVERAGE_DIR) not in sys.path:
        sys.path.insert(0, str(COVERAGE_DIR))
    spec = importlib.util.spec_from_file_location(name, COVERAGE_DIR / ("%s.py" % name))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tc():
    return _load("transcript_coords")


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


# ── coordinate mapping ───────────────────────────────────────────────────────

@pytest.mark.parametrize("tid", sorted(GEOMETRY))
def test_every_position_round_trips(tc, coords, tid):
    index = int(coords["transcripts"].index[
        coords["transcripts"]["transcript_id"] == tid][0])
    length = int(coords["transcripts"].at[index, "transcript_len"])
    positions = np.arange(length)
    genomic = tc.tx_to_genomic(coords, index, positions)
    assert (genomic >= 0).all()
    assert np.array_equal(tc.genomic_to_tx(coords, index, genomic), positions)


def test_plus_strand_mapping_is_explicit(tc, coords):
    index = 1                                              # ENSTP.1
    # tx 0 -> first base of exon 1; tx 11 -> last base of exon 1; tx 12 -> first of exon 2
    assert tc.tx_to_genomic(coords, index, [0, 11, 12, 19]).tolist() == [100, 111, 200, 207]


def test_minus_strand_mapping_mirrors_within_each_exon(tc, coords):
    index = 0                                              # ENSTM.2, chr2 '-'
    # 5' exon is [500, 510); its 5'-most base is the HIGHEST coordinate, 509
    assert tc.tx_to_genomic(coords, index, [0, 1, 9, 10, 14]).tolist() == [509, 508, 500, 304, 300]


def test_positions_outside_the_transcript_yield_minus_one(tc, coords):
    index = 1
    assert tc.tx_to_genomic(coords, index, [-1, 20, 10 ** 6]).tolist() == [-1, -1, -1]


def test_intronic_genomic_positions_yield_minus_one(tc, coords):
    """The gap between the two exons of ENSTP.1 is not part of the coordinate."""
    index = 1
    assert tc.genomic_to_tx(coords, index, [150, 112, 199]).tolist() == [-1, -1, -1]


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
