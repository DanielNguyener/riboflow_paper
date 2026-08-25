"""The assignment-policy exposure measurement, on hand-built overlap geometry.

Two rules pick one transcript for a read that falls on several: first CDS-exon overlap for
P-sites, maximum CDS-exon overlap for footprints. `validate_assignment_policy` counts how
often the rule rather than the data makes that choice. The counts are only meaningful if
"ambiguous" and "tied" mean what they say, so each case below constructs exactly one
situation and asserts the two numbers separately.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "code" / "coverage"))

import validate_assignment_policy as vap  # noqa: E402


def cds_pyranges(rows):
    """`rows` are (transcript_id, start, end) on one strand of one contig."""
    import pyranges as pr
    return pr.PyRanges(pd.DataFrame({
        "Chromosome": ["chr1"] * len(rows),
        "Start": [r[1] for r in rows],
        "End": [r[2] for r in rows],
        "Strand": ["+"] * len(rows),
        "transcript_id": [r[0] for r in rows],
        "cds_cum_start": [0] * len(rows)}))


def psites(positions):
    return (["chr1"] * len(positions),
            np.asarray(positions, dtype=np.int64),
            ["+"] * len(positions))


def blocks(intervals):
    """`intervals` is a list per read of [(start, end), ...] aligned blocks."""
    chroms, starts, ends, strands, read_ids = [], [], [], [], []
    for index, spans in enumerate(intervals):
        for start, end in spans:
            chroms.append("chr1")
            starts.append(start)
            ends.append(end)
            strands.append("+")
            read_ids.append(index)
    return (chroms, np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64),
            strands, np.asarray(read_ids, dtype=np.int64), len(intervals))


# ── P-sites: first overlap, so ambiguous == tied ─────────────────────────────

def test_one_candidate_transcript_is_never_ambiguous():
    exposure = vap.psite_exposure(*psites([150]), cds_pyranges([("TA", 100, 200)]))
    assert exposure == {"n_candidates_total": 1, "n_ambiguous": 0, "n_tied": 0}


def test_two_candidate_transcripts_make_the_read_ambiguous():
    overlapping = cds_pyranges([("TA", 100, 200), ("TB", 120, 260)])
    exposure = vap.psite_exposure(*psites([150]), overlapping)
    assert exposure["n_candidates_total"] == 1
    assert exposure["n_ambiguous"] == 1


def test_for_psites_every_ambiguous_read_is_a_tie():
    """`first_exon_overlap` scores nothing, so there is no criterion to break a tie with:
    whichever join row came first wins, and that is the rule deciding, not the data."""
    overlapping = cds_pyranges([("TA", 100, 200), ("TB", 120, 260), ("TC", 140, 300)])
    exposure = vap.psite_exposure(*psites([150, 151]), overlapping)
    assert exposure["n_ambiguous"] == exposure["n_tied"] == 2


def test_a_psite_outside_every_cds_exon_has_no_candidate():
    exposure = vap.psite_exposure(*psites([900]), cds_pyranges([("TA", 100, 200)]))
    assert exposure["n_candidates_total"] == 0


def test_a_mixed_population_counts_only_the_ambiguous_reads():
    overlapping = cds_pyranges([("TA", 100, 200), ("TB", 120, 260)])
    #                      unambiguous(TA only)  ambiguous  unambiguous(TB only)
    exposure = vap.psite_exposure(*psites([110, 150, 250]), overlapping)
    assert exposure["n_candidates_total"] == 3
    assert exposure["n_ambiguous"] == 1


# ── footprints: maximum overlap, so ambiguous > tied ─────────────────────────

def test_unequal_overlaps_are_ambiguous_but_not_tied():
    """The data decides: the read covers 40 nt of TA and 10 nt of TB, so `idxmax` is not
    breaking a tie and the assignment is not rule-dependent."""
    overlapping = cds_pyranges([("TA", 100, 200), ("TB", 190, 300)])
    exposure = vap.footprint_exposure(blocks([[(160, 200)]]), overlapping)
    assert exposure["n_ambiguous"] == 1
    assert exposure["n_tied"] == 0


def test_equal_overlaps_are_a_genuine_tie():
    """Both transcripts contain the read's whole span, so the summed overlaps are equal and
    `idxmax` takes the first -- the rule deciding."""
    overlapping = cds_pyranges([("TA", 100, 300), ("TB", 100, 300)])
    exposure = vap.footprint_exposure(blocks([[(150, 180)]]), overlapping)
    assert exposure["n_candidates_total"] == 1
    assert exposure["n_ambiguous"] == 1
    assert exposure["n_tied"] == 1


def test_overlap_is_summed_across_a_reads_blocks():
    """A spliced read's overlap is the total across its blocks, not the largest single one:
    TA takes 10+10, TB only 15, so TA wins and there is no tie."""
    overlapping = cds_pyranges([("TA", 100, 130), ("TA", 200, 230), ("TB", 205, 220)])
    exposure = vap.footprint_exposure(blocks([[(120, 130), (200, 210)]]), overlapping)
    assert exposure["n_ambiguous"] == 1
    assert exposure["n_tied"] == 0


def test_a_footprint_touching_one_transcript_is_not_ambiguous():
    exposure = vap.footprint_exposure(blocks([[(120, 150)]]),
                                      cds_pyranges([("TA", 100, 200)]))
    assert exposure == {"n_candidates_total": 1, "n_ambiguous": 0, "n_tied": 0}


def test_a_zero_length_overlap_is_not_a_candidate():
    """PyRanges reports a book-ended interval as a join hit; a 0-nt overlap is not evidence
    of anything and must not inflate the denominator."""
    exposure = vap.footprint_exposure(blocks([[(200, 230)]]),
                                      cds_pyranges([("TA", 100, 200), ("TB", 200, 300)]))
    assert exposure["n_candidates_total"] == 1, "only TB genuinely overlaps"
    assert exposure["n_ambiguous"] == 0


# ── determinism ──────────────────────────────────────────────────────────────

def test_repeated_execution_gives_identical_counts():
    """The rules depend on PyRanges join order, which is not a documented contract. That it
    is stable in the pinned environment is the reason the historical order was kept, so it
    is re-checked rather than assumed."""
    overlapping = cds_pyranges([("TA", 100, 200), ("TB", 120, 260), ("TC", 140, 300)])
    positions = list(range(140, 200))
    spans = [[(p, p + 30)] for p in positions]
    first_psite = vap.psite_exposure(*psites(positions), overlapping)
    first_footprint = vap.footprint_exposure(blocks(spans), overlapping)
    for _ in range(3):
        assert vap.psite_exposure(*psites(positions), overlapping) == first_psite
        assert vap.footprint_exposure(blocks(spans), overlapping) == first_footprint


def test_an_empty_population_reports_zeros_not_an_error():
    empty = cds_pyranges([("TA", 100, 200)])
    assert vap.psite_exposure([], np.empty(0, dtype=np.int64), [], empty)["n_ambiguous"] == 0
    assert vap.footprint_exposure(blocks([]), empty)["n_ambiguous"] == 0


# ── the program is a validator, not a pipeline stage ─────────────────────────

def test_it_is_not_a_make_tables_stage():
    text = (REPO / "code" / "make_tables.py").read_text()
    assert "validate_assignment_policy" not in text


def test_it_refuses_a_missing_input_by_name():
    import subprocess
    result = subprocess.run(
        [sys.executable, str(REPO / "code" / "coverage" / "validate_assignment_policy.py"),
         "--sample", "X", "--genome-bam", "/nope.bam", "--gtf", "/nope.gtf",
         "--appris", "/nope.tsv", "--qc-genome", "/nope.csv", "--output", "/tmp/x.json"],
        capture_output=True, text=True)
    assert result.returncode != 0
    assert "--genome-bam" in result.stderr
