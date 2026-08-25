"""CDS-boundary exclusion: what `--trim` drops, what it keeps, and why it must be a
multiple of three.

The trim exists because the first and last codons of a CDS carry initiation and
termination pile-ups that are not elongation signal. It is applied SYMMETRICALLY to both
routes over the same CDS window, so it cannot bias one against the other.
"""
from __future__ import annotations

import numpy as np
import pytest

from conftest import (CDS_START, GEOMETRY, TRIM, TX_MINUS, TX_PLUS, TX_SHORT,
                      build_config, cds_interior, expected_footprint_genome,
                      expected_psite)


def build(inputs, **overrides):
    import build_shared_coverage
    import coverage_schema
    path, _ = build_shared_coverage.build(build_config(inputs, **overrides))
    return coverage_schema.open_coverage(path)


# ── the geometry of the trimmed window ───────────────────────────────────────

@pytest.mark.parametrize("tid", [TX_PLUS, TX_MINUS])
def test_the_interior_is_the_cds_minus_twice_the_trim(coverage, tid):
    total = GEOMETRY[tid][3]
    assert cds_interior(coverage, tid, "genome_psite").size == total - 2 * TRIM


def test_the_trim_is_recorded_in_the_file(coverage):
    assert coverage.trim == TRIM
    assert int(coverage.handle.attrs["paper_cds_trim"]) == TRIM


def test_slice_region_uses_the_files_own_trim(coverage):
    """A consumer must not have to know the trim: it reads it off the file."""
    index = coverage.index_of_transcript(TX_PLUS)
    start, end = coverage.slice_region(index, "CDS")
    untrimmed_start, untrimmed_end = coverage.slice_region(index, "CDS", trim=0)
    assert start == untrimmed_start + TRIM
    assert end == untrimmed_end - TRIM


def test_the_cds_window_sits_where_the_regions_say(coverage):
    index = coverage.index_of_transcript(TX_PLUS)
    start, end = coverage.slice_region(index, "CDS", trim=0)
    assert (start, end) == (CDS_START[TX_PLUS],
                            CDS_START[TX_PLUS] + GEOMETRY[TX_PLUS][3])


# ── P-sites are dropped; footprints are partially kept ───────────────────────

def test_a_psite_inside_the_trim_zone_is_outside_the_interior(coverage):
    """The rel-5 read's P-site is inside the 5' trim zone. It is still STORED -- the file
    keeps the whole transcript -- but it is not part of the interior view."""
    index = coverage.index_of_transcript(TX_PLUS)
    stored = coverage.get_track(index, "genome_psite")
    assert stored[CDS_START[TX_PLUS] + 5] == 1, "the P-site itself must still be recorded"
    assert cds_interior(coverage, TX_PLUS, "genome_psite").sum() == expected_psite(
        TX_PLUS).sum()


def test_a_footprint_spanning_the_boundary_is_partially_counted(coverage):
    """A P-site is a point and is either in or out. A footprint is an interval, so a read
    straddling the trim boundary contributes only the bases inside it."""
    footprint = cds_interior(coverage, TX_PLUS, "genome_footprint")
    assert footprint[0] > 0, "the rel-5 read covers the first interior base"
    assert np.array_equal(footprint, expected_footprint_genome(TX_PLUS))


# ── the frame argument ───────────────────────────────────────────────────────

def test_the_trim_keeps_the_interior_in_frame(coverage):
    """15 is divisible by 3, so the interior still begins at a codon boundary. A trim that
    was not would silently rotate the reading frame of every downstream frame statistic."""
    assert TRIM % 3 == 0
    index = coverage.index_of_transcript(TX_PLUS)
    start, _end = coverage.slice_region(index, "CDS")
    cds_start, _ = coverage.slice_region(index, "CDS", trim=0)
    assert (start - cds_start) % 3 == 0


def test_a_trim_that_is_not_a_multiple_of_three_is_rejected():
    from build_shared_coverage import main
    with pytest.raises(SystemExit) as excinfo:
        main(["--sample", "X", "--genome-bam", "g", "--transcriptome-bam", "t",
              "--gtf", "g.gtf", "--appris", "a.tsv",
              "--qc-genome", "q.csv", "--qc-txome", "q.csv", "--trim", "10"])
    assert "multiple of 3" in str(excinfo.value)


# ── trim = 0 round trip ──────────────────────────────────────────────────────

def test_trim_zero_gives_the_whole_cds(inputs):
    with build(inputs, trim=0, output=inputs.root / "t0") as zero:
        assert zero.trim == 0
        index = zero.index_of_transcript(TX_PLUS)
        start, end = zero.slice_region(index, "CDS")
        assert end - start == GEOMETRY[TX_PLUS][3]


def test_the_trimmed_interior_is_a_slice_of_the_untrimmed_one(inputs):
    """The trim must not change any stored value -- only which range is looked at."""
    with build(inputs, trim=0, output=inputs.root / "t0") as zero:
        index = zero.index_of_transcript(TX_PLUS)
        whole = zero.get_track(index, "genome_psite")
    with build(inputs, trim=TRIM, output=inputs.root / "t15") as trimmed:
        index = trimmed.index_of_transcript(TX_PLUS)
        assert np.array_equal(trimmed.get_track(index, "genome_psite"), whole)


def test_a_larger_trim_shrinks_the_interior_symmetrically(inputs):
    with build(inputs, trim=21, output=inputs.root / "t21") as wide:
        index = wide.index_of_transcript(TX_PLUS)
        start, end = wide.slice_region(index, "CDS")
        cds_start, cds_end = wide.slice_region(index, "CDS", trim=0)
        assert start - cds_start == 21
        assert cds_end - end == 21


def test_a_trim_that_would_invert_the_window_collapses_to_empty(coverage):
    """TX_SHORT's CDS is 24 nt and the trim is 15 a side. The window must come back empty
    rather than inverted -- a negative-length slice would silently read backwards."""
    index = coverage.index_of_transcript(TX_SHORT)
    start, end = coverage.slice_region(index, "CDS")
    assert end == start
    assert cds_interior(coverage, TX_SHORT, "genome_psite").size == 0
