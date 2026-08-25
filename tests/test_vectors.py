"""Coverage-vector construction, against vectors written down by hand.

These run the real builder (`build_shared_coverage.build`) over the synthetic cohort, so
what is checked is the pipeline rather than a re-implementation of it. Every expected
value is derived in `conftest.py` from the geometry and the read plan, so a disagreement
here is a disagreement with arithmetic.
"""
from __future__ import annotations

import numpy as np
import pytest

from conftest import (CDS_START, GEOMETRY, JUNCTION_PSITE_CDS_REL, TRIM, TX_MINUS,
                      TX_PLUS, TX_SHORT, cds_interior, expected_footprint_genome,
                      expected_footprint_txome, expected_psite, expected_psite_txome)

SIGNALS = ("genome_psite", "txome_psite", "genome_footprint", "txome_footprint")


# ── the four vectors ─────────────────────────────────────────────────────────

def test_genome_psite_vector(coverage):
    assert np.array_equal(cds_interior(coverage, TX_PLUS, "genome_psite"),
                          expected_psite(TX_PLUS))


def test_txome_psite_vector(coverage):
    assert np.array_equal(cds_interior(coverage, TX_PLUS, "txome_psite"),
                          expected_psite_txome(TX_PLUS))


def test_genome_footprint_vector(coverage):
    assert np.array_equal(cds_interior(coverage, TX_PLUS, "genome_footprint"),
                          expected_footprint_genome(TX_PLUS))


def test_txome_footprint_vector(coverage):
    assert np.array_equal(cds_interior(coverage, TX_PLUS, "txome_footprint"),
                          expected_footprint_txome(TX_PLUS))


def test_the_two_routes_agree_on_psites_for_the_minus_strand_transcript(coverage):
    """TX_MINUS has no junction spanner and no boundary-straddling read, so the two routes
    must produce the SAME cds-relative vector -- which is the premise of comparing them.
    It is also the minus-strand case, where the coordinate runs against the genome."""
    genome = cds_interior(coverage, TX_MINUS, "genome_psite")
    txome = cds_interior(coverage, TX_MINUS, "txome_psite")
    assert np.array_equal(genome, txome)
    assert genome.sum() == 2


def test_minus_strand_footprints_agree_between_routes(coverage):
    assert np.array_equal(cds_interior(coverage, TX_MINUS, "genome_footprint"),
                          cds_interior(coverage, TX_MINUS, "txome_footprint"))


# ── the junction spanner: the case the placement rule exists for ─────────────

def test_a_junction_spanning_read_is_counted_at_the_base_it_covers(coverage):
    """Walking the offset along the REFERENCE would put this read's P-site at genomic
    1062, inside the intron, where the CDS clip discards it. Walking it along the READ
    puts it at cds_rel 62, which is where the read actually is."""
    genome = cds_interior(coverage, TX_PLUS, "genome_psite")
    assert genome[JUNCTION_PSITE_CDS_REL - TRIM] == 1


def test_the_junction_read_is_the_only_psite_difference_between_routes(coverage):
    """The genome and transcriptome P-site vectors for TX_PLUS differ by exactly one
    event at exactly one position: the junction spanner, which cannot exist on a spliced
    reference."""
    difference = (cds_interior(coverage, TX_PLUS, "genome_psite").astype(int)
                  - cds_interior(coverage, TX_PLUS, "txome_psite").astype(int))
    assert difference.sum() == 1
    assert list(np.nonzero(difference)[0]) == [JUNCTION_PSITE_CDS_REL - TRIM]


def test_the_junction_footprint_crosses_the_exon_boundary(coverage):
    """Its 30 bases are 10 in exon 1 and 20 in exon 2, contiguous in TRANSCRIPT
    coordinates -- cds_rel 50..79 with no gap -- even though 940 nt of genome separate
    them."""
    footprint = cds_interior(coverage, TX_PLUS, "genome_footprint")
    assert (footprint[50 - TRIM:80 - TRIM] >= 1).all()


# ── stacking, and the filters ────────────────────────────────────────────────

def test_psites_from_different_read_lengths_stack_on_one_base(coverage):
    """Two phase-1 lengths sharing an offset can put two P-sites on one base even in a
    position-deduplicated BAM. A vector that topped out at 1 would mean the second read
    was silently dropped."""
    assert cds_interior(coverage, TX_PLUS, "genome_psite")[20 - TRIM] == 2
    assert cds_interior(coverage, TX_PLUS, "txome_psite")[20 - TRIM] == 2


@pytest.mark.parametrize("rel,why", [
    (25, "MAPQ-3 multimapper (genome) / MAPQ-41 (transcriptome)"),
    (26, "read length outside the phase-1 window"),
    (27, "secondary alignment"),
])
def test_excluded_reads_never_produce_a_psite(coverage, rel, why):
    for signal in ("genome_psite", "txome_psite"):
        assert cds_interior(coverage, TX_PLUS, signal)[rel - TRIM] == 0, \
            "%s picked up a read it must exclude: %s" % (signal, why)


def test_supplementary_alignments_are_excluded(coverage):
    assert cds_interior(coverage, TX_PLUS, "genome_psite")[28 - TRIM] == 0


def test_intergenic_reads_reach_no_transcript(coverage):
    """8 P-sites total: the 6 planned, the junction spanner and the TX_SHORT read. The
    intergenic read, the multimapper, the out-of-window read, the secondary and the
    supplementary are all excluded."""
    total = sum(int(coverage.get_track(i, "genome_psite").sum())
                for i in range(coverage.n_transcripts))
    assert total == 8


def test_only_phase1_read_lengths_and_their_offsets_are_used(coverage):
    offsets = coverage.provenance["offsets"]["genome"]
    assert sorted(int(k) for k in offsets) == [28, 30]
    assert set(offsets.values()) == {12}


def test_the_other_sample_in_the_qc_master_is_not_used(coverage):
    """The QC master also describes OTHER_SAMPLE at a different offset. Picking that up
    would shift every P-site by one base and still look entirely plausible."""
    assert coverage.sample == "SYN"
    assert set(coverage.provenance["offsets"]["genome"].values()) == {12}


def test_changing_the_offset_shifts_every_psite(inputs, tmp_path):
    """Rebuild with a different offset and the P-sites must move by exactly that much."""
    import build_shared_coverage
    import coverage_schema

    from conftest import OFFSET, SAMPLE, build_config, build_qc_master

    shifted_qc = build_qc_master(tmp_path / "qc_shift.csv", sample=SAMPLE,
                                 offset=OFFSET + 3)
    path, _ = build_shared_coverage.build(build_config(
        inputs, qc_genome=shifted_qc, qc_txome=shifted_qc,
        output=tmp_path / "shifted"))
    with coverage_schema.open_coverage(path) as shifted:
        moved = cds_interior(shifted, TX_MINUS, "txome_psite")
    baseline = expected_psite_txome(TX_MINUS)
    # A larger offset walks further along the read, i.e. to a HIGHER cds_rel.
    assert list(np.nonzero(moved)[0]) == [i + 3 for i in np.nonzero(baseline)[0]]


# ── transcripts with no usable interior ──────────────────────────────────────

def test_a_transcript_shorter_than_twice_the_trim_has_an_empty_interior(coverage):
    """TX_SHORT has reads on both routes; it is excluded from interior comparisons for
    being too short, not for being empty -- and those are different states."""
    assert GEOMETRY[TX_SHORT][3] < 2 * TRIM
    assert cds_interior(coverage, TX_SHORT, "genome_psite").size == 0
    index = coverage.index_of_transcript(TX_SHORT)
    assert coverage.event_counts(index)["genome_psite"] == 1


def test_every_transcript_is_stored_even_when_it_has_no_reads(coverage):
    """Absence is recorded, not omitted: every transcript has a row and a full-length
    vector, so 'no reads' is distinguishable from 'not in the file'."""
    assert coverage.n_transcripts == 3
    for tid in (TX_PLUS, TX_MINUS, TX_SHORT):
        index = coverage.index_of_transcript(tid)
        assert coverage.get_track(index, "genome_psite").size == \
            int(coverage.transcript_len[index])


# ── the stored coordinate is the FULL transcript, not the CDS ────────────────

def test_the_stored_vector_spans_the_whole_transcript(coverage):
    assert coverage.get_track(coverage.index_of_transcript(TX_PLUS),
                              "genome_psite").size == 200


def test_utr_reads_are_kept_in_the_utr_rather_than_discarded(coverage):
    """Two transcriptome reads sit wholly inside the UTRs. A CDS-only pipeline threw them
    away; on the full transcript coordinate they are simply UTR coverage."""
    index = coverage.index_of_transcript(TX_PLUS)
    track = coverage.get_track(index, "txome_footprint")
    assert track[:CDS_START[TX_PLUS]].sum() > 0, "the wholly-5'UTR read is missing"
    assert track[CDS_START[TX_PLUS] + GEOMETRY[TX_PLUS][3]:].sum() > 0, \
        "the wholly-3'UTR read is missing"


def test_the_cds_view_is_a_slice_not_a_separate_computation(coverage):
    for tid in (TX_PLUS, TX_MINUS):
        assert np.array_equal(cds_interior(coverage, tid, "genome_footprint"),
                              expected_footprint_genome(tid))
        assert np.array_equal(cds_interior(coverage, tid, "txome_footprint"),
                              expected_footprint_txome(tid))


# ── per-transcript bookkeeping is derived, not stored ────────────────────────

@pytest.mark.parametrize("signal", SIGNALS)
def test_event_counts_are_the_sums_of_the_stored_vectors(coverage, signal):
    for index in range(coverage.n_transcripts):
        assert coverage.event_counts(index)[signal] == \
            int(coverage.get_track(index, signal).sum())


def test_cds_window_sums_agree_with_explicit_slicing(coverage):
    """`window_sums` is what the concordance keys are derived from; it must equal the
    obvious per-transcript slice for every trim."""
    for signal in SIGNALS:
        values = coverage.signal(signal)
        for trim in (0, TRIM, 99):
            sums = coverage.cds_window_sums(values, trim)
            for index in range(coverage.n_transcripts):
                start, end = coverage.slice_region(index, "CDS", trim=trim)
                assert int(sums[index]) == int(coverage.get_track(index, signal)[start:end].sum())


def test_the_file_stores_nothing_derivable(coverage):
    """Schema 3 carries the arrays and the CDS bounds and nothing else."""
    assert set(coverage.handle) == {"transcripts", "coverage"}
    assert set(coverage.handle["transcripts"]) == {
        "transcript_id", "gene_id", "gene_name", "transcript_len", "cds_start", "cds_end",
        "coverage_offset"}
    assert set(coverage.handle["coverage"]) == set(SIGNALS)
