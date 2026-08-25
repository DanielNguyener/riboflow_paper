"""The two full-coordinate summations, against a naive reference.

`per_transcript_sums` and `region_slice_sums` reduce a 70.5-million-position int32 array to
one number per transcript. Both sum each span separately rather than building a
full-coordinate int64 temporary, which on a human build is the difference between 0.8 MB and
564 MB (or 1.1 GB). They are checked here against the obvious implementation, on the cases
where a 64-bit accumulator actually matters.
"""
from __future__ import annotations

import numpy as np
import pytest


def naive_per_transcript(values, coverage_offset, transcript_len):
    return np.array([int(values[o:o + n].sum(dtype=np.int64))
                     for o, n in zip(coverage_offset, transcript_len)], dtype=np.int64)


def naive_region(values, coverage_offset, starts, ends):
    out = []
    for offset, start, end in zip(coverage_offset, starts, ends):
        lo = offset + max(int(start), 0)
        hi = offset + max(int(end), int(start))
        out.append(int(values[lo:hi].sum(dtype=np.int64)))
    return np.array(out, dtype=np.int64)


def geometry(lengths):
    lengths = np.asarray(lengths, dtype=np.int64)
    offsets = np.zeros(len(lengths), dtype=np.int64)
    offsets[1:] = np.cumsum(lengths)[:-1]
    return offsets, lengths


# ── per-transcript sums ──────────────────────────────────────────────────────

@pytest.mark.parametrize("lengths", [[5], [1, 1, 1], [3, 7, 2], [10, 1], [1, 10]])
def test_per_transcript_matches_the_naive_sum(bsc, lengths):
    offsets, lens = geometry(lengths)
    values = np.arange(int(lens.sum()), dtype=np.int32)
    assert np.array_equal(bsc.per_transcript_sums(values, offsets, lens),
                          naive_per_transcript(values, offsets, lens))


def test_per_transcript_on_an_all_zero_coordinate(bsc):
    offsets, lens = geometry([4, 6])
    values = np.zeros(10, dtype=np.int32)
    assert list(bsc.per_transcript_sums(values, offsets, lens)) == [0, 0]


def test_a_length_one_transcript_is_its_own_single_value(bsc):
    """A transcript of length 1 must sum to that one base, not to the rest of the array."""
    offsets, lens = geometry([1, 1, 3])
    values = np.array([7, 9, 1, 2, 3], dtype=np.int32)
    assert list(bsc.per_transcript_sums(values, offsets, lens)) == [7, 9, 6]


def test_the_result_is_int64_even_though_the_input_is_int32(bsc):
    offsets, lens = geometry([3])
    assert bsc.per_transcript_sums(np.ones(3, dtype=np.int32), offsets, lens).dtype == np.int64


def test_a_transcript_total_may_exceed_int32(bsc):
    """Each POSITION fits in int32, but a transcript's total need not: 200,000 bases at
    depth 20,000 overflows a 32-bit accumulator and must not wrap."""
    offsets, lens = geometry([200_000])
    values = np.full(200_000, 20_000, dtype=np.int32)
    total = bsc.per_transcript_sums(values, offsets, lens)
    assert total[0] == 200_000 * 20_000 == 4_000_000_000
    assert total[0] > np.iinfo(np.int32).max


def test_neither_summation_allocates_a_full_length_int64_copy(bsc):
    """A regression to `values.astype(np.int64)` or to a full `np.cumsum` prefix array is
    invisible in the result and costs hundreds of megabytes on a real build, so it is
    measured rather than reviewed."""
    import tracemalloc

    n = 2_000_000
    offsets, lens = geometry([n // 2, n // 2])
    values = np.ones(n, dtype=np.int32)
    int64_copy = n * 8

    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        bsc.per_transcript_sums(values, offsets, lens)
        per_transcript_peak = tracemalloc.get_traced_memory()[1]

        tracemalloc.reset_peak()
        bsc.region_slice_sums(values, offsets, np.array([0, 0]),
                              np.array([n // 2, n // 2]))
        region_peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()

    # Neither may allocate anything proportional to the coordinate: the result is one
    # number per transcript, so a full-length temporary is the whole cost being avoided.
    budget = int64_copy / 100
    assert per_transcript_peak < budget, (
        "per_transcript_sums allocated %.1f MB; a full-length int64 array is %.1f MB"
        % (per_transcript_peak / 1e6, int64_copy / 1e6))
    assert region_peak < budget, (
        "region_slice_sums allocated %.1f MB; a full-length int64 array is %.1f MB"
        % (region_peak / 1e6, int64_copy / 1e6))


# ── region slice sums ────────────────────────────────────────────────────────

def test_region_slices_match_the_naive_sum(bsc):
    offsets, lens = geometry([10, 10, 10])
    values = np.arange(30, dtype=np.int32)
    starts = np.array([0, 2, 5], dtype=np.int64)
    ends = np.array([10, 8, 6], dtype=np.int64)
    assert np.array_equal(bsc.region_slice_sums(values, offsets, starts, ends),
                          naive_region(values, offsets, starts, ends))


def test_a_window_at_each_transcript_boundary(bsc):
    """The first and last base of a transcript, and a window spanning the whole of it:
    an off-by-one in the window arithmetic shows up here and nowhere else."""
    offsets, lens = geometry([4, 4])
    values = np.array([1, 2, 3, 4, 10, 20, 30, 40], dtype=np.int32)
    first = bsc.region_slice_sums(values, offsets, np.array([0, 0]), np.array([1, 1]))
    last = bsc.region_slice_sums(values, offsets, np.array([3, 3]), np.array([4, 4]))
    whole = bsc.region_slice_sums(values, offsets, np.array([0, 0]), np.array([4, 4]))
    assert list(first) == [1, 10]
    assert list(last) == [4, 40]
    assert list(whole) == [10, 100]


def test_an_empty_window_contributes_zero(bsc):
    """A CDS shorter than twice the trim has end <= start and must sum to 0, not to a
    negative slice or a wrapped index."""
    offsets, lens = geometry([5, 5])
    values = np.ones(10, dtype=np.int32)
    sums = bsc.region_slice_sums(values, offsets, np.array([3, 1]), np.array([1, 4]))
    assert list(sums) == [0, 3]


def test_a_negative_start_is_clamped(bsc):
    offsets, lens = geometry([4])
    values = np.array([1, 2, 3, 4], dtype=np.int32)
    assert list(bsc.region_slice_sums(values, offsets, np.array([-5]), np.array([2]))) == [3]


def test_region_sums_on_an_all_zero_coordinate(bsc):
    offsets, lens = geometry([4, 4])
    values = np.zeros(8, dtype=np.int32)
    assert list(bsc.region_slice_sums(values, offsets, np.array([0, 0]),
                                      np.array([4, 4]))) == [0, 0]


def test_region_sums_survive_an_int32_overflow(bsc):
    offsets, lens = geometry([100_000, 100_000])
    values = np.full(200_000, 30_000, dtype=np.int32)
    sums = bsc.region_slice_sums(values, offsets, np.array([0, 0]),
                                 np.array([100_000, 100_000]))
    assert list(sums) == [3_000_000_000, 3_000_000_000]
    assert sums.dtype == np.int64


def test_an_empty_coordinate_does_not_raise(bsc):
    """A cohort with nothing to sum is a legitimate no-op, not a crash."""
    empty = np.zeros(0, dtype=np.int32)
    result = bsc.region_slice_sums(empty, np.array([], dtype=np.int64),
                                   np.array([], dtype=np.int64),
                                   np.array([], dtype=np.int64))
    assert result.size == 0
