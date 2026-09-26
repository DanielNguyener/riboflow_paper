"""Half-open [start, end) interval arithmetic shared by the read-ID analyses."""
from __future__ import annotations


def merge(intervals):
    """Union as sorted, non-overlapping (start, end) tuples."""
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [tuple(i) for i in out]


def subtract(intervals, holes):
    """Merged `intervals` minus merged `holes`."""
    out = []
    holes = merge(holes)
    for start, end in merge(intervals):
        cursor = start
        for h_start, h_end in holes:
            if h_end <= cursor or h_start >= end:
                continue
            if h_start > cursor:
                out.append((cursor, h_start))
            cursor = max(cursor, h_end)
        if cursor < end:
            out.append((cursor, end))
    return out


def merged_length(intervals):
    return sum(end - start for start, end in merge(intervals))
