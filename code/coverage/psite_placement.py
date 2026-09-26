#!/usr/bin/env python3
"""P-site placement: walk the offset along the READ, using the alignment's CIGAR."""
from __future__ import annotations

from pathlib import Path

PSITE_PLACEMENT = "cigar_aware"

class PlacementError(RuntimeError):
    pass

def place(read, offset: int):
    """Reference position `offset` aligned read bases from the read's 5' end.

    Returns a 0-based genomic position, or None when the read has fewer than
    `offset + 1` aligned bases.
    """
    if offset < 0:
        raise PlacementError("P-site offset must be >= 0, got %r" % (offset,))
    pairs = read.get_aligned_pairs(matches_only=True)
    if len(pairs) <= offset:
        return None
    return pairs[len(pairs) - 1 - offset][1] if read.is_reverse else pairs[offset][1]

def _selected_rows(qc_csv: Path, sample: str, columns):
    """The sample's SELECTED (`in_phase1`) rows of a `readlen_window_qc.csv`-shaped table.

    The one reader for the whole repository; a second loader could only disagree with it.
    """
    import pandas as pd

    frame = pd.read_csv(qc_csv)
    for column in ("sample", "in_phase1") + tuple(columns):
        if column not in frame.columns:
            raise PlacementError(
                "%s has no %r column; expected a readlen_window_qc.csv-shaped table with "
                "%s" % (qc_csv, column, "sample, read_length, in_phase1, psite_offset"))
    subset = frame[frame["sample"] == sample]
    if subset.empty:
        available = ", ".join(sorted(frame["sample"].astype(str).unique())[:8])
        raise PlacementError(
            "sample %r is not in %s. Samples present include: %s" % (sample, qc_csv, available))
    in_phase1 = subset["in_phase1"].map(
        lambda v: str(v).strip().lower() in ("true", "1"))
    subset = subset[in_phase1]
    if subset.empty:
        raise PlacementError("sample %r has no in_phase1 read lengths in %s" % (sample, qc_csv))
    return subset

def load_offsets(qc_csv: Path, sample: str) -> dict:
    """Selected read_length -> psite_offset for one sample, from a QC master table."""
    subset = _selected_rows(qc_csv, sample, ("read_length", "psite_offset"))
    return {int(r): int(o) for r, o in zip(subset["read_length"], subset["psite_offset"])}

def load_selected_lengths(qc_csv: Path, sample: str) -> list:
    """The sample's selected read lengths, ascending -- the window WITHOUT the offsets.

    Same table and filter as `load_offsets`, so Figure 3 (no offset applied) shares the population.
    """
    return sorted(int(r) for r in _selected_rows(qc_csv, sample, ("read_length",))["read_length"])
