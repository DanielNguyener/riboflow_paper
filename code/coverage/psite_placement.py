#!/usr/bin/env python3
"""P-site placement: walk the offset along the READ, using the alignment's CIGAR."""
from __future__ import annotations

import collections
import sys
from pathlib import Path

PSITE_PLACEMENT = "cigar_aware"

# CIGAR operations that consume the reference / the query, by pysam opcode.
_CONSUMES_REFERENCE = frozenset((0, 2, 3, 7, 8))
_CONSUMES_QUERY = frozenset((0, 1, 4, 7, 8))
_PURE_MATCH_OPS = frozenset((0, 7, 8))

CIGAR_CODES = "MIDNSHP=X"

class PlacementError(RuntimeError):
    pass

def is_pure_match(read) -> bool:
    """True when the alignment is a single run of match/mismatch operations."""
    cigar = read.cigartuples
    return bool(cigar) and all(op in _PURE_MATCH_OPS for op, _n in cigar)

def cigar_signature(read) -> str:
    """A compact op-set signature, e.g. 'M', 'MN', 'MDN' -- for grouping in reports."""
    return "".join(sorted({CIGAR_CODES[op] for op, _n in (read.cigartuples or [])}))

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

    Same table and filter as `load_offsets`, so Figure 4 (no offset applied) shares the population.
    """
    return sorted(int(r) for r in _selected_rows(qc_csv, sample, ("read_length",))["read_length"])

def summarize_placements(bam_path: Path, offsets: dict, limit: int = None) -> dict:
    """Stream a genome BAM and describe how its reads are placed.

    Read filtering matches the coverage builder exactly (same predicate, same window).
    """
    import pysam

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
    import bam_inputs

    counters = collections.Counter()
    by_cigar = collections.Counter()
    undefined_by_cigar = collections.Counter()

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if not bam_inputs.is_unique_genome_read(read):
                continue
            offset = offsets.get(read.query_length)
            if offset is None:
                continue

            counters["considered"] += 1
            signature = cigar_signature(read)
            by_cigar[signature] += 1
            counters["pure_match" if is_pure_match(read) else "non_pure_match"] += 1
            if "N" in signature:
                counters["spliced"] += 1
            if "I" in signature:
                counters["insertion"] += 1
            if "D" in signature:
                counters["deletion"] += 1
            if "S" in signature:
                counters["soft_clipped"] += 1

            if place(read, offset) is None:
                counters["undefined"] += 1
                undefined_by_cigar[signature] += 1
            else:
                counters["placed"] += 1

            if limit and counters["considered"] >= limit:
                break
    finally:
        bam.close()

    considered = counters["considered"]
    return {
        "bam": Path(bam_path).name,
        "psite_placement": PSITE_PLACEMENT,
        "genome_uniqueness": "NH==1",
        "read_lengths": sorted(offsets),
        "offsets": {str(k): v for k, v in sorted(offsets.items())},
        "counts": dict(counters),
        "pct": {
            key: (100.0 * counters[key] / considered) if considered else 0.0
            for key in ("pure_match", "non_pure_match", "spliced", "insertion", "deletion",
                        "soft_clipped", "placed", "undefined")
        },
        "cigar_signatures": dict(by_cigar.most_common(12)),
        "undefined_by_cigar_signature": dict(undefined_by_cigar.most_common(12)),
    }

