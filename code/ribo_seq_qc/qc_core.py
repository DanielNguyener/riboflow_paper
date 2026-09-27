#!/usr/bin/env python3
"""The route-independent half of the Ribo-seq QC steps."""
from __future__ import annotations

import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

_COMMON = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "common")
if _COMMON not in sys.path:
    sys.path.insert(0, _COMMON)

SELECT_MIN_LEN, SELECT_MAX_LEN = 21, 40
SELECT_CAPTURE = 0.85

def select_read_lengths(cds_length_counts,
                        min_len=SELECT_MIN_LEN, max_len=SELECT_MAX_LEN,
                        capture=SELECT_CAPTURE):
    """RiboBase's read-length interval (TE_model `src/utils.py::intevl`), ported EXACTLY.

    `cds_length_counts` must be the CDS-assigned histogram, NOT all accepted alignments.
    Returns `(lengths, lo, hi, captured)`; `captured / total` is the original's `read_pct`.
    """
    counts = {n: int(cds_length_counts.get(n, 0)) for n in range(min_len, max_len + 1)}
    total = sum(counts.values())
    if total == 0:
        raise ValueError("no CDS-assigned reads in %d-%d nt: cannot select a window"
                         % (min_len, max_len))

    threshold = total * capture
    peak = max(counts.values())
    lo = hi = min(n for n in counts if counts[n] == peak)
    captured = peak

    while captured <= threshold:
        if lo == min_len and hi == max_len:
            break
        if hi < max_len and lo > min_len:
            if counts[hi + 1] >= counts[lo - 1]:
                hi += 1
                captured += counts[hi]
            else:
                lo -= 1
                captured += counts[lo]
        elif hi == max_len:
            lo -= 1
            captured += counts[lo]
        else:
            hi += 1
            captured += counts[hi]

    return list(range(lo, hi + 1)), lo, hi, captured

# both histogram builders use raw 5' ends  selection runs BEFORE offset estimation

_CDS_HEADER_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def cds_length_hist_transcriptome(reads, refs, min_len, max_len):
    """`{read_length: n}` over transcriptome primaries whose 5' end lands in RiboPy's CDS.

    Computed from the pass-1 frame (`Chromosome`, `pos5`, `length`): pass 1 keeps exactly
    the MAPQ >= 42 primaries in 20-45 nt, a superset of this histogram's filter. bowtie2
    `--norc`: every read is forward, so `pos5` is `reference_start`.
    """
    import region_lib as rl

    core_lo, core_hi = {}, {}
    for ref in refs:
        match = _CDS_HEADER_RE.search(ref)
        if match:
            core_lo[ref] = int(match.group(1)) - 1 + rl.DEFAULT_RIGHT_SPAN + 1
            core_hi[ref] = int(match.group(2)) - rl.DEFAULT_LEFT_SPAN
    frame = reads[(reads["length"] >= min_len) & (reads["length"] <= max_len)]
    if frame.empty:
        return {}
    lo = frame["Chromosome"].map(core_lo)
    hi = frame["Chromosome"].map(core_hi)
    keep = lo.notna() & (frame["pos5"] >= lo) & (frame["pos5"] < hi)
    hits = frame.loc[keep, "length"].value_counts()
    return {int(k): int(v) for k, v in hits.items()}

def cds_length_hist_genome(reads, cds_intervals, min_len, max_len):
    """`{read_length: n}` over genome reads whose 5' end projects into RiboPy's CDS.

    A read is counted ONCE even where overlapping APPRIS CDS cores hit it several times.
    """
    import numpy as np
    import pyranges as pr

    frame = reads[(reads["length"] >= min_len) & (reads["length"] <= max_len)]
    if frame.empty:
        return {}
    query = pr.PyRanges(pd.DataFrame({
        "Chromosome": frame["Chromosome"].values,
        "Start": frame["pos5"].values.astype(np.int64),
        "End": frame["pos5"].values.astype(np.int64) + 1,
        "Strand": frame["Strand"].values,
        "length": frame["length"].values,
        "read_index": np.arange(len(frame), dtype=np.int64)}))
    hit = query.join(cds_intervals, strandedness="same")
    if len(hit) == 0:
        return {}
    hits = hit.df[["read_index", "length"]].drop_duplicates(subset="read_index")
    return {int(k): int(v) for k, v in hits["length"].value_counts().items()}

def genome_cds_core_intervals(left_span=None, right_span=None):
    """RiboPy's CDS core for every selected transcript, as genomic intervals.

    `cds_len_header` = GTF CDS length + 3: the reference header counts the stop codon.
    """
    import numpy as np
    import pyranges as pr
    import region_lib as rl

    left = rl.DEFAULT_LEFT_SPAN if left_span is None else left_span
    right = rl.DEFAULT_RIGHT_SPAN if right_span is None else right_span

    exons, _base2ver = rl.build_genome_exon_table()
    cds = exons[exons["region"] == "cds"].copy()
    total = cds.groupby("transcript_id")["exon_len"].transform("sum")
    core_lo = right + 1
    core_hi = total + 3 - left

    lo = np.maximum(cds["cum_offset"], core_lo)
    hi = np.minimum(cds["cum_offset"] + cds["exon_len"], core_hi)
    keep = lo < hi
    cds, lo, hi = cds[keep], lo[keep], hi[keep]
    if cds.empty:
        return pr.PyRanges(pd.DataFrame(
            {"Chromosome": [], "Start": [], "End": [], "Strand": []}))

    offset_lo = (lo - cds["cum_offset"]).astype(np.int64)
    offset_hi = (hi - cds["cum_offset"]).astype(np.int64)
    plus = cds["Strand"].values == "+"
    start = np.where(plus, cds["Start"] + offset_lo, cds["End"] - offset_hi)
    end = np.where(plus, cds["Start"] + offset_hi, cds["End"] - offset_lo)
    return pr.PyRanges(pd.DataFrame({
        "Chromosome": cds["Chromosome"].values,
        "Start": start.astype(np.int64),
        "End": end.astype(np.int64),
        "Strand": cds["Strand"].values})).merge(strand=True)

def metagene_counts(reads, up, down):
    """`{read_length: {position: count}}` over [-up, down), from a frame with `length` and
    `rel_pos` columns. Reads whose `rel_pos` is NaN (no coding reference) are excluded."""
    window = reads[reads["rel_pos"].between(-up, down - 1)].copy()
    window["rel_int"] = window["rel_pos"].astype(int)
    grouped = window.groupby(["length", "rel_int"]).size()
    counts = defaultdict(lambda: defaultdict(int))
    for (length, position), count in grouped.items():
        counts[int(length)][int(position)] = int(count)
    return counts

def detect_offsets(reads, phase1_lengths, pre_counts, offset_fn, up, down,
                   post_window, frame0_threshold):
    """Per selected read length: P-site offset, shifted first-10-codon frame %, and the
    periodicity pass/fail. Returns `phase2`; `offset_fn` is the detector from `psite_offset.py`."""
    phase2 = {}
    for length in phase1_lengths:
        counts = {p: pre_counts[length].get(p, 0) for p in range(-up, down)}
        offset = offset_fn(counts)

        relative = reads.loc[reads["length"].eq(length) & reads["rel_pos"].notna(),
                             "rel_pos"]
        shifted = relative.astype(int) + offset
        first10 = shifted[shifted.between(0, post_window - 1)]
        n = len(first10)
        if n > 0:
            frames = first10 % 3
            n0, n1, n2 = (int((frames == f).sum()) for f in range(3))
            f0, f1, f2 = (round(x / n * 100, 1) for x in (n0, n1, n2))
        else:
            n0 = n1 = n2 = 0
            f0 = f1 = f2 = 0.0

        periodic = f0 >= frame0_threshold
        print("  %d nt: offset=%+d  frame0=%.1f%%  (%s)"
              % (length, offset, f0, "PASS" if periodic else "FAIL"), flush=True)

        phase2[length] = {
            "psite_offset": offset,
            "frame0_pct": f0, "frame1_pct": f1, "frame2_pct": f2,
            "periodic": periodic,
            "n_first10": n,
            "n_first10_frame0": n0, "n_first10_frame1": n1, "n_first10_frame2": n2,
        }
    return phase2

P2_COLUMNS = ("psite_offset", "frame0_pct", "frame1_pct", "frame2_pct", "periodic",
              "n_first10", "n_first10_frame0", "n_first10_frame1", "n_first10_frame2")

def window_qc_table(length_counts, total_reads, phase2,
                    staging_dir, sample, frame0_threshold):
    """Write `<sample>_readlen_window_qc.csv`, the route's one QC master.

    Lengths outside the window carry nulls rather than being dropped.
    """
    rows = []
    for length in sorted(length_counts):
        row = {"read_length": length,
               "n_reads": length_counts[length],
               "pct_reads": round(length_counts[length] / total_reads * 100, 2),
               "in_phase1": length in phase2}
        selected = phase2.get(length)
        for column in P2_COLUMNS:
            row[column] = selected[column] if selected else (False if column == "periodic"
                                                             else None)
        rows.append(row)

    qc = pd.DataFrame(rows)
    qc_path = os.path.join(staging_dir, "%s_readlen_window_qc.csv" % sample)
    qc.to_csv(qc_path, index=False)
    print("  Saved: %s" % qc_path, flush=True)
    periodic = sorted(qc.loc[qc["periodic"].eq(True), "read_length"].tolist())
    print("  Periodic lengths (frame0 >= %.0f%%): %s" % (frame0_threshold, periodic),
          flush=True)
    return qc
