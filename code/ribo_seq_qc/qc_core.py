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

FRAME_COLORS = {0: "#2c7fb8", 1: "#7fcdbb", 2: "#edf8b1"}

def _fc():
    """`bam_inputs`, imported lazily so this module stays importable without pysam."""
    import bam_inputs
    return bam_inputs

def require(*packages):
    """Fail with an actionable message when a declared dependency is missing.

    Never installs anything: a script that pip-installs into whatever interpreter happens
    to be running can silently change the versions a published number was produced with.
    """
    missing = []
    for package in packages:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    if missing:
        raise SystemExit(
            "missing required package(s): %s\n"
            "Install the declared environment first:\n"
            "    python -m pip install -r requirements.txt"
            % ", ".join(sorted(missing)))

SELECT_MIN_LEN, SELECT_MAX_LEN = 21, 40
SELECT_CAPTURE = 0.85

def select_read_lengths(cds_length_counts,
                        min_len=SELECT_MIN_LEN, max_len=SELECT_MAX_LEN,
                        capture=SELECT_CAPTURE):
    """RiboBase's read-length interval: TE_model `src/utils.py::intevl`, ported.

    `cds_length_counts` is the CDS-assigned RPF length histogram -- RiboPy's
    `ribo_object.get_length_dist("CDS")`, NOT every accepted alignment. That distinction is
    the whole point: the two distributions have different modes and different tails, so
    selecting on all reads reproduces RiboBase's interval for only 15 of the 24 published
    libraries.

    Four details of the original that a paraphrase gets wrong, all load-bearing:

      * the loop runs `while value <= pct_85`, so a window capturing EXACTLY 85 % expands
        once more;
      * `pct_85` is 85 % of the CDS total WITHIN `min_len..max_len`, not of the library;
      * ties between the two candidate neighbours go to the LONGER length (`>=` on
        `count[mmax + 1]` vs `count[mmin - 1]`);
      * a tie for the mode resolves to the SHORTEST tied length, because the original takes
        `.values[0]` of a frame ordered by ascending `read_length`.

    Returns `(lengths, lo, hi, captured)` -- `lengths` is the complete inclusive interval,
    `captured` the count inside it. `captured / total` is the original's `read_pct`.
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

# Neither uses P-sites: selection runs BEFORE offset estimation, so the histogram is built
# from raw 5' ends.

_CDS_HEADER_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def cds_length_hist_transcriptome(bam, min_len, max_len):
    """`{read_length: n}` over transcriptome primaries whose 5' end lands in RiboPy's CDS.

    The transcript coordinate IS the alignment coordinate here, and bowtie2 `--norc` puts
    every read on the forward strand, so the 5' end is `reference_start`. One primary
    alignment per read, so a read is counted at most once.
    """
    import region_lib as rl

    bounds = {}
    for ref in bam.references:
        match = _CDS_HEADER_RE.search(ref)
        if match:
            bounds[ref] = (int(match.group(1)) - 1, int(match.group(2)))
    counts = defaultdict(int)
    for read in bam.fetch(until_eof=True):
        if not _fc().is_unique_txome_read(read):      # MAPQ >= 42
            continue
        length = read.query_length
        if not (min_len <= length <= max_len):
            continue
        bound = bounds.get(read.reference_name)
        if bound is None:
            continue
        start_site, stop_site = bound
        if (start_site + rl.DEFAULT_RIGHT_SPAN + 1
                <= read.reference_start
                < stop_site - rl.DEFAULT_LEFT_SPAN):
            counts[length] += 1
    return dict(counts)

def cds_length_hist_genome(reads, cds_intervals, min_len, max_len):
    """`{read_length: n}` over genome reads whose 5' end projects into RiboPy's CDS.

    `reads` is the `(Chromosome, pos5, Strand, length)` frame the genome step already
    builds; `cds_intervals` is `genome_cds_core_intervals()` below -- the RiboPy CDS of
    every selected transcript, expressed back in GENOMIC coordinates so no per-read
    projection is needed.

    A read is counted ONCE even where APPRIS isoforms overlap and its 5' end falls in
    several transcripts' CDS cores. RiboPy counts one alignment once; a per-transcript
    tally would inflate exactly the length bins that sit under dense annotation.
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

    Walks each transcript's CDS exons 5'->3', keeps the sub-interval whose CDS-relative
    offset lies in `[right_span + 1, cds_len_header - left_span)`, and maps it back to the
    genome. `cds_len_header` is the GTF CDS length PLUS 3: the transcriptome reference
    header counts the stop codon inside its CDS and the GTF does not, and that 3 nt shifts
    the downstream boundary. The resulting core is strictly inside the GTF CDS, so UTR
    exons never contribute and are not consulted.
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
    """Per selected read length: the P-site offset, the shifted first-10-codon frame
    percentages, and whether it clears the periodicity threshold.

    Returns `phase2`. `offset_fn` is the detector from `psite_offset.py`.
    """
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

    Every observed read length gets a row; lengths outside the window carry nulls rather
    than being dropped, so the table records what was seen as well as what was selected.
    The selected lengths and their P-site offsets are columns here, which is why there is
    no second `psite_shifts` table: it held the same two numbers and nothing read it.
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

def _panels(n_panels):
    import matplotlib.pyplot as plt
    ncols = 3
    nrows = max(1, (n_panels + ncols - 1) // ncols)
    figure, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows),
                               sharey=False)
    flat = np.array(axes).flat if n_panels > 1 else [axes]
    return figure, list(flat)

def plot_preshift(pre_counts, phase1_lengths, phase2, up, down, plots_dir, sample,
                  title):
    """The unshifted 5'-end metagene, one panel per selected length."""
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    positions = list(range(-up, down))
    figure, axes = _panels(len(phase1_lengths))
    for axis, length in zip(axes, phase1_lengths):
        raw = np.array([pre_counts[length].get(p, 0) for p in positions])
        # otherwise the CDS body dwarfs the upstream P-site peak.
        ys = raw / (raw.max() if raw.max() > 0 else 1)
        axis.axvspan(-up, 0, color="lightgrey", alpha=0.35, zorder=0)
        axis.axvline(0, color="black", linewidth=1.0, linestyle="--", zorder=2)
        axis.bar(positions, ys, width=1.0,
                 color=[FRAME_COLORS[p % 3] for p in positions], linewidth=0, zorder=1)
        if length in phase2 and phase2[length]["psite_offset"] is not None:
            offset = phase2[length]["psite_offset"]
            axis.axvline(-offset, color="red", linewidth=1.2, linestyle=":", zorder=3,
                         label="offset=%+d" % offset)
            axis.legend(fontsize=7, frameon=False)
        axis.set_title("%d nt" % length, fontsize=11)
        axis.set_xlabel("Position relative to CDS start (nt)", fontsize=8)
        axis.set_ylabel("Normalised coverage", fontsize=8)
        axis.set_xticks(range(-up, down + 1, 10))
        axis.set_ylim(0, 1.05)
        axis.tick_params(labelsize=7)
        axis.set_xlim(-up - 0.5, down - 0.5)
    for axis in axes[len(phase1_lengths):]:
        axis.set_visible(False)

    patches = [mpatches.Patch(color=FRAME_COLORS[f], label="Frame %d" % f)
               for f in range(3)]
    figure.legend(handles=patches, loc="lower right", fontsize=9, frameon=False)
    figure.suptitle(title, fontsize=12, y=1.01)
    figure.tight_layout()
    path = os.path.join(plots_dir, "%s_preshift.pdf" % sample)
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)
    print("  Saved: %s" % path, flush=True)

def plot_postshift(reads, phase1_lengths, phase2, length_counts, post_window,
                   plots_dir, sample, title):
    """The P-site-shifted metagene over the first `post_window` nt, in RPM."""
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    positions = list(range(post_window))
    post_counts = {}
    for length in phase1_lengths:
        selected = reads["length"].eq(length) & reads["rel_pos"].notna()
        shifted = reads.loc[selected, "rel_pos"].astype(int) + phase2[length]["psite_offset"]
        first10 = shifted[shifted.between(0, post_window - 1)]
        post_counts[length] = first10.groupby(first10).size().to_dict()

    figure, axes = _panels(len(phase1_lengths))
    for axis, length in zip(axes, phase1_lengths):
        n_reads = length_counts[length]
        rpm = 1e6 / n_reads if n_reads else 1
        summary = phase2[length]
        ys = np.array([post_counts[length].get(p, 0) * rpm for p in positions])
        axis.bar(positions, ys, width=1.0,
                 color=[FRAME_COLORS[p % 3] for p in positions], linewidth=0, zorder=1)
        for codon_start in range(0, post_window, 3):
            axis.axvline(codon_start - 0.5, color="grey", linewidth=0.4,
                         linestyle=":", zorder=0)
        axis.set_title("%d nt  %s  frame0=%.0f%%"
                       % (length, "PASS" if summary["periodic"] else "FAIL",
                          summary["frame0_pct"]), fontsize=10)
        axis.set_xlabel("P-site shifted position (nt)", fontsize=8)
        axis.set_ylabel("RPM", fontsize=8)
        axis.tick_params(labelsize=7)
        axis.set_xlim(-0.5, post_window - 0.5)
        axis.set_xticks(range(0, post_window, 3))
    for axis in axes[len(phase1_lengths):]:
        axis.set_visible(False)

    patches = [mpatches.Patch(color=FRAME_COLORS[f], label="Frame %d" % f)
               for f in range(3)]
    figure.legend(handles=patches, loc="lower right", fontsize=9, frameon=False)
    figure.suptitle(title, fontsize=12, y=1.01)
    figure.tight_layout()
    path = os.path.join(plots_dir, "%s_postshift.pdf" % sample)
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)
    print("  Saved: %s" % path, flush=True)

def cds_frame_table(per_length, qc_rows, staging_dir, sample):
    """Write `<sample>_cds_psite_frame.csv` from `{length: (n_total, n0, n1, n2)}`.

    `qc_rows` carries `in_phase1`, `periodic` and `psite_offset` per length, straight from
    step 01/01t, so the frame table records which window the counts came from.
    """
    rows = []
    for length in sorted(per_length):
        total, n0, n1, n2 = per_length[length]
        info = qc_rows.get(length, {})
        rows.append({
            "read_length": length,
            "in_phase1": info.get("in_phase1", True),
            "periodic": info.get("periodic", False),
            "psite_offset": info.get("psite_offset"),
            "n_psite_in_cds": total,
            "n_frame0": n0, "n_frame1": n1, "n_frame2": n2,
            "pct_frame0": round(n0 / total * 100, 2) if total else 0.0,
            "pct_frame1": round(n1 / total * 100, 2) if total else 0.0,
            "pct_frame2": round(n2 / total * 100, 2) if total else 0.0,
        })
    frame = pd.DataFrame(rows)
    path = os.path.join(staging_dir, "%s_cds_psite_frame.csv" % sample)
    frame.to_csv(path, index=False)
    print("  Saved: %s" % path, flush=True)
    return frame
