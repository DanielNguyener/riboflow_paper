#!/usr/bin/env python3
"""STEP 01t — Read-length selection + per-read-length P-site offset detection, run on TRANSCRIPTOME alignments (the genome-step twin)."""
import re

import os
import sys
import argparse

from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
import bam_inputs as fc          # the one uniqueness policy: fc.is_unique_txome_read
from psite_offset import ribotish_get_offset, get_offset_periodicity

import qc_core
qc_core.require("pysam", "matplotlib")

import pysam
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # backend for qc_core's lazy pyplot imports (--plots)

p = argparse.ArgumentParser(description="Per-sample read-length selection (transcriptome).")
p.add_argument("--sample",           required=True)
p.add_argument("--bam",              required=True)
p.add_argument("--out",              default=config.tx_out_dir())
p.add_argument("--frame0-threshold", type=float, default=50.0,
               help="Min frame0 %% after P-site shift to keep a length (default 50).")
p.add_argument("--offset-method", choices=["periodicity", "argmax"], default="periodicity",
               help="P-site offset frame selection: 'periodicity' (default; frame from "
                    "downstream 3-nt phasing, robust to bimodal start peaks) or 'argmax' "
                    "(ribotish's single-peak frame). Must match the genome step.")
p.add_argument("--plots", action="store_true",
               help="also write the pre- and post-shift metagene PDFs. Off by default: they\n"
                    "are diagnostics, and the read-length window and offsets they illustrate\n"
                    "are already in the tables this step writes.")
args = p.parse_args()

_OFFSET_FN = get_offset_periodicity if args.offset_method == "periodicity" else ribotish_get_offset

SAMPLE      = args.sample
BAM         = args.bam
OUT         = args.out
F0_THRESH   = args.frame0_threshold

MIN_LEN, MAX_LEN = config.MIN_LEN, config.MAX_LEN
FRAME_COLORS     = config.FRAME_COLORS

PRE_WIN_UP  = 50
PRE_WIN_DN  = 30
POST_WIN_DN = 30

dir_plots   = os.path.join(OUT, "plots", "metagene")
dir_staging = os.path.join(OUT, "tables", "_staging")
for d in ([dir_plots, dir_staging] if args.plots else [dir_staging]):
    os.makedirs(d, exist_ok=True)

print(f"=== [01t readlen_selection / transcriptome] sample={SAMPLE} ===", flush=True)

_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def _cds_start0_from_refname(ref):
    m = _CDS_RE.search(ref)
    if not m:
        return None
    return int(m.group(1)) - 1

print("Reading transcriptome BAM...", flush=True)
bam = pysam.AlignmentFile(BAM, "rb")

ref_cds0 = {ref: _cds_start0_from_refname(ref) for ref in bam.references}
n_refs_with_cds = sum(1 for v in ref_cds0.values() if v is not None)
print(f"  {len(ref_cds0):,} references; {n_refs_with_cds:,} carry a CDS region",
      flush=True)

rec = {"ref": [], "pos5": [], "length": []}
for read in bam.fetch(until_eof=True):
    if not fc.is_unique_txome_read(read):         # MAPQ >= 42
        continue
    rlen = read.query_length
    if not (MIN_LEN <= rlen <= MAX_LEN):
        continue
    # bowtie2 --norc: all reads forward on the transcript, 5' end = reference_start
    rec["ref"].append(read.reference_name)
    rec["pos5"].append(read.reference_start)
    rec["length"].append(rlen)
bam.close()

reads_df    = pd.DataFrame(rec)
total_reads = len(reads_df)
length_counts = Counter(reads_df["length"].tolist())
print(f"  {total_reads:,} reads in length range {MIN_LEN}-{MAX_LEN}", flush=True)

if total_reads == 0:
    print("  No reads — aborting.", flush=True)
    sys.exit(1)

print("Phase 1: CDS length distribution + 85 % expansion...", flush=True)
with pysam.AlignmentFile(BAM, "rb") as _cds_bam:
    cds_length_counts = qc_core.cds_length_hist_transcriptome(
        _cds_bam, qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN)
n_cds = sum(cds_length_counts.values())
print(f"  {n_cds:,} CDS-assigned reads in "
      f"{qc_core.SELECT_MIN_LEN}-{qc_core.SELECT_MAX_LEN} nt "
      f"({n_cds / total_reads * 100:.1f}% of accepted)", flush=True)
phase1_lengths, lo, hi, captured = qc_core.select_read_lengths(cds_length_counts)
print(f"  Peak window: {lo}-{hi} nt | captured {captured / n_cds * 100:.1f}% of CDS reads "
      f"| lengths: {phase1_lengths}", flush=True)

print("Computing rel_pos from transcript CDS starts...", flush=True)
cds0 = reads_df["ref"].map(ref_cds0)
reads_df["rel_pos"] = (reads_df["pos5"] - cds0).astype(float)
print(f"  {reads_df['rel_pos'].notna().sum():,} reads on coding refs", flush=True)

pre_counts = qc_core.metagene_counts(reads_df, PRE_WIN_UP, PRE_WIN_DN)

# ── 5. Per-length P-site offset and frame % ──────────────────────────────────
print("Phase 2: P-site detection and frame %...", flush=True)
phase2 = qc_core.detect_offsets(
    reads_df, phase1_lengths, pre_counts, _OFFSET_FN,
    PRE_WIN_UP, PRE_WIN_DN, POST_WIN_DN, F0_THRESH)

qc_df = qc_core.window_qc_table(
    length_counts, total_reads, phase2, dir_staging, SAMPLE, F0_THRESH)

if args.plots:
    print("Plotting metagenes...", flush=True)
    qc_core.plot_preshift(pre_counts, phase1_lengths, phase2, PRE_WIN_UP, PRE_WIN_DN,
                          dir_plots, SAMPLE,
                          f"{SAMPLE} (transcriptome) - 5' end metagene "
                          f"(unshifted; red = detected P-site offset)")
    qc_core.plot_postshift(reads_df, phase1_lengths, phase2, length_counts, POST_WIN_DN,
                           dir_plots, SAMPLE,
                           f"{SAMPLE} (transcriptome) - 5' end metagene "
                           f"(P-site shifted, first 10 codons)  [threshold={F0_THRESH:.0f}%]")
print("Done.", flush=True)
