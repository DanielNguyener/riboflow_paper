#!/usr/bin/env python3
"""STEP 01 — Read-length selection + per-read-length P-site offset detection (genome or transcriptome route)."""
import os
import re
import sys
import argparse

from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
import bam_inputs as fc          # the one uniqueness policy: fc.is_unique_{genome,txome}_read
from psite_offset import ribotish_get_offset, get_offset_periodicity

import qc_core

import pysam
import pyranges as pr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # backend for qc_core's lazy pyplot imports (--plots)

p = argparse.ArgumentParser(description="Per-sample read-length selection.")
p.add_argument("--sample",           required=True)
p.add_argument("--bam",              required=True)
p.add_argument("--route",            choices=["genome", "transcriptome"], default="genome",
               help="Which alignments --bam holds: genome (NH == 1 kept) or transcriptome "
                    "(MAPQ >= 42 kept; CDS bounds come from the reference names).")
p.add_argument("--gtf",              default=None, help="genome route only")
p.add_argument("--appris",           default=None, help="genome route only")
p.add_argument("--out",              default=None,
               help="Output root; default config.out_dir() / config.tx_out_dir() by route.")
p.add_argument("--plots", action="store_true",
               help="also write the pre- and post-shift metagene PDFs. Off by default: they\n"
                    "are diagnostics, and the read-length window and offsets they illustrate\n"
                    "are already in the tables this step writes.")
args = p.parse_args()

#: Offset frame from downstream 3-nt phasing (robust to bimodal start peaks); see psite_offset.py.
_OFFSET_FN = get_offset_periodicity

SAMPLE      = args.sample
BAM         = args.bam
TX          = args.route == "transcriptome"
OUT         = args.out or (config.tx_out_dir() if TX else config.out_dir())
F0_THRESH   = 50.0   # min frame0 % after the P-site shift to keep a length
TAG         = " (transcriptome)" if TX else ""

MIN_LEN, MAX_LEN = config.MIN_LEN, config.MAX_LEN
FRAME_COLORS     = config.FRAME_COLORS

PRE_WIN_UP  = 50   # nt upstream of start codon
PRE_WIN_DN  = 30   # nt downstream of start codon
POST_WIN_DN = 30

dir_plots   = os.path.join(OUT, "plots", "metagene")
dir_staging = os.path.join(OUT, "tables", "_staging")
for d in ([dir_plots, dir_staging] if args.plots else [dir_staging]):
    os.makedirs(d, exist_ok=True)

print(f"=== [01 readlen_selection{TAG}] sample={SAMPLE} ===", flush=True)

_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def _cds_start0_from_refname(ref):
    m = _CDS_RE.search(ref)
    return int(m.group(1)) - 1 if m else None

print("Reading BAM...", flush=True)
rec = {"Chromosome": [], "pos5": [], "Strand": [], "length": []}
bam = pysam.AlignmentFile(BAM, "rb")
if TX:
    # Reference names embed "|CDS:start-end|"; bowtie2 --norc, so the 5' end is reference_start.
    ref_cds0 = {ref: _cds_start0_from_refname(ref) for ref in bam.references}
    print(f"  {len(ref_cds0):,} references; "
          f"{sum(1 for v in ref_cds0.values() if v is not None):,} carry a CDS region", flush=True)
for read in bam.fetch(until_eof=TX):
    if TX:
        if not fc.is_unique_txome_read(read):     # MAPQ >= 42
            continue
    elif (read.is_unmapped or read.is_secondary or read.is_supplementary
          or not fc.is_unique_genome_read(read)):   # NH == 1
        continue
    rlen = read.query_length
    if not (MIN_LEN <= rlen <= MAX_LEN):
        continue
    if TX or not read.is_reverse:
        p5, strand = read.reference_start, "+"
    else:
        p5, strand = read.reference_end - 1, "-"
    rec["Chromosome"].append(read.reference_name)
    rec["pos5"].append(p5)
    rec["Strand"].append(strand)
    rec["length"].append(rlen)
bam.close()
reads_df    = pd.DataFrame(rec)
total_reads = len(reads_df)
length_counts = Counter(reads_df["length"].tolist())
print(f"  {total_reads:,} reads in length range {MIN_LEN}-{MAX_LEN}", flush=True)

if total_reads == 0:
    print("  No reads — aborting.", flush=True)
    sys.exit(1)

# CDS cores tested on raw 5' ends, no P-site shift.
print("Phase 1: CDS length distribution + 85 % expansion...", flush=True)
if TX:
    with pysam.AlignmentFile(BAM, "rb") as _cds_bam:
        cds_length_counts = qc_core.cds_length_hist_transcriptome(
            _cds_bam, qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN)
else:
    cds_length_counts = qc_core.cds_length_hist_genome(
        reads_df, qc_core.genome_cds_core_intervals(),
        qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN)
n_cds = sum(cds_length_counts.values())
print(f"  {n_cds:,} CDS-assigned reads in "
      f"{qc_core.SELECT_MIN_LEN}-{qc_core.SELECT_MAX_LEN} nt "
      f"({n_cds / total_reads * 100:.1f}% of accepted)", flush=True)
phase1_lengths, lo, hi, captured = qc_core.select_read_lengths(cds_length_counts)
print(f"  Peak window: {lo}-{hi} nt | captured {captured / n_cds * 100:.1f}% of CDS reads "
      f"| lengths: {phase1_lengths}", flush=True)

# rel_pos of each 5' end to its CDS start: from the reference name (transcriptome) or by
# joining to start-codon windows that extend upstream, where the P-site signal is (genome).
if TX:
    print("Computing rel_pos from transcript CDS starts...", flush=True)
    reads_df["rel_pos"] = (reads_df["pos5"] - reads_df["Chromosome"].map(ref_cds0)).astype(float)
else:
    print("Loading CDS annotation cache...", flush=True)
    ann = config.load_annotation(args.gtf, args.appris)
    tx_starts = (ann.groupby("transcript_id", sort=False)
                 .agg(Chromosome=("Chromosome", "first"),
                      Strand=("Strand", "first"),
                      cds_genomic_start=("cds_genomic_start", "first"))
                 .reset_index(drop=True))
    print(f"  {len(tx_starts):,} transcripts", flush=True)

    tx_p = tx_starts[tx_starts["Strand"] == "+"].copy()
    tx_m = tx_starts[tx_starts["Strand"] == "-"].copy()

    tx_p["Start"] = (tx_p["cds_genomic_start"] - PRE_WIN_UP).clip(lower=0)
    tx_p["End"]   =  tx_p["cds_genomic_start"] + PRE_WIN_DN

    tx_m["Start"] = (tx_m["cds_genomic_start"] - PRE_WIN_DN + 1).clip(lower=0)
    tx_m["End"]   =  tx_m["cds_genomic_start"] + PRE_WIN_UP + 1

    windows_df = pd.concat([tx_p, tx_m], ignore_index=True)

    print("Joining reads to CDS start windows...", flush=True)
    pos5_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": reads_df["Chromosome"].values,
        "Start":      reads_df["pos5"].values,
        "End":        reads_df["pos5"].values + 1,
        "Strand":     reads_df["Strand"].values,
        "read_idx":   reads_df.index.values,
    }))
    win_pr = pr.PyRanges(
        windows_df[["Chromosome", "Start", "End", "Strand", "cds_genomic_start"]]
    )
    joined = pos5_pr.join(win_pr, strandedness="same")

    reads_df["rel_pos"] = np.nan
    if not joined.df.empty:
        jdf  = joined.df.copy()
        plus = jdf["Strand"] == "+"
        jdf["rel_pos"] = np.where(
            plus,
            jdf["Start"] - jdf["cds_genomic_start"],
            jdf["cds_genomic_start"] - jdf["Start"],
        ).astype(float)
        jdf = jdf.drop_duplicates("read_idx", keep="first")
        reads_df.loc[jdf["read_idx"].values, "rel_pos"] = jdf["rel_pos"].values

print(f"  {reads_df['rel_pos'].notna().sum():,} reads mapped to CDS starts", flush=True)

pre_counts = qc_core.metagene_counts(reads_df, PRE_WIN_UP, PRE_WIN_DN)

# ── 7. Per-length P-site offset and frame % ──────────────────────────────────
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
                          f"{SAMPLE}{TAG} - 5' end metagene (unshifted; red = detected P-site offset)")
    qc_core.plot_postshift(reads_df, phase1_lengths, phase2, length_counts, POST_WIN_DN,
                           dir_plots, SAMPLE,
                           f"{SAMPLE}{TAG} - 5' end metagene (P-site shifted, first 10 codons)  "
                           f"[threshold={F0_THRESH:.0f}%]")
print("Done.", flush=True)
