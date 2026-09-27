#!/usr/bin/env python3
"""Per-sample Ribo-seq QC in one BAM traversal (genome or transcriptome route).

One pass collects the per-read state (reference, 5' end, strand, length); from it come
the read-length selection and P-site offsets (`<sample>_readlen_window_qc.csv`) and then,
applying the just-computed offsets to the same in-memory reads, the whole-CDS frame
table (`<sample>_cds_psite_frame.csv`).
"""
import os
import re
import sys
import argparse

from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
import bam_inputs as fc          # the one uniqueness policy  fc.is_unique_{genome,txome}_read
from psite_offset import ribotish_get_offset, get_offset_periodicity

import qc_core

import pysam
import pyranges as pr
import numpy as np
import pandas as pd
p = argparse.ArgumentParser(description="Per-sample read-length selection.")
p.add_argument("--sample",           required=True)
p.add_argument("--bam",              required=True)
p.add_argument("--route", choices=["genome", "transcriptome"], default="genome")
p.add_argument("--gtf", default=None)
p.add_argument("--appris", default=None)
p.add_argument("--out", default=None)
args = p.parse_args()

#: offset frame from downstream 3 nt phasing  robust to bimodal start peaks  see psite_offset.py
_OFFSET_FN = get_offset_periodicity

SAMPLE      = args.sample
BAM         = args.bam
TX          = args.route == "transcriptome"
OUT         = args.out or (config.tx_out_dir() if TX else config.out_dir())
F0_THRESH   = 50.0   # min frame0 % after P-site shift to keep a length
TAG         = " (transcriptome)" if TX else ""

MIN_LEN, MAX_LEN = config.MIN_LEN, config.MAX_LEN

PRE_WIN_UP  = 50   # nt upstream of start codon
PRE_WIN_DN  = 30   # nt downstream of start codon
POST_WIN_DN = 30

dir_staging = os.path.join(OUT, "tables", "_staging")
os.makedirs(dir_staging, exist_ok=True)

print(f"=== [sample_qc{TAG}] sample={SAMPLE} ===", flush=True)

_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def _cds_start0_from_refname(ref):
    m = _CDS_RE.search(ref)
    return int(m.group(1)) - 1 if m else None

def _cds_bounds_from_refname(ref):
    """(cds_start0, cds_len_nostop) or (None, None). Stop codon trimmed (-3)."""
    m = _CDS_RE.search(ref)
    if not m:
        return None, None
    start1, end1 = int(m.group(1)), int(m.group(2))
    cds_len_nostop = (end1 - start1 + 1) - 3  # drop the stop codon  CDS:end includes it
    if cds_len_nostop <= 0:
        return None, None
    return start1 - 1, cds_len_nostop

print("Reading BAM...", flush=True)
rec = {"Chromosome": [], "pos5": [], "Strand": [], "length": []}
bam = pysam.AlignmentFile(BAM, "rb")
if TX:
    # reference names embed "|CDS:start-end|"  bowtie2 --norc so the 5' end is reference_start
    tx_refs = list(bam.references)
    ref_cds0 = {ref: _cds_start0_from_refname(ref) for ref in tx_refs}
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

# CDS cores tested on raw 5' ends  no P-site shift
print("Phase 1: CDS length distribution + 85 % expansion...", flush=True)
if TX:
    cds_length_counts = qc_core.cds_length_hist_transcriptome(
        reads_df, tx_refs, qc_core.SELECT_MIN_LEN, qc_core.SELECT_MAX_LEN)
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

# rel_pos of each 5' end to its CDS start  from the reference name on transcriptome or by
# joining to start codon windows that extend upstream where the P-site signal is on genome
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

# ── per length P-site offset and frame % ──────────────────────────────────
print("Phase 2: P-site detection and frame %...", flush=True)
phase2 = qc_core.detect_offsets(
    reads_df, phase1_lengths, pre_counts, _OFFSET_FN,
    PRE_WIN_UP, PRE_WIN_DN, POST_WIN_DN, F0_THRESH)

qc_df = qc_core.window_qc_table(
    length_counts, total_reads, phase2, dir_staging, SAMPLE, F0_THRESH)

# ── P-site frame counts across whole CDS from the same in memory reads ──────────
print(f"\n=== [cds_frame{TAG}] sample={SAMPLE} ===", flush=True)

# inputs derived from qc_df with the exact expressions former step 03 used on the
# reread CSV  ints and bools survive the round trip identically
phase1_mask   = fc._as_bool(qc_df["in_phase1"])
periodic_mask = fc._as_bool(qc_df["periodic"])

phase1_rows      = qc_df[phase1_mask]
frame_lengths    = set(phase1_rows["read_length"].astype(int).tolist())
periodic_lengths = set(qc_df.loc[periodic_mask, "read_length"].astype(int).tolist())
psite_offsets    = dict(zip(phase1_rows["read_length"].astype(int),
                            phase1_rows["psite_offset"].astype(int)))

print(f"  Phase1 lengths:   {sorted(frame_lengths)}")
print(f"  Periodic lengths: {sorted(periodic_lengths)}")

# same read subset in same BAM order that 03's own pass loaded
sub = reads_df[reads_df["length"].isin(frame_lengths)].reset_index(drop=True)
offsets_s = sub["length"].map(psite_offsets)
n_loaded = len(sub)
print(f"  {n_loaded:,} phase1-length reads loaded")

if TX:
    ref_bounds = {ref: _cds_bounds_from_refname(ref) for ref in tx_refs}
    cds_start0 = sub["Chromosome"].map(
        {r: b[0] for r, b in ref_bounds.items() if b[0] is not None})
    cds_len_nostop = sub["Chromosome"].map(
        {r: b[1] for r, b in ref_bounds.items() if b[0] is not None})
    # bowtie2 --norc  all reads forward on the transcript  5' end = reference_start
    rel = sub["pos5"] + offsets_s - cds_start0
    keep = cds_start0.notna() & (rel >= 0) & (rel < cds_len_nostop)
    frame_reads = pd.DataFrame({
        "length":  sub.loc[keep, "length"].values,
        "frame":   (rel[keep] % 3).values,
        "rel_pos": rel[keep].values,
    })
    print(f"  {len(frame_reads):,} reads with P-site in CDS body "
          f"(0 <= rel_pos < CDS len, stop excluded)")
else:
    plus = sub["Strand"].values == "+"
    p_sites = np.where(plus, sub["pos5"].values + offsets_s.values,
                       sub["pos5"].values - offsets_s.values)

    print("Joining P-site positions to CDS exons...", flush=True)
    frame_reads = pd.DataFrame({
        "length":  sub["length"].values,
        "frame":   np.nan,
        "rel_pos": np.nan,
    })

    if n_loaded > 0 and len(ann) > 0:
        read_idx = np.arange(n_loaded)
        psite_pr = pr.PyRanges(pd.DataFrame({
            "Chromosome": sub["Chromosome"].values,
            "Start":      p_sites,
            "End":        p_sites + 1,
            "Strand":     sub["Strand"].values,
            "read_idx":   read_idx,
        }))
        cds_pr = pr.PyRanges(
            ann[["Chromosome", "Start", "End", "Strand", "Phase", "cds_genomic_start"]]
        )
        joined = psite_pr.join(cds_pr, strandedness="same")

        if not joined.df.empty:
            jdf  = joined.df.copy()
            plus = jdf["Strand"] == "+"

            # subtract Phase  using +Phase mislabels in frame P-sites in phase 1 and 2 exons
            jdf["frame"] = np.where(
                plus,
                (jdf["Start"] - jdf["Start_b"] - jdf["Phase"]) % 3,
                (jdf["End_b"] - 1 - jdf["Start"] - jdf["Phase"]) % 3,
            ).astype(float)

            jdf["rel_pos"] = np.where(
                plus,
                jdf["Start"] - jdf["cds_genomic_start"],
                jdf["cds_genomic_start"] - jdf["Start"],
            ).astype(float)

            jdf = jdf.drop_duplicates("read_idx", keep="first")
            frame_reads.loc[jdf["read_idx"].values, "frame"]   = jdf["frame"].values
            frame_reads.loc[jdf["read_idx"].values, "rel_pos"] = jdf["rel_pos"].values

            n_in_cds = int((frame_reads["rel_pos"] >= 0).sum())
            print(f"  {n_in_cds:,} reads with P-site in CDS body (rel_pos >= 0)")
        else:
            print("  No P-site positions overlapped any CDS exon.")

print("Aggregating frame counts per read length...", flush=True)
out_rows = []
for rlen in sorted(frame_lengths):
    mask = (
        frame_reads["length"].eq(rlen)
        & frame_reads["rel_pos"].notna()
        & (frame_reads["rel_pos"] >= 0)
    )
    frames   = frame_reads.loc[mask, "frame"]
    n_psite  = len(frames)
    n_f0 = int((frames == 0).sum())
    n_f1 = int((frames == 1).sum())
    n_f2 = int((frames == 2).sum())
    pct_f0 = round(n_f0 / n_psite * 100, 1) if n_psite else 0.0
    pct_f1 = round(n_f1 / n_psite * 100, 1) if n_psite else 0.0
    pct_f2 = round(n_f2 / n_psite * 100, 1) if n_psite else 0.0
    out_rows.append({
        "read_length":   rlen,
        "in_phase1":     True,
        "periodic":      rlen in periodic_lengths,
        "psite_offset":  psite_offsets[rlen],
        "n_psite_in_cds": n_psite,
        "n_frame0":      n_f0,
        "n_frame1":      n_f1,
        "n_frame2":      n_f2,
        "pct_frame0":    pct_f0,
        "pct_frame1":    pct_f1,
        "pct_frame2":    pct_f2,
    })

frame_df = pd.DataFrame(out_rows)

frame_path = os.path.join(dir_staging, f"{SAMPLE}_cds_psite_frame.csv")
frame_df.to_csv(frame_path, index=False)

total_in_cds = frame_df["n_psite_in_cds"].sum()
total_f0     = frame_df["n_frame0"].sum()
pct_f0_total = round(total_f0 / total_in_cds * 100, 1) if total_in_cds else 0.0
print(f"  P-sites in CDS body: {total_in_cds:,} | frame 0: {total_f0:,} ({pct_f0_total:.1f}%)")
print(f"  Saved: {frame_path}")
print("Done.", flush=True)
