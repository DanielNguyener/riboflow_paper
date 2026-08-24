#!/usr/bin/env python3
"""STEP 05 — P-site frame counts across APPRIS CDS (phase1 lengths)."""

import sys
import os
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
import bam_inputs as fc          # the one uniqueness policy: fc.is_unique_genome_read

import qc_core
qc_core.require("pysam", "pyranges")

import pysam
import pyranges as pr
import pandas as pd
import numpy as np

p = argparse.ArgumentParser(description="P-site frame counts across APPRIS CDS.")
p.add_argument("--sample", required=True)
p.add_argument("--bam",    required=True)
p.add_argument("--gtf",    default=config.gtf_path())
p.add_argument("--appris", default=config.appris_path())
p.add_argument("--out",    default=config.out_dir())
args = p.parse_args()

SAMPLE = args.sample
BAM    = args.bam
GTF    = args.gtf
APPRIS = args.appris
OUT    = args.out

dir_tables  = os.path.join(OUT, "tables")
dir_staging = os.path.join(dir_tables, "_staging")
for d in [dir_tables, dir_staging]:
    os.makedirs(d, exist_ok=True)

QC_CSV = os.path.join(dir_staging, f"{SAMPLE}_readlen_window_qc.csv")
print(f"=== [05 cds_frame] sample={SAMPLE} ===", flush=True)

# ── 1. Load phase1 lengths, P-site offsets ────────────────────────────────────
print("Loading phase1 lengths and P-site offsets from step 01...", flush=True)
qc_df = pd.read_csv(QC_CSV)
_as_bool = lambda s: s.map(lambda x: str(x).strip().lower() in ("true", "1"))

phase1_mask  = _as_bool(qc_df["in_phase1"])
periodic_mask = _as_bool(qc_df["periodic"])

phase1_rows = qc_df[phase1_mask]
phase1_lengths   = set(phase1_rows["read_length"].astype(int).tolist())
periodic_lengths = set(qc_df.loc[periodic_mask, "read_length"].astype(int).tolist())
psite_offsets    = dict(zip(phase1_rows["read_length"].astype(int),
                            phase1_rows["psite_offset"].astype(int)))

print(f"  Phase1 lengths:   {sorted(phase1_lengths)}")
print(f"  Periodic lengths: {sorted(periodic_lengths)}")

print("Loading CDS annotation cache...", flush=True)
cds_df = config.load_annotation()
print(f"  {len(cds_df):,} CDS exon records, "
      f"{cds_df['transcript_id'].nunique():,} transcripts")

# ── 4. Load BAM reads (phase1 lengths), compute P-site positions ──────────────
print("Reading BAM (phase1 lengths only)...", flush=True)
chrom_col  = []
psite_col  = []
strand_col = []
length_col = []

bam = pysam.AlignmentFile(BAM, "rb")
for read in bam.fetch():
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    if not fc.is_unique_genome_read(read):        # NH == 1
        continue
    rlen = read.query_length
    if rlen not in phase1_lengths:
        continue
    offset = psite_offsets.get(rlen)
    if offset is None:
        continue
    if read.is_reverse:
        pos5   = read.reference_end - 1
        p_site = pos5 - offset
        strand = "-"
    else:
        pos5   = read.reference_start
        p_site = pos5 + offset
        strand = "+"
    chrom_col.append(read.reference_name)
    psite_col.append(p_site)
    strand_col.append(strand)
    length_col.append(rlen)
bam.close()

n_loaded = len(chrom_col)
print(f"  {n_loaded:,} phase1-length reads loaded")

# ── 5. PyRanges join: P-sites vs CDS exons (with Phase + cds_genomic_start) ───
print("Joining P-site positions to CDS exons...", flush=True)

reads_df = pd.DataFrame({
    "length": length_col,
    "frame":   np.nan,
    "rel_pos": np.nan,
})

if n_loaded > 0 and len(cds_df) > 0:
    read_idx = np.arange(n_loaded)
    psite_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": chrom_col,
        "Start":      psite_col,
        "End":        np.array(psite_col) + 1,
        "Strand":     strand_col,
        "read_idx":   read_idx,
    }))
    cds_pr = pr.PyRanges(
        cds_df[["Chromosome", "Start", "End", "Strand", "Phase", "cds_genomic_start"]]
    )
    joined = psite_pr.join(cds_pr, strandedness="same")

    if not joined.df.empty:
        jdf  = joined.df.copy()
        plus = jdf["Strand"] == "+"

        # subtract Phase: using +Phase mislabels in-frame P-sites in phase-1/2 exons.
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
        reads_df.loc[jdf["read_idx"].values, "frame"]   = jdf["frame"].values
        reads_df.loc[jdf["read_idx"].values, "rel_pos"] = jdf["rel_pos"].values

        n_in_cds = int((reads_df["rel_pos"] >= 0).sum())
        print(f"  {n_in_cds:,} reads with P-site in CDS body (rel_pos >= 0)")
    else:
        print("  No P-site positions overlapped any CDS exon.")

print("Aggregating frame counts per read length...", flush=True)
out_rows = []
for rlen in sorted(phase1_lengths):
    mask = (
        reads_df["length"].eq(rlen)
        & reads_df["rel_pos"].notna()
        & (reads_df["rel_pos"] >= 0)
    )
    frames   = reads_df.loc[mask, "frame"]
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

out_path = os.path.join(dir_staging, f"{SAMPLE}_cds_psite_frame.csv")
frame_df.to_csv(out_path, index=False)

total_in_cds = frame_df["n_psite_in_cds"].sum()
total_f0     = frame_df["n_frame0"].sum()
total_f1     = frame_df["n_frame1"].sum()
total_f2     = frame_df["n_frame2"].sum()
pct_f0_total = round(total_f0 / total_in_cds * 100, 1) if total_in_cds else 0.0

print(f"\n  Per-length breakdown (periodic lengths):")
sel_df = frame_df[frame_df["periodic"]]
if len(sel_df):
    print(sel_df[["read_length", "psite_offset", "n_psite_in_cds",
                  "pct_frame0", "pct_frame1", "pct_frame2"]].to_string(index=False))

print(f"\n  Grand total (all phase1 lengths):")
print(f"    P-sites in CDS body: {total_in_cds:,}")
print(f"    Frame 0 (in-frame):  {total_f0:,}  ({pct_f0_total:.1f}%)")
print(f"    Frame 1:             {total_f1:,}  ({round(total_f1/total_in_cds*100,1) if total_in_cds else 0:.1f}%)")
print(f"    Frame 2:             {total_f2:,}  ({round(total_f2/total_in_cds*100,1) if total_in_cds else 0:.1f}%)")
print(f"  Saved: {out_path}")
print("Done.")
