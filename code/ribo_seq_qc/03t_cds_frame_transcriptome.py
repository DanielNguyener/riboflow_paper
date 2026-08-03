#!/usr/bin/env python3
"""STEP 03t — Whole-CDS P-site frame counts, run on TRANSCRIPTOME alignments."""
import re

import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
import bam_inputs as fc          # the one uniqueness policy: fc.is_unique_txome_read

import qc_core
qc_core.require("pysam")

import pysam
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description="Whole-CDS P-site frame counts (transcriptome).")
p.add_argument("--sample", required=True)
p.add_argument("--bam",    required=True)
p.add_argument("--out",    default=config.tx_out_dir())
args = p.parse_args()

SAMPLE = args.sample
BAM    = args.bam
OUT    = args.out

MIN_LEN, MAX_LEN = config.MIN_LEN, config.MAX_LEN

dir_tables  = os.path.join(OUT, "tables")
dir_staging = os.path.join(dir_tables, "_staging")
for d in [dir_tables, dir_staging]:
    os.makedirs(d, exist_ok=True)

QC_CSV = os.path.join(dir_staging, f"{SAMPLE}_readlen_window_qc.csv")
print(f"=== [03t cds_frame / transcriptome] sample={SAMPLE} ===", flush=True)

# ── 1. Load phase-1 lengths, periodic flags, P-site offsets (from step 01t) ───
print("Loading phase-1 lengths and P-site offsets from step 01t...", flush=True)
qc_df = pd.read_csv(QC_CSV)
_as_bool = lambda s: s.map(lambda x: str(x).strip().lower() in ("true", "1"))

phase1_rows      = qc_df[_as_bool(qc_df["in_phase1"])]
phase1_lengths   = set(phase1_rows["read_length"].astype(int).tolist())
periodic_lengths = set(qc_df.loc[_as_bool(qc_df["periodic"]), "read_length"].astype(int).tolist())
psite_offsets    = dict(zip(phase1_rows["read_length"].astype(int),
                            phase1_rows["psite_offset"].astype(int)))
print(f"  Phase-1 lengths:  {sorted(phase1_lengths)}")
print(f"  Periodic lengths: {sorted(periodic_lengths)}")

# Reference names embed "|CDS:start-end|" (1-based, inclusive, stop codon INCLUDED).
_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def _cds_bounds_from_refname(ref):
    """(cds_start0, cds_len_nostop) or (None, None). Stop codon trimmed (-3)."""
    m = _CDS_RE.search(ref)
    if not m:
        return None, None
    start1, end1 = int(m.group(1)), int(m.group(2))
    cds_start0 = start1 - 1                 # 1-based -> 0-based start codon position
    cds_len_nostop = (end1 - start1 + 1) - 3  # drop the stop codon (CDS:end includes it)
    if cds_len_nostop <= 0:
        return None, None
    return cds_start0, cds_len_nostop

print("Reading transcriptome BAM (phase-1 lengths only)...", flush=True)
bam = pysam.AlignmentFile(BAM, "rb")
ref_bounds = {ref: _cds_bounds_from_refname(ref) for ref in bam.references}
n_refs_with_cds = sum(1 for v in ref_bounds.values() if v[0] is not None)
print(f"  {len(ref_bounds):,} references; {n_refs_with_cds:,} carry a CDS region",
      flush=True)

length_col = []
frame_col  = []
n_seen = 0
for read in bam.fetch(until_eof=True):
    if not fc.is_unique_txome_read(read):         # MAPQ >= 42
        continue
    # NO MAPQ filter: transcriptome BAM used as delivered (see module docstring).
    rlen = read.query_length
    if rlen not in phase1_lengths:
        continue
    cds_start0, cds_len_nostop = ref_bounds.get(read.reference_name, (None, None))
    if cds_start0 is None:
        continue
    offset = psite_offsets.get(rlen)
    if offset is None:
        continue
    # bowtie2 --norc: all reads forward on the transcript, 5' end = reference_start
    p_site  = read.reference_start + offset
    rel_pos = p_site - cds_start0
    if rel_pos < 0 or rel_pos >= cds_len_nostop:
        continue
    n_seen += 1
    length_col.append(rlen)
    frame_col.append(rel_pos % 3)
bam.close()
print(f"  {n_seen:,} reads with P-site in CDS body (0 <= rel_pos < CDS len, stop excluded)")

reads_df = pd.DataFrame({"length": length_col, "frame": frame_col})

print("Aggregating frame counts per read length...", flush=True)
out_rows = []
for rlen in sorted(phase1_lengths):
    frames  = reads_df.loc[reads_df["length"].eq(rlen), "frame"]
    n_psite = len(frames)
    n_f0 = int((frames == 0).sum())
    n_f1 = int((frames == 1).sum())
    n_f2 = int((frames == 2).sum())
    pct_f0 = round(n_f0 / n_psite * 100, 1) if n_psite else 0.0
    pct_f1 = round(n_f1 / n_psite * 100, 1) if n_psite else 0.0
    pct_f2 = round(n_f2 / n_psite * 100, 1) if n_psite else 0.0
    out_rows.append({
        "read_length":    rlen,
        "in_phase1":      True,
        "periodic":       rlen in periodic_lengths,
        "psite_offset":   psite_offsets[rlen],
        "n_psite_in_cds": n_psite,
        "n_frame0":       n_f0,
        "n_frame1":       n_f1,
        "n_frame2":       n_f2,
        "pct_frame0":     pct_f0,
        "pct_frame1":     pct_f1,
        "pct_frame2":     pct_f2,
    })

frame_df = pd.DataFrame(out_rows)

out_path = os.path.join(dir_staging, f"{SAMPLE}_cds_psite_frame.csv")
frame_df.to_csv(out_path, index=False)

total_in_cds = frame_df["n_psite_in_cds"].sum()
total_f0     = frame_df["n_frame0"].sum()
pct_f0_total = round(total_f0 / total_in_cds * 100, 1) if total_in_cds else 0.0
sel_df = frame_df[frame_df["periodic"]]
if len(sel_df):
    print("\n  Per-length breakdown (periodic lengths):")
    print(sel_df[["read_length", "psite_offset", "n_psite_in_cds",
                  "pct_frame0", "pct_frame1", "pct_frame2"]].to_string(index=False))
print(f"\n  Grand total (all phase-1 lengths): P-sites in CDS body {total_in_cds:,}, "
      f"frame0 {pct_f0_total:.1f}%")
print(f"  Saved: {out_path}")
print("Done.")
