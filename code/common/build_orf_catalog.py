#!/usr/bin/env python3
"""Step 0 — build the matched ORF set + transcript<->genome coordinate map (once)."""
from __future__ import annotations

import argparse
import re
import sys

import numpy as np
import pandas as pd
import pysam

import bam_inputs as fc

_CDS_RE = re.compile(r"\|CDS:(\d+)-(\d+)\|")

def txome_ref_transcripts(bam_path) -> set:
    """base ENST present as @SQ in the transcriptome reference (with a CDS header)."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    refs = set()
    for ref in bam.references:
        if _CDS_RE.search(ref):
            refs.add(ref.split("|", 1)[0].split(".", 1)[0])
    bam.close()
    return refs

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--txome-bam", default=None,
                    help="Any transcriptome BAM (for the @SQ reference set). "
                         "Default: first discovered sample.")
    ap.add_argument("--junc-win", type=int, default=3,
                    help="±nt around an internal exon boundary counted as junction-proximal.")
    ap.add_argument("--out-dir", default=str(fc.output_root()))
    args = ap.parse_args()

    out = fc.output_root() if args.out_dir is None else __import__("pathlib").Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    txbam = args.txome_bam
    if txbam is None:
        samples = fc.discover_samples()
        if not samples:
            sys.exit("no samples with both genome+transcriptome BAMs found")
        txbam = fc.txome_bam(samples[0])
    print(f"Reading transcriptome @SQ reference set from:\n  {txbam}", flush=True)
    tx_refs = txome_ref_transcripts(txbam)
    print(f"  {len(tx_refs):,} transcripts in transcriptome reference", flush=True)

    print("Building CDS-exon coordinate table from the QC annotation cache …", flush=True)
    tbl = fc.build_cds_table()
    cds = tbl["cds"]
    meta = fc.config.load_appris_meta().set_index("transcript_id")
    lf = {tid.split(".", 1)[0]: bool(v)
          for tid, v in meta["length_filtered"].items()}

    transcripts = sorted(tbl["cds_total"])
    rows = []
    junctions = {}
    for tid in transcripts:
        base = tid.split(".", 1)[0]
        cds_len = int(tbl["cds_total"][tid])
        n_exons = int(tbl["n_exons"][tid])
        internal = tbl["junctions"][tid]
        junctions[tid] = internal
        rows.append({
            "transcript_id": tid,
            "base_enst": base,
            "gene_id": tbl["gene_id"][tid],
            "gene_name": tbl["gene_name"][tid],
            "chrom": tbl["chrom"][tid],
            "strand": tbl["strand"][tid],
            "cds_len": cds_len,
            "n_codons": cds_len // 3,
            "n_cds_exons": n_exons,
            "n_internal_junctions": int(internal.size),
            "div_by_3": (cds_len % 3 == 0),
            "length_filtered": lf.get(base, None),
            "in_txome_ref": base in tx_refs,
        })

    cat = pd.DataFrame(rows)

    n = len(cat)
    n_txref = int(cat["in_txome_ref"].sum())
    n_div3 = int(cat["div_by_3"].sum())
    n_lf = int((cat["length_filtered"] == True).sum())
    tot_pos = int(cat["cds_len"].sum())
    jp = 0
    for tid in transcripts:
        L = int(tbl["cds_total"][tid])
        if L <= 0:
            continue
        mask = np.zeros(L, dtype=bool)
        for b in junctions[tid]:
            lo = max(0, b - args.junc_win)
            hi = min(L, b + args.junc_win)
            mask[lo:hi] = True
        jp += int(mask.sum())

    print("\n=== ORF catalog accounting ===")
    print(f"transcripts (APPRIS principal, CDS in cache): {n:,}")
    print(f"  in transcriptome reference (@SQ):           {n_txref:,}  ({100*n_txref/n:.1f}%)")
    print(f"  CDS length divisible by 3:                  {n_div3:,}  ({100*n_div3/n:.1f}%)  "
          f"-> {n-n_div3} flagged not-divisible")
    print(f"  length_filtered==True (carried, not dropped): {n_lf:,}")
    print(f"junction-proximal CDS positions (±{args.junc_win} nt): "
          f"{jp:,} / {tot_pos:,}  ({100*jp/tot_pos:.1f}%)")

    cat_path = out / "orf_catalog.tsv"
    cat.to_csv(cat_path, sep="\t", index=False)
    print(f"\nWrote:\n  {cat_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
