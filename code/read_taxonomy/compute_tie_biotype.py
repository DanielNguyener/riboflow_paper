#!/usr/bin/env python3
"""Per-sample count of shared genome-multimapper reads whose primary alignment is score-tied with a protein_coding / processed_pseudogene secondary."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tie_biotype_lib as tie
cl = tie.cl
fc = tie.fc

STAGING = tie.OUTDIR / "_staging_tie"
TAXONOMY = fc.output_root() / "read_taxonomy" / "taxonomy" / "taxonomy_all.tsv"
CATS = ["cross_pc_pp", "cross_pp_pc", "same_pc_pc", "same_pp_pp"]

def compute_sample(sample, log=print):
    exon_pr = cl.load_exon_gene_pr()
    gene_pr = tie.bl.gene_body_pr()

    log(f"[{sample}] reading txome BAM (present qname set)...")
    t_all, _ = tie.tl.status_sets(fc.txome_bam(sample), "txome")
    log(f"[{sample}] n_txome_present={len(t_all):,}; enumerating genome multimapper loci...")
    recs = tie.read_genome_multi_records_flagged(fc.genome_bam(sample), t_all)
    n_reads = len(recs)

    if TAXONOMY.exists():
        tr = pd.read_csv(TAXONOMY, sep="\t")
        tr = tr[tr["sample"] == sample]
        if len(tr):
            exp = int(tr.iloc[0]["n_gM_tU"]) + int(tr.iloc[0]["n_gM_tM"])
            log(f"[{sample}] taxonomy check dark-green reads {n_reads} vs {exp} -> "
                f"{'OK' if n_reads == exp else 'MISMATCH'}")
            assert n_reads == exp, f"[{sample}] dark-green population disagrees with taxonomy_all.tsv"

    counts, n_reads = tie.categorize(recs, exon_pr, gene_pr)
    row = {"sample": sample, "n_reads": n_reads, **{c: counts[c] for c in CATS},
           "n_qualifying": counts["n_qualifying"]}
    for c in CATS + ["n_qualifying"]:
        row[f"pct_{c}"] = 100.0 * counts.get(c, row.get(c)) / n_reads if n_reads else float("nan")
    row["pct_n_qualifying"] = 100.0 * counts["n_qualifying"] / n_reads if n_reads else float("nan")

    log(f"[{sample}] n_reads={n_reads:,} qualifying={counts['n_qualifying']:,} "
        f"({row['pct_n_qualifying']:.1f}%): "
        f"cross_pc_pp {counts['cross_pc_pp']:,}, cross_pp_pc {counts['cross_pp_pc']:,}, "
        f"same_pc_pc {counts['same_pc_pc']:,}, same_pp_pp {counts['same_pp_pp']:,}")

    STAGING.mkdir(parents=True, exist_ok=True)
    dest = STAGING / f"{sample}.tsv"
    pd.DataFrame([row]).to_csv(dest, sep="\t", index=False)
    log(f"[{sample}] wrote {dest}")
    return row

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    args = ap.parse_args()
    compute_sample(args.sample)

if __name__ == "__main__":
    main()
