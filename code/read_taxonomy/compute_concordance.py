#!/usr/bin/env python3
"""Per-sample read-ID-level alignment concordance."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import concordance_lib as cl
fc = cl.fc

STAGING = cl.OUTDIR / "_staging"

def compute_sample(sample, log=print):
    txome_payload = cl.build_transcript_table()
    base2ver = txome_payload["base2ver"]

    log(f"[{sample}] reading txome BAM...")
    # One pass keeps EVERY primary (any MAPQ) + the unique subset, so the gU_tU and
    txome_present, txome_uniq_q, txome_all_q = cl.read_txome_primary(fc.txome_bam(sample), base2ver)
    txome_dict = {q: txome_present[q] for q in txome_uniq_q}
    log(f"[{sample}] reading genome BAM...")
    genome_dict = cl.read_genome_unique(fc.genome_bam(sample))

    genome_q = set(genome_dict)
    shared = set(txome_dict) & genome_q
    log(f"[{sample}] n_txome_unique={len(txome_dict)} n_txome_present={len(txome_all_q)} "
        f"n_genome_unique={len(genome_dict)} n_gUtU={len(shared)}")

    exon_gene_pr = cl.load_exon_gene_pr()
    result = cl.classify_all(shared, txome_dict, genome_dict, txome_payload, exon_gene_pr)

    counts = result["label"].value_counts().reindex(cl.CATEGORIES, fill_value=0)
    n = int(counts.sum())
    row = {"sample": sample, "n_gUtU": n, "n_genome_unique": len(genome_dict),
           "n_txome_unique": len(txome_dict)}
    for cat in cl.CATEGORIES:
        row[f"n_{cat}"] = int(counts[cat])
        row[f"pct_{cat}"] = 100.0 * counts[cat] / n if n else float("nan")

    # ── broader population: genome-unique ∩ txome-PRESENT (any MAPQ) ──────────────
    shared_present = genome_q & txome_all_q
    classifiable = genome_q & set(txome_present)
    n_present = len(shared_present)
    n_unresolved = n_present - len(classifiable)
    result_p = cl.classify_all(classifiable, txome_present, genome_dict, txome_payload, exon_gene_pr)
    counts_p = result_p["label"].value_counts().reindex(cl.CATEGORIES, fill_value=0)
    row["n_gU_txPresent"] = n_present
    row["n_unresolved_present"] = n_unresolved
    row["pct_unresolved_present"] = 100.0 * n_unresolved / n_present if n_present else float("nan")
    for cat in cl.CATEGORIES:
        row[f"n_{cat}_present"] = int(counts_p[cat])
        row[f"pct_{cat}_present"] = 100.0 * counts_p[cat] / n_present if n_present else float("nan")
    log(f"[{sample}] n_gU_txPresent={n_present} n_concordant_present={int(counts_p['concordant'])} "
        f"pct_concordant_present={row['pct_concordant_present']:.4f} n_unresolved={n_unresolved}")

    STAGING.mkdir(parents=True, exist_ok=True)
    out = STAGING / f"{sample}.tsv"
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False)
    log(f"[{sample}] wrote {out}")
    return row

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    args = ap.parse_args()
    compute_sample(args.sample)

if __name__ == "__main__":
    main()
