#!/usr/bin/env python3
"""Per-sample genome-anchored reach: where the genome-unique reads that the transcriptome
route lacks fall relative to their gene's selected transcript."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reach_lib as rl
import reference_lib as cl
fc = rl.fc

STAGING = rl.OUTDIR / "_staging"
TAXONOMY_TSV = fc.output_root() / "read_taxonomy" / "taxonomy" / "taxonomy_all.tsv"
def compute_sample(sample, dump_reads=False, log=print):
    tax_row = pd.read_csv(TAXONOMY_TSV, sep="\t").set_index("sample").loc[sample]

    n_genome_unique = int(tax_row["n_genome_unique"])
    n_gUtP = int(tax_row["n_gU_tP"])
    n_gUtA = int(tax_row["n_gU_tA"])
    assert n_gUtP + n_gUtA == n_genome_unique, "gU partition mismatch vs taxonomy_all.tsv"
    log(f"[{sample}] N_G={n_genome_unique} gU_tP={n_gUtP} gU_tA={n_gUtA}")

    log(f"[{sample}] recomputing gU_tA qname set...")
    genome_all, genome_uniq = rl.tl.genome_status_sets(fc.genome_bam(sample))
    txome_all = rl.tl.txome_present_qnames(fc.txome_bam(sample))
    gUtA_qnames = genome_uniq - txome_all
    assert len(gUtA_qnames) == n_gUtA, (
        f"recomputed gU_tA ({len(gUtA_qnames)}) != taxonomy_all.tsv ({n_gUtA})")

    log(f"[{sample}] reading genome blocks for {len(gUtA_qnames)} gU_tA reads...")
    genome_blocks = rl.read_genome_blocks(fc.genome_bam(sample), gUtA_qnames)

    transcript_payload = cl.build_transcript_table()
    table = transcript_payload["table"]
    selected_genes = set(v["gene_id"] for v in table.values())
    gene2tid = rl.gene_to_transcript_map(table)
    exon_gene_df = cl.build_exon_gene_table()
    exon_gene_pr = cl.load_exon_gene_pr()
    all_gene_body_pr = _all_gene_body_pr()
    omitted_genes = rl.omitted_pc_genes(exon_gene_df, selected_genes)

    log(f"[{sample}] classifying gU_tA reads...")
    labels = rl.classify_gU_tA(gUtA_qnames, genome_blocks, exon_gene_pr, exon_gene_df,
                               all_gene_body_pr, table, gene2tid, omitted_genes)
    counts = labels.value_counts()

    log(f"[{sample}] testing direct overlap with omitted exonic sequence...")
    omitted_index = rl.build_omitted_exon_index(exon_gene_df, table)
    overlap_qnames = rl.omitted_exon_overlap_qnames(gUtA_qnames, genome_blocks, omitted_index)

    row = {"sample": sample, "n_genome_unique": n_genome_unique,
           "n_gU_tP": n_gUtP, "n_gU_tA": n_gUtA}
    for cat in rl.REACH_CATEGORIES:
        row[f"n_{cat}"] = int(counts.get(cat, 0))
    check_sum = n_gUtP + sum(row[f"n_{c}"] for c in rl.REACH_CATEGORIES)
    assert check_sum == n_genome_unique, f"partition does not sum to N_G ({check_sum} != {n_genome_unique})"
    for cat in rl.REACH_CATEGORIES:
        row[f"pct_{cat}"] = 100.0 * row[f"n_{cat}"] / n_gUtA if n_gUtA else float("nan")
    # Figure 4D: NOT a partition slice. The gene-level alternative-exon test of Figure 5A
    # applied to the whole gU_tA population; each read ID counts at most once.
    row["n_omitted_exon_overlap"] = len(overlap_qnames)
    row["pct_omitted_exon_overlap"] = (
        100.0 * len(overlap_qnames) / n_gUtA if n_gUtA else float("nan"))
    log(f"[{sample}] omitted-exon overlap {len(overlap_qnames)} / {n_gUtA} "
        f"({row['pct_omitted_exon_overlap']:.2f}%)")

    STAGING.mkdir(parents=True, exist_ok=True)
    out = STAGING / f"{sample}.tsv"
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False)
    log(f"[{sample}] wrote {out}")

    if dump_reads:
        _dump_read_table(sample, gUtA_qnames, genome_blocks, table, gene2tid, exon_gene_pr,
                         labels, overlap_qnames)

    return row

def _all_gene_body_pr():
    import pyranges as pr
    df = fc.config.load_all_gene_bodies()
    return pr.PyRanges(df.reset_index(drop=True))

def _dump_read_table(sample, gUtA_qnames, genome_blocks, table, gene2tid, exon_gene_pr,
                     labels, overlap_qnames):
    """Debug dump: one row per gU_tA read, per the user's requested ~20-column schema.
    Single-sample only: `gU_tA` alone is ~1M rows per sample."""
    import pyranges as pr

    rows = []
    for q in gUtA_qnames:
        rec = genome_blocks.get(q)
        if rec is None:
            continue
        chrom, strand, blocks = rec
        rows.append({
            "sample_id": sample, "read_id": q,
            "genome_chromosome": chrom,
            "genome_start": min(b[0] for b in blocks),
            "genome_end": max(b[1] for b in blocks),
            "genome_strand": strand,
            "genome_n_blocks": len(blocks),
            "comparison_category": labels.get(q, "other_unclassified"),
            "overlaps_omitted_exonic_sequence": q in overlap_qnames,
            "present_in_deduplicated_transcriptome_bam": False,
            "transcriptome_mapping_status": "absent",
        })
    df = pd.DataFrame(rows)
    out = rl.OUTDIR / f"{sample}_gUtA_read_dump.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[{sample}] wrote debug read dump {out} ({len(df)} rows)")

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--dump-reads", action="store_true",
                    help="also write a per-read debug TSV for this sample only")
    args = ap.parse_args()
    compute_sample(args.sample, dump_reads=args.dump_reads)

if __name__ == "__main__":
    main()
