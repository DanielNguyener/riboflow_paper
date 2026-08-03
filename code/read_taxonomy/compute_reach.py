#!/usr/bin/env python3
"""Per-sample genome-anchored reach + concordance."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reach_lib as rl
fc = rl.fc
cl = rl.cl

STAGING = rl.OUTDIR / "_staging"
TAXONOMY_TSV = fc.output_root() / "read_taxonomy" / "taxonomy" / "taxonomy_all.tsv"
CONCORDANCE_TSV = (fc.output_root() / ".cache" / "read_taxonomy"
                   / "alignment_concordance" / "alignment_concordance_all.tsv")

def compute_sample(sample, dump_reads=False, log=print):
    tax_row = pd.read_csv(TAXONOMY_TSV, sep="\t").set_index("sample").loc[sample]
    conc_row = pd.read_csv(CONCORDANCE_TSV, sep="\t").set_index("sample").loc[sample]

    n_genome_unique = int(tax_row["n_genome_unique"])
    n_gUtU = int(tax_row["n_gU_tU"])
    n_gUtM = int(tax_row["n_gU_tM"])
    n_gUtA = int(tax_row["n_gU_tA"])
    assert n_gUtU + n_gUtM + n_gUtA == n_genome_unique, "gU partition mismatch vs taxonomy_all.tsv"
    assert int(conc_row["n_gUtU"]) == n_gUtU, "concordance n_gUtU mismatch vs taxonomy_all.tsv"

    n_concordant = int(conc_row["n_concordant"])
    n_discordant = n_gUtU - n_concordant

    log(f"[{sample}] N_G={n_genome_unique} gU_tU={n_gUtU} gU_tM={n_gUtM} gU_tA={n_gUtA}")

    log(f"[{sample}] recomputing gU_tA qname set...")
    genome_all, genome_uniq = rl.genome_status_sets(fc.genome_bam(sample))
    txome_all = rl.txome_all_qnames(fc.txome_bam(sample))
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

    row = {"sample": sample, "n_genome_unique": n_genome_unique,
           "n_gUtU": n_gUtU, "n_gUtM": n_gUtM, "n_gUtA": n_gUtA}
    row["n_shared_unique_concordant"] = n_concordant
    row["n_shared_unique_discordant"] = n_discordant
    row["n_genome_unique_transcriptome_multimapped"] = n_gUtM
    for cat in rl.UNREACHABLE_CATEGORIES + ["representable_not_present_in_dedup_bam"]:
        row[f"n_{cat}"] = int(counts.get(cat, 0))

    check_sum = (row["n_shared_unique_concordant"] + row["n_shared_unique_discordant"] +
                 row["n_genome_unique_transcriptome_multimapped"] +
                 sum(row[f"n_{c}"] for c in rl.UNREACHABLE_CATEGORIES + ["representable_not_present_in_dedup_bam"]))
    assert check_sum == n_genome_unique, f"partition does not sum to N_G ({check_sum} != {n_genome_unique})"

    for cat in rl.REACH_CATEGORIES:
        row[f"pct_{cat}"] = 100.0 * row[f"n_{cat}"] / n_genome_unique

    c_shared = 100.0 * n_concordant / n_gUtU if n_gUtU else float("nan")
    r_shared = 100.0 * n_gUtU / n_genome_unique
    c_overall = 100.0 * n_concordant / n_genome_unique
    row["C_shared_pct"] = c_shared
    row["R_shared_pct"] = r_shared
    row["C_overall_pct"] = c_overall
    assert abs(c_overall - r_shared * c_shared / 100.0) < 1e-6, "C_overall != R_shared*C_shared"

    STAGING.mkdir(parents=True, exist_ok=True)
    out = STAGING / f"{sample}.tsv"
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False)
    log(f"[{sample}] wrote {out}")

    if dump_reads:
        _dump_read_table(sample, gUtA_qnames, genome_blocks, table, gene2tid, exon_gene_pr, labels)

    return row

def _all_gene_body_pr():
    import pyranges as pr
    df = fc.config.load_all_gene_bodies()
    return pr.PyRanges(df.reset_index(drop=True))

def _dump_read_table(sample, gUtA_qnames, genome_blocks, table, gene2tid, exon_gene_pr, labels):
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
