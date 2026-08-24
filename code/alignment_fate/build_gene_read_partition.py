#!/usr/bin/env python3
"""Gene-anchored read-ID partition: every read at a gene, on either route, in one chain.

`--dump-reads` writes the per-read table Figure 6A is folded from.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def log(message):
    print("[partition] %s" % message, flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome-bam", required=True, type=Path,
                        help="coordinate-sorted and INDEXED: the gene side is a region fetch")
    parser.add_argument("--transcriptome-bam", required=True, type=Path)
    parser.add_argument("--gene-id", default="", help="comma-separated gene IDs")
    parser.add_argument("--transcript-id", default="", help="comma-separated transcript IDs")
    parser.add_argument("--coverage", type=Path,
                        help="a shared_coverage.h5, used only to resolve gene IDs and names")
    parser.add_argument("--output", type=Path, default=Path("results/alignment_fate"))
    parser.add_argument("--dump-reads", action="store_true",
                        help="also write one row per read id, with its uncollapsed reach "
                             "category (a few thousand rows per gene, not a cohort dump)")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import gene_read_partition_lib as lib

    genes = [g.strip() for g in args.gene_id.split(",") if g.strip()]
    transcripts = [t.strip() for t in args.transcript_id.split(",") if t.strip()]
    if not genes and not transcripts:
        raise SystemExit("give --gene-id and/or --transcript-id")

    coverage = None
    if args.coverage:
        sys.path.insert(0, str(REPO / "code" / "coverage"))
        import coverage_schema
        coverage = coverage_schema.open_coverage(args.coverage)
    try:
        log("%s: resolving %d gene(s) and %d transcript(s)"
            % (args.sample, len(genes), len(transcripts)))
        wide, tidy, dump = lib.compute_partition(
            args.sample, args.genome_bam, args.transcriptome_bam,
            gene_ids=genes, transcript_ids=transcripts, coverage=coverage, log=log)
    finally:
        if coverage is not None:
            coverage.close()

    args.output.mkdir(parents=True, exist_ok=True)
    tidy_path = args.output / ("%s.gene_read_partition.tsv" % args.sample)
    tidy.to_csv(tidy_path, sep="\t", index=False, lineterminator="\n")
    wide_path = args.output / ("%s.gene_read_partition_summary.tsv" % args.sample)
    wide.to_csv(wide_path, sep="\t", index=False, lineterminator="\n")

    for _, row in wide.iterrows():
        log("  %-20s union %6d = genome-side %6d + transcriptome-only %5d "
            "(shared %6d, genome unique %6d / multi %6d)"
            % (row["gene_name"] or row["transcript_id"], row["n_union"],
               row["n_genome_side"], row["n_txome_only"], row["n_shared"],
               row["n_genome_unique"], row["n_genome_multi"]))

    # Categories must sum to the union and to 100%; re-checked here so a hand-edited
    # table cannot be plotted.
    for tid, group in tidy.groupby("transcript_id", sort=False):
        n_union = int(wide.loc[wide["transcript_id"] == tid, "n_union"].iat[0])
        total = int(group["n_reads"].sum())
        if total != n_union:
            raise SystemExit("%s: categories sum to %d, not the union %d"
                             % (tid, total, n_union))
        if abs(float(group["pct_of_union"].sum()) - 100.0) > 1e-4 and n_union:
            raise SystemExit("%s: percentages sum to %.6f, not 100"
                             % (tid, group["pct_of_union"].sum()))

    log("wrote %s" % tidy_path)
    log("wrote %s" % wide_path)
    if args.dump_reads:
        dump_path = args.output / ("%s.gene_read_partition_reads.tsv" % args.sample)
        dump.to_csv(dump_path, sep="\t", index=False, lineterminator="\n")
        log("wrote %s (%d rows)" % (dump_path, len(dump)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
