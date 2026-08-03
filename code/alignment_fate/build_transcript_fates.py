#!/usr/bin/env python3
"""Per-transcript alignment fates for arbitrary samples and gene/transcript IDs."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

def log(message):
    print("[fates] %s" % message, flush=True)

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome-bam", required=True, type=Path)
    parser.add_argument("--transcriptome-bam", required=True, type=Path)
    parser.add_argument("--gene-id", default="", help="comma-separated gene IDs")
    parser.add_argument("--transcript-id", default="", help="comma-separated transcript IDs")
    parser.add_argument("--coverage", type=Path,
                        help="a shared_coverage.h5, used only to resolve gene IDs")
    parser.add_argument("--output", type=Path, default=Path("results/alignment_fate"))
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import transcript_fate_lib as lib

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
        wide, tidy = lib.compute_fates(
            args.sample, args.genome_bam, args.transcriptome_bam,
            gene_ids=genes, transcript_ids=transcripts, coverage=coverage)
    finally:
        if coverage is not None:
            coverage.close()

    args.output.mkdir(parents=True, exist_ok=True)
    tidy_path = args.output / ("%s.transcript_alignment_fates.tsv" % args.sample)
    tidy.to_csv(tidy_path, sep="\t", index=False, lineterminator="\n")

    for _, row in wide.iterrows():
        log("  %-20s assigned %7d   unique %7d   pseudogene-tie %7d   multi %7d"
            % (row["tid"], row["n_txome_assigned"], row["n_genome_unique"],
               row["n_top_on_target_and_pp_tie"], row["n_genome_multi"]))
    # The three categories must sum to the transcriptome-assigned denominator, or a stacked
    for tid, group in tidy.groupby("transcript_id", sort=False):
        total = int(group["n_reads"].sum())
        denominator = int(wide.loc[wide["tid"] == tid, "n_txome_assigned"].iloc[0])
        assert total == denominator, "%s: %d != %d" % (tid, total, denominator)
        log("  %-20s categories sum to %d == n_txome_assigned" % (tid, total))
    log("wrote %s" % tidy_path.name)
    return 0

if __name__ == "__main__":
    sys.exit(main())
