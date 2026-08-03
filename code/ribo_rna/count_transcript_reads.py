#!/usr/bin/env python3
"""CDS-assigned read counts, and the ribo-vs-RNA route comparison, for one sample."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

COUNT_COLUMNS = ("genome_ribo_reads", "genome_rna_reads",
                 "txome_ribo_reads", "txome_rna_reads")

ROUTE_COLUMNS = ("sample", "route", "region", "n_transcripts",
                 "spearman_rho", "pearson_log2_raw",
                 "ribo_reads", "rna_reads", "ribo_library", "rna_library",
                 "ribo_cds_frac", "rna_cds_frac",
                 "ribo_ambiguous_excluded", "rna_ambiguous_excluded")

REGION = "cds"

class CountError(RuntimeError):
    pass

def log(message):
    print("[counts] %s" % message, flush=True)

def count_sample(sample, ribo_genome_bam, ribo_txome_bam, rna_genome_bam, rna_txome_bam,
                 gtf, appris, qc_genome, qc_txome, regions=None, annotation_cache=None):
    """One pass over each of the four BAMs -> (counts frame, route frame, universe report).

    The counts frame has one row per transcript in the shared universe; the route frame has
    one row per route.
    """
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str(HERE.parent / "coverage"))
    import psite_placement
    import ribo_rna_lib as rrl

    bundle, reused = rrl.load_bundle(gtf, appris, regions, annotation_cache)
    log("%s: annotation bundle %s" % (sample, "reused" if reused else "built"))

    universe, report = rrl.build_universe(bundle, ribo_txome_bam)
    log("%s: shared reference set %d APPRIS transcripts "
        "(%d annotated but unnamed by the reference, %d named but unannotated)"
        % (sample, report["n_universe"], report["n_annotated_not_in_reference"],
           report["n_reference_not_annotated"]))
    if not universe:
        raise CountError(
            "the annotation and the transcriptome reference share no transcript. Check "
            "that --appris and --ribo-txome-bam describe the same reference build.")

    cds_spans = rrl.transcript_cds_spans(bundle, universe)
    cds_pr = rrl.genome_cds_intervals(bundle, universe)

    genome_lengths = set(psite_placement.load_selected_lengths(qc_genome, sample))
    txome_lengths = set(psite_placement.load_selected_lengths(qc_txome, sample))
    log("%s: selected read lengths  genome %s  transcriptome %s"
        % (sample, sorted(genome_lengths), sorted(txome_lengths)))

    for path in (ribo_genome_bam, ribo_txome_bam, rna_genome_bam, rna_txome_bam):
        rrl.assert_single_end(path)

    counts = {}
    counts["genome_ribo_reads"], n_rb_g, amb_rb_g, lib_rb_g = rrl.count_genome_cds(
        ribo_genome_bam, cds_pr, stranded=True, read_lengths=genome_lengths)
    counts["genome_rna_reads"], n_rn_g, amb_rn_g, lib_rn_g = rrl.count_genome_cds(
        rna_genome_bam, cds_pr, stranded=False, read_lengths=None)
    counts["txome_ribo_reads"], n_rb_t, lib_rb_t = rrl.count_txome_cds(
        ribo_txome_bam, cds_spans, read_lengths=txome_lengths)
    counts["txome_rna_reads"], n_rn_t, lib_rn_t = rrl.count_txome_cds(
        rna_txome_bam, cds_spans, read_lengths=None)

    frame = pd.DataFrame({"transcript_id": universe})
    for column in COUNT_COLUMNS:
        mapping = counts[column]
        frame[column] = np.array([mapping.get(k, 0) for k in universe], dtype=int)

    rows = []
    for route, ribo, rna, n_rb, n_rn, lib_rb, lib_rn, amb_rb, amb_rn in (
            ("genome", counts["genome_ribo_reads"], counts["genome_rna_reads"],
             n_rb_g, n_rn_g, lib_rb_g, lib_rn_g, amb_rb_g, amb_rn_g),
            ("transcriptome", counts["txome_ribo_reads"], counts["txome_rna_reads"],
             n_rb_t, n_rn_t, lib_rb_t, lib_rn_t, 0, 0)):
        row = rrl.correlate(ribo, rna, universe)
        row.update({
            "sample": sample, "route": route, "region": REGION,
            "ribo_reads": n_rb, "rna_reads": n_rn,
            "ribo_library": lib_rb, "rna_library": lib_rn,
            "ribo_cds_frac": n_rb / lib_rb if lib_rb else float("nan"),
            "rna_cds_frac": n_rn / lib_rn if lib_rn else float("nan"),
            "ribo_ambiguous_excluded": amb_rb, "rna_ambiguous_excluded": amb_rn,
        })
        rows.append(row)
    route_frame = pd.DataFrame(rows)[list(ROUTE_COLUMNS)]
    return frame, route_frame, report

def write_table(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
    log("wrote %s (%d rows)" % (path, len(frame)))
    return path

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--ribo-genome-bam", required=True, type=Path)
    parser.add_argument("--ribo-txome-bam", required=True, type=Path)
    parser.add_argument("--rna-genome-bam", required=True, type=Path)
    parser.add_argument("--rna-txome-bam", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path,
                        help="GENCODE GTF -- the source of the canonical CDS")
    parser.add_argument("--appris", required=True, type=Path,
                        help="APPRIS transcript-lengths table (the reference headers)")
    parser.add_argument("--regions", type=Path,
                        help="optional actual-regions BED, as for the coverage build")
    parser.add_argument("--annotation-cache", type=Path,
                        help="the shared coordinate bundle; rebuilt if absent or stale")
    parser.add_argument("--qc-genome", required=True, type=Path,
                        help="genome readlen_window_qc.csv (the selected read lengths)")
    parser.add_argument("--qc-txome", required=True, type=Path,
                        help="transcriptome readlen_window_qc.csv")
    parser.add_argument("--route-output", type=Path,
                        help="route-comparison TSV (default results/ribo_rna/"
                             "_staging_cds/<sample>.tsv)")
    parser.add_argument("--counts-output", type=Path,
                        help="per-transcript count TSV; omit to skip writing it")
    return parser

def main(argv=None):
    args = _build_parser().parse_args(argv)
    missing = ["  %-20s %s" % (flag, getattr(args, attr))
               for flag, attr in (("--ribo-genome-bam", "ribo_genome_bam"),
                                  ("--ribo-txome-bam", "ribo_txome_bam"),
                                  ("--rna-genome-bam", "rna_genome_bam"),
                                  ("--rna-txome-bam", "rna_txome_bam"),
                                  ("--gtf", "gtf"), ("--appris", "appris"),
                                  ("--qc-genome", "qc_genome"), ("--qc-txome", "qc_txome"))
               if not getattr(args, attr).exists()]
    if missing:
        raise SystemExit("these required inputs do not exist:\n" + "\n".join(missing))
    frame, route_frame, _report = count_sample(
        args.sample, args.ribo_genome_bam, args.ribo_txome_bam,
        args.rna_genome_bam, args.rna_txome_bam,
        args.gtf, args.appris, args.qc_genome, args.qc_txome,
        args.regions, args.annotation_cache)

    route_output = args.route_output or (
        Path("results/ribo_rna/_staging_%s" % REGION) / ("%s.tsv" % args.sample))
    write_table(route_frame, route_output)
    for row in route_frame.to_dict("records"):
        log("  %-13s Spearman rho=%.4f  Pearson log2(raw+1)=%.4f  "
            "(ribo CDS %d/%d = %.1f%%, rna %.1f%%)"
            % (row["route"], row["spearman_rho"], row["pearson_log2_raw"],
               row["ribo_reads"], row["ribo_library"], 100 * row["ribo_cds_frac"],
               100 * row["rna_cds_frac"]))
    excluded = int(route_frame["ribo_ambiguous_excluded"].sum()
                   + route_frame["rna_ambiguous_excluded"].sum())
    log("  cross-gene CDS overlaps excluded (genome route): %d reads" % excluded)

    if args.counts_output:
        write_table(frame, args.counts_output)
        for column in COUNT_COLUMNS:
            log("  %-20s %12d reads" % (column, int(frame[column].sum())))

    return 0

if __name__ == "__main__":
    sys.exit(main())
