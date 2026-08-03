#!/usr/bin/env python3
"""Measure how often the assignment RULE, rather than the data, decides a read's transcript."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

def log(message):
    print("[assignment] %s" % message, flush=True)

def psite_exposure(chroms, positions, strands, cds_pr):
    """Candidate counts for the P-site rule, over the stage-1 CDS join.

    `first_exon_overlap` applies no criterion at all -- `drop_duplicates(keep="first")`
    takes whichever row the join emitted first -- so every read with two or more candidate
    transcripts is decided by the rule, and ambiguous == tied by construction.
    """
    import pyranges as pr

    if not len(chroms):
        return {"n_candidates_total": 0, "n_ambiguous": 0, "n_tied": 0}
    reads = pr.PyRanges(pd.DataFrame({
        "Chromosome": chroms, "Start": positions, "End": positions + 1,
        "Strand": strands, "read_idx": np.arange(len(chroms), dtype=np.int64)}))
    joined = reads.join(cds_pr, strandedness="same").df
    if joined.empty:
        return {"n_candidates_total": 0, "n_ambiguous": 0, "n_tied": 0}

    per_read = joined.groupby("read_idx", sort=False)["transcript_id"].nunique()
    n_ambiguous = int((per_read >= 2).sum())
    return {"n_candidates_total": int(len(per_read)),
            "n_ambiguous": n_ambiguous,
            "n_tied": n_ambiguous}

def footprint_exposure(blocks, cds_pr):
    """Candidate counts for the footprint rule, over the stage-1 CDS join.

    `max_exon_overlap` sums the overlap per (read, transcript) and takes the argmax, so a
    read is TIED only when two or more transcripts share that maximum. That is strictly
    fewer reads than are ambiguous.
    """
    import pyranges as pr

    chroms, starts, ends, strands, read_ids, _n_reads = blocks
    if not len(chroms):
        return {"n_candidates_total": 0, "n_ambiguous": 0, "n_tied": 0}
    block_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": chroms, "Start": starts, "End": ends,
        "Strand": strands, "read_idx": read_ids}))
    joined = block_pr.join(cds_pr, strandedness="same").df
    if joined.empty:
        return {"n_candidates_total": 0, "n_ambiguous": 0, "n_tied": 0}

    overlap = (np.minimum(joined["End"].to_numpy(), joined["End_b"].to_numpy())
               - np.maximum(joined["Start"].to_numpy(), joined["Start_b"].to_numpy()))
    joined = joined.assign(olen=overlap)
    joined = joined[joined["olen"] > 0]
    if joined.empty:
        return {"n_candidates_total": 0, "n_ambiguous": 0, "n_tied": 0}

    totals = joined.groupby(["read_idx", "transcript_id"], sort=False)["olen"].sum()
    totals = totals.reset_index()
    per_read = totals.groupby("read_idx", sort=False)
    n_candidates = per_read.size()
    best = per_read["olen"].transform("max")
    at_best = totals.assign(at_best=totals["olen"].eq(best)).groupby(
        "read_idx", sort=False)["at_best"].sum()
    return {"n_candidates_total": int(len(n_candidates)),
            "n_ambiguous": int((n_candidates >= 2).sum()),
            "n_tied": int((at_best >= 2).sum())}

def measure(sample, genome_bam, gtf, appris, regions, qc_genome, annotation_cache=None):
    """Run both rules over one sample and return the exposure report."""
    import annotation_cache as ac
    import build_shared_coverage as bsc
    import psite_placement

    started = time.time()
    bundle, reused = ac.load_or_build(annotation_cache, gtf, appris, regions, 35, 10)
    transcripts = bundle["transcripts"]
    index_of_id = bundle["index_of_id"]
    cds_table = bundle["cds_table"]
    coverage_offset = transcripts["coverage_offset"].to_numpy()
    cds_total_by_id = dict(zip(transcripts["transcript_id"], transcripts["cds_len_gtf"]))

    exon_pr = bsc._exon_pyranges(bundle["exons"], transcripts)
    cds_pr = bsc._cds_pyranges(cds_table)
    offsets = psite_placement.load_offsets(qc_genome, sample)
    log("offsets: %s" % offsets)

    log("P-sites: streaming the genome BAM")
    chroms, positions, strands = bsc.read_genome_psites(genome_bam, offsets)
    psite_counts = psite_exposure(chroms, positions, strands, cds_pr)
    _indices, psite_stats = bsc.project_genome_psites(
        chroms, positions, strands, exon_pr, cds_pr, cds_total_by_id,
        index_of_id, coverage_offset)
    del chroms, positions, strands

    log("footprints: streaming the genome BAM")
    blocks = bsc.read_genome_blocks(genome_bam, set(offsets))
    footprint_counts = footprint_exposure(blocks, cds_pr)
    _starts, _ends, footprint_stats = bsc.project_genome_footprints(
        blocks, exon_pr, cds_pr, index_of_id, coverage_offset)
    del blocks

    def merge(exposure, stats):
        merged = dict(exposure)
        merged.update(stats)
        total = merged["n_candidates_total"]
        merged["pct_ambiguous"] = round(100.0 * merged["n_ambiguous"] / total, 3) if total else 0.0
        merged["pct_tied"] = round(100.0 * merged["n_tied"] / total, 3) if total else 0.0
        return merged

    return {
        "sample": sample,
        "genome_bam": Path(genome_bam).name,
        "genome_uniqueness": "NH==1",
        "read_lengths": sorted(offsets),
        "annotation_cache_reused": reused,
        "policies": {
            "psite": dict(rule="first_exon_overlap", stage1_exon_set="cds_only",
                          tie_break="historical_join_order",
                          utr_fallback="first_exon_overlap",
                          **merge(psite_counts, psite_stats)),
            "footprint": dict(rule="max_exon_overlap", stage1_exon_set="cds_only",
                              tie_break="historical_join_order",
                              utr_fallback="max_exon_overlap",
                              **merge(footprint_counts, footprint_stats)),
        },
        "elapsed_minutes": round((time.time() - started) / 60, 2),
    }

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome-bam", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--appris", required=True, type=Path)
    parser.add_argument("--regions", type=Path, default=None)
    parser.add_argument("--qc-genome", required=True, type=Path,
                        help="genome readlen_window_qc.csv")
    parser.add_argument("--annotation-cache", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None,
                        help="optionally write the full JSON report here. The numbers are "
                             "printed either way; this is a validation program and its "
                             "report is not a pipeline artifact.")
    return parser

def main(argv=None):
    args = _build_parser().parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    for label, path in (("--genome-bam", args.genome_bam), ("--gtf", args.gtf),
                        ("--appris", args.appris), ("--qc-genome", args.qc_genome)):
        if not path.exists():
            raise SystemExit("%s does not exist: %s" % (label, path))

    report = measure(args.sample, args.genome_bam, args.gtf, args.appris, args.regions,
                     args.qc_genome, args.annotation_cache)

    for name, policy in report["policies"].items():
        print("  %-10s %s" % (name, policy["rule"]))
        print("      reads with a candidate   %12d" % policy["n_candidates_total"])
        print("      ambiguous (>= 2)         %12d  %6.3f %%"
              % (policy["n_ambiguous"], policy["pct_ambiguous"]))
        print("      tied under the criterion %12d  %6.3f %%"
              % (policy["n_tied"], policy["pct_tied"]))

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        log("wrote %s" % args.output)
    return 0

if __name__ == "__main__":
    sys.exit(main())
