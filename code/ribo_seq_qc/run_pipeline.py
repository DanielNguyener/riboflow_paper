#!/usr/bin/env python3
"""Cohort driver for the Ribo-seq QC, on the genome route (default) or the transcriptome route."""

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

_HERE = os.path.dirname(os.path.abspath(__file__))
_COMMON = os.path.join(os.path.dirname(_HERE), "common")
for _entry in (_HERE, _COMMON, os.path.join(_COMMON, "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import config

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_WORKERS = 10

SAMPLE_SCRIPT = "sample_qc.py"   # one program  one BAM traversal  both staging CSVs

BAM_GLOB = {
    "genome": "*/genome/alignment_ribo/merged/*.post_dedup.bam",
    "transcriptome": "*/transcriptome/alignment_ribo/merged/*.transcriptome.post_dedup.bam",
}
TX_SUFFIX = ".transcriptome.post_dedup.bam"

MASTER_TABLES = {
    "readlen_window_qc.csv": "readlen_window_qc",
    "cds_psite_frame.csv": "cds_psite_frame",
}

def sample_from_bam(path):
    """A2780.transcriptome.post_dedup.bam -> A2780; otherwise `config.sample_from_bam`."""
    base = os.path.basename(path)
    if base.endswith(TX_SUFFIX):
        return base[: -len(TX_SUFFIX)]
    return config.sample_from_bam(path)

def discover_samples(bam_dir, pattern):
    """Sorted [(sample, bam_path)] for every BAM matching the `bam_dir`-relative glob `pattern`."""
    bams = sorted(glob.glob(os.path.join(bam_dir, pattern)))
    return [(sample_from_bam(b), b) for b in bams]

def out_dir(route):
    """The route's output root, kept distinct so the genome masters are never clobbered."""
    return config.tx_out_dir() if route == "transcriptome" else config.out_dir()

def staging_path(sample, suffix, route):
    return os.path.join(out_dir(route), "tables", "_staging", "%s_%s.csv" % (sample, suffix))

def run_sample(sample, bam, route):
    """Run the per-sample QC program -> list of failed samples."""
    print("\n%s\nSAMPLE: %s\n  BAM: %s\n%s" % ("=" * 70, sample, bam, "=" * 70), flush=True)
    command = [
        sys.executable, os.path.join(HERE, SAMPLE_SCRIPT),
        "--sample", sample,
        "--bam", bam,
        "--route", route,
        "--out", out_dir(route),
    ]
    if route == "genome":
        command += ["--appris", config.appris_path(), "--gtf", config.gtf_path()]
    print("\n$ %s" % " ".join(command), flush=True)
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        print("  !! %s FAILED (exit %d); continuing with the next sample"
              % (sample, exc.returncode), file=sys.stderr, flush=True)
        return [sample]
    return []

def aggregate(samples_done, route):
    """Concatenate the staging CSVs into the two master tables."""
    import pandas as pd

    tables_dir = os.path.join(out_dir(route), "tables")
    os.makedirs(tables_dir, exist_ok=True)
    print("\n=== Aggregating master tables ===", flush=True)
    for filename, suffix in MASTER_TABLES.items():
        rows = []
        for sample in samples_done:
            path = staging_path(sample, suffix, route)
            if not os.path.exists(path):
                continue
            frame = pd.read_csv(path)
            frame.insert(0, "sample", sample)
            rows.append(frame)
        if rows:
            master = pd.concat(rows, ignore_index=True)
            master.to_csv(os.path.join(tables_dir, filename), index=False)
            print("  %s: %d rows, %d samples"
                  % (filename, len(master), master["sample"].nunique()))
        else:
            print("  %s: no staging inputs found -- skipped" % filename)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--route", choices=["genome", "transcriptome"], default="genome")
    parser.add_argument("--samples", default=None)
    args = parser.parse_args()

    bam_glob = BAM_GLOB[args.route]
    samples = discover_samples(args.bam_dir, bam_glob)
    if not samples:
        parser.error("no BAMs matching %r found in %s" % (bam_glob, args.bam_dir))
    if args.samples:
        wanted = {x.strip() for x in args.samples.split(",")}
        samples = [(s, b) for (s, b) in samples if s in wanted]
        if not samples:
            parser.error("none of --samples %s matched BAMs in %s"
                         % (sorted(wanted), args.bam_dir))

    os.makedirs(os.path.join(out_dir(args.route), "tables", "_staging"), exist_ok=True)

    print("Discovered %d %s sample(s): %s" % (len(samples), args.route, [s for s, _ in samples]))

    if args.route == "genome":
        print("\nBuilding/loading annotation cache...", flush=True)
        config.load_annotation()

    failures = []
    n_workers = min(MAX_WORKERS, len(samples))
    if n_workers > 1:
        print("\nRunning %d sample(s) in parallel (workers=%d)..." % (len(samples), n_workers),
              flush=True)
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(run_sample, s, b, args.route): s
                       for s, b in samples}
            for future in as_completed(futures):
                failures += future.result()
    else:
        for sample, bam in samples:
            failures += run_sample(sample, bam, args.route)

    aggregate([s for s, _ in samples], args.route)

    print("\nDone. %d sample(s) processed." % len(samples))
    if failures:
        print("FAILURES (%d): %s" % (len(failures), failures), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
