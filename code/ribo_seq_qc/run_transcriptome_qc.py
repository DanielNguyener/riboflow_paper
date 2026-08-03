#!/usr/bin/env python3
"""Cohort driver for the TRANSCRIPTOME-route Ribo-seq QC -- the twin of run_pipeline.py."""

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

#: step name -> (script, staging suffix). Order is load-bearing.
STEP_SCRIPTS = {
    "qc":        ("01t_readlen_psite_qc_transcriptome.py", "readlen_window_qc"),
    "cds_frame": ("03t_cds_frame_transcriptome.py", "cds_psite_frame"),
}
STEP_ORDER = ["qc", "cds_frame"]

MASTER_TABLES = {
    "readlen_window_qc.csv": "readlen_window_qc",
    "cds_psite_frame.csv": "cds_psite_frame",
}

DEFAULT_BAM_GLOB = "*/transcriptome/alignment_ribo/merged/*.transcriptome.post_dedup.bam"
TX_SUFFIX = ".transcriptome.post_dedup.bam"

def sample_from_tx_bam(path):
    """A2780.transcriptome.post_dedup.bam -> A2780 (also strips a flat .bam)."""
    base = os.path.basename(path)
    if base.endswith(TX_SUFFIX):
        return base[: -len(TX_SUFFIX)]
    return config.sample_from_bam(path)

def discover_samples(bam_dir, pattern):
    bams = sorted(glob.glob(os.path.join(bam_dir, pattern)))
    return [(sample_from_tx_bam(b), b) for b in bams]

def staging_dir():
    return os.path.join(config.tx_out_dir(), "tables", "_staging")

def staging_path(sample, suffix):
    return os.path.join(staging_dir(), "%s_%s.csv" % (sample, suffix))

def run_step(step, sample, bam, frame0_threshold, plots=False):
    script = os.path.join(HERE, STEP_SCRIPTS[step][0])
    command = [sys.executable, script, "--sample", sample, "--bam", bam,
               "--out", config.tx_out_dir()]
    if step == "qc":
        command += ["--frame0-threshold", str(frame0_threshold)]
        if plots:
            command.append("--plots")
    print("\n$ %s" % " ".join(command), flush=True)
    subprocess.run(command, check=True)

def run_sample(sample, bam, steps, skip_existing, frame0_threshold, plots=False):
    """Run the requested steps for one sample -> list of (sample, step) failures."""
    failures = []
    print("\n%s\nSAMPLE: %s\n  BAM: %s\n%s" % ("=" * 70, sample, bam, "=" * 70), flush=True)
    for step in steps:
        suffix = STEP_SCRIPTS[step][1]
        if skip_existing and os.path.exists(staging_path(sample, suffix)):
            print("  [%s] [skip-existing] %s already staged" % (sample, step), flush=True)
            continue
        try:
            run_step(step, sample, bam, frame0_threshold, plots)
        except subprocess.CalledProcessError as exc:
            print("  !! %s/%s FAILED (exit %d)" % (sample, step, exc.returncode),
                  file=sys.stderr, flush=True)
            failures.append((sample, step))
            if step == "qc":
                break
    return failures

def aggregate(samples_done):
    """Concatenate the staging CSVs into the two master tables."""
    import pandas as pd

    tables_dir = os.path.join(config.tx_out_dir(), "tables")
    os.makedirs(tables_dir, exist_ok=True)
    print("\n=== Aggregating transcriptome master tables ===", flush=True)
    for filename, suffix in MASTER_TABLES.items():
        rows = []
        for sample in samples_done:
            path = staging_path(sample, suffix)
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
    parser.add_argument("--bam-dir", required=True,
                        help="Folder containing the per-sample BAM trees. Required.")
    parser.add_argument("--bam-glob", default=DEFAULT_BAM_GLOB,
                        help="Glob relative to --bam-dir for transcriptome ribo BAMs.")
    parser.add_argument("--samples", default=None,
                        help="Comma-separated subset of sample names to process.")
    parser.add_argument("--steps", default=",".join(STEP_ORDER),
                        help="Comma-separated steps to run (qc,cds_frame).")
    parser.add_argument("--frame0-threshold", type=float, default=50.0)
    parser.add_argument("--plots", action="store_true",
                        help="also write the per-sample metagene PDFs. Off by default: "
                             "diagnostics, not results.")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip a step if its staging output already exists.")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="Skip the per-sample steps; rebuild the masters from staging.")
    args = parser.parse_args()

    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    unknown = [s for s in steps if s not in STEP_SCRIPTS]
    if unknown:
        parser.error("unknown step(s): %s; valid: %s" % (unknown, list(STEP_SCRIPTS)))
    steps = [s for s in STEP_ORDER if s in steps]

    samples = discover_samples(args.bam_dir, args.bam_glob)
    if not samples:
        parser.error("no BAMs matching %r found in %s" % (args.bam_glob, args.bam_dir))
    if args.samples:
        wanted = {x.strip() for x in args.samples.split(",")}
        samples = [(s, b) for (s, b) in samples if s in wanted]
        if not samples:
            parser.error("none of --samples %s matched BAMs in %s"
                         % (sorted(wanted), args.bam_dir))

    os.makedirs(staging_dir(), exist_ok=True)

    if args.aggregate_only:
        aggregate([s for s, _ in samples])
        print("\nDone. Aggregation only -- no per-sample steps run.")
        return

    print("Discovered %d transcriptome sample(s): %s"
          % (len(samples), [s for s, _ in samples]))
    print("Steps: %s" % steps)

    failures = []
    n_workers = min(MAX_WORKERS, len(samples))
    if n_workers > 1:
        print("\nRunning %d sample(s) in parallel (workers=%d)..."
              % (len(samples), n_workers), flush=True)
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(run_sample, s, b, steps, args.skip_existing,
                                   args.frame0_threshold, args.plots): s
                       for s, b in samples}
            for future in as_completed(futures):
                failures += future.result()
    else:
        for sample, bam in samples:
            failures += run_sample(sample, bam, steps, args.skip_existing,
                                   args.frame0_threshold, args.plots)

    aggregate([s for s, _ in samples])

    print("\nDone. %d sample(s) processed." % len(samples))
    if failures:
        print("FAILURES (%d): %s" % (len(failures), failures), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
