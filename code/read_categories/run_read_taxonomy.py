#!/usr/bin/env python3
"""Cohort driver for every read-taxonomy analysis."""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import taxonomy_lib as tl
fc = tl.fc

DEFAULT_WORKERS = 2

def _out(*parts):
    return fc.output_root().joinpath("read_taxonomy", *parts)

ANALYSES = {
    "taxonomy": {
        "worker": "compute_taxonomy.py",
        "staging": _out("taxonomy", "_staging"),
        "master": _out("taxonomy", "taxonomy_all.tsv"),
    },
    "reach": {
        "worker": "compute_reach.py",
        "staging": _out("reach", "_staging"),
        "master": _out("reach", "genome_anchored_reach_all.tsv"),
    },
    "tie_biotype": {
        "worker": "compute_tie_biotype.py",
        "staging": _out("multimap_biotype", "_staging_tie"),
        "master": _out("multimap_biotype", "multimap_tie_biotype_all.tsv"),
    },
}


def run_sample(analysis, sample):
    spec = ANALYSES[analysis]
    staged = spec["staging"] / ("%s.tsv" % sample)
    command = [sys.executable, str(HERE / spec["worker"]), "--sample", sample]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        print("  !! %s FAILED (exit %d)" % (sample, exc.returncode),
              file=sys.stderr, flush=True)
        return sample
    return None

def aggregate(analysis, samples):
    """Concatenate the per-sample staging TSVs into the analysis's master table."""
    import pandas as pd

    spec = ANALYSES[analysis]
    staging, master = spec["staging"], spec["master"]
    frames = [pd.read_csv(staging / ("%s.tsv" % s), sep="\t")
              for s in samples if (staging / ("%s.tsv" % s)).exists()]
    if not frames:
        print("no staging rows in %s -- nothing aggregated" % staging, file=sys.stderr)
        return None
    table = pd.concat(frames, ignore_index=True).sort_values("sample")
    master.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(master, sep="\t", index=False, lineterminator="\n")
    print("\nwrote %s (%d samples)" % (master, table["sample"].nunique()))
    _summarise(analysis, table)
    return master

def _summarise(analysis, table):
    """A short console summary. Console only -- nothing downstream parses this."""
    if analysis == "taxonomy":
        cells = ["pct_" + tl.cell_key(g, t) for g, t in tl.CELLS]
        if all(c in table.columns for c in cells):
            median = table[cells].median()
            print("core-cell medians: " + "  ".join(
                "%s=%.2f%%" % (c.replace("pct_", ""), median[c]) for c in cells))

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("analysis", choices=sorted(ANALYSES),
                        help="which read-taxonomy analysis to run over the cohort")
    parser.add_argument("--samples", default=None, help="comma-separated subset")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args(argv)

    samples = fc.discover_samples()
    if args.samples:
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        samples = [s for s in samples if s in wanted]
    if not samples:
        parser.error("no samples discovered")

    spec = ANALYSES[args.analysis]
    spec["staging"].mkdir(parents=True, exist_ok=True)

    print("[%s] %d sample(s), %d worker(s)"
          % (args.analysis, len(samples), args.workers), flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=min(args.workers, len(samples))) as pool:
        futures = {pool.submit(run_sample, args.analysis, s): s
                   for s in samples}
        for future in as_completed(futures):
            failed = future.result()
            if failed:
                failures.append(failed)

    aggregate(args.analysis, samples)
    if failures:
        print("FAILURES (%d): %s" % (len(failures), failures), file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
