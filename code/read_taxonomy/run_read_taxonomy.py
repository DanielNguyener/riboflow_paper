#!/usr/bin/env python3
"""Cohort driver for every read-taxonomy analysis."""
from __future__ import annotations

import argparse
import glob
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import taxonomy_lib as tl
fc = tl.fc

DEFAULT_WORKERS = 2
TX_GLOB = "*/transcriptome/alignment_ribo/merged/*.transcriptome.post_dedup.bam"

def _out(*parts):
    return fc.output_root().joinpath("read_taxonomy", *parts)

def _cache(*parts):
    return fc.output_root().joinpath(".cache", "read_taxonomy", *parts)

ANALYSES = {
    "taxonomy": {
        "worker": "compute_taxonomy.py",
        "staging": _out("taxonomy", "_staging"),
        "master": _out("taxonomy", "taxonomy_all.tsv"),
    },
    "alignment_concordance": {
        "worker": "compute_concordance.py",
        "staging": _cache("alignment_concordance", "_staging"),
        "master": _cache("alignment_concordance", "alignment_concordance_all.tsv"),
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

def discover_samples():
    """Samples with BOTH a transcriptome and a genome ribo BAM."""
    found = []
    for path in sorted(glob.glob(str(fc.bams_root() / TX_GLOB))):
        sample = Path(path).name[: -len(".transcriptome.post_dedup.bam")]
        if fc.genome_bam(sample).exists():
            found.append(sample)
    return found

def run_sample(analysis, sample, skip_existing):
    spec = ANALYSES[analysis]
    staged = spec["staging"] / ("%s.tsv" % sample)
    if skip_existing and staged.exists():
        print("  [%s] [skip] staged" % sample, flush=True)
        return None
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
        cells = ["pct_gU_tU", "pct_gU_tM", "pct_gM_tU", "pct_gM_tM"]
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
    parser.add_argument("--skip-existing", action="store_true",
                        help="skip samples that already have a staged TSV")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="rebuild the master from existing staging, run nothing")
    args = parser.parse_args(argv)

    samples = discover_samples()
    if args.samples:
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        samples = [s for s in samples if s in wanted]
    if not samples:
        parser.error("no samples discovered")

    spec = ANALYSES[args.analysis]
    spec["staging"].mkdir(parents=True, exist_ok=True)

    if args.aggregate_only:
        aggregate(args.analysis, samples)
        return 0

    print("[%s] %d sample(s), %d worker(s)"
          % (args.analysis, len(samples), args.workers), flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=min(args.workers, len(samples))) as pool:
        futures = {pool.submit(run_sample, args.analysis, s, args.skip_existing): s
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
