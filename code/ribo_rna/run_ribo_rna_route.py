#!/usr/bin/env python3
"""Driver: the CDS ribo-vs-RNA route comparison over every sample in the cohort."""
import argparse
import glob
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import count_transcript_reads as ctr
import ribo_rna_lib as rrl
fc = rrl.fc

HERE = Path(__file__).resolve().parent
MAX_WORKERS = 6
OUT = fc.output_root() / "ribo_rna"
TX_GLOB = "*/transcriptome/alignment_ribo/merged/*.transcriptome.post_dedup.bam"

STAGING = OUT / ("_staging_%s" % ctr.REGION)
MASTER = OUT / ("ribo_rna_route_%s.tsv" % ctr.REGION)

def counts_path(sample):
    return OUT / ("ribo_rna_counts_raw_%s_%s.tsv" % (sample, ctr.REGION))

def discover_samples():
    bams = sorted(glob.glob(str(fc.bams_root() / TX_GLOB)))
    out = []
    for b in bams:
        s = Path(b).name[: -len(".transcriptome.post_dedup.bam")]
        if (fc.genome_bam(s).exists() and rrl.rna_genome_bam(s).exists()
                and rrl.rna_txome_bam(s).exists()):
            out.append(s)
    return out

def run_sample(sample, args):
    if args.skip_existing and (STAGING / ("%s.tsv" % sample)).exists():
        print("  [%s] [skip] staged" % sample, flush=True)
        return None
    cmd = [sys.executable, str(HERE / "count_transcript_reads.py"),
           "--sample", sample,
           "--ribo-genome-bam", str(fc.genome_bam(sample)),
           "--ribo-txome-bam", str(fc.txome_bam(sample)),
           "--rna-genome-bam", str(rrl.rna_genome_bam(sample)),
           "--rna-txome-bam", str(rrl.rna_txome_bam(sample)),
           "--gtf", str(args.gtf), "--appris", str(args.appris),
           "--qc-genome", str(args.qc_genome), "--qc-txome", str(args.qc_txome),
           "--route-output", str(STAGING / ("%s.tsv" % sample))]
    if args.regions:
        cmd += ["--regions", str(args.regions)]
    if args.annotation_cache:
        cmd += ["--annotation-cache", str(args.annotation_cache)]
    if sample == args.counts_sample:
        cmd += ["--counts-output", str(counts_path(sample))]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print("  !! %s FAILED (exit %d)" % (sample, exc.returncode),
              file=sys.stderr, flush=True)
        return sample
    return None

def aggregate(samples):
    import pandas as pd

    rows = [pd.read_csv(STAGING / ("%s.tsv" % s), sep="\t")
            for s in samples if (STAGING / ("%s.tsv" % s)).exists()]
    if not rows:
        print("no staging rows -- nothing aggregated", file=sys.stderr)
        return
    master = pd.concat(rows, ignore_index=True).sort_values(["sample", "route"])
    MASTER.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(MASTER, sep="\t", index=False, lineterminator="\n")
    print("\nwrote %s: %d samples" % (MASTER, master["sample"].nunique()))
    print(master.to_string(index=False))

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default=None)
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--gtf", type=Path)
    ap.add_argument("--appris", type=Path)
    ap.add_argument("--regions", type=Path)
    ap.add_argument("--annotation-cache", type=Path)
    ap.add_argument("--qc-genome", type=Path)
    ap.add_argument("--qc-txome", type=Path)
    ap.add_argument("--counts-sample", default="HeLa",
                    help="the one sample that also gets a per-transcript count table")
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args()

    samples = discover_samples()
    if args.samples:
        want = {x.strip() for x in args.samples.split(",")}
        samples = [s for s in samples if s in want]
    if not samples:
        ap.error("no samples discovered")

    if args.aggregate_only:
        aggregate(samples)
        return

    for flag, value in (("--gtf", args.gtf), ("--appris", args.appris),
                        ("--qc-genome", args.qc_genome), ("--qc-txome", args.qc_txome)):
        if value is None:
            ap.error("%s is required (only --aggregate-only can do without it)" % flag)
        if not value.exists():
            ap.error("%s does not exist: %s" % (flag, value))

    STAGING.mkdir(parents=True, exist_ok=True)
    print("Discovered %d samples: %s" % (len(samples), samples), flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=min(args.workers, len(samples))) as pool:
        futures = {pool.submit(run_sample, s, args): s for s in samples}
        for future in as_completed(futures):
            failed = future.result()
            if failed:
                failures.append(failed)

    aggregate(samples)
    if failures:
        print("FAILURES (%d): %s" % (len(failures), failures), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
