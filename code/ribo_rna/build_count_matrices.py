#!/usr/bin/env python3
"""Raw CDS counts for the cohort: runs `count_transcript_reads.py` per sample and pivots the
four count columns into the transcripts x samples matrices `code/te_route/normalization.R`
consumes; `--check` compares with data/ribo_rna/counts/. Run with `python` (3.9)."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
COUNTER = HERE / "count_transcript_reads.py"
MANIFEST = REPO / "config" / "cohort_manifest.tsv"

#: count_transcript_reads.py column per (assay, route) pair, each normalized independently:
#: ribo/RNA are never pooled, and neither are the two routes.
MATRICES = {
    ("ribo", "genome"): "genome_ribo_reads",
    ("ribo", "txome"): "txome_ribo_reads",
    ("rna", "genome"): "genome_rna_reads",
    ("rna", "txome"): "txome_rna_reads",
}

#: The shipped matrices, used by --check.
REFERENCE_COUNTS = REPO / "data" / "ribo_rna" / "counts"

class BuildError(RuntimeError):
    pass

def log(message):
    print("[counts] %s" % message, flush=True)

# ── the cohort ───────────────────────────────────────────────────────────────

def load_manifest(path=MANIFEST):
    """Sample ids and their four BAM paths, in manifest order.

    Sample ORDER comes from the manifest so all four matrices share column order.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str)
    required = ["sample_id", "ribo_genome_bam", "ribo_txome_bam",
                "rna_genome_bam", "rna_txome_bam"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise BuildError("%s has no %s column(s)" % (path, ", ".join(missing)))
    return frame[required].copy()

def resolve_bams(row, bams_root):
    """Manifest-relative BAM paths against the tree root."""
    return {flag: bams_root / row[column] for flag, column in (
        ("--ribo-genome-bam", "ribo_genome_bam"),
        ("--ribo-txome-bam", "ribo_txome_bam"),
        ("--rna-genome-bam", "rna_genome_bam"),
        ("--rna-txome-bam", "rna_txome_bam"))}

def require_inputs(manifest, args):
    """Fail before any BAM is opened, naming every missing input at once."""
    missing = []
    for path, label in ((args.gtf, "--gtf"), (args.appris, "--appris"),
                        (args.qc_genome, "--qc-genome"), (args.qc_txome, "--qc-txome")):
        if not Path(path).exists():
            missing.append("  %-18s %s" % (label, path))
    for _, row in manifest.iterrows():
        for flag, path in resolve_bams(row, args.bams).items():
            if not path.exists():
                missing.append("  %-18s %s" % (flag, path))
    if missing:
        raise BuildError("these inputs do not exist:\n" + "\n".join(missing))

# ── running the counter ──────────────────────────────────────────────────────

def counts_path(sample, output):
    return Path(output) / "per_sample" / ("%s_cds_counts.tsv" % sample)

def run_sample(row, args):
    """One sample through count_transcript_reads.py. Returns the sample id on failure."""
    sample = row["sample_id"]
    target = counts_path(sample, args.output)
    if args.skip_existing and target.exists():
        log("[skip] %s (already counted)" % sample)
        return None

    command = [sys.executable, str(COUNTER), "--sample", sample,
               "--gtf", str(args.gtf), "--appris", str(args.appris),
               "--qc-genome", str(args.qc_genome), "--qc-txome", str(args.qc_txome),
               "--counts-output", str(target),
               # route table is an unread by-product here; scratch it.
               "--route-output", str(Path(args.output) / "_route_scratch" / ("%s.tsv" % sample)),
               "--annotation-cache", str(Path(args.output) / ".cache" / "annotation.pkl")]
    for flag, path in resolve_bams(row, args.bams).items():
        command += [flag, str(path)]
    if args.regions:
        command += ["--regions", str(args.regions)]

    environment = dict(os.environ)
    # imported modules resolve indirect BAM paths through this variable.
    environment["RIBOFLOW_PAPER_BAMS"] = str(Path(args.bams).resolve())
    environment.setdefault("MPLBACKEND", "Agg")

    log("%s: counting" % sample)
    try:
        subprocess.run(command, check=True, env=environment)
    except subprocess.CalledProcessError as exc:
        print("  !! %s FAILED (exit %d)" % (sample, exc.returncode), file=sys.stderr,
              flush=True)
        return sample
    return None

# ── pivoting ─────────────────────────────────────────────────────────────────

def build_matrices(samples, output):
    """The four per-sample count columns, pivoted to transcripts x samples.

    All matrices are asserted to share the same transcript index and order.
    """
    per_sample = {}
    for sample in samples:
        path = counts_path(sample, output)
        if not path.exists():
            raise BuildError("%s was never written; rerun without --skip-existing" % path)
        frame = pd.read_csv(path, sep="\t")
        frame = frame.set_index("transcript_id")
        per_sample[sample] = frame

    index = per_sample[samples[0]].index
    for sample, frame in per_sample.items():
        if not frame.index.equals(index):
            raise BuildError(
                "%s has a different transcript set or order than %s. The two runs saw "
                "different annotations; do not normalize across them."
                % (sample, samples[0]))

    written = []
    for (assay, route), column in sorted(MATRICES.items()):
        matrix = pd.DataFrame({s: per_sample[s][column] for s in samples}, index=index)
        if matrix.isna().any().any():
            raise BuildError("%s/%s contains NA" % (assay, route))
        if (matrix < 0).any().any():
            raise BuildError("%s/%s contains a negative count" % (assay, route))
        matrix = matrix.astype("int64")
        target = Path(output) / "counts" / ("%s_counts_%s.csv" % (assay, route))
        target.parent.mkdir(parents=True, exist_ok=True)
        matrix.to_csv(target, index_label="transcript_id", lineterminator="\n")
        log("wrote %s (%d transcripts x %d samples, %d reads)"
            % (target, matrix.shape[0], matrix.shape[1], int(matrix.to_numpy().sum())))
        written.append(target)
    return written

# ── checks against the shipped tables ────────────────────────────────────────

def check_against_reference(samples, output):
    """Every regenerated count must equal the shipped matrix, column by column."""
    failures = []
    compared = 0
    for (assay, route), _column in sorted(MATRICES.items()):
        name = "%s_counts_%s.csv" % (assay, route)
        shipped_path = REFERENCE_COUNTS / name
        if not shipped_path.exists():
            log("check: %s is not shipped, skipping" % name)
            continue
        shipped = pd.read_csv(shipped_path, index_col="transcript_id")
        mine = pd.read_csv(Path(output) / "counts" / name, index_col="transcript_id")
        if not mine.index.equals(shipped.index):
            failures.append("%s: transcript index differs from the shipped matrix" % name)
            continue
        for sample in samples:
            if sample not in shipped.columns:
                log("check: %s has no shipped column %s, skipping" % (name, sample))
                continue
            compared += 1
            if not mine[sample].equals(shipped[sample]):
                n = int((mine[sample] != shipped[sample]).sum())
                failures.append("%s: %s differs in %d transcript(s)" % (name, sample, n))
    if failures:
        raise BuildError("the regenerated counts do not match the shipped matrices:\n  "
                         + "\n  ".join(failures))
    log("check: %d (matrix, sample) column(s) reproduce the shipped matrices exactly"
        % compared)

# ── cli ──────────────────────────────────────────────────────────────────────

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bams", required=True, type=Path,
                        help="RiboFlow output tree; the manifest's paths are relative to it")
    parser.add_argument("--gtf", required=True, type=Path, help="GENCODE v34 GTF")
    parser.add_argument("--appris", required=True, type=Path,
                        help="APPRIS transcript-lengths table matching the transcriptome BAMs")
    parser.add_argument("--qc-genome", type=Path, required=True,
                        help="the genome route's readlen_window_qc.csv")
    parser.add_argument("--qc-txome", type=Path, required=True,
                        help="the transcriptome route's readlen_window_qc.csv")
    parser.add_argument("--regions", type=Path, help="optional actual-regions BED")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--output", type=Path, default=REPO / "results" / "ribo_rna")
    parser.add_argument("--samples", default=None,
                        help="comma-separated subset; default is every manifest sample")
    parser.add_argument("--workers", type=int, default=2,
                        help="concurrent samples; each holds one sample's annotation and "
                             "four open BAMs, so 2 is the tested setting")
    parser.add_argument("--skip-existing", action="store_true",
                        help="leave already-written per-sample count tables alone")
    parser.add_argument("--no-check", action="store_true",
                        help="skip the comparison against the shipped matrices")
    parser.add_argument("--pivot-only", action="store_true",
                        help="rebuild the matrices from existing per-sample tables")
    return parser

def main(argv=None):
    args = _build_parser().parse_args(argv)

    manifest = load_manifest(args.manifest)
    if args.samples:
        wanted = [s.strip() for s in args.samples.split(",") if s.strip()]
        unknown = [s for s in wanted if s not in set(manifest["sample_id"])]
        if unknown:
            raise BuildError("not in %s: %s" % (args.manifest, ", ".join(unknown)))
        manifest = manifest[manifest["sample_id"].isin(wanted)]
    samples = list(manifest["sample_id"])
    if not samples:
        raise BuildError("no samples selected")

    if not args.pivot_only:
        require_inputs(manifest, args)
        log("counting %d sample(s) with %d worker(s): %s"
            % (len(samples), args.workers, ", ".join(samples)))
        rows = [row for _, row in manifest.iterrows()]
        with ThreadPoolExecutor(max_workers=max(1, min(args.workers, len(rows)))) as pool:
            failed = [s for s in pool.map(lambda r: run_sample(r, args), rows) if s]
        if failed:
            raise BuildError("these samples failed: %s" % ", ".join(failed))

    build_matrices(samples, args.output)
    if not args.no_check:
        check_against_reference(samples, args.output)
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except BuildError as error:
        raise SystemExit("[counts] %s" % error)
