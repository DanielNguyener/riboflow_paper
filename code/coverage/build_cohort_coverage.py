#!/usr/bin/env python3
"""Build shared-coverage HDF5 files for a set of samples from the cohort manifest."""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCHEMA_VERSION = "riboflow_paper/cohort-manifest/1"

BAM_COLUMNS = ("ribo_genome_bam", "ribo_genome_bai", "ribo_txome_bam", "ribo_txome_bai",
               "rna_genome_bam", "rna_genome_bai", "rna_txome_bam", "rna_txome_bai")
COVERAGE_COLUMNS = ("ribo_genome_bam", "ribo_genome_bai", "ribo_txome_bam", "ribo_txome_bai")

class CohortError(RuntimeError):
    pass

def log(message):
    print("[cohort] %s" % message, flush=True)

def read_manifest(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise CohortError("%s has no rows" % path)
    problems = []
    seen = set()
    for number, row in enumerate(rows, 2):
        version = (row.get("schema_version") or "").strip()
        if version != SCHEMA_VERSION:
            problems.append("line %d: schema_version %r, expected %r"
                            % (number, version, SCHEMA_VERSION))
        sample = (row.get("sample_id") or "").strip()
        if not sample:
            problems.append("line %d: empty sample_id" % number)
        elif sample in seen:
            problems.append("line %d: duplicate sample_id %r" % (number, sample))
        seen.add(sample)
    if problems:
        raise CohortError("%s is not a valid cohort manifest:\n%s"
                          % (path, "\n".join("  - %s" % p for p in problems)))
    return rows

def resolve(row, column, bams_root):
    value = (row.get(column) or "").strip()
    if not value:
        return None
    path = Path(value)
    return path if path.is_absolute() else (Path(bams_root) / path if bams_root else path)

def validate(rows, bams_root, columns=BAM_COLUMNS):
    """Report EVERY missing or empty file at once, before any compute starts."""
    problems = []
    for row in rows:
        for column in columns:
            path = resolve(row, column, bams_root)
            if path is None:
                problems.append("%s: no %s in the manifest" % (row["sample_id"], column))
            elif not path.exists():
                problems.append("%s: %s does not exist -- %s"
                                % (row["sample_id"], column, path))
            elif path.stat().st_size == 0:
                problems.append("%s: %s is empty -- %s" % (row["sample_id"], column, path))
    return problems

def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def build_one(row, args):
    sample = row["sample_id"]
    command = [
        sys.executable, str(HERE / "build_shared_coverage.py"),
        "--sample", sample,
        "--genome-bam", str(resolve(row, "ribo_genome_bam", args.bams)),
        "--transcriptome-bam", str(resolve(row, "ribo_txome_bam", args.bams)),
        "--gtf", str(args.gtf), "--appris", str(args.appris),
        "--qc-genome", str(args.qc_genome), "--qc-txome", str(args.qc_txome),
        "--output", str(args.output),
        "--trim", str(args.trim), "--assay", args.assay,
        "--annotation-cache", str(args.annotation_cache),
        "--gzip-level", str(args.gzip_level), "--chunk", str(args.chunk),
    ]
    # by content digest, and a second copy beside it can only go stale.
    if args.regions:
        command += ["--regions", str(args.regions)]
    if args.hash_bams:
        command.append("--hash-bams")
    if args.record_input_paths:
        command.append("--record-input-paths")

    started = time.time()
    completed = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.time() - started
    if completed.returncode != 0:
        tail = "\n".join((completed.stderr or completed.stdout).strip().splitlines()[-12:])
        return {"sample": sample, "ok": False, "seconds": round(elapsed, 1), "detail": tail}
    return {"sample": sample, "ok": True, "seconds": round(elapsed, 1)}

def write_checksums(output_dir, samples):
    import h5py

    rows = []
    for sample in samples:
        path = Path(output_dir) / ("%s.shared_coverage.h5" % sample)
        if not path.exists():
            continue
        with h5py.File(path, "r") as handle:
            provenance = handle["provenance"].attrs.get("json", "")
            record = {
                "sample_id": sample,
                "filename": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "schema_version": handle.attrs.get("schema", ""),
                "provenance_sha256": hashlib.sha256(
                    provenance.encode("utf-8")).hexdigest(),
                "n_transcripts": int(handle.attrs["n_transcripts"]),
                "n_positions": int(handle.attrs["n_positions"]),
                "created_utc": handle.attrs.get("created_utc", ""),
            }
        rows.append(record)
    if not rows:
        return None
    destination = Path(output_dir) / "coverage_checksums.tsv"
    with open(destination, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return destination

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path,
                        default=Path("config/cohort_manifest.tsv"))
    parser.add_argument("--bams", type=Path, help="root for relative manifest paths")
    parser.add_argument("--samples", help="comma-separated sample_id list")
    parser.add_argument("--all", action="store_true",
                        help="process every row -- must be given explicitly")
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--appris", type=Path)
    parser.add_argument("--regions", type=Path)
    parser.add_argument("--qc-genome", type=Path)
    parser.add_argument("--qc-txome", type=Path)
    parser.add_argument("--output", type=Path, default=Path("results/coverage"))
    parser.add_argument("--trim", type=int, default=15)
    parser.add_argument("--left-span", type=int, default=35)
    parser.add_argument("--right-span", type=int, default=10)
    parser.add_argument("--annotation-cache", type=Path, default=None,
                        help="the shared annotation bundle; built once and reused by every "
                             "sample (default <output>/../.cache/annotation/"
                             "coverage_annotation.pkl)")
    parser.add_argument("--assay", default="ribo", choices=("ribo", "rna"))
    parser.add_argument("--gzip-level", type=int, default=9)
    parser.add_argument("--chunk", type=int, default=1 << 16)
    parser.add_argument("--workers", type=int, default=1,
                        help="concurrent samples. A build peaks near 5 GB resident.")
    parser.add_argument("--hash-bams", action="store_true")
    parser.add_argument("--record-input-paths", action="store_true",
                        help="store full filesystem paths in each file's provenance")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--validate", action="store_true",
                        help="check every BAM and index, then exit")
    return parser

def main(argv=None):
    args = _build_parser().parse_args(argv)
    rows = read_manifest(args.manifest)
    log("manifest: %d samples, schema %s" % (len(rows), SCHEMA_VERSION))

    if args.validate:
        problems = validate(rows, args.bams)
        if problems:
            print("INVALID -- %d problem(s):" % len(problems))
            for problem in problems:
                print("  - %s" % problem)
            return 1
        print("VALID -- all %d samples have all %d alignment files and indexes."
              % (len(rows), len(BAM_COLUMNS)))
        return 0

    if args.all and args.samples:
        raise SystemExit("give --samples or --all, not both")
    if not args.all and not args.samples:
        raise SystemExit(
            "nothing selected. Pass --samples A,B or --all.\n"
            "There is no implicit whole-cohort run: the product is ~1 GB and hours of "
            "I/O, so it has to be asked for.")

    for required in ("gtf", "appris", "qc_genome", "qc_txome"):
        if getattr(args, required) is None:
            raise SystemExit("--%s is required to build" % required.replace("_", "-"))

    wanted = [r["sample_id"] for r in rows] if args.all else \
        [s.strip() for s in args.samples.split(",") if s.strip()]
    known = {r["sample_id"]: r for r in rows}
    unknown = [s for s in wanted if s not in known]
    if unknown:
        raise SystemExit("not in the manifest: %s" % ", ".join(unknown))

    selected = [known[s] for s in wanted]
    problems = validate(selected, args.bams, COVERAGE_COLUMNS)
    if problems:
        print("INVALID -- %d problem(s) with the selected samples:" % len(problems))
        for problem in problems:
            print("  - %s" % problem)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)
    if args.skip_existing:
        before = len(selected)
        selected = [r for r in selected
                    if not (args.output / ("%s.shared_coverage.h5" % r["sample_id"])).exists()]
        if before != len(selected):
            log("skipping %d sample(s) already built" % (before - len(selected)))

    if args.annotation_cache is None:
        args.annotation_cache = (args.output.parent / ".cache" / "annotation"
                                 / "coverage_annotation.pkl")
    sys.path.insert(0, str(HERE))
    import annotation_cache as ac
    _bundle, reused = ac.load_or_build(args.annotation_cache, args.gtf, args.appris,
                                       args.regions, args.left_span, args.right_span)
    del _bundle
    log("annotation cache %s: %s" % ("reused" if reused else "built", args.annotation_cache))

    log("building %d sample(s) with %d worker(s)" % (len(selected), args.workers))
    started = time.time()
    results = []
    if args.workers <= 1:
        for row in selected:
            log("  %s ..." % row["sample_id"])
            outcome = build_one(row, args)
            results.append(outcome)
            log("  %s %s in %.1f min" % (row["sample_id"],
                                         "OK" if outcome["ok"] else "FAILED",
                                         outcome["seconds"] / 60.0))
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(build_one, row, args): row for row in selected}
            for future in concurrent.futures.as_completed(futures):
                outcome = future.result()
                results.append(outcome)
                log("  %s %s in %.1f min" % (outcome["sample"],
                                             "OK" if outcome["ok"] else "FAILED",
                                             outcome["seconds"] / 60.0))

    failed = [r for r in results if not r["ok"]]
    checksums = write_checksums(args.output, [r["sample"] for r in results if r["ok"]])

    print()
    log("%d succeeded, %d failed, %.1f min total"
        % (len(results) - len(failed), len(failed), (time.time() - started) / 60.0))
    for outcome in failed:
        print("  FAILED %s:\n%s" % (outcome["sample"],
                                    "\n".join("      " + line
                                              for line in outcome["detail"].splitlines())))
    if checksums:
        log("wrote %s" % checksums)
    return 1 if failed else 0

if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    sys.exit(main())
