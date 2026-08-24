#!/usr/bin/env python3
"""Rebuild every analysis table in this repository from indexed BAMs.

Stages write under --output (default results/); `--into-data` copies them over data/.
The two `te_*` R stages need `Rscript` (base R, no packages).
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"

SAMPLES_CSV = REPO / "supporting_information" / "S1_Table" / "samples.csv"
DEFAULT_MANIFEST = REPO / "config" / "cohort_manifest.tsv"

EXAMPLE_SAMPLE = "HeLa"
EXAMPLE_GSM = "GSM2100602"
#: The genes Figure 6A partitions (gene IDs resolve through the annotation cache).
PARTITION_GENES = ("ENSG00000093010", "ENSG00000111640", "ENSG00000124831")   # COMT, GAPDH, LRRFIP1
LOCUS_GENE = "LRRFIP1"

def shipped_for(relative):
    return REPO / "data" / relative

BAM_TEMPLATES = {
    "ribo_genome_bam": "{s}/genome/alignment_ribo/merged/{s}.post_dedup.bam",
    "ribo_txome_bam": "{s}/transcriptome/alignment_ribo/merged/{s}.transcriptome.post_dedup.bam",
    "rna_genome_bam": "{s}/rnaseq/genome/alignment_ribo/merged/{s}.rnaseq.post_dedup.bam",
    "rna_txome_bam": ("{s}/rnaseq/transcriptome/alignment_ribo/merged/"
                      "{s}.rnaseq.transcriptome.post_dedup.bam"),
}

def log(msg):
    print("[make_tables] %s" % msg, flush=True)

def prepare_environment(args, create_dirs=True):
    """Point the drivers at the output root. `create_dirs=False` for read-only modes."""
    out = Path(args.output).resolve()
    os.environ["RIBOFLOW_PAPER_BAMS"] = str(args.bams)
    os.environ["RIBOFLOW_PAPER_OUT"] = str(out)
    os.environ["RIBOFLOW_PAPER_QC_OUT"] = str(out / "ribo_seq_qc" / "genome")
    os.environ["RIBOFLOW_PAPER_QC_TX_OUT"] = str(out / "ribo_seq_qc" / "transcriptome")
    if args.gtf:
        os.environ["RIBOFLOW_PAPER_GTF"] = str(Path(args.gtf).resolve())
    if args.appris:
        os.environ["RIBOFLOW_PAPER_APPRIS"] = str(Path(args.appris).resolve())
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("RIBOFLOW_PAPER_SAMPLES_CSV", str(SAMPLES_CSV))

    if create_dirs:
        for directory in (out / "ribo_seq_qc" / "genome",
                          out / "ribo_seq_qc" / "transcriptome",
                          out / "annotation"):
            directory.mkdir(parents=True, exist_ok=True)

    shared = [str(CODE / "common"), str(CODE / "common" / "ribo_seq_qc")]
    existing = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = os.pathsep.join(shared + ([existing] if existing else []))
    for entry in shared:
        if entry in sys.path:
            sys.path.remove(entry)
        sys.path.insert(0, entry)
    return out

def read_manifest(path):
    """sample_id -> row; {} when there is no manifest (the BAM templates suffice)."""
    path = Path(path)
    if not path.exists():
        return {}
    with open(path) as handle:
        return {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

def resolve_bam(manifest, sample, column, bams_root):
    """A BAM path: the manifest's entry when it has one, else the declared template.
    Relative manifest paths resolve against --bams."""
    row = manifest.get(sample) or {}
    value = row.get(column)
    if value:
        path = Path(value)
        return path if path.is_absolute() else Path(bams_root) / path
    return Path(bams_root) / BAM_TEMPLATES[column].format(s=sample)

def sh(cmd, cwd=None):
    log("  $ " + " ".join(str(c) for c in cmd))
    started = time.time()
    result = subprocess.run([str(c) for c in cmd], cwd=str(cwd) if cwd else None)
    log("    -> exit %d in %.1f min" % (result.returncode, (time.time() - started) / 60))
    return result.returncode

def stage_annotation(samples, args):
    """Build the shared annotation caches once; a missing GTF/APPRIS fails here, clearly."""
    import config
    try:
        gtf, appris = config.gtf_path(), config.appris_path()
    except config.AnnotationError as exc:
        log("  " + str(exc).replace("\n", "\n  "))
        return 1
    for label, path in (("GTF", gtf), ("APPRIS", appris)):
        if not os.path.exists(path):
            log("  MISSING %s: %s" % (label, path))
            return 1
    log("  GTF    %s" % gtf)
    log("  APPRIS %s" % appris)
    config.load_annotation()
    config.load_appris_meta()
    log("  annotation cache ready at %s" % config.cache_dir())
    return 0

def stage_qc(samples, args):
    qc = CODE / "ribo_seq_qc"
    selection = ["--samples", ",".join(samples)] if samples else []
    plots = ["--plots"] if args.qc_plots else []
    code = sh([sys.executable, qc / "run_pipeline.py",
               "--bam-dir", args.bams, "--steps", "qc,cds_frame",
               "--bam-glob", "*/genome/alignment_ribo/merged/*.post_dedup.bam"]
              + selection + plots)
    code |= sh([sys.executable, qc / "run_transcriptome_qc.py",
                "--bam-dir", args.bams] + selection + plots)
    return code

def stage_offsets(samples, args):
    """Record which detector branch produced each `psite_offset`; consumed by no stage."""
    return sh([sys.executable, CODE / "ribo_seq_qc" / "determine_offset_method.py",
               "--workers", str(args.workers)])

def stage_orf_catalog(samples, args):
    reference = samples[0] if samples else EXAMPLE_SAMPLE
    return sh([sys.executable, CODE / "common" / "build_orf_catalog.py",
               "--txome-bam", args.bam_for(reference, "ribo_txome_bam"),
               "--out-dir", args.out / "annotation"])

def stage_coverage(samples, args):
    """The shared-coordinate coverage HDF5, one per sample; a durable product, never deleted."""
    if not (args.gtf and args.appris):
        log("  the coverage stage needs --gtf and --appris")
        return 1
    command = [sys.executable, CODE / "coverage" / "build_cohort_coverage.py",
               "--manifest", args.manifest, "--bams", args.bams,
               "--gtf", args.gtf, "--appris", args.appris,
               "--qc-genome",
               args.out / "ribo_seq_qc" / "genome" / "tables" / "readlen_window_qc.csv",
               "--qc-txome",
               args.out / "ribo_seq_qc" / "transcriptome" / "tables" / "readlen_window_qc.csv",
               "--output", args.out / "coverage",
               "--workers", str(min(args.workers, 2))]
    command += ["--samples", ",".join(samples)] if samples else ["--all"]
    if args.regions:
        command += ["--regions", args.regions]
    if args.skip_existing:
        command.append("--skip-existing")
    return sh(command)

def stage_concordance(samples, args):
    """The four concordance tables, computed from the HDF5 cohort. No BAM is opened."""
    command = [sys.executable, CODE / "coverage" / "compute_coverage_concordance.py",
               "--coverage", args.out / "coverage",
               "--output", args.out / "coverage" / "concordance",
               "--gzip-per-transcript"]
    if samples:
        command += ["--samples", ",".join(samples)]
    return sh(command)

def _rscript():
    import shutil
    path = shutil.which("Rscript")
    if not path:
        log("  Rscript is not on PATH; the te_normalize and te_stats stages need base R")
    return path

def _qc_tables(args):
    qc = args.out / "ribo_seq_qc"
    return (qc / "genome" / "tables" / "readlen_window_qc.csv",
            qc / "transcriptome" / "tables" / "readlen_window_qc.csv")

def stage_te_counts(samples, args):
    """Four transcripts x samples CDS count matrices, one count program per sample."""
    if not (args.gtf and args.appris):
        log("  the te_counts stage needs --gtf and --appris")
        return 1
    qc_genome, qc_txome = _qc_tables(args)
    command = [sys.executable, CODE / "ribo_rna" / "build_count_matrices.py",
               "--bams", args.bams, "--gtf", args.gtf, "--appris", args.appris,
               "--qc-genome", qc_genome, "--qc-txome", qc_txome,
               "--manifest", args.manifest, "--output", args.out / "ribo_rna",
               "--workers", str(min(args.workers, 2))]
    if args.regions:
        command += ["--regions", args.regions]
    if samples:
        command += ["--samples", ",".join(samples)]
    return sh(command)

def stage_te_normalize(samples, args):
    rscript = _rscript()
    if not rscript:
        return 1
    return sh([rscript, CODE / "te_route" / "normalization.R",
               "--counts", args.out / "ribo_rna" / "counts",
               "--output", args.out / "te_route" / "normalized"])

def stage_te_stats(samples, args):
    rscript = _rscript()
    if not rscript:
        return 1
    return sh([rscript, CODE / "te_route" / "te_statistics.R",
               "--normalized", args.out / "te_route" / "normalized",
               "--orf-catalog", args.out / "annotation" / "orf_catalog.tsv",
               "--output", args.out / "te_route" / "tables"])

def stage_gene_partition(samples, args):
    """Figure 6A: the per-read gene partition dump, folded to the seven-segment table."""
    coverage = args.out / "coverage" / ("%s.shared_coverage.h5" % EXAMPLE_SAMPLE)
    command = [sys.executable, CODE / "alignment_fate" / "build_gene_read_partition.py",
               "--sample", EXAMPLE_SAMPLE,
               "--genome-bam", args.bam_for(EXAMPLE_SAMPLE, "ribo_genome_bam"),
               "--transcriptome-bam", args.bam_for(EXAMPLE_SAMPLE, "ribo_txome_bam"),
               "--gene-id", ",".join(PARTITION_GENES),
               "--output", args.out / "alignment_fate", "--dump-reads"]
    if coverage.exists():
        command += ["--coverage", coverage]
    code = sh(command)
    if code:
        return code
    return sh([sys.executable, CODE / "alignment_fate" / "build_gene_partition_data.py",
               "--reads", args.out / "alignment_fate"
               / ("%s.gene_read_partition_reads.tsv" % EXAMPLE_SAMPLE),
               "--sample", EXAMPLE_SAMPLE, "--gsm", EXAMPLE_GSM,
               "--output", args.out / "alignment_fate" / "gene_partition_route7", "--force"])

def stage_locus(samples, args):
    """Figure 6B: the LRRFIP1 locus coverage artifact."""
    if not (args.gtf and args.appris):
        log("  the locus stage needs --gtf and --appris")
        return 1
    qc_genome, qc_txome = _qc_tables(args)
    return sh([sys.executable, CODE / "alignment_fate" / "build_locus_data.py",
               "--gene", LOCUS_GENE, "--sample", EXAMPLE_SAMPLE, "--gsm", EXAMPLE_GSM,
               "--bams", args.bams, "--gtf", args.gtf, "--appris", args.appris,
               "--qc-genome", qc_genome, "--qc-txome", qc_txome,
               "--output", args.out / "alignment_fate" / ("locus_%s" % LOCUS_GENE),
               "--force"])

def _taxonomy_driver(analysis, samples, args):
    """One cohort driver, told which analysis to run.

    Capped at 2 workers regardless of `--workers`: each subprocess peaks near 5 GB.
    """
    selection = ["--samples", ",".join(samples)] if samples else []
    return sh([sys.executable, CODE / "read_taxonomy" / "run_read_taxonomy.py", analysis,
               "--workers", str(min(args.workers, 2))] + selection)

def stage_taxonomy(samples, args):
    return _taxonomy_driver("taxonomy", samples, args)

def stage_alignment_concordance(samples, args):
    return _taxonomy_driver("alignment_concordance", samples, args)

def stage_reach(samples, args):
    return _taxonomy_driver("reach", samples, args)

def stage_multimap_biotype(samples, args):
    return _taxonomy_driver("tie_biotype", samples, args)

STAGES = [
    ("annotation",   stage_annotation,   (),                       True,  ()),
    ("qc",           stage_qc,           ("annotation",),          True,
     ("ribo_seq_qc/genome/tables/readlen_window_qc.csv",
      "ribo_seq_qc/genome/tables/cds_psite_frame.csv",
      "ribo_seq_qc/transcriptome/tables/readlen_window_qc.csv",
      "ribo_seq_qc/transcriptome/tables/cds_psite_frame.csv")),
    ("offsets",      stage_offsets,      ("qc",),                  True,
     ("ribo_seq_qc/offsets/offset_method_per_length.tsv",)),
    ("orf_catalog",  stage_orf_catalog,  ("annotation",),          True,
     ("annotation/orf_catalog.tsv",)),
    ("coverage",     stage_coverage,     ("annotation", "qc"),     True,  ()),
    ("concordance",  stage_concordance,  ("coverage",),            False,
     ("coverage/concordance/region_concordance_per_sample.tsv",
      "coverage/concordance/region_coverage_per_sample.tsv",
      "coverage/concordance/region_concordance_per_transcript.tsv.gz",
      "coverage/concordance/region_coverage_per_transcript.tsv.gz")),
    ("te_counts",    stage_te_counts,    ("annotation", "qc"),     True,
     ("ribo_rna/counts/ribo_counts_genome.csv",
      "ribo_rna/counts/ribo_counts_txome.csv",
      "ribo_rna/counts/rna_counts_genome.csv",
      "ribo_rna/counts/rna_counts_txome.csv")),
    ("te_normalize", stage_te_normalize, ("te_counts",),           False,
     ("te_route/normalized/size_factors.csv",)),
    ("te_stats",     stage_te_stats,     ("te_normalize", "orf_catalog"), False,
     ("te_route/tables/per_gene_delta.tsv",
      "te_route/tables/route_correlation.tsv")),
    ("gene_partition", stage_gene_partition, ("annotation",),     True,
     ("alignment_fate/gene_partition_route7.tsv",
      "alignment_fate/gene_partition_route7.json")),
    ("locus",        stage_locus,        ("annotation", "qc"),     True,
     ("alignment_fate/locus_LRRFIP1.npz",
      "alignment_fate/locus_LRRFIP1.json")),
    ("taxonomy",     stage_taxonomy,     ("annotation",),          True,
     ("read_taxonomy/taxonomy/taxonomy_all.tsv",)),
    ("alignment_concordance", stage_alignment_concordance, ("annotation",), True, ()),
    ("reach",        stage_reach,        ("taxonomy", "alignment_concordance"), True,
     ("read_taxonomy/reach/genome_anchored_reach_all.tsv",)),
    ("multimap_biotype", stage_multimap_biotype, ("annotation",),  True,
     ("read_taxonomy/multimap_biotype/multimap_tie_biotype_all.tsv",)),
]

STAGE_STAGING = {
    "qc": ("ribo_seq_qc/genome/tables/_staging",
           "ribo_seq_qc/transcriptome/tables/_staging"),
    "te_counts": ("ribo_rna/_route_scratch",),
    "taxonomy": ("read_taxonomy/taxonomy/_staging",),
    "alignment_concordance": ("read_taxonomy/alignment_concordance/_staging",),
    "reach": ("read_taxonomy/reach/_staging",),
    "multimap_biotype": ("read_taxonomy/multimap_biotype/_staging_tie",),
}

#: Shipped under data/ but built by no stage: third-party inputs, recorded with their source.
EXTERNAL_INPUTS = {
    "te_route/housekeeping/Housekeeping_GenesHuman.csv": "HRT Atlas v1.0",
    "te_route/housekeeping/Housekeeping_TranscriptsHuman.csv": "HRT Atlas v1.0",
}

STAGE_ORDER = [name for name, _run, _needs, _anno, _out in STAGES]
STAGE_RUN = {name: run for name, run, _needs, _anno, _out in STAGES}
NEEDS_ANNOTATION = {name for name, _r, _n, anno, _o in STAGES if anno}
STAGE_OUTPUTS = {name: outs for name, _r, _n, _a, outs in STAGES}
OUTPUTS = [rel for _n, _r, _nd, _a, outs in STAGES for rel in outs]

def prune_staging(stage, out):
    """Remove the per-sample staging directories `stage` owns, once it has succeeded."""
    import shutil
    for relative in STAGE_STAGING.get(stage, ()):
        directory = out / relative
        if not directory.is_dir():
            continue
        n = sum(1 for p in directory.rglob("*") if p.is_file())
        shutil.rmtree(directory)
        log("  pruned %s (%d intermediate file(s))" % (relative, n))

def required_stages(requested):
    """`requested` plus everything it depends on, in dependency order."""
    needed, pending = set(), list(requested)
    depends = {name: set(deps) for name, _r, deps, _a, _o in STAGES}
    while pending:
        name = pending.pop()
        if name in needed:
            continue
        needed.add(name)
        pending.extend(depends.get(name, ()))
    return [name for name in STAGE_ORDER if name in needed]

def do_into_data(out: Path):
    """Copy the regenerated artifacts over the shipped ones, naming each one that changed."""
    import shutil
    same = replaced = 0
    for rel in OUTPUTS:
        source = out / rel
        if not source.exists():
            continue
        destination = shipped_for(rel)
        if destination.exists() and destination.read_bytes() == source.read_bytes():
            same += 1
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        replaced += 1
        log("  REPLACED data/%s (the regenerated bytes differ from the published ones)"
            % rel)
    log("--into-data: %d already identical, %d replaced" % (same, replaced))
    return 0

def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bams", help="RiboFlow output tree")
    parser.add_argument("--gtf", default=None, help="GENCODE annotation GTF")
    parser.add_argument("--appris", default=None, help="APPRIS transcript-lengths TSV")
    parser.add_argument("--regions", default=None,
                        help="APPRIS actual-regions BED; a cross-check for the coverage stage")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST),
                        help="sample manifest (default config/cohort_manifest.tsv)")
    parser.add_argument("--output", default=str(REPO / "results"),
                        help="output root (default results/)")
    parser.add_argument("--samples", default=None,
                        help="comma-separated subset (default: every discovered sample)")
    parser.add_argument("--stages", default=None,
                        help="comma-separated subset of: " + ", ".join(STAGE_ORDER))
    parser.add_argument("--all", action="store_true", help="run every stage")
    parser.add_argument("--workers", type=int, default=2,
                        help="parallel samples; memory-heavy stages cap at 2 regardless")
    parser.add_argument("--skip-existing", action="store_true",
                        help="the coverage stage skips samples whose HDF5 already exists")
    parser.add_argument("--qc-plots", action="store_true",
                        help="also write the per-sample metagene PDFs from the QC stage. "
                             "Off by default: they are diagnostics, and the window and "
                             "offsets they illustrate are in the QC tables.")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="keep the per-sample _staging/ directories. They are merged "
                             "into the master tables and deleted once their stage succeeds.")
    parser.add_argument("--validate", action="store_true",
                        help="report the discovered samples and exit without computing")
    parser.add_argument("--into-data", action="store_true",
                        help="OVERWRITE data/ with the regenerated tables (off by default)")
    return parser

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.bams:
        parser.error("--bams is required")
    bams = Path(args.bams).resolve()
    if not bams.is_dir():
        parser.error("--bams is not a directory: %s" % bams)
    args.bams = str(bams)

    unknown = [s.strip() for s in (args.stages or "").split(",")
               if s.strip() and s.strip() not in STAGE_RUN]
    if unknown:
        parser.error("unknown stage(s): %s -- choose from %s"
                     % (", ".join(unknown), ", ".join(STAGE_ORDER)))
    if not args.stages and not args.all and not args.validate:
        parser.error("choose --all or --stages (see --help)")

    if sys.version_info[:2] != (3, 9):
        log("WARNING: running Python %d.%d; this pipeline is developed on 3.9"
            % sys.version_info[:2])

    args.out = prepare_environment(args, create_dirs=not args.validate)
    manifest = read_manifest(args.manifest)
    args.bam_for = lambda sample, column: resolve_bam(manifest, sample, column, args.bams)

    import bam_inputs

    found = bam_inputs.discover_samples()
    log("BAMs   = %s" % bams)
    log("output = %s" % args.out)
    log("discovered %d sample(s) with both a genome and a transcriptome BAM" % len(found))

    if args.samples:
        wanted = [s.strip() for s in args.samples.split(",") if s.strip()]
        missing = [s for s in wanted if s not in found]
        if missing:
            log("NOT FOUND in the BAM tree: %s" % ", ".join(missing))
            return 1
        samples = wanted
    else:
        samples = found

    if args.all or not args.stages:
        selected = STAGE_ORDER
    else:
        requested = {x.strip() for x in args.stages.split(",") if x.strip()}
        selected = []
        for stage in required_stages(requested):
            if stage in requested:
                selected.append(stage)
                continue
            outputs = STAGE_OUTPUTS[stage]
            if outputs and all((args.out / rel).exists() for rel in outputs):
                log("dependency %s satisfied by existing output" % stage)
            else:
                log("adding dependency %s (%s)"
                    % (stage, "no declared output to check" if not outputs
                       else "its output is missing"))
                selected.append(stage)

    if args.validate:
        for sample in samples:
            print("  %-24s genome=%s txome=%s"
                  % (sample,
                     bam_inputs.genome_bam(sample).exists(),
                     bam_inputs.txome_bam(sample).exists()))
        print()
        print("%d sample(s) usable. Stages that would run: %s"
              % (len(samples), ", ".join(selected)))
        print("Outputs would go to %s (data/ untouched)." % args.out)
        if set(selected) & NEEDS_ANNOTATION and not (args.gtf and args.appris):
            print("NOTE: those stages need --gtf and --appris, which were not given.")
        return 0

    log("samples: %s" % (", ".join(samples) if len(samples) <= 6
                         else "%d samples" % len(samples)))
    log("stages : %s" % ", ".join(selected))
    results = []
    for stage in selected:
        log("-- stage %s" % stage)
        started = time.time()
        code = STAGE_RUN[stage](samples, args)
        results.append((stage, code, (time.time() - started) / 60))
        if code == 0 and not args.keep_intermediates:
            prune_staging(stage, args.out)
        if code:
            log("stage %s FAILED (exit %d) -- stopping; later stages depend on it"
                % (stage, code))
            break

    print()
    for stage, code, minutes in results:
        print("%-22s %s  %.1f min" % (stage, "OK  " if code == 0 else "FAIL", minutes))
    failed = [s for s, code, _ in results if code]

    if args.into_data and not failed:
        do_into_data(args.out)
    elif args.into_data:
        log("refusing --into-data: %d stage(s) failed" % len(failed))

    print()
    print("Regenerated tables are in %s. The shipped tables in data/ were NOT modified%s."
          % (args.out, " (use --into-data to replace them)" if not args.into_data else ""))
    return 1 if failed else 0

if __name__ == "__main__":
    sys.exit(main())
