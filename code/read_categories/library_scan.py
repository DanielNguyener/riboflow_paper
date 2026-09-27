#!/usr/bin/env python3
"""The per-library category tables of Figure 4, over the whole cohort.

    python code/read_categories/library_scan.py taxonomy     # 4A/4B: the five categories
    python code/read_categories/library_scan.py tie_biotype  # 4C: pseudogene ties
    python code/read_categories/library_scan.py reach        # 4D: omitted-exon overlap

Each analysis scans the two post-dedup BAMs per library (transcriptome presence = any
primary; see categories.py), stages one TSV per library, and concatenates them into its
master table. `--sample X` runs one library in this process; the cohort driver spawns one
such worker per library.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reach_lib as rl  # noqa: E402
import reference_lib as cl  # noqa: E402
import taxonomy_lib as tl  # noqa: E402
import tie_biotype_lib as tie  # noqa: E402
fc = tl.fc

DEFAULT_WORKERS = 2

def _out(*parts):
    return fc.output_root().joinpath("read_taxonomy", *parts)

ANALYSES = {
    "taxonomy": {
        "staging": _out("taxonomy", "_staging"),
        "master": _out("taxonomy", "taxonomy_all.tsv"),
    },
    "reach": {
        "staging": _out("reach", "_staging"),
        "master": _out("reach", "genome_anchored_reach_all.tsv"),
    },
    "tie_biotype": {
        "staging": _out("multimap_biotype", "_staging_tie"),
        "master": _out("multimap_biotype", "multimap_tie_biotype_all.tsv"),
    },
}

TAXONOMY_TSV = ANALYSES["taxonomy"]["master"]
CATS = ["cross_pc_pp", "cross_pp_pc", "same_pc_pc", "same_pp_pp"]


# ── taxonomy: one row of the five categories per library (4A, 4B) ─────────────

def build_row(sample, counts, n_universe):
    row = {"sample": sample, "n_universe": n_universe}
    for g, t in tl.CELLS:
        n = counts[(g, t)]
        key = tl.cell_key(g, t)
        row[f"n_{key}"] = n
        row[f"pct_{key}"] = 100.0 * n / n_universe if n_universe else float("nan")
    for g in tl.GENOME_STATES:
        row[f"n_genome_{g}"] = sum(counts[(g, t)] for t in tl.TXOME_STATES if (g, t) in counts)
    for t in tl.TXOME_STATES:
        row[f"n_txome_{t}"] = sum(counts[(g, t)] for g in tl.GENOME_STATES if (g, t) in counts)
    return row

def taxonomy_sample(sample, log=print):
    counts, n_universe = tl.classify_sample(sample)
    row = build_row(sample, counts, n_universe)
    _stage("taxonomy", sample, row)
    log(f"[{sample}] n_universe={n_universe:,}  " + "  ".join(
        "%s=%.2f%%" % (tl.cell_key(g, t), row["pct_" + tl.cell_key(g, t)])
        for g, t in tl.CELLS))
    return row


# ── tie_biotype: protein-coding-pseudogene ties per library (4C) ──────────────

def tie_sample(sample, log=print):
    exon_pr = cl.load_exon_gene_pr()
    gene_pr = cl.gene_body_pr()

    log(f"[{sample}] reading txome BAM (present qname set)...")
    t_all = tie.tl.txome_present_qnames(fc.txome_bam(sample))
    log(f"[{sample}] n_txome_present={len(t_all):,}; enumerating genome multimapper loci...")
    recs = tie.read_genome_multi_records_flagged(fc.genome_bam(sample), t_all)
    n_reads = len(recs)

    if TAXONOMY_TSV.exists():
        tr = pd.read_csv(TAXONOMY_TSV, sep="\t")
        tr = tr[tr["sample"] == sample]
        if len(tr):
            exp = int(tr.iloc[0]["n_gM_tP"])
            log(f"[{sample}] taxonomy check dark-green reads {n_reads} vs {exp} -> "
                f"{'OK' if n_reads == exp else 'MISMATCH'}")
            assert n_reads == exp, f"[{sample}] dark-green population disagrees with taxonomy_all.tsv"

    counts, n_reads = tie.categorize(recs, exon_pr, gene_pr)
    row = {"sample": sample, "n_reads": n_reads, **{c: counts[c] for c in CATS},
           "n_qualifying": counts["n_qualifying"]}
    for c in CATS + ["n_qualifying"]:
        row[f"pct_{c}"] = 100.0 * counts.get(c, row.get(c)) / n_reads if n_reads else float("nan")
    row["pct_n_qualifying"] = 100.0 * counts["n_qualifying"] / n_reads if n_reads else float("nan")

    log(f"[{sample}] n_reads={n_reads:,} qualifying={counts['n_qualifying']:,} "
        f"({row['pct_n_qualifying']:.1f}%): "
        f"cross_pc_pp {counts['cross_pc_pp']:,}, cross_pp_pc {counts['cross_pp_pc']:,}, "
        f"same_pc_pc {counts['same_pc_pc']:,}, same_pp_pp {counts['same_pp_pp']:,}")

    _stage("tie_biotype", sample, row)
    return row


# ── reach: omitted alternative-exon overlap per library (4D) ──────────────────

def reach_sample(sample, log=print):
    tax_row = pd.read_csv(TAXONOMY_TSV, sep="\t").set_index("sample").loc[sample]

    n_genome_unique = int(tax_row["n_genome_unique"])
    n_gUtP = int(tax_row["n_gU_tP"])
    n_gUtA = int(tax_row["n_gU_tA"])
    assert n_gUtP + n_gUtA == n_genome_unique, "gU partition mismatch vs taxonomy_all.tsv"
    log(f"[{sample}] N_G={n_genome_unique} gU_tP={n_gUtP} gU_tA={n_gUtA}")

    log(f"[{sample}] recomputing gU_tA qname set...")
    genome_all, genome_uniq = rl.tl.genome_status_sets(fc.genome_bam(sample))
    txome_all = rl.tl.txome_present_qnames(fc.txome_bam(sample))
    gUtA_qnames = genome_uniq - txome_all
    assert len(gUtA_qnames) == n_gUtA, (
        f"recomputed gU_tA ({len(gUtA_qnames)}) != taxonomy_all.tsv ({n_gUtA})")

    log(f"[{sample}] reading genome blocks for {len(gUtA_qnames)} gU_tA reads...")
    genome_blocks = rl.read_genome_blocks(fc.genome_bam(sample), gUtA_qnames)

    transcript_payload = cl.build_transcript_table()
    table = transcript_payload["table"]
    selected_genes = set(v["gene_id"] for v in table.values())
    gene2tid = rl.gene_to_transcript_map(table)
    exon_gene_df = cl.build_exon_gene_table()
    exon_gene_pr = cl.load_exon_gene_pr()
    all_gene_body_pr = _all_gene_body_pr()
    omitted_genes = rl.omitted_pc_genes(exon_gene_df, selected_genes)

    log(f"[{sample}] classifying gU_tA reads...")
    labels = rl.classify_gU_tA(gUtA_qnames, genome_blocks, exon_gene_pr, exon_gene_df,
                               all_gene_body_pr, table, gene2tid, omitted_genes)
    counts = labels.value_counts()

    log(f"[{sample}] testing direct overlap with omitted exonic sequence...")
    omitted_index = rl.build_omitted_exon_index(exon_gene_df, table)
    overlap_qnames = rl.omitted_exon_overlap_qnames(gUtA_qnames, genome_blocks, omitted_index)

    row = {"sample": sample, "n_genome_unique": n_genome_unique,
           "n_gU_tP": n_gUtP, "n_gU_tA": n_gUtA}
    for cat in rl.REACH_CATEGORIES:
        row[f"n_{cat}"] = int(counts.get(cat, 0))
    check_sum = n_gUtP + sum(row[f"n_{c}"] for c in rl.REACH_CATEGORIES)
    assert check_sum == n_genome_unique, f"partition does not sum to N_G ({check_sum} != {n_genome_unique})"
    for cat in rl.REACH_CATEGORIES:
        row[f"pct_{cat}"] = 100.0 * row[f"n_{cat}"] / n_gUtA if n_gUtA else float("nan")
    # Figure 4D: NOT a partition slice. The gene-level alternative-exon test of Figure 5A
    # applied to the whole gU_tA population; each read ID counts at most once.
    row["n_omitted_exon_overlap"] = len(overlap_qnames)
    row["pct_omitted_exon_overlap"] = (
        100.0 * len(overlap_qnames) / n_gUtA if n_gUtA else float("nan"))
    log(f"[{sample}] omitted-exon overlap {len(overlap_qnames)} / {n_gUtA} "
        f"({row['pct_omitted_exon_overlap']:.2f}%)")

    _stage("reach", sample, row)
    return row

def _all_gene_body_pr():
    import pyranges as pr
    df = fc.config.load_all_gene_bodies()
    return pr.PyRanges(df.reset_index(drop=True))

def _all_gene_body_pr():
    import pyranges as pr
    df = fc.config.load_all_gene_bodies()
    return pr.PyRanges(df.reset_index(drop=True))


WORKERS = {"taxonomy": taxonomy_sample, "tie_biotype": tie_sample, "reach": reach_sample}


def _stage(analysis, sample, row):
    staging = ANALYSES[analysis]["staging"]
    staging.mkdir(parents=True, exist_ok=True)
    dest = staging / f"{sample}.tsv"
    pd.DataFrame([row]).to_csv(dest, sep="\t", index=False)
    print(f"[{sample}] wrote {dest}", flush=True)


# ── the cohort driver ──────────────────────────────────────────────────────────

def run_sample(analysis, sample):
    command = [sys.executable, str(HERE / "library_scan.py"), analysis, "--sample", sample]
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
    parser.add_argument("--sample", default=None,
                        help="run ONE library in this process (the worker mode)")
    parser.add_argument("--samples", default=None, help="comma-separated subset")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args(argv)

    if args.sample:
        WORKERS[args.analysis](args.sample)
        return 0

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

