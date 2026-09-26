#!/usr/bin/env python3
"""Step 1 of the clustering: the gene x five-metric read-count matrix for ONE sample.

Figure 5A already answers "what happens to every read at a gene" for three hand-picked
genes. This runs the SAME classifier over every APPRIS-selected gene and reduces each
gene to five read counts -- the seven route-7 segments of the Figure 5A bar, folded:

    n_shared_genome_unique        r7_shared_unique
    n_shared_genome_multimapped   r7_shared_multi_pp + r7_shared_multi_other
    n_genome_only_unique          r7_gonly_unique_omit + r7_gonly_unique_other
    n_genome_only_multimapped     r7_gonly_multi
    n_transcriptome_only          r7_txonly

The fold keeps WHERE a read went (both routes or one; genome-unique or multimapped) and
drops the mechanism the seven-way split attributes (pseudogene tie, omitted exon). The
count table carries only the five; the seven are checked (they partition the union, and
the fold uses each exactly once) and compared to the shipped Figure 5A counts. The two
mechanism segments the cluster characterization needs are written to a second table,
`<stem>.gene_mechanisms.tsv`:

    n_alt_exon_genome_unique      r7_gonly_unique_omit   genome-only unique reads on an
                                                         exon the selected transcript omits
    n_pseudogene_tie_shared_multi r7_shared_multi_pp     shared multimapped reads with a
                                                         top-score placement on this gene's
                                                         exon and one on a processed
                                                         pseudogene

Nothing here reimplements the classification. The ten-category chain is
`gene_read_partition_lib.classify_union` and the ten-to-seven fold is
`plot_gene_read_partition._route7_segment`, both imported and called verbatim. What IS
replaced is the I/O around them: `compute_partition` reads both BAMs end to end on every
call, which is affordable for three genes and not for twenty thousand. The two passes are
done once by `read_state.py` into an HDF5 store, and each gene is then a slice of it.

The gene-assignment rule is Figure 5A's own: a read joins a gene's union if a top-score
genome placement (its primary, or a secondary with the primary's AS) overlaps the gene's
full multi-isoform span, or if its transcriptome primary is the gene's selected APPRIS
transcript. Genome unique/multimapped is the primary's NH (Figures 3-4). Unions overlap
between genes; the columns do not sum to the library.

    python code/clustering/read_state.py --sample HeLa
    python code/clustering/build_gene_counts.py --sample HeLa

`--verify N` proves the swap: it runs the untouched `compute_partition` on N genes from
the same BAMs and asserts every one of the seven segment counts matches.

Writes `results/clustering/<sample>.<label>.gene_counts.tsv`, one row per gene, sorted by
`n_union` descending then gene symbol -- the row order every later step inherits. The
HeLa table is shipped as `data/clustering/HeLa.post_dedup.gene_counts.tsv` (Figure 6).
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "code" / "common"))
import inputs  # noqa: E402

log = inputs.make_log("clustering/counts")

#: The seven route-7 segment keys, in the order `ROUTE7_SEGMENTS` declares them. Checked
#: against the panel module at run time.
SEGMENT_KEYS = ("r7_shared_unique", "r7_shared_multi_pp", "r7_shared_multi_other",
                "r7_gonly_unique_omit", "r7_gonly_unique_other", "r7_gonly_multi",
                "r7_txonly")

#: The five metrics, each the SUM of named route-7 segments. Named explicitly, never
#: derived by subtraction: `fold_route5` asserts every segment is used exactly once and
#: that the five sum to `n_union`, which is what gives a typo here something to fail on.
ROUTE5_FOLD = (
    ("n_shared_genome_unique", ("r7_shared_unique",)),
    ("n_shared_genome_multimapped", ("r7_shared_multi_pp", "r7_shared_multi_other")),
    ("n_genome_only_unique", ("r7_gonly_unique_omit", "r7_gonly_unique_other")),
    ("n_genome_only_multimapped", ("r7_gonly_multi",)),
    ("n_transcriptome_only", ("r7_txonly",)),
)
METRIC_COLUMNS = tuple(name for name, _segments in ROUTE5_FOLD)

#: The shipped Figure 5A counts. Reproducible only from the post-dedup BAMs they were
#: built from, so the comparison is reported, not enforced.
EXPECT_TABLE = REPO / "data" / "alignment_fate" / "gene_partition_route7.tsv"
CHECK_GENES = ("COMT", "GAPDH", "LRRFIP1")

COLUMNS = ["gene", "gene_id", "transcript_id", "status", "n_union"] + list(METRIC_COLUMNS)

#: The mechanism table: two route-7 segments under their own names, as subsets of the
#: metrics they belong to (genome_only_unique and shared_genome_multimapped).
MECHANISMS = (("n_alt_exon_genome_unique", "r7_gonly_unique_omit"),
              ("n_pseudogene_tie_shared_multi", "r7_shared_multi_pp"))
MECHANISM_COLUMNS = (["gene", "gene_id", "transcript_id", "status", "n_union"]
                     + [name for name, _segment in MECHANISMS])


# ── the chain, imported by path ──────────────────────────────────────────────

def load_partition_lib():
    """`alignment_fate.gene_read_partition_lib` -- the ten-category chain."""
    return inputs.import_from(REPO / "code" / "read_categories", "gene_read_partition_lib",
                              extra=(REPO / "code" / "common",
                                     REPO / "code" / "common" / "ribo_seq_qc"))


def load_panel_fold():
    """`panels.plot_gene_read_partition` -- the ten-to-seven fold."""
    return inputs.import_from(REPO / "code" / "panels", "plot_gene_read_partition")


# ── the gene universe ────────────────────────────────────────────────────────

def gene_universe(annotation, spans, log):
    """(tids, gene_of, dropped) over every gene with a selected APPRIS transcript."""
    gene2tid = annotation["gene2tid"]
    if len(set(gene2tid.values())) != len(gene2tid):
        raise SystemExit("gene2tid is not 1:1; the gene universe would double-count reads")

    tids, gene_of, dropped = [], {}, []
    for gene_id, tid in sorted(gene2tid.items()):
        base = str(gene_id).split(".", 1)[0]
        n_chrom = int(spans["n_chrom"].get(base, 0))
        if n_chrom == 0:
            dropped.append((gene_id, tid, "no exon in the annotation"))
        elif n_chrom != 1:
            # PAR genes: the version-stripped id collapses the X and Y copies, so
            # `gene_locus` raises. Excluded up front rather than mid-run.
            dropped.append((gene_id, tid, "gene_locus spans %d chromosomes" % n_chrom))
        else:
            tids.append(tid)
            gene_of[tid] = gene_id
    log("gene universe: %d genes, %d excluded" % (len(tids), len(dropped)))
    return tids, gene_of, dropped


def transcript_names(log):
    """transcript id -> gene symbol, for display only; never used to group reads."""
    sys.path.insert(0, str(REPO / "code" / "common" / "ribo_seq_qc"))
    import config
    frame = config.load_annotation()
    if "gene_name" not in frame.columns:
        log("annotation carries no gene_name column; falling back to transcript ids")
        return {}
    return (frame[["transcript_id", "gene_name"]].drop_duplicates("transcript_id")
            .set_index("transcript_id")["gene_name"].astype(str).to_dict())


# ── the core: one gene, sliced out of the store, classified by the chain ─────

def classify_gene(partition_lib, fold, libs, annotation, state, tid):
    """(n_union, {segment: count}) for one gene. The chain and the fold, unchanged."""
    locus = partition_lib.gene_locus(annotation, annotation["table"][tid]["gene_id"])
    genome_side = state.gene_side(*locus)
    txome_side = state.transcript_side(tid)

    labels, _detail = partition_lib.classify_union(
        libs, annotation, tid, locus, genome_side, txome_side,
        state.txome_present, state.genome_present,
        state.genome_unique, state.primary, state.records)

    counts = dict.fromkeys(SEGMENT_KEYS, 0)
    for read_id, category in labels.items():
        # The seventh dimension: a genome multimapper is "shared" only if the read itself
        # has a transcriptome primary. Same test the panel makes on the per-read dump.
        segment = fold._route7_segment(category, read_id in state.txome_present)
        counts[segment] += 1
    return int(len(labels)), counts


def classify_all(partition_lib, fold, libs, annotation, state, tids, log, every=500):
    rows, started = [], time.time()
    for index, tid in enumerate(tids, start=1):
        n_union, counts = classify_gene(partition_lib, fold, libs, annotation, state, tid)
        row = {"transcript_id": tid, "n_union": n_union}
        row.update(counts)
        rows.append(row)
        if index % every == 0 or index == len(tids):
            elapsed = time.time() - started
            log("  %d/%d genes, %.0f s (%.1f genes/s)"
                % (index, len(tids), elapsed, index / max(elapsed, 1e-9)))
    return pd.DataFrame(rows, columns=["transcript_id", "n_union"] + list(SEGMENT_KEYS))


# ── checks and the fold ──────────────────────────────────────────────────────

def check_counts(table, label):
    """The assertion with content: the seven segments partition the union.

    `n_union` is the size of the union the chain labelled; the seven counts come from
    folding those labels. Equality means every union read was seen once and assigned
    exactly one segment -- it catches a category falling through the fold, a read counted
    twice, or one gene's reads leaking into another's.
    """
    segments = table[list(SEGMENT_KEYS)]
    if (segments < 0).any().any():
        raise SystemExit("%s: negative segment count" % label)
    bad = table[segments.sum(axis=1) != table["n_union"]]
    if not bad.empty:
        raise SystemExit("%s: %d gene(s) whose segments do not sum to n_union, e.g. %s"
                         % (label, len(bad), bad.head(3)["transcript_id"].tolist()))


def fold_route5(table):
    """Add the five metric columns, each an integer SUM of named route-7 segments.

    Two assertions restate the partition in the arithmetic the output uses: every route-7
    segment is folded into exactly one metric, and the five metrics sum to `n_union` on
    every row. A typo in `ROUTE5_FOLD` (a segment folded twice, or into nothing) fails
    here even though `check_counts` still holds.
    """
    used = [segment for _name, segments in ROUTE5_FOLD for segment in segments]
    if sorted(used) != sorted(SEGMENT_KEYS):
        raise SystemExit("ROUTE5_FOLD does not use every route-7 segment exactly once: %s"
                         % sorted(used))
    out = table.copy()
    for name, segments in ROUTE5_FOLD:
        out[name] = out[list(segments)].sum(axis=1).astype(int)
    bad = out[out[list(METRIC_COLUMNS)].sum(axis=1) != out["n_union"]]
    if not bad.empty:
        raise SystemExit("the five metrics do not sum to n_union for %d gene(s), e.g. %s"
                         % (len(bad), bad.head(3)["transcript_id"].tolist()))
    return out


def verify_against_compute_partition(partition_lib, fold, state, tids, gene_of,
                                     genome_bam, txome_bam, sample, table, log):
    """Run the untouched `compute_partition` on the same BAMs and compare every segment.

    This is the proof that replacing the I/O did not change the answer: same inputs, same
    genes, one path through the store and one through the original reader.
    """
    if not tids:
        return
    log("verifying %d gene(s) against compute_partition on the same BAMs" % len(tids))
    wide, _tidy, dump = partition_lib.compute_partition(
        sample, genome_bam, txome_bam, transcript_ids=list(tids), coverage=None,
        log=lambda m: log("    " + m))

    present = dump["txome_primary_transcript"].fillna("").astype(str) != ""
    segments = [fold._route7_segment(c, p) for c, p in zip(dump["category"], present)]
    reference = (dump.assign(_seg=segments).groupby("transcript_id")["_seg"]
                 .value_counts().unstack(fill_value=0)
                 .reindex(columns=list(SEGMENT_KEYS), fill_value=0))
    unions = wide.set_index("transcript_id")["n_union"]

    failures = []
    indexed = table.set_index("transcript_id")
    for tid in tids:
        if int(indexed.loc[tid, "n_union"]) != int(unions.loc[tid]):
            failures.append("%s n_union %d vs %d"
                            % (tid, indexed.loc[tid, "n_union"], unions.loc[tid]))
        for key in SEGMENT_KEYS:
            want = int(reference.loc[tid, key]) if tid in reference.index else 0
            if int(indexed.loc[tid, key]) != want:
                failures.append("%s %s %d vs %d" % (tid, key, indexed.loc[tid, key], want))
    if failures:
        for line in failures[:20]:
            log("  MISMATCH %s" % line)
        raise SystemExit("the store-backed path disagrees with compute_partition")
    log("all %d gene(s) match compute_partition exactly, on every segment" % len(tids))


# ── reporting ────────────────────────────────────────────────────────────────

def expected_counts(path, sample):
    """The shipped Figure 5A counts: {gene: (n_union, {segment: n_reads})}."""
    frame = pd.read_csv(path, sep="\t")
    frame = frame[frame["sample"].astype(str) == str(sample)]
    return {str(gene): (int(rows["n_union"].iloc[0]),
                        dict(zip(rows["segment_key"], rows["n_reads"].astype(int))))
            for gene, rows in frame.groupby("gene_name")}


def named_gene_report(table, expect, log):
    """Print the five counts for the check genes, with the shipped counts as a reference.

    Advisory, not a gate: the shipped counts come from the post-dedup BAMs, so a run over
    any other BAM set is expected to differ, and differing is the finding, not a failure.
    """
    for gene in CHECK_GENES:
        rows = table[table["gene"] == gene]
        if rows.empty:
            log("%-8s not in this run" % gene)
            continue
        row = rows.iloc[0]
        log("%-8s n_union %7d   %s"
            % (gene, row["n_union"],
               "  ".join("%s %d" % (name.replace("n_", "", 1), row[name])
                         for name in METRIC_COLUMNS)))
        if gene not in expect:
            continue
        want_union, want_counts = expect[gene]
        deltas = [("n_union", int(row["n_union"]), want_union)]
        deltas += [(key, int(row[key]), int(want_counts.get(key, 0)))
                   for key in SEGMENT_KEYS]
        differing = [(name, got, want) for name, got, want in deltas if got != want]
        if not differing:
            log("         matches the shipped Figure 5A counts exactly")
        else:
            log("         differs from the shipped Figure 5A counts (post-dedup) in %d of "
                "%d values:" % (len(differing), len(deltas)))
            for name, got, want in differing:
                log("           %-22s this run %8d   Figure 5A %8d" % (name, got, want))


def summarise(table, log):
    covered = table[table["n_union"] > 0]
    log("genes: %d total, %d with reads, %d with an empty union, %d excluded"
        % (len(table), len(covered), int((table["status"] == "zero_union").sum()),
           int((table["status"] == "excluded").sum())))
    if not len(covered):
        return
    for name in METRIC_COLUMNS:
        share = 100.0 * covered[name].sum() / covered["n_union"].sum()
        zero = 100.0 * (covered[name] == 0).mean()
        log("  %-30s %5.2f %% of all union reads; zero in %5.2f %% of genes"
            % (name, share, zero))
    log("n_union deciles (genes with reads):")
    for quantile in np.arange(0.1, 1.001, 0.1):
        log("  %3d%%  %10.0f" % (round(quantile * 100), covered["n_union"].quantile(quantile)))


# ── CLI ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True, help="RiboFlow sample directory name")
    parser.add_argument("--label", default="post_dedup",
                        help="which read-state store to use; the output stem is "
                             "<sample>.<label>")
    parser.add_argument("--state", help="path to a read_state HDF5 (overrides --label)")
    parser.add_argument("--genome-bam", help="only needed by --verify")
    parser.add_argument("--transcriptome-bam")
    parser.add_argument("--bams")
    parser.add_argument("--gtf")
    parser.add_argument("--appris")
    parser.add_argument("--output", help="default: results/clustering")
    parser.add_argument("--limit", type=int, help="first N genes of the universe")
    parser.add_argument("--genes", default="",
                        help="comma-separated symbols, gene ids or transcript ids")
    parser.add_argument("--verify", type=int, default=0,
                        help="genes to cross-check against compute_partition (slow: it "
                             "re-reads both BAMs once)")
    parser.add_argument("--expect", help="default: data/alignment_fate/gene_partition_route7.tsv")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import read_state as read_state_module

    paths = inputs.resolve_external_inputs(
        args.bams, args.gtf, args.appris, sample=args.sample)
    os.environ["RIBOFLOW_PAPER_GTF"] = str(Path(paths["gtf"]).resolve())
    os.environ["RIBOFLOW_PAPER_APPRIS"] = str(Path(paths["appris"]).resolve())

    root = REPO / "results" / "clustering"
    output = Path(args.output) if args.output else root
    output.mkdir(parents=True, exist_ok=True)
    state_path = (Path(args.state) if args.state
                  else read_state_module.default_path(root, args.sample, args.label))
    if not state_path.exists():
        raise SystemExit(
            "no read-state store at %s. Build it first:\n"
            "    python code/clustering/read_state.py --sample %s --label %s"
            % (state_path, args.sample, args.label))

    partition_lib = load_partition_lib()
    fold = load_panel_fold()
    if [key for key, *_rest in fold.ROUTE7_SEGMENTS] != list(SEGMENT_KEYS):
        raise SystemExit("ROUTE7_SEGMENTS changed; this script's segment map is stale")

    log("loading the annotation the chain uses")
    libs = partition_lib.load_libraries()
    annotation = partition_lib.load_annotation(libs)
    names = transcript_names(log)
    tids, gene_of, dropped = gene_universe(annotation, annotation["spans"], log)

    wanted = [g.strip() for g in args.genes.split(",") if g.strip()]
    if wanted:
        by_key = {}
        for tid in tids:
            by_key.setdefault(names.get(tid, ""), tid)
            by_key.setdefault(tid, tid)
            by_key.setdefault(str(gene_of[tid]).split(".", 1)[0], tid)
        chosen = [by_key.get(item) for item in wanted]
        unknown = [item for item, tid in zip(wanted, chosen) if tid is None]
        if unknown:
            raise SystemExit("not in the gene universe: %s" % ", ".join(unknown))
        selected = chosen
        log("restricted to %d requested gene(s)" % len(selected))
    else:
        selected = tids[:args.limit] if args.limit else tids

    with read_state_module.ReadState(state_path) as state:
        log("read state %s: %d reads, %d alignments"
            % (state_path.name, state.n_reads,
               int(state._handle.attrs["n_alignments"])))
        counts = classify_all(partition_lib, fold, libs, annotation, state,
                              selected, log)
        check_counts(counts, "whole run")
        table = fold_route5(counts)
        table["gene"] = [names.get(t, t) for t in table["transcript_id"]]
        table["gene_id"] = [gene_of[t] for t in table["transcript_id"]]
        table["status"] = np.where(table["n_union"] > 0, "ok", "zero_union")

        if args.verify:
            genome_bam = Path(args.genome_bam or paths["ribo_genome"])
            txome_bam = Path(args.transcriptome_bam or paths["ribo_txome"])
            verify_against_compute_partition(
                partition_lib, fold, state, selected[:args.verify], gene_of,
                genome_bam, txome_bam, args.sample, table, log)

    if dropped and not wanted:
        table = pd.concat([table, pd.DataFrame(
            [{"gene": names.get(tid, tid), "gene_id": gene_id, "transcript_id": tid,
              "status": "excluded", "n_union": 0} for gene_id, tid, _r in dropped])],
            ignore_index=True, sort=False)

    count_columns = ["n_union"] + list(METRIC_COLUMNS) + list(SEGMENT_KEYS)
    table[count_columns] = table[count_columns].fillna(0).astype(int)
    table = table.sort_values(["n_union", "gene"], ascending=[False, True])
    stem = "%s.%s" % (args.sample, args.label)
    destination = output / ("%s.gene_counts.tsv" % stem)
    table[COLUMNS].to_csv(destination, sep="\t", index=False, lineterminator="\n")
    log("wrote %s (%d genes)" % (destination, len(table)))
    for name, segment in MECHANISMS:
        table[name] = table[segment]
    mechanisms = output / ("%s.gene_mechanisms.tsv" % stem)
    table[MECHANISM_COLUMNS].to_csv(mechanisms, sep="\t", index=False, lineterminator="\n")
    log("wrote %s" % mechanisms)

    expect_path = Path(args.expect) if args.expect else EXPECT_TABLE
    named_gene_report(
        table, expected_counts(expect_path, args.sample) if expect_path.exists() else {},
        log)
    summarise(table, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
