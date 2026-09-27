#!/usr/bin/env python3
"""The three per-gene annotation tables behind Figure 6 E-G, from annotation alone.

    python code/read_categories/cluster_validation.py pseudogene_counts \\
        --clusters results/clustering/HeLa.post_dedup.clusters_k4.tsv \\
        --pseudogenes data/clustering/gencode.v34.2wayconspseudos.gtf.gz \\
        --output results/clustering [--gtf G]
    python code/read_categories/cluster_validation.py omitted_sequence   ... [--gtf G]
    python code/read_categories/cluster_validation.py reference_duplication ... [--appris A]

Nothing read-derived is opened: each table annotates the clustered genes from the GTF
(and, for `reference_duplication`, the APPRIS reference itself).

pseudogene_counts (Figure 6E: genes with >= 1 pseudogene). The main GENCODE GTF carries
no parent attribute on a pseudogene; GENCODE's Yale-UCSC 2-way consensus pseudogene GTF
does (`parent_id`, an unversioned ENSG). The consensus file types every record only as
`pseudogene`, so the processed / unprocessed split comes from the main v34 GTF: each
consensus record takes the biotype of the v34 pseudogene gene on the same chromosome and
strand that overlaps it most (`unprocessed` tested before `processed`; anything else
`other`; no overlap `unclassified`). Writes <stem>.pseudogene_counts_genes.tsv.

omitted_sequence (Figure 6F: exonic sequence absent from the transcriptome). Per gene,
omitted_exon_fraction = 1 - exon_nt_selected / exon_nt_gene, where exon_nt_gene is the
merged length of every transcript's exons and exon_nt_selected the selected transcript's;
likewise omitted_cds_fraction over CDS features (NA when the gene has no CDS). Writes
<stem>.omitted_sequence_genes.tsv.

reference_duplication (Figure 6G: duplicated transcriptome exonic sequence). Per
reference entry, the fraction of its exonic bases whose genomic coordinates overlap an
exon of ANOTHER entry on the same chromosome and strand (and, as a companion column, on
either strand). Writes <stem>.reference_duplication_entries.tsv.

`panels/plot_read_category_clusters.py` draws Figure 6 E-G from these tables.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "code" / "common"))
import inputs  # noqa: E402
from intervals import merge, merged_length  # noqa: E402

KEY = ["gene", "gene_id", "transcript_id"]
ATTR = re.compile(r'(\S+) "([^"]*)"')

log = inputs.make_log("clustering/validation")


# ── shared ────────────────────────────────────────────────────────────────────

def gtf_records(gtf, features):
    """The tab-split fields of every `features` line of a (gzipped) GTF."""
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] in features:
                yield f


def resolve_gtf(a):
    gtf = a.gtf or os.environ.get("RIBOFLOW_PAPER_GTF") or inputs._local_config().get("gtf")
    if not gtf or not Path(gtf).exists():
        inputs.die("GTF not configured or missing: --gtf / RIBOFLOW_PAPER_GTF / local.yaml:gtf")
    return gtf


def load_clusters(clusters_path):
    """(stem, cluster table) from a `<stem>.clusters_k<K>.tsv` path."""
    return clusters_path.name.split(".clusters_k")[0], pd.read_csv(clusters_path, sep="\t")


# ── pseudogene_counts (6E) ────────────────────────────────────────────────────

PSEUDOGENE_CATEGORIES = ("processed", "unprocessed", "other", "unclassified")


def _pseudogene_category(gene_type):
    if "unprocessed" in gene_type:
        return "unprocessed"
    if "processed" in gene_type:
        return "processed"
    return "other"


def parse_v34_pseudogenes(gtf):
    """Per (chrom, strand): (starts, ends, categories) of every v34 pseudogene gene."""
    by = defaultdict(list)
    for f in gtf_records(gtf, ("gene",)):
        gene_type = dict(ATTR.findall(f[8]))["gene_type"]
        if "pseudogene" in gene_type:
            by[(f[0], f[6])].append((int(f[3]) - 1, int(f[4]), _pseudogene_category(gene_type)))
    return {key: (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]),
                  [r[2] for r in rows]) for key, rows in by.items()}


def parse_consensus(path):
    """One row per 2-way consensus pseudogene: chrom, strand, start, end, id, parent ENSG."""
    rows = []
    for f in gtf_records(path, ("transcript",)):
        attrs = dict(ATTR.findall(f[8]))
        if "parent_id" not in attrs:
            raise SystemExit("%s: %s has no parent_id" % (path, attrs.get("gene_id")))
        rows.append({"chrom": f[0], "strand": f[6], "start": int(f[3]) - 1, "end": int(f[4]),
                     "pseudogene_id": attrs["gene_id"], "parent_id": attrs["parent_id"]})
    return pd.DataFrame(rows)


def classify_consensus(consensus, v34):
    """The category of each consensus record: biotype of the v34 pseudogene gene on the same
    chromosome and strand with the largest overlap, or `unclassified`."""
    out = []
    for r in consensus.itertuples(index=False):
        if (r.chrom, r.strand) not in v34:
            out.append("unclassified")
            continue
        starts, ends, cats = v34[(r.chrom, r.strand)]
        overlap = np.minimum(ends, r.end) - np.maximum(starts, r.start)
        if overlap.max() <= 0:
            out.append("unclassified")
        else:
            out.append(cats[int(overlap.argmax())])
    return out


def pseudogene_counts(a, stem, clusters, gtf):
    v34 = parse_v34_pseudogenes(gtf)
    log("GTF: %d pseudogene genes" % sum(len(v[0]) for v in v34.values()))
    consensus = parse_consensus(a.pseudogenes)
    consensus["category"] = classify_consensus(consensus, v34)
    tally = consensus["category"].value_counts()
    log("consensus: %d pseudogenes, %d parents; %s" % (
        len(consensus), consensus["parent_id"].nunique(),
        ", ".join("%s %d" % (c, tally.get(c, 0)) for c in PSEUDOGENE_CATEGORIES)))

    counts = Counter(zip(consensus["parent_id"], consensus["category"]))
    rows = []
    for gene, gene_id, transcript_id, cluster in clusters[KEY + ["cluster"]].itertuples(index=False):
        parent = str(gene_id).split(".")[0]
        row = {"gene": gene, "gene_id": gene_id, "transcript_id": transcript_id,
               "cluster": int(cluster)}
        for cat in PSEUDOGENE_CATEGORIES:
            row["n_%s" % cat] = counts.get((parent, cat), 0)
        row["n_total"] = sum(row["n_%s" % cat] for cat in PSEUDOGENE_CATEGORIES)
        row["has_processed"] = row["n_processed"] > 0
        row["has_unprocessed"] = row["n_unprocessed"] > 0
        row["has_any"] = row["n_total"] > 0
        rows.append(row)
    table = pd.DataFrame(rows)

    genes_tsv = a.output / ("%s.pseudogene_counts_genes.tsv" % stem)
    table.to_csv(genes_tsv, sep="\t", index=False)
    log("wrote %s" % genes_tsv)
    log("genes with >= 1 pseudogene, by cluster: %s" % "  ".join(
        "%s=%.1f%%" % (c, 100.0 * table.loc[table["cluster"] == c, "has_any"].mean())
        for c in sorted(table["cluster"].unique())))


# ── omitted_sequence (6F) ─────────────────────────────────────────────────────

OMITTED_FEATURES = ("exon", "CDS")


def parse_gene_features(gtf, wanted_genes):
    """Per gene: {feature: {transcript_id: [(start, end), ...]}} and transcript types.

    Coordinates become 0-based half-open. Only the wanted genes are kept.
    """
    per_gene = {g: {f: defaultdict(list) for f in OMITTED_FEATURES} for g in wanted_genes}
    tx_type = {}
    for fields in gtf_records(gtf, OMITTED_FEATURES):
        attrs = dict(ATTR.findall(fields[8]))
        gene_id = attrs["gene_id"]
        if gene_id not in per_gene:
            continue
        tx = attrs["transcript_id"]
        tx_type[tx] = attrs.get("transcript_type", "")
        per_gene[gene_id][fields[2]][tx].append((int(fields[3]) - 1, int(fields[4])))
    return per_gene, tx_type


def omitted_sequence(a, stem, clusters, gtf):
    per_gene, tx_type = parse_gene_features(gtf, set(clusters["gene_id"]))
    log("GTF: %d transcripts of the clustered genes" % len(tx_type))
    rows = []
    for gene, gene_id, transcript_id, cluster in clusters[KEY + ["cluster"]].itertuples(index=False):
        feats = per_gene[gene_id]
        row = {"gene": gene, "gene_id": gene_id, "transcript_id": transcript_id,
               "cluster": int(cluster),
               "selected_transcript_type": tx_type.get(transcript_id, "")}
        if transcript_id not in feats["exon"]:
            raise SystemExit("%s: selected transcript %s has no exon in the GTF" % (gene, transcript_id))
        for feature, tag in (("exon", "exon"), ("CDS", "cds")):
            by_tx = feats[feature]
            gene_union = [iv for ivs in by_tx.values() for iv in ivs]
            n_gene = merged_length(gene_union)
            n_sel = merged_length(by_tx.get(transcript_id, []))
            if n_sel > n_gene:
                raise SystemExit("%s: selected %s length exceeds the gene union" % (gene, feature))
            row["%s_nt_gene" % tag] = n_gene
            row["%s_nt_selected" % tag] = n_sel
            row["n_transcripts_%s" % tag] = len(by_tx)
            row["omitted_%s_fraction" % tag] = (1.0 - n_sel / n_gene) if n_gene else np.nan
        rows.append(row)
    table = pd.DataFrame(rows)

    genes_tsv = a.output / ("%s.omitted_sequence_genes.tsv" % stem)
    table.to_csv(genes_tsv, sep="\t", index=False, float_format="%.6f", na_rep="NA")
    log("wrote %s" % genes_tsv)
    log("median omitted exon fraction, by cluster: %s" % "  ".join(
        "%s=%.2f" % (c, table.loc[table["cluster"] == c, "omitted_exon_fraction"].median())
        for c in sorted(table["cluster"].unique())))


# ── reference_duplication (6G) ────────────────────────────────────────────────

def reference_entries(appris_path):
    """transcript id -> gene name, from the APPRIS lengths table's headers."""
    entries = {}
    with open(appris_path) as handle:
        for line in handle:
            fields = line.split("\t")[0].split("|")
            if len(fields) > 6 and fields[0].startswith("ENST"):
                entries[fields[0]] = fields[5]
    if not entries:
        raise SystemExit("%s holds no reference entries" % appris_path)
    return entries


def parse_reference_exons(gtf, tids):
    exons = defaultdict(list)
    where = {}
    for f in gtf_records(gtf, ("exon",)):
        tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        if tid not in tids:
            continue
        exons[tid].append((int(f[3]) - 1, int(f[4])))
        where[tid] = (f[0], f[6])
    return {t: merge(v) for t, v in exons.items()}, where


def duplication(exons, where, strand_aware=True):
    """Per entry: duplicated_nt and the set of other entries sharing a base."""
    groups = defaultdict(list)
    for tid, (chrom, strand) in where.items():
        groups[(chrom, strand if strand_aware else ".")].append(tid)
    dup_nt, partners = {}, defaultdict(set)
    for key, tids in groups.items():
        # sweep  coverage count over sorted boundaries
        events = []
        for tid in tids:
            for s, e in exons[tid]:
                events.append((s, 1, tid))
                events.append((e, -1, tid))
        events.sort(key=lambda x: (x[0], x[1]))
        active = set()
        last = None
        # per entry duplicated length  accumulate segments where >= 2 active
        for pos, delta, tid in events:
            if last is not None and pos > last and len(active) >= 2:
                length = pos - last
                for t in active:
                    dup_nt[t] = dup_nt.get(t, 0) + length
                    partners[t].update(active - {t})
            if delta == 1:
                active.add(tid)
            else:
                active.discard(tid)
            last = pos
    return dup_nt, partners


def reference_duplication(a, stem, clusters, gtf):
    appris = a.appris or os.environ.get("RIBOFLOW_PAPER_APPRIS") or inputs._local_config().get("appris")
    if not appris or not Path(appris).exists():
        inputs.die("appris not configured or missing")

    entries = reference_entries(appris)
    log("%d reference entries; %d clustered genes" % (len(entries), len(clusters)))
    exons, where = parse_reference_exons(gtf, set(entries))
    missing = [t for t in entries if t not in exons]
    if missing:
        log("note: %d reference entries have no exon in the GTF and are skipped: %s"
            % (len(missing), ", ".join(missing[:5])))
    dup_same, partners_same = duplication(exons, where, strand_aware=True)
    dup_any, _partners_any = duplication(exons, where, strand_aware=False)

    cluster_of = dict(zip(clusters["transcript_id"], clusters["cluster"]))
    rows = []
    for tid in exons:
        total = sum(e - s for s, e in exons[tid])
        d = dup_same.get(tid, 0)
        d_any = dup_any.get(tid, 0)
        rows.append({"gene": entries[tid], "transcript_id": tid, "chrom": where[tid][0],
                     "strand": where[tid][1], "cluster": cluster_of.get(tid, np.nan),
                     "exonic_nt": total, "duplicated_nt": d, "duplicated_fraction": d / total,
                     "n_duplicating_entries": len(partners_same.get(tid, ())),
                     "duplicating_entries": ";".join(sorted(partners_same.get(tid, ()))),
                     "duplicating_genes": ";".join(entries[t] for t in sorted(partners_same.get(tid, ()))),
                     "duplicated_nt_any_strand": d_any, "duplicated_fraction_any_strand": d_any / total})
    table = pd.DataFrame(rows).sort_values(["duplicated_fraction", "gene"], ascending=[False, True])
    not_found = [t for t in clusters["transcript_id"] if t not in exons]
    if not_found:
        log("note: %d clustered transcripts are not reference entries with exons" % len(not_found))

    entries_tsv = a.output / ("%s.reference_duplication_entries.tsv" % stem)
    table.to_csv(entries_tsv, sep="\t", index=False, float_format="%.6f", na_rep="NA")
    log("wrote %s" % entries_tsv)


# ── main ──────────────────────────────────────────────────────────────────────

TABLES = {"pseudogene_counts": pseudogene_counts,
          "omitted_sequence": omitted_sequence,
          "reference_duplication": reference_duplication}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("table", choices=sorted(TABLES), help="which per-gene table to build")
    p.add_argument("--clusters", required=True, type=Path, help="<stem>.clusters_k<K>.tsv")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--gtf", type=Path, help="GENCODE v34 GTF (default: env / config/local.yaml)")
    p.add_argument("--pseudogenes", type=Path,
                   help="GENCODE v34 2-way consensus pseudogene GTF (pseudogene_counts)")
    p.add_argument("--appris", type=Path, help="reference lengths TSV (reference_duplication)")
    a = p.parse_args(argv)
    gtf = resolve_gtf(a)
    if a.table == "pseudogene_counts":
        if not a.pseudogenes or not a.pseudogenes.exists():
            inputs.die("--pseudogenes not found: %s" % a.pseudogenes)

    a.output.mkdir(parents=True, exist_ok=True)
    stem, clusters = load_clusters(a.clusters)
    log("%s: %d clustered genes" % (stem, len(clusters)))
    TABLES[a.table](a, stem, clusters, gtf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
