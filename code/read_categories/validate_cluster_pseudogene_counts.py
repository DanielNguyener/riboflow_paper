#!/usr/bin/env python3
"""GENCODE v34 pseudogenes per clustered gene, processed and unprocessed, from annotation alone.

    python code/clustering/validate_cluster_pseudogene_counts.py \\
        --clusters results/clustering/HeLa.post_dedup.clusters_k4.tsv \\
        --pseudogenes data/clustering/gencode.v34.2wayconspseudos.gtf.gz \\
        --output results/clustering [--gtf G]

The main GENCODE GTF carries no parent attribute on a pseudogene. GENCODE's companion
Yale-UCSC 2-way consensus pseudogene GTF does: every `transcript` record names the
protein-coding gene it derives from (`parent_id`, an unversioned ENSG). That is the
association used here. The consensus file types every record only as `pseudogene`, so the
processed / unprocessed split comes from the main v34 GTF: each consensus record is matched
to the v34 pseudogene gene on the same chromosome and strand that overlaps it most, and takes
that gene's biotype:

    unprocessed    gene_type contains "unprocessed" (tested first: it also contains "processed")
    processed      gene_type contains "processed"
    other          any remaining *pseudogene* biotype (unitary, polymorphic, IG/TR, rRNA, ...)
    unclassified   no v34 pseudogene gene overlaps the record on its strand

For every clustered gene: `n_processed`, `n_unprocessed`, `n_other`, `n_unclassified`,
`n_total` (their sum) and the three booleans `has_processed`, `has_unprocessed`, `has_any`.
Counts are zero, not NA, for a gene that no consensus record names as parent.

Outputs, under <output> with the cluster table's stem:

    <stem>.pseudogene_counts_genes.tsv        per gene: cluster and the values above
    <stem>.pseudogene_counts_by_cluster.tsv   per cluster (and all): n, n / % with >= 1 of each
                                              category, odds ratios vs the other clusters
                                              (Woolf CI, Fisher p) for the three booleans, and
                                              the quartiles of each count over genes with >= 1
    <stem>.pseudogene_counts_tests.tsv        Kruskal-Wallis and Holm-adjusted pairwise
                                              Mann-Whitney U on the processed and unprocessed
                                              counts (genes with >= 1)
    <stem>.pseudogene_counts_sources.json     the inputs by path, size and sha256

The genes table is Figure 6E's input (`panels/plot_read_category_clusters.py`).
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import validate_common as vc
from validate_common import ATTR, KEY, inputs

CATEGORIES = ("processed", "unprocessed", "other", "unclassified")
BOOLEANS = (("has_processed", "n_processed"), ("has_unprocessed", "n_unprocessed"),
            ("has_any", "n_total"))
TESTED = ("n_processed", "n_unprocessed")

log = inputs.make_log("clustering/pseudogene_counts")


def category(gene_type):
    if "unprocessed" in gene_type:
        return "unprocessed"
    if "processed" in gene_type:
        return "processed"
    return "other"


# ── inputs ───────────────────────────────────────────────────────────────────────────────
def parse_v34_pseudogenes(gtf):
    """Per (chrom, strand): (starts, ends, categories) of every v34 pseudogene gene."""
    by = defaultdict(list)
    for f in vc.gtf_records(gtf, ("gene",)):
        gene_type = dict(ATTR.findall(f[8]))["gene_type"]
        if "pseudogene" in gene_type:
            by[(f[0], f[6])].append((int(f[3]) - 1, int(f[4]), category(gene_type)))
    return {key: (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]),
                  [r[2] for r in rows]) for key, rows in by.items()}


def parse_consensus(path):
    """One row per 2-way consensus pseudogene: chrom, strand, start, end, id, parent ENSG."""
    rows = []
    for f in vc.gtf_records(path, ("transcript",)):
        attrs = dict(ATTR.findall(f[8]))
        if "parent_id" not in attrs:
            raise SystemExit("%s: %s has no parent_id" % (path, attrs.get("gene_id")))
        rows.append({"chrom": f[0], "strand": f[6], "start": int(f[3]) - 1, "end": int(f[4]),
                     "pseudogene_id": attrs["gene_id"], "parent_id": attrs["parent_id"]})
    return pd.DataFrame(rows)


def classify(consensus, v34):
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


# ── the per-gene table ───────────────────────────────────────────────────────────────────
def annotate(clusters, consensus):
    counts = Counter(zip(consensus["parent_id"], consensus["category"]))
    rows = []
    for gene, gene_id, transcript_id, cluster in clusters[KEY + ["cluster"]].itertuples(index=False):
        parent = str(gene_id).split(".")[0]
        row = {"gene": gene, "gene_id": gene_id, "transcript_id": transcript_id,
               "cluster": int(cluster)}
        for cat in CATEGORIES:
            row["n_%s" % cat] = counts.get((parent, cat), 0)
        row["n_total"] = sum(row["n_%s" % cat] for cat in CATEGORIES)
        for flag, column in BOOLEANS:
            row[flag] = row[column] > 0
        rows.append(row)
    return pd.DataFrame(rows)


# ── summary ──────────────────────────────────────────────────────────────────────────────
def summarize(table):
    clusters = sorted(table["cluster"].unique())
    rows, tests = [], []
    for label in ["all"] + clusters:
        g = table if label == "all" else table[table["cluster"] == label]
        row = {"cluster": str(label), "n_genes": len(g)}
        for flag, column in BOOLEANS:
            row["n_%s" % flag[4:]] = int(g[flag].sum())
            row["pct_%s" % flag[4:]] = 100.0 * g[flag].mean()
            if label != "all":
                orr, lo, hi, p, _cells = vc.odds_ratio(g[flag], table.loc[table["cluster"] != label, flag])
                row.update({"%s_or" % flag: orr, "%s_or_ci_low" % flag: lo,
                            "%s_or_ci_high" % flag: hi, "%s_fisher_p" % flag: p})
        for cat in CATEGORIES:
            column = "n_%s" % cat
            values = g.loc[g[column] > 0, column].to_numpy(float)
            row.update({"%s_%s" % (column, k): v for k, v in vc.quartiles(values).items()})
        rows.append(row)
    for column in TESTED:
        groups = {c: table.loc[(table["cluster"] == c) & (table[column] > 0), column].to_numpy(float)
                  for c in clusters}
        tests += [{"metric": column, **t} for t in vc.cluster_tests(groups, clusters)]
    return pd.DataFrame(rows), pd.DataFrame(tests)


# ── main ─────────────────────────────────────────────────────────────────────────────────
def main(argv=None):
    p = vc.parser(__doc__)
    p.add_argument("--pseudogenes", required=True, type=Path,
                   help="GENCODE v34 2-way consensus pseudogene GTF (gencode.v34.2wayconspseudos.gtf.gz)")
    a = p.parse_args(argv)
    gtf = vc.resolve_gtf(a)
    if not a.pseudogenes.exists():
        inputs.die("--pseudogenes not found: %s" % a.pseudogenes)

    a.output.mkdir(parents=True, exist_ok=True)
    stem, clusters = vc.load_clusters(a.clusters)
    log("%s: %d clustered genes" % (stem, len(clusters)))
    v34 = parse_v34_pseudogenes(gtf)
    log("GTF: %d pseudogene genes" % sum(len(v[0]) for v in v34.values()))
    consensus = parse_consensus(a.pseudogenes)
    consensus["category"] = classify(consensus, v34)
    tally = consensus["category"].value_counts()
    log("consensus: %d pseudogenes, %d parents; %s" % (
        len(consensus), consensus["parent_id"].nunique(),
        ", ".join("%s %d" % (c, tally.get(c, 0)) for c in CATEGORIES)))

    table = annotate(clusters, consensus)
    summary, tests = summarize(table)

    genes_tsv = a.output / ("%s.pseudogene_counts_genes.tsv" % stem)
    by_cluster_tsv = a.output / ("%s.pseudogene_counts_by_cluster.tsv" % stem)
    tests_tsv = a.output / ("%s.pseudogene_counts_tests.tsv" % stem)
    sources_json = a.output / ("%s.pseudogene_counts_sources.json" % stem)
    table.to_csv(genes_tsv, sep="\t", index=False)
    summary.to_csv(by_cluster_tsv, sep="\t", index=False, float_format="%.4g", na_rep="NA")
    tests.to_csv(tests_tsv, sep="\t", index=False, float_format="%.4g")
    sources_json.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "clusters": vc.describe(a.clusters), "gtf": vc.describe(gtf),
        "pseudogenes": vc.describe(a.pseudogenes),
        "consensus_categories": {c: int(tally.get(c, 0)) for c in CATEGORIES}}, indent=2) + "\n")
    for path in (genes_tsv, by_cluster_tsv, tests_tsv, sources_json):
        log("wrote %s" % path)
    show = ["cluster", "n_genes"] + ["%s_%s" % (k, f[4:]) for f, _c in BOOLEANS for k in ("n", "pct")] \
        + ["n_processed_median", "n_unprocessed_median"]
    print(summary[show].to_string(index=False))
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
