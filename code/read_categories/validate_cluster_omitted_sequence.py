#!/usr/bin/env python3
"""Fraction of each clustered gene's annotated sequence that the transcriptome reference omits.

    python code/clustering/validate_cluster_omitted_sequence.py \\
        --clusters results/clustering/HeLa.post_dedup.clusters_k4.tsv \\
        --output results/clustering [--gtf G]

For every clustered gene, from the GENCODE v34 GTF alone (nothing read-derived is opened):

    omitted_exon_fraction = 1 - exon_nt_selected / exon_nt_gene
    omitted_cds_fraction  = 1 - cds_nt_selected  / cds_nt_gene

where `exon_nt_gene` is the length of the union of every `exon` feature of every
transcript of the gene (overlapping intervals merged, so a nucleotide counts once whatever
number of isoforms carry it), `exon_nt_selected` is the merged exon length of the
APPRIS-selected transcript (the one model of the gene the transcriptome reference holds),
and the CDS pair is the same on `CDS` features. The selected transcript's exons are a
subset of the gene's, so each fraction lies in [0, 1]; it is the share of the gene's
annotated exonic (or coding) nucleotides that no read can reach on the transcriptome
route. A gene whose union of CDS features is empty has an undefined CDS fraction (NA).
`--transcript-types` restricts the union to transcripts of the listed GENCODE
transcript_type(s); the default is every transcript of the gene.

Outputs, under <output> with the cluster table's stem:

    <stem>.omitted_sequence_genes.tsv       per gene: cluster, the four lengths, the two
                                            fractions, number of transcripts in the union
    <stem>.omitted_sequence_by_cluster.tsv  per cluster (and all): n, n defined, quartiles,
                                            mean, % of genes with fraction > 0; plus the
                                            Kruskal-Wallis test across clusters and the
                                            pairwise Mann-Whitney U p-values (Holm) for
                                            each metric
    <stem>.omitted_sequence_sources.json    the inputs by path, size and sha256

The genes table (exon fraction) is Figure 6F's input (`panels/plot_read_fate_clusters.py`).
"""
from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import validate_common as vc
from validate_common import ATTR, KEY, inputs, merged_length

FEATURES = ("exon", "CDS")
METRICS = (("omitted_exon_fraction", "exon", "annotated exonic sequence"),
           ("omitted_cds_fraction", "cds", "annotated CDS sequence"))

log = inputs.make_log("clustering/omitted")


def parse_gtf(gtf, wanted_genes, transcript_types):
    """Per gene: {feature: {transcript_id: [(start, end), ...]}} and transcript types.

    Coordinates become 0-based half-open. Only the wanted genes are kept.
    """
    per_gene = {g: {f: defaultdict(list) for f in FEATURES} for g in wanted_genes}
    tx_type = {}
    for fields in vc.gtf_records(gtf, FEATURES):
        attrs = dict(ATTR.findall(fields[8]))
        gene_id = attrs["gene_id"]
        if gene_id not in per_gene:
            continue
        tx = attrs["transcript_id"]
        tx_type[tx] = attrs.get("transcript_type", "")
        if transcript_types and tx_type[tx] not in transcript_types:
            continue
        per_gene[gene_id][fields[2]][tx].append((int(fields[3]) - 1, int(fields[4])))
    return per_gene, tx_type


def annotate(clusters, per_gene, tx_type, transcript_types):
    rows = []
    for gene, gene_id, transcript_id, cluster in clusters[KEY + ["cluster"]].itertuples(index=False):
        feats = per_gene[gene_id]
        row = {"gene": gene, "gene_id": gene_id, "transcript_id": transcript_id,
               "cluster": int(cluster),
               "selected_transcript_type": tx_type.get(transcript_id, "")}
        if transcript_id not in feats["exon"]:
            if transcript_types and tx_type.get(transcript_id, "") not in transcript_types:
                raise SystemExit("%s: selected transcript %s is of type %r, outside "
                                 "--transcript-types" % (gene, transcript_id, tx_type.get(transcript_id)))
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
    return pd.DataFrame(rows)


# ── summary ──────────────────────────────────────────────────────────────────────────────
def summarize(table):
    clusters = sorted(table["cluster"].unique())
    rows, tests = [], []
    for metric, _tag, _label in METRICS:
        groups = {c: table.loc[table["cluster"] == c, metric].dropna().to_numpy(float)
                  for c in clusters}
        for label, values in [("all", table[metric].dropna().to_numpy(float))] + list(groups.items()):
            n_total = len(table) if label == "all" else int((table["cluster"] == label).sum())
            rows.append({"metric": metric, "cluster": str(label), "n_genes": n_total,
                         "n_defined": len(values),
                         "pct_omitting_any": 100.0 * np.mean(values > 0) if len(values) else np.nan,
                         **vc.quartiles(values)})
        tests += [{"metric": metric, **t} for t in vc.cluster_tests(groups, clusters)]
    return pd.DataFrame(rows), pd.DataFrame(tests)


# ── figure ───────────────────────────────────────────────────────────────────────────────
# ── main ─────────────────────────────────────────────────────────────────────────────────
def main(argv=None):
    p = vc.parser(__doc__)
    p.add_argument("--transcript-types", default="",
                   help="comma-separated GENCODE transcript_type(s) to include in the gene "
                        "union; default: every transcript of the gene")
    a = p.parse_args(argv)
    gtf = vc.resolve_gtf(a)
    types = {t.strip() for t in a.transcript_types.split(",") if t.strip()}

    a.output.mkdir(parents=True, exist_ok=True)
    stem, clusters = vc.load_clusters(a.clusters)
    log("%s: %d clustered genes" % (stem, len(clusters)))
    per_gene, tx_type = parse_gtf(gtf, set(clusters["gene_id"]), types)
    log("GTF: %d transcripts of the clustered genes%s" % (
        len(tx_type), " (union restricted to %s)" % ", ".join(sorted(types)) if types else ""))

    table = annotate(clusters, per_gene, tx_type, types)
    summary, tests = summarize(table)

    genes_tsv = a.output / ("%s.omitted_sequence_genes.tsv" % stem)
    by_cluster_tsv = a.output / ("%s.omitted_sequence_by_cluster.tsv" % stem)
    tests_tsv = a.output / ("%s.omitted_sequence_tests.tsv" % stem)
    sources_json = a.output / ("%s.omitted_sequence_sources.json" % stem)
    table.to_csv(genes_tsv, sep="\t", index=False, float_format="%.6f", na_rep="NA")
    summary.to_csv(by_cluster_tsv, sep="\t", index=False, float_format="%.4f", na_rep="NA")
    tests.to_csv(tests_tsv, sep="\t", index=False, float_format="%.4g")
    sources_json.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "clusters": vc.describe(a.clusters), "gtf": vc.describe(gtf),
        "transcript_types": sorted(types) or "all"}, indent=2) + "\n")
    for path in (genes_tsv, by_cluster_tsv, tests_tsv, sources_json):
        log("wrote %s" % path)
    print(summary.to_string(index=False))
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
