#!/usr/bin/env python3
"""Step 2 of the clustering: keep the genes whose composition is worth clustering.

    python code/clustering/filter_genes.py --input results/clustering/HeLa.post_dedup.gene_counts.tsv

Two conditions, both on columns `build_gene_counts.py` wrote:

    status == "ok"          the gene has a non-empty read-ID union (PAR genes and genes
                            with no exon in the annotation are `excluded`; genes with no
                            reads are `zero_union` -- a composition of an empty union is
                            undefined, not zero)
    n_union >= MIN_UNION    enough reads that the five proportions are estimates rather
                            than coin flips; 100 reads resolve a component down to 1 %

`MIN_UNION` is the ONE cutoff the whole directory uses. Columns and row order pass
through unchanged, so the filtered table is the count table minus rows.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "code" / "common"))
import inputs  # noqa: E402

#: The single gene-depth cutoff, in union reads. Every later step inherits it through
#: the filtered table; nothing downstream re-filters.
MIN_UNION = 100

REQUIRED = ("gene", "gene_id", "transcript_id", "status", "n_union")


log = inputs.make_log("clustering/filter")


def default_output(input_path):
    name = input_path.name
    if name.endswith(".gene_counts.tsv"):
        return input_path.with_name(name[:-len(".gene_counts.tsv")]
                                    + ".gene_counts_filtered.tsv")
    return input_path.with_name(input_path.stem + "_filtered.tsv")


def filter_genes(table, min_union):
    missing = [c for c in REQUIRED if c not in table.columns]
    if missing:
        raise SystemExit("input lacks column(s): %s" % ", ".join(missing))
    ok = table["status"].astype(str) == "ok"
    deep = table["n_union"] >= min_union
    kept = table[ok & deep]
    return kept, int((~ok).sum()), int((ok & ~deep).sum())


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=Path,
                        help="<stem>.gene_counts.tsv from build_gene_counts.py")
    parser.add_argument("--output", type=Path,
                        help="default: the input with .gene_counts_filtered.tsv")
    parser.add_argument("--min-union", type=int, default=MIN_UNION,
                        help="minimum union reads per gene (default %d)" % MIN_UNION)
    args = parser.parse_args(argv)

    table = pd.read_csv(args.input, sep="\t")
    log("%s: %d genes" % (args.input.name, len(table)))
    kept, n_not_ok, n_shallow = filter_genes(table, args.min_union)
    log("dropped %d genes whose status is not ok, %d ok genes with n_union < %d"
        % (n_not_ok, n_shallow, args.min_union))
    if kept.empty:
        raise SystemExit("no gene survives the filter; nothing written")

    output = args.output or default_output(args.input)
    kept.to_csv(output, sep="\t", index=False, lineterminator="\n")
    log("wrote %s (%d genes, n_union %d..%d)"
        % (output, len(kept), int(kept["n_union"].min()), int(kept["n_union"].max())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
