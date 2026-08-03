#!/usr/bin/env python3
"""Figure 4 A and B -- Ribo-seq versus RNA-seq per transcript, one alignment route."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROUTE_COLUMNS = {"genome": ("genome_ribo_reads", "genome_rna_reads"),
                 "transcriptome": ("txome_ribo_reads", "txome_rna_reads")}
ROUTE_COLOUR = {"genome": "#3a923a", "transcriptome": "#cc3d3d"}

def compute_axis_max(counts_path):
    """One log2 maximum over ALL four count vectors, so both routes share a scale."""
    frame = pd.read_csv(counts_path, sep="\t")
    peak = max(float(frame[c].max()) for pair in ROUTE_COLUMNS.values() for c in pair)
    return float(np.log2(peak + 1)) * 1.03

def library_size(frame, ribo_column, rna_column):
    """N_R, N_M -- total CDS-assigned reads for this route, summed over EVERY transcript.

    The denominator a raw `Ribo / RNA` ratio omits: Ribo-seq and RNA-seq are separately
    normalised libraries, so two transcripts with identical raw counts do not carry the
    same signal unless their libraries are also the same size. Summed over the same table
    `prepare()` already loaded, before any gene is picked out, so it never depends on
    which genes are marked.
    """
    return float(frame[ribo_column].sum()), float(frame[rna_column].sum())

def translation_efficiency(ribo, rna, ribo_total, rna_total):
    """Library-normalised TE = (Ribo/N_R) / (RNA/N_M). No log, no pseudocount.

    `ribo_total`/`rna_total` (N_R/N_M, from `library_size`) rescale each raw count by its
    own route's library size before dividing, so this is comparable ACROSS transcripts on
    one panel -- a within-transcript ratio, not an absolute rate. Still computed on raw
    counts, never on the plotted log2(count + 1) values: the axes are log-scaled for
    display only, and taking a ratio of those transformed values would be a different
    quantity (the pseudocount does not cancel).

    `nan` where RNA is zero -- 7,643 of the 19,736 transcripts in the HeLa table have no
    RNA-seq CDS read at all, and a division there is undefined rather than infinite in any
    useful sense. Callers render it as "n/a" instead of printing `inf`.
    """
    ribo = np.asarray(ribo, dtype=float)
    rna = np.asarray(rna, dtype=float)
    numerator = ribo / ribo_total
    denominator = np.divide(rna, rna_total, out=np.full(rna.shape, np.nan), where=rna > 0)
    return np.divide(numerator, denominator,
                     out=np.full(ribo.shape, np.nan), where=rna > 0)

def format_te(value):
    """`TE_GENE = 0.10`-style text, or `n/a` when RNA was zero."""
    return "n/a" if not np.isfinite(value) else "%.2f" % value

def prepare(counts_path, route, orf_catalog=None, mark_genes=(), axis_max=None):
    sys.path.insert(0, str(HERE))
    import panel_style as ps

    ribo_column, rna_column = ROUTE_COLUMNS[route]
    frame = pd.read_csv(counts_path, sep="\t")
    ps.require_columns(frame, ("transcript_id", ribo_column, rna_column), str(counts_path))
    ribo = frame[ribo_column].to_numpy(float)
    rna = frame[rna_column].to_numpy(float)

    from scipy import stats
    spearman = float(stats.spearmanr(ribo, rna).correlation)
    pearson = float(stats.pearsonr(np.log2(ribo + 1.0), np.log2(rna + 1.0))[0])

    ribo_total, rna_total = library_size(frame, ribo_column, rna_column)

    marks = {"x": np.array([]), "y": np.array([]), "labels": [], "te": np.array([])}
    if mark_genes and orf_catalog:
        catalog = pd.read_csv(orf_catalog, sep="\t", usecols=["base_enst", "gene_name"])
        order = [g.strip().upper() for g in mark_genes]
        lookup = catalog[catalog["gene_name"].astype(str).str.upper().isin(set(order))]
        base = frame["transcript_id"].astype(str).str.split(".").str[0]
        joined = frame.assign(base_enst=base).merge(lookup, on="base_enst", how="inner")
        joined = joined.assign(
            _rank=joined["gene_name"].astype(str).str.upper().map(order.index)
        ).sort_values("_rank")
        ribo_marked = joined[ribo_column].to_numpy(float)
        rna_marked = joined[rna_column].to_numpy(float)
        marks = {"x": ribo_marked,
                 "y": rna_marked,
                 "labels": joined["gene_name"].tolist(),
                 "transcript_ids": joined["transcript_id"].tolist(),
                 "te": translation_efficiency(ribo_marked, rna_marked,
                                              ribo_total, rna_total)}

    if axis_max is None:
        # four-column maximum, because scaling a paired figure to one route's own data is
        axis_max = float(np.log2(max(ribo.max(), rna.max()) + 1)) * 1.03
        axis_provenance = "derived from this route alone"
    else:
        axis_provenance = "shared across both routes"
    return {"route": route, "ribo": ribo, "rna": rna, "n_transcripts": int(len(frame)),
            "spearman": spearman, "pearson": pearson, "marks": marks,
            "ribo_library": ribo_total, "rna_library": rna_total,
            "axis_max": float(axis_max), "axis_provenance": axis_provenance,
            "source": str(counts_path)}

def draw(prepared, metric="spearman", show_ylabel=True, figsize=(4.6, 4.6)):
    import matplotlib.pyplot as plt
    import panel_style as ps

    ps.apply_rcparams()
    hi = prepared["axis_max"]
    figure, axis = plt.subplots(figsize=figsize)
    axis.plot([0, hi], [0, hi], ls="--", lw=1.0, color="#999999", zorder=1)
    # The point CLOUD is rasterised on purpose: ~20,000 individual vector circles would make
    # axes, identity line, ringed genes, all text -- stays vector and editable. Contrast the
    # Figure-2 heatmaps, whose flat cells SHOULD be vector; here the data is a dense cloud
    axis.scatter(np.log2(prepared["ribo"] + 1), np.log2(prepared["rna"] + 1), s=5,
                 alpha=0.40, linewidths=0, color=ROUTE_COLOUR[prepared["route"]],
                 zorder=2, rasterized=True)

    marks = prepared["marks"]
    if len(marks["x"]):
        mx, my = np.log2(marks["x"] + 1), np.log2(marks["y"] + 1)
        axis.scatter(mx, my, s=46, facecolors="white", edgecolors="black", linewidths=1.6,
                     zorder=6)
        for gx, gy, name in zip(mx, my, marks["labels"]):
            offset, align = (9, "left") if gx < 0.25 * hi else (-10, "right")
            axis.annotate(name, (gx, gy), textcoords="offset points", xytext=(offset, 6),
                          ha=align, fontsize=ps.FONT_ANNOTATION, color="black", fontweight="bold",
                          zorder=7)

    ticks = [v for v in range(0, 20, 2) if v <= hi]
    axis.set_xticks(ticks)
    axis.set_yticks(ticks)
    axis.set_xlim(0, hi)
    axis.set_ylim(0, hi)
    axis.set_aspect("equal")
    axis.set_xlabel("Ribo-seq reads, log$_2$(count + 1)", fontsize=ps.FONT_LABEL)
    axis.set_ylabel("RNA-seq reads, log$_2$(count + 1)" if show_ylabel else "",
                    fontsize=ps.FONT_LABEL)
    if not show_ylabel:
        axis.tick_params(axis="y", labelleft=False)
    symbol, value = (("Spearman $\\rho$", prepared["spearman"]) if metric == "spearman"
                     else ("Pearson $r$", prepared["pearson"]))
    lines = ["n = %s transcripts" % format(prepared["n_transcripts"], ","),
             "%s = %.3f" % (symbol, value)]
    for name, te in zip(marks["labels"], marks.get("te", [])):
        lines.append("TE$_{\\mathrm{%s}}$ = %s" % (name, format_te(te)))
    axis.text(0.03, 0.97, "\n".join(lines),
              transform=axis.transAxes, va="top", ha="left", fontsize=ps.FONT_ANNOTATION)
    for side in ("top", "right", "left", "bottom"):
        axis.spines[side].set_visible(True)
        axis.spines[side].set_color("black")
        axis.spines[side].set_linewidth(0.8)
    figure.tight_layout()
    return figure, axis

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--route", choices=tuple(ROUTE_COLUMNS), default="genome")
    parser.add_argument("--orf-catalog", type=Path)
    parser.add_argument("--mark-genes", default="COMT,GAPDH")
    parser.add_argument("--metric", choices=("spearman", "pearson"), default="spearman")
    parser.add_argument("--axis-max", type=float, default=None,
                        help="override the shared log2 axis maximum. Omit it and both "
                             "routes derive the same value from all four count columns.")
    parser.add_argument("--hide-ylabel", action="store_true")
    parser.add_argument("--figsize", nargs=2, type=float, default=(4.6, 4.6))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    axis_max = args.axis_max if args.axis_max is not None else compute_axis_max(args.counts)
    if not args.output:
        raise SystemExit("--output is required")

    genes = [g.strip() for g in args.mark_genes.split(",") if g.strip()]
    prepared = prepare(args.counts, args.route, args.orf_catalog, genes, axis_max)
    print("[panel] %s route: %d transcripts, Spearman %.4f, Pearson(log2) %.4f"
          % (prepared["route"], prepared["n_transcripts"], prepared["spearman"],
             prepared["pearson"]))
    print("[panel] axis max %.4f (%s); marked %s"
          % (prepared["axis_max"], prepared["axis_provenance"],
             ", ".join(prepared["marks"]["labels"]) or "nothing"))
    print("[panel] TE library normalisation: N_R=%d Ribo-seq / N_M=%d RNA-seq CDS reads, "
          "summed over %d transcripts"
          % (prepared["ribo_library"], prepared["rna_library"], prepared["n_transcripts"]))

    figure, _axis = draw(prepared, args.metric, not args.hide_ylabel, tuple(args.figsize))
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
