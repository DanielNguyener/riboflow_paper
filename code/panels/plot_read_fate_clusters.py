#!/usr/bin/env python3
"""Figure 6: genes clustered by how their reads split between the two alignment routes.

    python code/panels/plot_read_fate_clusters.py \\
        --clusters data/clustering/HeLa.post_dedup.clusters_k4.tsv \\
        --centroids data/clustering/HeLa.post_dedup.cluster_centroids.tsv \\
        --tree data/clustering/HeLa.post_dedup.tree_merge.tsv \\
        --pseudogene-counts data/clustering/HeLa.post_dedup.pseudogene_counts_genes.tsv \\
        --omitted-sequence data/clustering/HeLa.post_dedup.omitted_sequence_genes.tsv \\
        --reference-duplication data/clustering/HeLa.post_dedup.reference_duplication_entries.tsv \\
        --output results/panels/fig06_read_fate_clusters

Every table comes from `code/clustering/` (the `clustering` stage of make_tables.py): the
Ward tree (`ward_cluster.R`, R's merge matrix and heights), its cut at k with each gene's
composition in %, the centroid of each cluster, and the three per-gene validation tables.
Seven panels, lettered here because the whole page is one panel to the assembler:

    A  dendrogram   every subtree below the cut in the colour of the fate that dominates its
                    cluster's centroid, grey above; the cut as a dotted line; cluster sizes
    B  heatmap      the five shares of every gene, genes in leaf order under the dendrogram,
                    with a numbered cluster strip between
    C  composition  per fate, the % of each gene's union it holds, one box per cluster
    D  centroids    one stacked bar per cluster (member mean composition)
    E  pseudogenes  % of each cluster's genes with >= 1 GENCODE consensus pseudogene
    F  omitted      % of each gene's annotated exonic sequence absent from its reference transcript
    G  duplication  % of each gene's exonic sequence duplicated by another reference entry

The R merge matrix is turned into a scipy linkage matrix, and the scipy cut at k must
reproduce the clusters table before anything is drawn.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import panel_style as ps  # noqa: E402

#: The five fates in table order; index 0 is the concordant one (cluster 1 by construction).
COMPONENTS = ("shared_genome_unique", "shared_genome_multimapped", "genome_only_unique",
              "genome_only_multimapped", "transcriptome_only")
#: The Figure 5A fate abbreviations (`plot_gene_read_partition.ROUTE7_KEY`), in
#: COMPONENTS order: SH shared, GO genome-only, TO transcriptome-only; U unique, M multi.
SHORT = ("SH-U", "SH-M", "GO-U", "GO-M", "TO")
#: The Fig 4B fate colours, in COMPONENTS order (`plot_read_id_union.SEGMENTS`).
COLOURS = ("#a6d96a", "#1a7d1a", "#7fb9da", "#0d57a1", "#cc3d3d")
EDGE = "#333333"
ABOVE_CUT = "#8c8c8c"
CUT = "#4d4d4d"
HEAT_CMAP = "Blues"
FIGSIZE = (6.9, 8.9)   # tight bbox lands inside PLOS's 7.5 x 8.75 in page
#: Genes pointed out under the heatmap: the three Figure 5A genes.
LABEL_GENES = ("COMT", "GAPDH", "LRRFIP1")


# ── the tree ─────────────────────────────────────────────────────────────────────────────
def r_merge_to_linkage(merge, height):
    """scipy linkage matrix from R's `hclust$merge` (negative = leaf, positive = merge)."""
    merge = np.asarray(merge, dtype=int)
    height = np.asarray(height, dtype=float)
    n = len(height) + 1
    if merge.shape != (n - 1, 2):
        raise ValueError("merge must be (n - 1) x 2")
    size = np.ones(2 * n - 1)
    Z = np.empty((n - 1, 4))
    for i, ((a, b), h) in enumerate(zip(merge, height)):
        ia = -a - 1 if a < 0 else n + a - 1
        ib = -b - 1 if b < 0 else n + b - 1
        size[n + i] = size[ia] + size[ib]
        Z[i] = (min(ia, ib), max(ia, ib), h, size[n + i])
    return Z


def linkage_to_r_merge(Z):
    """The inverse: R's merge matrix and heights from a scipy linkage matrix."""
    n = len(Z) + 1
    ids = Z[:, :2].astype(int)
    merge = np.where(ids < n, -(ids + 1), ids - n + 1)
    return merge, Z[:, 2].copy()


def same_partition(a, b):
    """True when two labellings induce the same partition (labels may differ)."""
    table = pd.crosstab(np.asarray(a), np.asarray(b))
    return ((table > 0).sum(axis=1) == 1).all() and ((table > 0).sum(axis=0) == 1).all()


def load(args):
    clusters = pd.read_csv(args.clusters, sep="\t")
    ps.require_columns(clusters, ["gene_id", "transcript_id", "cluster"]
                       + ["pct_" + c for c in COMPONENTS], args.clusters.name)
    centroids = pd.read_csv(args.centroids, sep="\t").sort_values("cluster")
    ps.require_columns(centroids, ["cluster", "n"] + ["pct_" + c for c in COMPONENTS],
                       args.centroids.name)
    k = int(centroids["cluster"].max())
    if int(clusters["cluster"].max()) != k or len(centroids) != k:
        ps.die("cluster count disagrees between %s and %s" % (args.clusters.name,
                                                               args.centroids.name))
    tree = pd.read_csv(args.tree, sep="\t")
    ps.require_columns(tree, ["merge1", "merge2", "height"], args.tree.name)
    Z = r_merge_to_linkage(tree[["merge1", "merge2"]].to_numpy(), tree["height"].to_numpy())
    if len(Z) + 1 != len(clusters):
        ps.die("tree and clusters table differ in size")
    if not same_partition(hierarchy.fcluster(Z, k, "maxclust"), clusters["cluster"]):
        ps.die("the clusters table is not the k = %d cut of this tree" % k)
    X = clusters[["pct_" + c for c in COMPONENTS]].to_numpy(float) / 100.0
    validation = {}
    for name, key in (("pseudogene_counts", "gene_id"), ("omitted_sequence", "gene_id"),
                      ("reference_duplication", "transcript_id")):
        path = getattr(args, name)
        if path is None:
            continue
        table = pd.read_csv(path, sep="\t")
        ours = clusters.set_index(key)["cluster"]
        theirs = table.dropna(subset=["cluster"]).set_index(key)["cluster"].astype(int)
        if len(theirs) != len(ours) or not theirs.reindex(ours.index).eq(ours).all():
            ps.die("%s: cluster labels do not match %s" % (path.name, args.clusters.name))
        validation[name] = table
    return k, clusters, centroids, X, Z, validation


# ── drawing helpers ──────────────────────────────────────────────────────────────────────
def cluster_colours(centroids):
    """One colour per cluster (index cluster - 1): the fate that dominates its centroid."""
    shares = centroids.sort_values("cluster")[["pct_" + c for c in COMPONENTS]].to_numpy(float)
    return [COLOURS[j] for j in shares.argmax(axis=1)]


def is_dark(hex_colour):
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255.0 for i in (1, 3, 5))
    return 0.299 * r + 0.587 * g + 0.114 * b < 0.5


def box(axis, groups, x, colours):
    """One box per cluster (5-95 % whiskers, no fliers), in the cluster colour."""
    bp = axis.boxplot(groups, positions=x, widths=0.6, whis=(5, 95), showfliers=False,
                      patch_artist=True, zorder=3,
                      whiskerprops=dict(color=EDGE, linewidth=0.8),
                      capprops=dict(color=EDGE, linewidth=0.8),
                      medianprops=dict(color="black", linewidth=1.3))
    for patch, median, colour in zip(bp["boxes"], bp["medians"], colours):
        patch.set(facecolor=colour, edgecolor=EDGE, linewidth=0.8)
        if is_dark(colour):     # a black median vanishes on the dark fills
            median.set_color("white")


def _cluster_axis(axis, clusters):
    x = list(range(1, len(clusters) + 1))
    axis.set_xticks(x)
    axis.set_xticklabels([str(c) for c in clusters])
    axis.set_xlim(0.4, len(clusters) + 0.6)
    axis.set_xlabel("cluster", fontsize=ps.FONT_TICK)
    axis.tick_params(labelsize=ps.FONT_TICK)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)
    return x


def draw_pseudogene_bars(axis, table, colours):
    """E: % of each cluster's genes with >= 1 consensus pseudogene of any category."""
    clusters = sorted(table["cluster"].unique())
    x = _cluster_axis(axis, clusters)
    pct = [100.0 * table.loc[table["cluster"] == c, "has_any"].mean() for c in clusters]
    axis.bar(x, pct, width=0.68, color=colours, edgecolor=EDGE, linewidth=0.6)
    for xi, value in zip(x, pct):
        axis.text(xi, value + max(pct) * 0.02, "%.0f%%" % value, ha="center", va="bottom",
                  fontsize=ps.FONT_INSET)
    axis.set_ylim(0, max(pct) * 1.18)
    axis.set_ylabel("genes with ≥1\npseudogene (%)", fontsize=ps.FONT_TICK)


def draw_omitted_box(axis, table, colours):
    """F: per cluster, a box of the % of each gene's exonic sequence its reference transcript
    omits, the median above each."""
    clusters = sorted(table["cluster"].unique())
    x = _cluster_axis(axis, clusters)
    groups = [100.0 * table.loc[table["cluster"] == c, "omitted_exon_fraction"]
              .dropna().to_numpy(float) for c in clusters]
    box(axis, groups, x, colours)
    for xi, g in zip(x, groups):
        axis.text(xi, 103, "%.0f%%" % np.median(g), ha="center", va="bottom",
                  fontsize=ps.FONT_INSET)
    axis.set_ylim(-3, 112)
    axis.set_yticks([0, 20, 40, 60, 80, 100])
    axis.set_ylabel("exonic sequence absent\nfrom transcriptome (%)", fontsize=ps.FONT_TICK)


def draw_duplication_box(axis, table, colours):
    """G: box of the same-strand duplicated fraction per cluster (clustered entries only),
    the median above each."""
    clustered = table[table["cluster"].notna()]
    clusters = sorted(clustered["cluster"].astype(int).unique())
    x = _cluster_axis(axis, clusters)
    groups = [100.0 * clustered.loc[clustered["cluster"] == c, "duplicated_fraction"]
              .to_numpy(float) for c in clusters]
    drawn = axis.boxplot(groups, positions=x, widths=0.68, patch_artist=True, showfliers=False,
                         medianprops={"color": EDGE, "linewidth": 1.0})
    for patch, colour in zip(drawn["boxes"], colours):
        patch.set(facecolor=colour, edgecolor=EDGE, linewidth=0.6)
    for xi, g in zip(x, groups):
        axis.text(xi, 103, "%.0f%%" % np.median(g), ha="center", va="bottom",
                  fontsize=ps.FONT_INSET)
    axis.set_ylim(-3, 112)
    axis.set_yticks([0, 20, 40, 60, 80, 100])
    axis.set_ylabel("duplicated transcriptome\nexonic sequence (%)", fontsize=ps.FONT_TICK)


# ── the figure ───────────────────────────────────────────────────────────────────────────
def draw(k, clusters, centroids, X, Z, validation, figsize=FIGSIZE):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    n = len(clusters)
    x = np.arange(1, k + 1)
    labels = clusters["cluster"].to_numpy()
    colours = cluster_colours(centroids)            # indexed by cluster - 1
    # Any height strictly between the merge that leaves k clusters and the one that leaves
    # k - 1 cuts the tree into the k clusters of the table.
    cut_h = float(Z[n - k - 1:n - k + 1, 2].mean())

    # Each link's colour: the cluster of the leaves under it when below the cut.
    link_cluster = np.zeros(2 * n - 1, dtype=int)
    link_cluster[:n] = labels
    for i, (a, b, h, _) in enumerate(Z):
        link_cluster[n + i] = link_cluster[int(a)] if h < cut_h else 0

    def link_colour(link_id):
        c = link_cluster[link_id]
        return colours[c - 1] if c else ABOVE_CUT

    extra = [name for name in ("pseudogene_counts", "omitted_sequence", "reference_duplication")
             if name in validation]
    fig = plt.figure(figsize=figsize)
    # an empty spacer row holds the fate legend between the composition and validation rows
    # (hspace is a fraction of the mean row height, so the spacer row needs a larger one)
    outer = fig.add_gridspec(4 if extra else 2, 1,
                             height_ratios=(4.0, 2.2, 0.25, 2.0) if extra else (4.0, 2.2),
                             hspace=0.44 if extra else 0.30)
    grid = outer[0].subgridspec(3, 1, height_ratios=(1.6, 0.2, 2.2), hspace=0.08)
    ax_tree = fig.add_subplot(grid[0])
    ax_strip = fig.add_subplot(grid[1], sharex=ax_tree)
    ax_heat = fig.add_subplot(grid[2], sharex=ax_tree)

    # scipy's dendrogram recurses once per level, so the limit must scale with the leaves
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 4 * n + 1000))
    dend = hierarchy.dendrogram(Z, ax=ax_tree, no_labels=True, link_color_func=link_colour,
                                color_threshold=cut_h)
    leaves = np.asarray(dend["leaves"])
    for line in ax_tree.collections + ax_tree.lines:
        line.set_linewidth(0.6)
    ax_tree.axhline(cut_h, linestyle=":", color=CUT, linewidth=0.9)
    ax_tree.set_ylabel("height", fontsize=ps.FONT_TICK)
    ax_tree.set_ylim(0, Z[-1, 2] * 1.04)
    ax_tree.yaxis.set_major_locator(MaxNLocator(nbins=8, steps=[1, 2, 5, 10]))
    ax_tree.tick_params(axis="x", bottom=False, labelbottom=False)
    for side in ("top", "right"):
        ax_tree.spines[side].set_visible(False)
    sizes = centroids.set_index("cluster")["n"]
    # cluster sizes as a boxed two-column table under the total (columns align, unlike
    # spaces in a proportional font)
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea, VPacker
    props = dict(fontsize=ps.FONT_ANNOTATION)
    cell = lambda c: TextArea("%d: n = %s" % (c, format(int(sizes[c]), ",")), textprops=props)  # noqa: E731
    columns = [VPacker(children=[cell(c) for c in x[j::2]], align="left", pad=0, sep=3)
               for j in range(2)]
    table = VPacker(children=[TextArea("n = %s" % format(n, ","), textprops=props),
                              HPacker(children=columns, align="top", pad=0, sep=10)],
                    align="left", pad=0, sep=3)
    sizes_box = AnchoredOffsetbox(loc="upper left", child=table, pad=0.35, borderpad=0.3,
                                  frameon=True, bbox_to_anchor=(0.005, 0.98),
                                  bbox_transform=ax_tree.transAxes)
    sizes_box.patch.set(edgecolor="black", linewidth=0.7)
    sizes_box.set_zorder(5)
    ax_tree.add_artist(sizes_box)

    # Leaves sit at x = 5, 15, 25, ... in dendrogram units; the strip and heatmap share
    # that axis so one gene is one column of both.
    x_edges = np.arange(n + 1) * 10.0
    leaf_labels = labels[leaves]
    ax_strip.pcolormesh(x_edges, [0, 1], leaf_labels[np.newaxis, :] - 1,
                        cmap=ListedColormap(colours), vmin=-0.5, vmax=k - 0.5,
                        rasterized=True)
    ax_strip.set_yticks([])
    ax_strip.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_strip.set_ylabel("cluster", rotation=0, ha="right", va="center", fontsize=ps.FONT_TICK)
    boundaries = np.flatnonzero(np.diff(leaf_labels)) + 1
    for start, stop in zip(np.r_[0, boundaries], np.r_[boundaries, n]):
        if stop - start >= n * 0.03:     # a run of leaves wide enough to hold its number
            ax_strip.text((start + stop) * 5.0, 0.5, str(leaf_labels[start]), ha="center",
                          va="center", fontsize=ps.FONT_INSET, color="white")

    # 5 x n: one row per fate, genes in leaf order; rasterized, since 12,000 columns cannot
    # be vector cells.
    mesh = ax_heat.pcolormesh(x_edges, np.arange(len(COMPONENTS) + 1), X[leaves].T,
                              cmap=HEAT_CMAP, vmin=0, vmax=1, rasterized=True)
    ax_heat.set_ylim(len(COMPONENTS), 0)
    ax_heat.set_yticks(np.arange(len(COMPONENTS)) + 0.5)
    ax_heat.set_yticklabels(SHORT, fontsize=ps.FONT_TICK)
    ax_heat.set_xlim(0, n * 10.0)
    ax_heat.tick_params(axis="x", bottom=False, labelbottom=False)
    # each cluster's run of columns boxed (the cut is contiguous in leaf order)
    from matplotlib.patches import Rectangle
    for start, stop in zip(np.r_[0, boundaries], np.r_[boundaries, n]):
        ax_heat.add_patch(Rectangle((start * 10.0, 0), (stop - start) * 10.0, len(COMPONENTS),
                                    fill=False, edgecolor="black", linewidth=0.7, zorder=5,
                                    clip_on=False))
    # the Figure 5A genes, pointed out at their column
    position = {gene: i for i, gene in enumerate(clusters["gene"])}
    column = np.empty(n, dtype=int)
    column[leaves] = np.arange(n)
    # Labels whose columns sit closer than their text is wide are spread apart, left to
    # right; each leader still ends at its gene's column.
    heat_box = ax_heat.get_position()
    data_per_pt = n * 10.0 / (heat_box.width * fig.get_figwidth() * 72.0)
    drop = -11.0 / (heat_box.height * fig.get_figheight() * 72.0)       # 11 pt, in axes fraction
    placed = []
    for gene in sorted((g for g in LABEL_GENES if g in position),
                       key=lambda g: column[position[g]]):
        x_gene = column[position[gene]] * 10.0 + 5.0
        half = 0.5 * len(gene) * 0.6 * ps.FONT_INSET * data_per_pt
        x_text = x_gene
        if placed:
            x_prev, half_prev = placed[-1]
            x_text = max(x_gene, x_prev + half_prev + half + 6.0 * data_per_pt)
        placed.append((x_text, half))
        ax_heat.annotate(gene, xy=(x_gene, len(COMPONENTS)), xycoords="data",
                         xytext=(x_text, drop), textcoords=("data", "axes fraction"),
                         ha="center", va="top",
                         fontsize=ps.FONT_INSET, annotation_clip=False,
                         arrowprops=dict(arrowstyle="-", color=EDGE, linewidth=0.7,
                                         shrinkA=0, shrinkB=1))
    # The colour bar lives in an inset outside the heatmap: `fig.colorbar(ax=ax_heat)` would
    # shrink ax_heat alone, and its columns would no longer sit under the tree's leaves.
    bar = fig.colorbar(mesh, cax=ax_heat.inset_axes([1.012, 0.0, 0.018, 1.0]))
    bar.set_label("fraction of union", fontsize=ps.FONT_TICK)
    bar.ax.tick_params(labelsize=ps.FONT_TICK)

    # composition: one axis per fate, x = cluster, y = % of the gene's union; then the
    # centroids as stacked bars
    bottom_row = outer[1].subgridspec(1, 2, width_ratios=(len(COMPONENTS), 1.3), wspace=0.28)
    inner = bottom_row[0].subgridspec(1, len(COMPONENTS), wspace=0.2)
    axes = []
    for j, (component, short) in enumerate(zip(COMPONENTS, SHORT)):
        ax = fig.add_subplot(inner[j], sharey=axes[0] if axes else None)
        axes.append(ax)
        pct = clusters["pct_" + component].to_numpy(float)
        box(ax, [pct[labels == c] for c in x], x, colours)
        ax.set_title(short, fontsize=ps.FONT_TICK, pad=3)
        _cluster_axis(ax, list(x))
        for side in ("top", "right"):
            ax.spines[side].set_visible(True)
        if j:
            ax.tick_params(axis="y", labelleft=False)
        else:
            ax.set_ylabel("fraction of union (%)", fontsize=ps.FONT_TICK)
    axes[0].set_ylim(-3, 103)

    ax_bar = fig.add_subplot(bottom_row[1])
    shares = centroids[["pct_" + c for c in COMPONENTS]].to_numpy(float)
    bottom = np.zeros(k)
    for j, colour in enumerate(COLOURS):
        ax_bar.bar(x, shares[:, j], bottom=bottom, color=colour, width=0.72,
                   edgecolor="white", linewidth=0.3)
        bottom += shares[:, j]
    for i in range(k):
        for j in range(len(COMPONENTS)):
            if shares[i, j] >= 8:
                ax_bar.text(x[i], shares[i, :j].sum() + shares[i, j] / 2, "%.0f" % shares[i, j],
                            ha="center", va="center", fontsize=ps.FONT_INSET,
                            color="white" if is_dark(COLOURS[j]) else "black")
    _cluster_axis(ax_bar, list(x))
    for side in ("top", "right"):
        ax_bar.spines[side].set_visible(True)
    ax_bar.set_ylim(0, 100)
    ax_bar.set_yticks([0, 25, 50, 75, 100])
    ax_bar.set_ylabel("centroid composition (%)", fontsize=ps.FONT_TICK)
    ax_bar.set_title("centroids", fontsize=ps.FONT_INSET, pad=3)

    # The fate colours, once, under the whole composition row; anchored to the middle
    # composition axis so it centres on the figure.
    handles = [Patch(facecolor=c, edgecolor=EDGE, label=s) for s, c in zip(SHORT, COLOURS)]
    legend = ps.legend_below(axes[len(axes) // 2], handles=handles, ncol=len(handles),
                             fontsize=ps.FONT_INSET, handlelength=1.0, columnspacing=1.2,
                             x=1.15, pad_pt=ps.LEGEND_PAD_PT + 4)

    # (axis, x pad, y pad) in points; E-G sit higher, clear of the panel F medians
    lettered = [(ax_tree, 70, 4), (ax_heat, 70, 4), (axes[0], 42, 4), (ax_bar, 14, 4)]
    if extra:
        drawers = {"pseudogene_counts": draw_pseudogene_bars,
                   "omitted_sequence": draw_omitted_box,
                   "reference_duplication": draw_duplication_box}
        row = outer[3].subgridspec(1, len(extra), wspace=0.7)
        for slot, name in zip(row, extra):
            axis = fig.add_subplot(slot)
            drawers[name](axis, validation[name], colours)
            lettered.append((axis, 42, 14))
    # bold letters just outside each panel's top-left corner, clear of its y label
    for (axis, pad, lift), letter in zip(lettered, "ABCDEFG"):
        axis.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-pad, lift),
                      textcoords="offset points", ha="left", va="bottom",
                      fontsize=ps.FONT_TITLE + 1, fontweight="bold")   # 12 pt: PLOS's maximum
    return fig, legend


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clusters", required=True, type=Path, help="<stem>.clusters_k<K>.tsv")
    parser.add_argument("--centroids", required=True, type=Path, help="<stem>.cluster_centroids.tsv")
    parser.add_argument("--tree", required=True, type=Path, help="<stem>.tree_merge.tsv")
    parser.add_argument("--pseudogene-counts", type=Path,
                        help="<stem>.pseudogene_counts_genes.tsv: panel E")
    parser.add_argument("--omitted-sequence", type=Path,
                        help="<stem>.omitted_sequence_genes.tsv: panel F")
    parser.add_argument("--reference-duplication", type=Path,
                        help="<stem>.reference_duplication_entries.tsv: panel G")
    parser.add_argument("--figsize", nargs=2, type=float, default=FIGSIZE)
    parser.add_argument("--output", required=True, type=Path, help="stem, no extension")
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    k, clusters, centroids, X, Z, validation = load(args)
    print("[panel] %s: k = %d, %d genes, cut verified against the tree"
          % (args.clusters.name, k, len(clusters)))
    ps.apply_rcparams()
    ps.resolve_font()
    fig, _legend = draw(k, clusters, centroids, X, Z, validation, tuple(args.figsize))
    for path in ps.save(fig, args.output, ps.resolve_formats(args.formats), args.force):
        print("[panel] wrote %s" % path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
