#!/usr/bin/env python3
"""The three-panel alignment-route TE figure: one combined page, or one panel on its own.

Panels A (route agreement), B (ranked delta TE + CI), C (delta RNA vs delta Ribo plane);
one solved `Geometry` shared by all pages. Mathematics: docs/methods_te_route.md.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

sys.path.insert(0, str(HERE))
import te_panel_style as ps  # noqa: E402

PANELS = ("combined", "A", "B", "C")

#: White at exactly zero, then viridis; the white stop needs the vmin solve below.
WHITE_VIRIDIS = [(0.0, "#ffffff"), (1e-20, "#440053"), (0.2, "#404388"), (0.4, "#2a788e"),
                 (0.6, "#21a784"), (0.8, "#78d151"), (1.0, "#fde624")]
BAND, MEAN = "#b9c9d8", "#2c5f8a"
#: INK is text (black, print-safe); ZERO is a reference rule and stays grey.
IDENTITY, ZERO, INK, HILITE = "#1a1a19", "#9b9b97", "#000000", "#e34948"

LINTHRESH = 0.02
TICKS = (-4, -2, -1, -0.5, -0.2, -0.05, 0, 0.05, 0.2, 0.5, 1)
PADJ, LFC, STATIC = 0.05, 1.0, 0.5
EXAMPLES = ("GAPDH", "COMT", "LRRFIP1")

#: Panel A matches fig03D_pooled_concordance's 2.082 x 3.669 in axes-box proportions.
A_RATIO = 2.082 / 3.669
#: ...but never narrower than this: two 10-pt tick labels need ~2 in between centres.
A_MIN_WIDTH = 2.1


def point_density(x, y, bins=200, sigma=4.0):
    """Local density at each point: a smoothed 2-D histogram read back at the points."""
    counts, xe, ye = np.histogram2d(x, y, bins=bins)
    try:
        from scipy.ndimage import gaussian_filter
        counts = gaussian_filter(counts, sigma=sigma, mode="nearest")
    except ImportError:
        pass
    ix = np.clip(np.digitize(x, xe) - 1, 0, counts.shape[0] - 1)
    iy = np.clip(np.digitize(y, ye) - 1, 0, counts.shape[1] - 1)
    return counts[ix, iy]


def housekeeping_ids(paths):
    """Ensembl ids from the HRT Atlas tables (semicolon-separated, `Ensembl` column)."""
    ids = set()
    for path in paths:
        if path and Path(path).exists():
            ids |= set(pd.read_csv(path, sep=";", dtype=str)["Ensembl"])
    return ids


class Geometry:
    """Every length in the figure, in INCHES, solved once and shared by both layouts.

    A single-panel page reuses the identical solved axes box, so a panel on its own is
    byte-for-byte the same plot area; tight_layout would equalise pages instead.
    """

    left, right, top = 0.72, 0.10, 0.32
    #: Under row 1: B's x label, then C's letter.
    row_gap, bottom = 0.95, 0.62
    gap = 0.85
    cbar_pad, cbar_w, cbar_label = 0.10, 0.16, 0.55

    def __init__(self, width=ps.PAGE_WIDTH_MAX, height=ps.PAGE_HEIGHT_MAX):
        for value, cap, what in ((width, ps.PAGE_WIDTH_MAX, "width"),
                                 (height, ps.PAGE_HEIGHT_MAX, "height")):
            if value > cap:
                raise SystemExit("page %s %.2f in exceeds PLOS maximum %.2f in"
                                 % (what, value, cap))
        # PLOS measures the file including the 2-pt TIFF border, so draw the page LESS it.
        border = 2 * ps.BORDER_PT / 72.0
        by_width = (width - border - self.left - self.gap - self.right) / (1.0 + A_RATIO)
        by_height = (height - border - self.top - self.row_gap - self.bottom) / 2.0
        self.height = min(by_width, by_height)
        self.boxes = (max(self.height * A_RATIO, A_MIN_WIDTH), self.height, self.height)
        row1 = self.left + self.boxes[0] + self.gap + self.boxes[1] + self.right
        row2 = self.left + self.boxes[2] + self.colorbar_width() + self.right
        self.width = max(row1, row2)

    @property
    def page_height(self):
        return self.top + 2 * self.height + self.row_gap + self.bottom

    #: The colorbar rides with panel C, so C's page carries it and the others do not.
    def colorbar_width(self):
        return self.cbar_pad + self.cbar_w + self.cbar_label

    def _adder(self, figure, width):
        def add(x0, y0, w):
            return figure.add_axes([x0 / width, y0 / self.page_height,
                                    w / width, self.height / self.page_height])
        return add

    def combined_page(self, figure_factory):
        figure = figure_factory((self.width, self.page_height))
        add = self._adder(figure, self.width)
        y1 = self.bottom + self.height + self.row_gap
        axes = [add(self.left, y1, self.boxes[0]),
                add(self.left + self.boxes[0] + self.gap, y1, self.boxes[1]),
                add(self.left, self.bottom, self.boxes[2])]
        cax = add(self.left + self.boxes[2] + self.cbar_pad, self.bottom, self.cbar_w)
        return figure, axes, cax

    def single_page(self, figure_factory, index):
        """One panel on its own page, at exactly its combined-page axes box."""
        box = self.boxes[index]
        bar = self.colorbar_width() if index == 2 else 0.0
        width = self.left + box + bar + self.right
        figure = figure_factory((width, self.page_height))
        add = self._adder(figure, width)
        axis = add(self.left, self.bottom, box)
        cax = add(self.left + box + self.cbar_pad, self.bottom, self.cbar_w) if index == 2 else None
        return figure, axis, cax


def dress(ax, letter):
    ax.tick_params(labelsize=ps.FONT_TICK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_axisbelow(True)
    ax.text(0.0, 1.02, letter, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=ps.FONT_PANEL_LETTER, fontweight="bold")


def draw_a(ax, corr):
    """A: route agreement per cell line."""
    from matplotlib.colors import to_rgba

    metrics = (("spearman_rho", "Spearman $\\rho$", ps.SPEARMAN_FILL, ps.SPEARMAN_LINE),
               ("pearson_r", "Pearson $r$", ps.PEARSON_FILL, ps.PEARSON_LINE))
    series = [corr[c].to_numpy(float) for c, _l, _f, _e in metrics]
    positions = np.arange(len(metrics))
    bp = ax.boxplot(series, positions=positions, widths=0.42, showfliers=False,
                    showcaps=False, patch_artist=True, zorder=10)
    for k, (_c, _l, fill, line) in enumerate(metrics):
        bp["boxes"][k].set(facecolor=to_rgba(fill, 0.35), edgecolor=line,
                           linewidth=ps.lw(1.4), zorder=11)
        bp["medians"][k].set(color=line, linewidth=ps.lw(2.0), zorder=13)
        for stroke in bp["whiskers"][2 * k:2 * k + 2]:
            stroke.set(color=line, linewidth=ps.lw(1.4), zorder=12)
        ax.annotate("%.3f" % np.nanmedian(series[k]), (positions[k], np.nanmax(series[k])),
                    textcoords="offset points", xytext=(0, 7), ha="center",
                    fontsize=ps.FONT_ANNOTATION, color=INK, zorder=15)
    ax.set_xticks(positions)
    ax.set_xticklabels([l for _c, l, _f, _e in metrics], fontsize=ps.FONT_TICK)
    ax.set_xlim(-0.6, len(metrics) - 0.4)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("correlation", fontsize=ps.FONT_LABEL)
    ax.grid(axis="y", alpha=0.13, zorder=0)
    dress(ax, "A")


def draw_b(ax, genes):
    """B: ranked delta TE with its confidence band."""
    mean = genes.dte_mean.to_numpy(float)
    lo, hi = genes.dte_ci_low.to_numpy(float), genes.dte_ci_high.to_numpy(float)
    order = np.argsort(mean)
    rank = np.arange(len(order))
    ax.axhline(0.0, color=ZERO, linewidth=ps.lw(0.9), zorder=2)
    ax.fill_between(rank, lo[order], hi[order], color=BAND, linewidth=0, alpha=0.8,
                    zorder=3, rasterized=True, label="95% CI")
    ax.plot(rank, mean[order], color=MEAN, linewidth=ps.lw(1.4), zorder=4,
            label="mean over cell lines")
    # symlog crossover sits BELOW the IQR so the tail lives in the logarithmic part.
    ax.set_yscale("symlog", linthresh=LINTHRESH, linscale=0.4)
    ax.set_yticks(list(TICKS))
    ax.set_yticklabels(["%g" % t for t in TICKS])
    ax.set_ylim(-6, 6)
    ax.set_xlim(0, len(order))
    ax.set_xlabel("transcript, ranked by $\\Delta$TE", fontsize=ps.FONT_LABEL)
    ax.set_ylabel("$\\Delta$TE (log$_2$)", fontsize=ps.FONT_LABEL)
    ax.grid(axis="y", alpha=0.13, zorder=0)
    leg = ax.legend(fontsize=ps.FONT_ANNOTATION, loc="upper left", **ps.KEY_FRAME)
    leg.get_frame().set_linewidth(ps.lw(0.8))
    for text in leg.get_texts():
        text.set_color(INK)
    dress(ax, "B")


def draw_c(figure, ax, cax, genes, marker, housekeeping=()):
    """C: the assay plane, with its density colorbar. Returns (labelled, dropped)."""
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    x, y = genes.drna_mean.to_numpy(float), genes.dribo_mean.to_numpy(float)
    te, padj = genes.dte_mean.to_numpy(float), genes.dte_padj.to_numpy(float)
    half = float(np.nanmax(np.abs(np.concatenate([x, y]))))
    pad = 0.03 * half

    density = point_density(x, y)
    dorder = np.argsort(density)
    cmap = LinearSegmentedColormap.from_list("white_viridis", WHITE_VIRIDIS, N=256)
    dlo, dhi, floor = float(density.min()), float(density.max()), 2.0 / 255.0
    vmin = max(float(np.exp((np.log(dlo) - floor * np.log(dhi)) / (1.0 - floor))), 1e-12)

    ax.axhline(0.0, color=ZERO, linewidth=ps.lw(0.8), zorder=2)
    ax.axvline(0.0, color=ZERO, linewidth=ps.lw(0.8), zorder=2)
    edge = np.array([-half - pad, half + pad])
    ax.plot(edge, edge, color=IDENTITY, linewidth=ps.lw(1.0), zorder=2.5,
            label="$y = x$  ($\\Delta$TE $= 0$)")
    dots = ax.scatter(x[dorder], y[dorder], c=density[dorder], s=7 * marker, cmap=cmap,
                      norm=LogNorm(vmin=vmin, vmax=dhi), linewidths=0, zorder=3,
                      rasterized=True)
    hi_mask = np.isfinite(padj) & (padj < PADJ) & (np.abs(te) > LFC)
    sel = hi_mask[dorder]
    ax.scatter(x[dorder][sel], y[dorder][sel], c=density[dorder][sel], cmap=cmap,
               norm=dots.norm, s=21 * marker, edgecolors=HILITE,
               linewidths=ps.lw(0.7), zorder=4, rasterized=True,
               label="$p_{adj}$ (BH) < %g, |$\\Delta$TE| > %g" % (PADJ, LFC))
    ax.set_xlim(*edge)
    ax.set_ylim(*edge)
    ax.set_aspect("equal")
    from matplotlib.ticker import MultipleLocator
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.set_xlabel("mean $\\Delta$RNA (log$_2$)", fontsize=ps.FONT_LABEL)
    ax.set_ylabel("mean $\\Delta$Ribo (log$_2$)", fontsize=ps.FONT_LABEL)
    ax.grid(alpha=0.13, zorder=0)
    dress(ax, "C")
    # Framed key inside the plane, upper left (empty corner); highlighted entry first.
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(handles)), key=lambda i: 0 if "p_{adj}" in labels[i] else 1)
    leg = ax.legend([handles[i] for i in order], [labels[i] for i in order],
                    fontsize=ps.FONT_ANNOTATION, loc="upper left", ncol=1, **ps.KEY_FRAME)
    leg.get_frame().set_linewidth(ps.lw(0.8))
    for text in leg.get_texts():
        text.set_color(INK)

    # Labels chosen by ROLE, not rank; ranking alone would show none of the atlas genes.
    hk = housekeeping_ids(housekeeping)
    base = np.array([i.split(".")[0] for i in genes.transcript_id])
    is_hk = np.isin(base, list(hk)) if hk else np.zeros(len(base), bool)
    idx = np.where(hi_mask)[0]

    def extreme(mask):
        cand = idx[mask[idx]]
        return cand[np.argmax(np.abs(te[cand]))] if len(cand) else None

    picked = [k for k in (extreme((np.abs(x) <= STATIC) & (np.abs(y) > STATIC)),
                          extreme(te > 0)) if k is not None]
    pool = idx[(np.abs(y[idx]) <= STATIC) & (np.abs(x[idx]) > STATIC) & is_hk[idx]]
    for k in pool[np.argsort(-np.abs(te[pool]))]:
        if len(picked) >= 5:
            break
        if k not in picked:
            picked.append(k)

    names = genes.gene_name.to_numpy()
    colour_of = {k: IDENTITY for k in picked}
    for gene in EXAMPLES:
        hits = np.where(names == gene)[0]          # exact match: the catalog holds GAPDHS too
        if len(hits) and int(hits[0]) not in colour_of:
            colour_of[int(hits[0])] = HILITE
            picked.append(int(hits[0]))

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    placed = [t.get_window_extent(renderer) for t in ax.texts]
    placed.append(leg.get_window_extent(renderer))
    disp = ax.transData.transform(np.column_stack([x, y]))
    px = figure.dpi / 72.0
    radius = np.where(hi_mask, np.sqrt(21 * marker / np.pi), np.sqrt(7 * marker / np.pi)) * px

    def hits_a_dot(b):
        return bool(np.any((disp[:, 0] >= b.x0 - radius) & (disp[:, 0] <= b.x1 + radius)
                           & (disp[:, 1] >= b.y0 - radius) & (disp[:, 1] <= b.y1 + radius)))

    inv = ax.transData.inverted()

    def crosses_identity(b, k):
        """A label on the far side of y = x from its point reads as the wrong gene.

        Only for points clearly OFF the line: an on-line point has no far side.
        """
        if abs(y[k] - x[k]) <= LFC:
            return False
        corners = inv.transform([(b.x0, b.y0), (b.x1, b.y0), (b.x0, b.y1), (b.x1, b.y1)])
        side = np.sign(corners[:, 1] - corners[:, 0])
        return bool(np.any(side != np.sign(y[k] - x[k])))

    HORIZONTAL = [(1, 0, "left", "center"), (-1, 0, "right", "center"),
                  (1, 1, "left", "bottom"), (1, -1, "left", "top"),
                  (-1, 1, "right", "bottom"), (-1, -1, "right", "top"),
                  (0, 1, "center", "bottom"), (0, -1, "center", "top")]
    # A point on the delta-Ribo ~ 0 band is offered vertical positions first.
    VERTICAL = [(0, -1, "center", "top"), (0, 1, "center", "bottom"),
                (1, -1, "left", "top"), (1, 1, "left", "bottom"),
                (-1, -1, "right", "top"), (-1, 1, "right", "bottom"),
                (1, 0, "left", "center"), (-1, 0, "right", "center")]
    dropped = []
    for k in picked:
        name = names[k]
        if not isinstance(name, str) or not name.strip():
            continue
        order_k = VERTICAL if abs(y[k]) <= STATIC else HORIZONTAL
        done = False
        for step in (0.04, 0.08, 0.14, 0.22, 0.32, 0.45):
            for sx, sy, ha, va in order_k:
                d = step * half
                txt = ax.text(x[k] + sx * d, y[k] + sy * d, name, ha=ha, va=va,
                              fontsize=ps.FONT_ANNOTATION,
                              color=colour_of.get(k, IDENTITY), zorder=8)
                bbox = txt.get_window_extent(renderer)
                if (any(bbox.overlaps(b) for b in placed) or hits_a_dot(bbox)
                        or crosses_identity(bbox, k)):
                    txt.remove()
                    continue
                placed.append(bbox)
                ax.plot([x[k], x[k] + sx * d * 0.9], [y[k], y[k] + sy * d * 0.9],
                        color=colour_of.get(k, IDENTITY), linewidth=ps.lw(0.6), zorder=7)
                done = True
                break
            if done:
                break
        if not done:
            dropped.append(name)

    bar = figure.colorbar(dots, cax=cax)
    # Not "transcripts per dot": one dot IS one transcript; the bar counts its neighbours.
    bar.set_label("local transcript density", fontsize=ps.FONT_ANNOTATION, color=INK)
    bar.ax.tick_params(labelsize=ps.FONT_TICK, length=0, labelcolor=INK)
    bar.outline.set_visible(False)
    return [str(names[k]) for k in picked], dropped


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--per-gene-delta", required=True, type=Path,
                        help="per_gene_delta.tsv from te_statistics.R")
    parser.add_argument("--route-correlation", required=True, type=Path,
                        help="route_correlation.tsv from te_statistics.R")
    parser.add_argument("--housekeeping-genes", type=Path,
                        help="HRT Atlas Housekeeping_GenesHuman.csv (panel C labels)")
    parser.add_argument("--housekeeping-transcripts", type=Path,
                        help="HRT Atlas Housekeeping_TranscriptsHuman.csv (panel C labels)")
    parser.add_argument("--panel", choices=PANELS, default="combined",
                        help="the combined page, or one panel on its own page")
    parser.add_argument("--width", type=float, default=ps.PAGE_WIDTH_MAX,
                        help="page width cap in inches (PLOS maximum %.2f)" % ps.PAGE_WIDTH_MAX)
    parser.add_argument("--height", type=float, default=ps.PAGE_HEIGHT_MAX,
                        help="page height cap in inches (PLOS maximum %.2f)" % ps.PAGE_HEIGHT_MAX)
    parser.add_argument("--output", required=True, type=Path, help="path stem, no extension")
    parser.add_argument("--format", dest="formats", default="pdf", help="pdf,png,tif")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    for path in (args.per_gene_delta, args.route_correlation):
        if not path.exists():
            raise SystemExit("%s does not exist; run code/te_route/normalization.R then "
                             "te_statistics.R (or use the shipped data/te_route/tables)" % path)
    genes = pd.read_csv(args.per_gene_delta, sep="\t")
    corr = pd.read_csv(args.route_correlation, sep="\t")
    housekeeping = (args.housekeeping_genes, args.housekeeping_transcripts)

    import matplotlib.pyplot as plt

    ps.apply_rcparams()
    formats = ps.resolve_formats(args.formats)
    marker = ps.MARKER_AREA
    geom = Geometry(args.width, args.height)
    labelled, dropped = None, None

    def new_figure(size):
        return plt.figure(figsize=size)

    if args.panel == "combined":
        figure, axes, cax = geom.combined_page(new_figure)
        draw_a(axes[0], corr)
        draw_b(axes[1], genes)
        labelled, dropped = draw_c(figure, axes[2], cax, genes, marker, housekeeping)
        written = ps.save(figure, args.output, formats, args.force, tight=False)
    else:
        index = PANELS.index(args.panel) - 1
        figure, axis, cax = geom.single_page(new_figure, index)
        if index == 0:
            draw_a(axis, corr)
        elif index == 1:
            draw_b(axis, genes)
        else:
            labelled, dropped = draw_c(figure, axis, cax, genes, marker, housekeeping)
        written = ps.save(figure, args.output, formats, args.force, tight=True)
    plt.close(figure)

    if labelled:
        print("labelled: %s" % ", ".join(labelled))
    if dropped:
        print("NO ROOM for: %s" % ", ".join(dropped))
    print("axes boxes A %.2f  B %.2f  C %.2f in; page %.2f x %.2f in"
          % (geom.boxes[0], geom.boxes[1], geom.boxes[2], geom.width, geom.page_height))
    for path in written:
        print("wrote %s" % path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
