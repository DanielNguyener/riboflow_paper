#!/usr/bin/env python3
"""The exon-compressed locus panel, drawn from the compact locus artifact.

Reads ONLY `locus_<GENE>.npz` + `.json`; exons to scale, introns collapsed. Blocks absent
from the selected reference are shaded, never drawn as zero. Run with `python` (3.9).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
from intervals import merge  # noqa: E402
from panel_style import die  # noqa: E402

#: Type scale, module constants so one caller can rebind them (PLOS allows 8-12 pt).
FONT_TITLE = 11.0
FONT_LABEL = 9.0
FONT_XLABEL = 10.0
FONT_MARK = 9.0
FONT_TICK = 10.0
FONT_MODEL = 7.5

#: The coverage is stacked by the Figure 5A read populations and painted in 5A's key
#: colours (`plot_gene_read_partition.ROUTE7_KEY`), so one colour means one population in
#: both panels. Route identity is carried by the mirrored layout and the boxed marks, not
#: by colour. Stacking order: baseline first (the values are the npz array suffixes).
GENOME_LAYERS = (("shared_unique", "#a6d96a", "SH-U"),
                 ("shared_multi", "#1a7d1a", "SH-M"),
                 ("genome_only", "#7fb9da", "GO-U"),
                 ("genome_only_multi", "#0d57a1", "GO-M"))
TXOME_LAYERS = (("shared_unique", "#a6d96a", "SH-U"),
                ("shared_multi", "#1a7d1a", "SH-M"),
                ("txome_only", "#cc3d3d", "TO"))
#: Neutral grey for the region the reference omits: a statement about the annotation, not
#: a measurement, so it carries no signal colour.
ABSENT_COLOUR = "#9a9a9a"
#: Both models in black: they are structure, not signal; the labels say which is which.
SELECTED_COLOUR = "#000000"
ALT_COLOUR = "#000000"




class SplicedAxis:
    """Genomic positions -> plot coordinates, exons to scale, introns a constant gap."""

    def __init__(self, blocks, strand, gap):
        self.blocks = list(blocks) if strand == "+" else list(reversed(blocks))
        self.strand = strand
        self.gap = gap
        self.offsets = []
        cursor = 0.0
        for start, end in self.blocks:
            self.offsets.append(cursor)
            cursor += (end - start) + gap
        self.width = cursor - gap if self.blocks else 0.0

    def span(self, start, end):
        out = []
        for (b_start, b_end), offset in zip(self.blocks, self.offsets):
            lo, hi = max(start, b_start), min(end, b_end)
            if hi <= lo:
                continue
            if self.strand == "+":
                out.append((offset + lo - b_start, offset + hi - b_start))
            else:
                out.append((offset + b_end - hi, offset + b_end - lo))
        return out

    def base_axis(self):
        xs, gs = [], []
        for (start, end), offset in zip(self.blocks, self.offsets):
            n = end - start
            xs.append(offset + np.arange(n))
            gs.append(np.arange(start, end) if self.strand == "+"
                      else np.arange(end - 1, start - 1, -1))
        if not xs:
            return np.array([]), np.array([])
        return np.concatenate(xs), np.concatenate(gs)


def load_locus(npz_path, meta_path):
    with np.load(npz_path) as data:
        arrays = {key: data[key] for key in data.files}
    with open(meta_path) as handle:
        meta = json.load(handle)
    return arrays, meta


def apply_rcparams():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "font.size": FONT_TICK,
                         "xtick.labelsize": FONT_TICK,
                         "ytick.labelsize": FONT_TICK})


def render_locus(ax_ribo, ax_model, arrays, meta, title=True, xlabel=True,
                 model_labels=True, title_inside=False, ylabel=None,
                 model_names=None, no_xticks=False, marks_at_axis=False):
    """Draw one locus into axes the CALLER owns."""
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    strand = meta["strand"]
    sel_exons = [tuple(int(v) for v in row) for row in arrays["sel_exons"]]
    alt_exons = [tuple(int(v) for v in row) for row in arrays["alt_exons"]]
    absent_blocks = [tuple(int(v) for v in row) for row in arrays["absent_blocks"]]
    union = merge(list(sel_exons) + list(alt_exons))
    axis = SplicedAxis(union, strand, meta["intron_gap_plot_units"])
    xs, gs = axis.base_axis()
    if not np.array_equal(gs, arrays["genomic_position"]):
        die("the derived coverage is not on this locus' exonic bases")
    genome_cov = np.asarray(arrays["genome_cov"], dtype=float)
    txome_cov = np.asarray(arrays["txome_cov"], dtype=float)

    axes = [ax_ribo, ax_model]

    # Mirrored: genome above the line, transcriptome below, one shared magnitude scale; the
    # y tick labels are absolute values. Each half is a stack of its 5A populations.
    for sign, layers, track in ((1.0, GENOME_LAYERS, "genome_cov"),
                                (-1.0, TXOME_LAYERS, "txome_cov")):
        base = np.zeros_like(genome_cov)
        for key, colour, label in layers:
            name = "%s_%s" % (track, key)
            if name not in arrays:
                die("%s lacks the %r layer; rebuild it with build_locus_data.py"
                    % (meta.get("gene", "the locus artifact"), name))
            layer = np.asarray(arrays[name], dtype=float)
            top = base + layer
            ax_ribo.fill_between(xs, sign * base, sign * top, step="mid", color=colour,
                                 linewidth=0, zorder=1, label=label)
            # A one-base P-site spike is narrower than a device pixel at page width, so
            # the fill alone can vanish: stem each layer over ITS OWN extent (never the
            # whole envelope, which would paint an empty layer's colour everywhere).
            covered = layer > 0
            if covered.any():
                ax_ribo.vlines(xs[covered], sign * base[covered], sign * top[covered],
                               color=colour, lw=0.9, zorder=3)
            base = top
        if not np.array_equal(base, np.asarray(arrays[track], dtype=float)):
            die("%s layers do not sum to the track" % track)
    ax_ribo.axhline(0, color="black", lw=0.6, zorder=4)

    # Each half cropped to its own peak; on ONE axes units-per-pixel is shared, so both
    # halves render at exactly the same scale.
    g_peak = float(genome_cov.max()) if genome_cov.size else 0.0
    t_peak = float(txome_cov.max()) if txome_cov.size else 0.0
    unit = max(g_peak, t_peak, 1.0) * 0.06
    ax_ribo.set_ylim(-(t_peak + unit), g_peak + unit)
    ax_ribo.yaxis.set_major_locator(MaxNLocator(nbins=7, integer=True))
    ax_ribo.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: "%g" % abs(v)))
    ax_ribo.set_ylabel(ylabel or ("%s coverage" % ("P-site" if meta["signal"] == "psite"
                                                   else "footprint")),
                       fontsize=FONT_LABEL)
    ax_ribo.grid(axis="y", alpha=0.15)

    label_box = dict(boxstyle="round,pad=0.35", fc="white", ec="#555", lw=0.9)
    mark_top, mark_bottom = "genome", "transcriptome"
    if marks_at_axis:
        # Unboxed, at the right edge, hugging the mirror line: genome sits just above it and
        # transcriptome just below, so each name is read against its own half.
        from matplotlib.transforms import blended_transform_factory
        trans = blended_transform_factory(ax_ribo.transAxes, ax_ribo.transData)
        for dy, va, text in ((2.0, "bottom", mark_top), (-2.0, "top", mark_bottom)):
            ax_ribo.annotate(text, xy=(0.995, 0.0), xycoords=trans, xytext=(0.0, dy),
                             textcoords="offset points", ha="right", va=va,
                             fontsize=FONT_MARK, zorder=8)
    else:
        # Boxed, top and bottom left inside the axes.
        ax_ribo.text(0.012, 0.94, mark_top, transform=ax_ribo.transAxes, ha="left",
                     va="top", fontsize=FONT_MARK, bbox=label_box, zorder=8)
        ax_ribo.text(0.012, 0.06, mark_bottom, transform=ax_ribo.transAxes, ha="left",
                     va="bottom", fontsize=FONT_MARK, bbox=label_box, zorder=8)
    if title:
        heading = ("%s  %s" % (meta["gene"], meta["locus"]["label"])
                   if title is True else title)
        if title_inside:
            ax_ribo.text(0.995, 0.96, heading, transform=ax_ribo.transAxes,
                         ha="right", va="top", fontsize=FONT_TITLE, zorder=9)
        else:
            ax_ribo.set_title(heading, fontsize=FONT_TITLE, loc="left", pad=2.0)

    # Isoform models; names are y-tick labels so the shared left margin makes room.
    models = [(meta["selected_transcript"], sel_exons, SELECTED_COLOUR,
               "selected\n(txome reference)")]
    if alt_exons:
        models.append((meta["alternative_transcript"], alt_exons, ALT_COLOUR,
                       "alternative\nGENCODE isoform"))
    ticks, tick_labels = [], []
    for row, (tid, exons, colour, label) in enumerate(models):
        y = (0.5 if len(models) == 1 else 0.68) - row * 0.40
        spans = []
        for start, end in exons:
            for x0, x1 in axis.span(start, end):
                ax_model.add_patch(Rectangle((x0, y - 0.10), max(x1 - x0, 1.0), 0.20,
                                             facecolor=colour, edgecolor="none"))
                spans.append((x0, x1))
        spans.sort()
        for (_a, a_end), (b_start, _b) in zip(spans, spans[1:]):
            if b_start > a_end:
                ax_model.plot([a_end, b_start], [y, y], linestyle=(0, (4, 3)),
                              color=colour, linewidth=0.9)
        ticks.append(y)
        if model_names is not None:
            tick_labels.append(model_names[row])
        else:
            tick_labels.append(("%s\n%s" % (label, tid)) if model_labels
                               else label.split("\n")[0])
    ax_model.set_ylim(0, 1)
    ax_model.set_yticks(ticks)
    ax_model.set_yticklabels(tick_labels, fontsize=FONT_MODEL)
    ax_model.tick_params(axis="y", length=0)
    if xlabel:
        ax_model.set_xlabel(xlabel if isinstance(xlabel, str) else
                            "spliced locus position, 5' to 3'   "
                            "(exons to scale, introns collapsed)",
                            fontsize=FONT_XLABEL, labelpad=2.0)
    if no_xticks:
        # The coordinate is exon-compressed; a number on that axis is neither a genomic
        # nor a transcript position, so ticks would assert a scale the axis does not have.
        for ax in axes:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)

    # The shading that ties the tracks together: sequence the selected reference lacks.
    for start, end in absent_blocks:
        for x0, x1 in axis.span(start, end):
            for ax in axes:
                ax.axvspan(x0, x1, color=ABSENT_COLOUR, alpha=0.22, linewidth=0, zorder=0)

    for ax in axes:
        ax.set_xlim(-axis.width * 0.01, axis.width * 1.01)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    ax_model.spines["left"].set_visible(False)
    return axis


def draw(args, arrays, meta):
    # Size overrides BEFORE `apply_rcparams`, which copies FONT_TICK into the rcParams.
    if args.font_size:
        for name in ("FONT_TITLE", "FONT_LABEL", "FONT_XLABEL", "FONT_MARK",
                     "FONT_MODEL", "FONT_TICK"):
            globals()[name] = args.font_size
    if args.title_size:
        globals()["FONT_TITLE"] = args.title_size

    apply_rcparams()
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(
        2, 1, figsize=(args.figsize[0], args.figsize[1]), sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.18})

    heading = ("%s  %s   %s   spliced locus, introns not to scale"
               % (meta["gene"], meta["locus"]["label"], meta["sample"]))
    if args.title is not None:
        heading = args.title or False
    render_locus(axes[0], axes[-1], arrays, meta,
                 title=heading, model_labels=not args.compact,
                 title_inside=args.compact and not args.title_outside,
                 ylabel=args.ylabel, xlabel=args.xlabel or True,
                 model_names=tuple(args.model_names.split(",")) if args.model_names else None,
                 no_xticks=args.no_xticks, marks_at_axis=args.marks_at_axis)

    # Declared margins, not `tight_layout`: the model track's tick labels are two lines of
    # text on a patch-only axes, which tight_layout measures badly.
    figure.subplots_adjust(left=0.155, right=0.985, top=0.93, bottom=0.11)
    outputs = ["%s.%s" % (args.output, suffix.strip()) for suffix in args.format.split(",")]
    existing = [o for o in outputs if os.path.exists(o)]
    if existing and not args.force:
        die("refusing to overwrite %s; pass --force" % ", ".join(existing))
    for out in outputs:
        figure.savefig(out, bbox_inches="tight")
        print("wrote %s" % out)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--derived-npz", required=True, help="locus_<GENE>.npz")
    parser.add_argument("--derived-meta", required=True, help="locus_<GENE>.json")
    parser.add_argument("--ylabel")
    parser.add_argument("--xlabel")
    parser.add_argument("--model-names", help="model-track names, comma-separated, selected first")
    parser.add_argument("--marks-at-axis", action="store_true",
                        help="unboxed route marks at the right edge, hugging the mirror "
                             "line: the top mark just above it, the bottom mark just below")
    parser.add_argument("--title-outside", action="store_true")
    parser.add_argument("--no-xticks", action="store_true")
    parser.add_argument("--font-size", type=float)
    parser.add_argument("--title-size", type=float)
    parser.add_argument("--title", help="replace the default heading; empty string draws none")
    parser.add_argument("--compact", action="store_true")
    parser.add_argument("--figsize", nargs=2, type=float, default=(13.0, 7.0))
    parser.add_argument("--output", required=True, help="path stem, no extension")
    parser.add_argument("--format", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    arrays, meta = load_locus(args.derived_npz, args.derived_meta)
    draw(args, arrays, meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
