#!/usr/bin/env python3
"""Figure 2 A -- read-length selection and P-site offsets, genome versus transcriptome."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

BLANK_C, AGREE_C, DIFFER_C, GENOME_C, TXOME_C = (
    "#eeeeee", "#3182bd", "#984ea3", "#2ca25f", "#d62728")
STATUS_LABELS = ("both selected, equal offset", "both selected, different offsets",
                 "genome-only", "transcriptome-only", "not selected")

def prepare(qc_genome, qc_txome, samples_csv):
    """The status grid, its cell labels, and the agreement counts the caption quotes."""
    sys.path.insert(0, str(HERE))
    import _fig02_common as common

    samples, lengths, genome, txome = common.load_pair(qc_genome, qc_txome, "psite_offset")
    status = np.zeros((len(samples), len(lengths)), dtype=int)
    labels = np.empty((len(samples), len(lengths)), dtype=object)
    n_agree = n_disagree = identical_sets = 0

    for i in range(len(samples)):
        genome_selected, txome_selected = set(), set()
        for j in range(len(lengths)):
            gv, tv = genome[i, j], txome[i, j]
            has_g, has_t = not np.isnan(gv), not np.isnan(tv)
            if has_g:
                genome_selected.add(lengths[j])
            if has_t:
                txome_selected.add(lengths[j])
            if not has_g and not has_t:
                labels[i, j] = ""
                continue
            if has_g and has_t and int(gv) == int(tv):
                status[i, j], labels[i, j] = 1, "%d" % int(gv)
                n_agree += 1
            else:
                n_disagree += 1
                if has_g and has_t:
                    status[i, j], labels[i, j] = 2, "%d|%d" % (int(gv), int(tv))
                elif has_g:
                    status[i, j], labels[i, j] = 3, "%d" % int(gv)
                else:
                    status[i, j], labels[i, j] = 4, "%d" % int(tv)
        identical_sets += genome_selected == txome_selected

    return {"samples": samples, "lengths": lengths, "status": status, "labels": labels,
            "n_agree": n_agree, "n_disagree": n_disagree,
            "identical_sets": identical_sets,
            "gsm": common.gsm_labels(samples_csv, samples),
            "sources": {"qc_genome": str(qc_genome), "qc_txome": str(qc_txome),
                        "samples_csv": str(samples_csv)}}

def draw(prepared, axes_size=None, margins=None, type_scale="large", legend_ncol=2):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    sys.path.insert(0, str(HERE))
    import _fig02_common as common
    import panel_style as ps

    ps.apply_rcparams()
    sizes = common.grid_type(type_scale)
    samples, lengths = prepared["samples"], prepared["lengths"]
    # Shared with fig02B (`_fig02_common.MARGINS`); right gutter deliberately blank here.
    width, height = axes_size or common.AXES_SIZE
    left, bottom, right, top = margins or common.MARGINS
    fig_w, fig_h = left + width + right, bottom + height + top
    figure, axis = plt.subplots(figsize=(fig_w, fig_h))
    axis.set_position([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])
    from matplotlib.colors import BoundaryNorm
    colormap = ListedColormap([BLANK_C, AGREE_C, DIFFER_C, GENOME_C, TXOME_C])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1.0), colormap.N)
    common.draw_grid(axis, prepared["status"].astype(float), colormap, norm)
    for i in range(len(samples)):
        for j in range(len(lengths)):
            if prepared["labels"][i, j]:
                axis.text(common.cell_centre(j), common.cell_centre(i),
                          prepared["labels"][i, j], ha="center", va="center",
                          fontsize=sizes["annotation"],
                          color="white" if prepared["status"][i, j] else "black",
                          linespacing=0.9)
    common.style_grid(axis, samples, lengths, prepared["gsm"], sizes)
    # Only the statuses that occur in the grid; the caption names the full scheme.
    present = set(np.unique(prepared["status"]))
    entries = [(AGREE_C, STATUS_LABELS[0], 1), (DIFFER_C, STATUS_LABELS[1], 2),
               (GENOME_C, STATUS_LABELS[2], 3), (TXOME_C, STATUS_LABELS[3], 4),
               (BLANK_C, STATUS_LABELS[4], 0)]
    entries = [(c, label) for c, label, code in entries if code in present and code != 0]
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c, _ in entries]
    labels = [label for _, label in entries]
    # Two columns: ncol=3 overruns the 7 in axes. Anchor is MEASURED, not an axes fraction.
    legend = ps.legend_below(axis, handles, labels, ncol=legend_ncol,
                             fontsize=sizes["tick"], frameon=False, handlelength=1.2,
                             handletextpad=0.5)
    return figure, axis, legend

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--qc-genome", required=True, type=Path)
    parser.add_argument("--qc-txome", required=True, type=Path)
    parser.add_argument("--samples-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--axes-size", nargs=2, type=float, metavar=("W", "H"),
                        help="grid size in inches (default: _fig02_common.AXES_SIZE)")
    parser.add_argument("--margins", nargs=4, type=float, metavar=("L", "B", "R", "T"),
                        help="margins in inches (default: _fig02_common.MARGINS)")
    parser.add_argument("--type-scale", choices=("large", "base"), default="large",
                        help="large: standalone panel type; base: journal-page type (8-12 pt)")
    parser.add_argument("--legend-ncol", type=int, default=2)
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps

    prepared = prepare(args.qc_genome, args.qc_txome, args.samples_csv)
    total = prepared["n_agree"] + prepared["n_disagree"]
    print("[panel] %d libraries x %d read lengths"
          % (len(prepared["samples"]), len(prepared["lengths"])))
    print("[panel] phase-1 sets identical in %d/%d libraries"
          % (prepared["identical_sets"], len(prepared["samples"])))
    print("[panel] %d cells selected by either method; agree %d (%.1f%%), disagree %d"
          % (total, prepared["n_agree"], 100.0 * prepared["n_agree"] / total,
             prepared["n_disagree"]))

    figure, _axis, legend = draw(
        prepared, tuple(args.axes_size) if args.axes_size else None,
        tuple(args.margins) if args.margins else None, args.type_scale, args.legend_ncol)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force,
                      extra_artists=[legend], tight=False)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
