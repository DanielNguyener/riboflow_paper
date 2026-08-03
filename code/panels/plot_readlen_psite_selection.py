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

def draw(prepared, axes_size=(7.0, 7.8)):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    sys.path.insert(0, str(HERE))
    import _fig02_common as common
    import panel_style as ps

    ps.apply_rcparams()
    samples, lengths = prepared["samples"], prepared["lengths"]
    width, height = axes_size
    left, bottom, right, top = 1.25, 1.55, 0.25, 0.15
    figure, axis = plt.subplots(figsize=(left + width + right, bottom + height + top))
    axis.set_position([left / (left + width + right), bottom / (bottom + height + top),
                       width / (left + width + right), height / (bottom + height + top)])
    from matplotlib.colors import BoundaryNorm
    colormap = ListedColormap([BLANK_C, AGREE_C, DIFFER_C, GENOME_C, TXOME_C])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1.0), colormap.N)
    common.draw_grid(axis, prepared["status"].astype(float), colormap, norm)
    for i in range(len(samples)):
        for j in range(len(lengths)):
            if prepared["labels"][i, j]:
                axis.text(common.cell_centre(j), common.cell_centre(i),
                          prepared["labels"][i, j], ha="center", va="center",
                          fontsize=ps.FONT_ANNOTATION,
                          color="white" if prepared["status"][i, j] else "black",
                          linespacing=0.9)
    common.style_grid(axis, samples, lengths, prepared["gsm"])
    handles = [plt.Rectangle((0, 0), 1, 1, color=c)
               for c in (AGREE_C, DIFFER_C, GENOME_C, TXOME_C, BLANK_C)]
    axis.legend(handles, STATUS_LABELS, loc="upper center", bbox_to_anchor=(0.5, -0.13),
                ncol=3, fontsize=ps.FONT_TICK, frameon=False, handlelength=1.2,
                handletextpad=0.5)
    return figure, axis

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--qc-genome", required=True, type=Path)
    parser.add_argument("--qc-txome", required=True, type=Path)
    parser.add_argument("--samples-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
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

    figure, _axis = draw(prepared)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
