#!/usr/bin/env python3
"""Figure 5 A-D at PLOS page size: the four cohort panels re-rendered with 8 pt type.

Rebinds `fig05_common` box constants and `panel_style.FONT_*` in-process BEFORE importing
the generator modules; one panel per process (`--panel A`).
"""
from __future__ import annotations

import argparse
import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

#: The cohort panels' shared plot box, shrunk from the published 470 / 40 pt.
BOX_HEIGHT_PT = 240.0
BOX_TOP_OFFSET_PT = 22.0

#: `plot_read_id_union.SEGMENTS` colours with the verified route-status names.
UNION_KEY = (
    ("Shared, unique", "#a6d96a"),
    ("Shared, multimapping", "#1a7d1a"),
    ("Genome only, unique", "#7fb9da"),
    ("Genome only, multimapping", "#0d57a1"),
    ("Transcriptome only", "#cc3d3d"),
)

MODULES = {"A": "plot_route_read_counts", "B": "plot_read_id_union",
           "C": "plot_multimap_biotype", "D": "plot_nonselected_isoform_reach"}


def die(message):
    raise SystemExit("error: %s" % message)


def patched_panels(font_pt):
    common = importlib.import_module("fig05_common")
    common.STACK_AXES_HEIGHT_PT = BOX_HEIGHT_PT
    common.STACK_AXES_TOP_OFFSET_PT = BOX_TOP_OFFSET_PT
    style = importlib.import_module("panel_style")
    style.FONT_TITLE = style.FONT_LABEL = style.FONT_TICK = font_pt
    style.FONT_ANNOTATION = style.FONT_INSET = font_pt
    return common


def render_panel(letter, width_pt, height_in, inputs, font_pt, output, formats=("pdf",),
                 force=False):
    """One cohort panel to `output.<fmt>`, final labels applied in-process."""
    patched_panels(font_pt)
    style = importlib.import_module("panel_style")
    side = importlib.import_module("_fig05_side_panel")
    common = importlib.import_module("fig05_common")
    common.median_iqr = lambda values: ""       # summary statistics belong in the caption
    figsize = (width_pt / 72.0, height_in)
    extras = []
    module = importlib.import_module(MODULES[letter])

    def rp(key):
        if not inputs.get(key):
            die("panel %s needs --%s" % (letter, key.replace("_", "-")))
        return inputs[key]

    if letter == "A":
        prepared = module.prepare(rp("taxonomy"), rp("samples_csv"))
        values = prepared["genome_m"] - prepared["txome_m"]
        figure, axis = side.draw_side_panel(
            values, prepared["labels"], style.GENOME,
            "Difference in aligned\nread IDs", "", figsize, True)
        for artist in list(axis.texts):          # per-bar values: a table, not a panel
            artist.remove()
        # Largest gap is 1.286 M; a final tick at 1.5 keeps every bar short of it.
        axis.set_xlim(0, 1.55)
        axis.set_xticks([0.0, 0.5, 1.0, 1.5])
        axis.set_xlabel("Genome −\ntranscriptome\n(millions)", fontsize=font_pt)
        figure.subplots_adjust(left=0.47, right=0.96)

    elif letter == "B":
        prepared = module.prepare(rp("taxonomy"), rp("samples_csv"))
        figure, axis, _ = module.draw(prepared, figsize, show_labels=False, show_key=False)
        axis.set_title("Composition of the\nread-ID union", fontsize=font_pt, loc="left",
                       fontweight="normal")
        axis.set_xlabel("Read IDs in\nalignment-route union (%)", fontsize=font_pt)
        axis.set_xlim(0, 128)
        axis.set_xticks([0, 20, 40, 60, 80, 100])
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=colour, edgecolor="white", linewidth=0.4, label=label)
                   for label, colour in UNION_KEY]
        extras.append(style.legend_below(axis, handles=handles, ncol=1,
                                         fontsize=font_pt, handlelength=1.0,
                                         labelspacing=0.25, borderpad=0.0, pad_pt=4.0))

    elif letter in ("C", "D"):
        if letter == "C":
            prepared = module.prepare(rp("tie_master"), rp("taxonomy"))
            colour, title = "#1a7d1a", "Protein-coding–\npseudogene ties"
            xlabel = "Shared genome-\nmultimapping\nreads (%)"
        else:
            prepared = module.prepare(rp("reach_master"), rp("taxonomy"))
            colour, title = "#7fb9da", "Overlap with omitted\nalternative exons"
            xlabel = "Uniquely mapped\ngenome-only\nreads (%)"
        figure, axis = side.draw_side_panel(
            prepared["values"], prepared["labels"], colour, title, xlabel, figsize, False)
        if letter == "C":
            axis.set_xticks(range(0, 81, 20))
            axis.set_xticks(range(10, 80, 20), minor=True)
    else:
        die("unknown cohort panel %s" % letter)

    # tight=False, as the generators themselves save: a crop would undo the row alignment.
    return style.save(figure, output, formats, force=force,
                      extra_artists=extras or None, tight=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", required=True, choices=sorted(MODULES))
    parser.add_argument("--taxonomy", help="taxonomy_all.tsv (A, B, C, D)")
    parser.add_argument("--samples-csv", help="S1 Table samples.csv (A, B)")
    parser.add_argument("--tie-master", help="multimap_tie_biotype_all.tsv (C)")
    parser.add_argument("--reach-master", help="genome_anchored_reach_all.tsv (D)")
    parser.add_argument("--width-pt", type=float, required=True, help="page width")
    parser.add_argument("--height-in", type=float, required=True, help="page height")
    parser.add_argument("--font-pt", type=float, default=8.0)
    parser.add_argument("--output", required=True, help="path stem, no extension")
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    inputs = {"taxonomy": args.taxonomy, "samples_csv": args.samples_csv,
              "tie_master": args.tie_master, "reach_master": args.reach_master}
    style = importlib.import_module("panel_style")
    written = render_panel(args.panel, args.width_pt, args.height_in, inputs, args.font_pt,
                           args.output, style.resolve_formats(args.formats), args.force)
    for path in written:
        print("wrote %s" % path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
