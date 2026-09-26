#!/usr/bin/env python3
"""Figure 4 C and D -- one cohort share per cell line, drawn as a side panel.

C: protein-coding to pseudogene tie share of genome multimapper reads (`--tie-master`).
D: omitted alternative-exon overlap of genome-only unique reads (`--reach-master`):
   the Figure 5A gene-level test, cohort-wide -- a read counts once when an aligned block
   of its primary genomic alignment overlaps exonic sequence absent from the overlapped
   gene's selected transcript.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

#: Per panel: the master's flag, its required columns, the share, and the panel's dress.
PANELS = {
    "C": {"master": "tie_master",
          "required": ("sample", "pct_cross_pp_pc", "pct_cross_pc_pp"),
          "share": lambda f: f["pct_cross_pp_pc"] + f["pct_cross_pc_pp"],
          "stat": "cross tie", "colour": "#1a7d1a",
          "title": "protein-coding ↔\npseudogene tie",
          "xlabel": "% of multimapper\nreads (primary)"},
    "D": {"master": "reach_master",
          "required": ("sample", "n_omitted_exon_overlap", "n_gU_tA"),
          "share": lambda f: 100.0 * f["n_omitted_exon_overlap"] / f["n_gU_tA"],
          "stat": "omitted-exon overlap", "colour": "#7fb9da",
          "title": "omitted\nalternative exon",
          "xlabel": "% of genome-only\nunique reads"},
}

def prepare(letter, master, taxonomy, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import cohort_common as common
    import panel_style as ps

    spec = PANELS[letter]
    frame = pd.read_csv(master, sep="\t")
    ps.require_columns(frame, spec["required"], str(master))
    frame = frame.copy()
    frame["share"] = spec["share"](frame)

    order = common.sample_order(common.load_taxonomy(taxonomy))
    lookup = frame.set_index("sample")["share"]
    values = np.array([lookup.get(s, np.nan) for s in order], dtype=float)
    labels_map = common.load_labels(samples_csv)
    return {"order": order, "values": values,
            "order_provenance": "shared cohort order, derived from %s" % taxonomy,
            "labels": [labels_map.get(s, s) for s in order]}

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", required=True, choices=sorted(PANELS))
    parser.add_argument("--tie-master", type=Path, help="multimap_tie_biotype_all.tsv (C)")
    parser.add_argument("--reach-master", type=Path, help="genome_anchored_reach_all.tsv (D)")
    parser.add_argument("--taxonomy", required=True, type=Path,
                        help="the taxonomy master, for the shared Figure-5 cohort order")
    parser.add_argument("--samples-csv", type=Path)
    parser.add_argument("--show-labels", action="store_true")
    parser.add_argument("--figsize", nargs=2, type=float, default=(3.0, 8.0))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--format", dest="formats", default="pdf")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(HERE))
    import panel_style as ps
    from _cohort_side_panel import draw_side_panel

    spec = PANELS[args.panel]
    master = getattr(args, spec["master"])
    if not master:
        parser.error("panel %s needs --%s" % (args.panel, spec["master"].replace("_", "-")))
    prepared = prepare(args.panel, master, args.taxonomy, args.samples_csv)
    values = prepared["values"]
    print("[panel] %d cell lines, order %s"
          % (len(prepared["order"]), prepared["order_provenance"]))
    print("[panel] %s %% median %.2f, range [%.2f, %.2f]"
          % (spec["stat"], np.nanmedian(values), np.nanmin(values), np.nanmax(values)))

    figure, _axis = draw_side_panel(
        values, prepared["labels"], spec["colour"], spec["title"], spec["xlabel"],
        tuple(args.figsize), args.show_labels)
    # tight=False: a tight crop would undo the shared-box row alignment.
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats),
                      args.force, tight=False)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
