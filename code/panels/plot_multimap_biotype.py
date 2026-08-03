#!/usr/bin/env python3
"""Figure 5 C -- protein-coding to pseudogene tie share of genome multimapper reads."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REQUIRED = ("sample", "pct_cross_pp_pc", "pct_cross_pc_pp")

def prepare(tie_master, taxonomy=None, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import fig05_common as common
    import panel_style as ps

    frame = pd.read_csv(tie_master, sep="\t")
    ps.require_columns(frame, REQUIRED, str(tie_master))
    frame = frame.copy()
    frame["cross_tie_pct"] = frame["pct_cross_pp_pc"] + frame["pct_cross_pc_pp"]

    if taxonomy:
        order = common.sample_order(common.load_taxonomy(taxonomy))
        provenance = "shared Figure-5 order, derived from %s" % taxonomy
    else:
        order = frame.sort_values("cross_tie_pct")["sample"].tolist()
        provenance = "this panel alone, by ascending cross_tie_pct"
    lookup = frame.set_index("sample")["cross_tie_pct"]
    values = np.array([lookup.get(s, np.nan) for s in order], dtype=float)
    labels_map = common.load_labels(samples_csv)
    return {"order": order, "values": values, "order_provenance": provenance,
            "labels": [labels_map.get(s, s) for s in order], "source": str(tie_master)}

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tie-master", required=True, type=Path)
    parser.add_argument("--taxonomy", type=Path,
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
    from _fig05_side_panel import draw_side_panel

    prepared = prepare(args.tie_master, args.taxonomy, args.samples_csv)
    values = prepared["values"]
    print("[panel] %d cell lines, order %s"
          % (len(prepared["order"]), prepared["order_provenance"]))
    print("[panel] cross tie %% median %.2f, range [%.2f, %.2f]"
          % (np.nanmedian(values), np.nanmin(values), np.nanmax(values)))

    figure, _axis = draw_side_panel(
        values, prepared["labels"], "#1a7d1a",
        "protein-coding ↔\npseudogene tie", "% of multimapper\nreads (primary)",
        tuple(args.figsize), args.show_labels)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
