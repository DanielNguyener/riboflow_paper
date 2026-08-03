#!/usr/bin/env python3
"""Figure 5 D -- alternative-isoform exon share of genome-only unique reads."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REQUIRED = ("sample", "n_nonselected_isoform_exon", "n_gUtA")

def prepare(reach_master, taxonomy=None, samples_csv=None):
    sys.path.insert(0, str(HERE))
    import fig05_common as common
    import panel_style as ps

    frame = pd.read_csv(reach_master, sep="\t")
    ps.require_columns(frame, REQUIRED, str(reach_master))
    frame = frame.copy()
    frame["nonselected_pct"] = 100.0 * frame["n_nonselected_isoform_exon"] / frame["n_gUtA"]

    if taxonomy:
        order = common.sample_order(common.load_taxonomy(taxonomy))
        provenance = "shared Figure-5 order, derived from %s" % taxonomy
    else:
        order = frame.sort_values("nonselected_pct")["sample"].tolist()
        provenance = "this panel alone, by ascending nonselected_pct"
    lookup = frame.set_index("sample")["nonselected_pct"]
    values = np.array([lookup.get(s, np.nan) for s in order], dtype=float)
    labels_map = common.load_labels(samples_csv)
    return {"order": order, "values": values, "order_provenance": provenance,
            "labels": [labels_map.get(s, s) for s in order], "source": str(reach_master)}

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reach-master", required=True, type=Path)
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

    prepared = prepare(args.reach_master, args.taxonomy, args.samples_csv)
    values = prepared["values"]
    print("[panel] %d cell lines, order %s"
          % (len(prepared["order"]), prepared["order_provenance"]))
    print("[panel] alternative-isoform exon %% median %.2f, range [%.2f, %.2f]"
          % (np.nanmedian(values), np.nanmin(values), np.nanmax(values)))

    figure, _axis = draw_side_panel(
        values, prepared["labels"], "#7fb9da",
        "alternative\nisoform exon", "% of genome-only\nunique reads",
        tuple(args.figsize), args.show_labels)
    written = ps.save(figure, args.output, ps.resolve_formats(args.formats), args.force)
    for path in written:
        print("[panel] wrote %s" % path)
    return 0

if __name__ == "__main__":
    sys.exit(main())
