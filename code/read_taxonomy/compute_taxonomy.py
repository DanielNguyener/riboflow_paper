#!/usr/bin/env python3
"""Multi-mapping read taxonomy for ONE cell line."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import taxonomy_lib as tl
fc = tl.fc

OUT = fc.output_root() / "read_taxonomy" / "taxonomy"
STAGING = OUT / "_staging"

_ABBR = {"unique": "U", "multi": "M", "absent": "A"}
CELLS = [(g, t) for g in tl.STATES for t in tl.STATES if not (g == "absent" and t == "absent")]

def build_row(sample, counts, n_universe):
    row = {"sample": sample, "n_universe": n_universe}
    for g, t in CELLS:
        n = counts[(g, t)]
        key = f"g{_ABBR[g]}_t{_ABBR[t]}"
        row[f"n_{key}"] = n
        row[f"pct_{key}"] = 100.0 * n / n_universe if n_universe else float("nan")
    for g in tl.STATES:
        row[f"n_genome_{g}"] = sum(counts[(g, t)] for t in tl.STATES if (g, t) in counts)
    for t in tl.STATES:
        row[f"n_txome_{t}"] = sum(counts[(g, t)] for g in tl.STATES if (g, t) in counts)
    return row

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    args = ap.parse_args()

    counts, n_universe = tl.classify_sample(args.sample)
    row = build_row(args.sample, counts, n_universe)

    STAGING.mkdir(parents=True, exist_ok=True)
    out = STAGING / f"{args.sample}.tsv"
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False)

    core = {k: row[f"pct_g{_ABBR[g]}_t{_ABBR[t]}"]
            for k, (g, t) in {"gU_tU": ("unique", "unique"), "gU_tM": ("unique", "multi"),
                              "gM_tU": ("multi", "unique"), "gM_tM": ("multi", "multi")}.items()}
    print(f"[{args.sample}] n_universe={n_universe:,}  "
          f"gU_tU={core['gU_tU']:.2f}%  gU_tM={core['gU_tM']:.2f}%  "
          f"gM_tU={core['gM_tU']:.2f}%  gM_tM={core['gM_tM']:.2f}%", flush=True)
    print(f"wrote {out}", flush=True)

if __name__ == "__main__":
    main()
