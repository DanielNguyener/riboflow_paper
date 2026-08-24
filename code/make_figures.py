#!/usr/bin/env python3
"""End to end: panels -> the five published figures, optionally from BAMs first.

Figures 3A/3B additionally need results/coverage/HeLa.shared_coverage.h5 (see README).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"

#: The `make_tables.py` stages each figure's inputs come from.
STAGES = {2: ["qc"], 3: ["coverage", "concordance"],
          4: ["te_counts", "te_normalize", "te_stats"],
          5: ["taxonomy", "reach", "multimap_biotype"],
          6: ["gene_partition", "locus"]}


def log(message):
    print("[make_figures] %s" % message, flush=True)


def sh(cmd):
    log("$ " + " ".join(str(c) for c in cmd))
    return subprocess.run([str(c) for c in cmd]).returncode


def panel_ids(spec):
    rows = spec.get("rows") or [spec.get("panels", [])]
    return [p for row in rows for p in row]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--figure", type=int, action="append", help="2, 3, 4, 5 or 6")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--check", action="store_true", help="verify every figure against the spec")
    parser.add_argument("--verify", action="store_true",
                        help="compare every panel with figures/panel_references/")
    parser.add_argument("--bams", help="rebuild the tables from this RiboFlow_v2 tree first")
    parser.add_argument("--gtf")
    parser.add_argument("--appris")
    parser.add_argument("--into-data", action="store_true",
                        help="with --bams: replace data/ with the regenerated tables")
    args = parser.parse_args(argv)
    if not args.all and not args.figure:
        parser.error("name a figure (--figure 5) or pass --all")
    sys.path.insert(0, str(CODE))
    import make_panels
    document, panels = make_panels.load_manifest(REPO / "config" / "panel_manifest.yaml")
    figures = sorted(document["figures"]) if args.all else sorted(set(args.figure))

    if args.bams:
        stages = [s for n in figures for s in STAGES[n]]
        cmd = [sys.executable, CODE / "make_tables.py", "--bams", args.bams,
               "--stages", ",".join(stages)]
        if args.gtf:
            cmd += ["--gtf", args.gtf]
        if args.appris:
            cmd += ["--appris", args.appris]
        if args.into_data:
            cmd.append("--into-data")
        if sh(cmd):
            return 1

    # Figures 2-4 re-render their panels at page size inside the assembler; 5 and 6 place
    # the manifest's panel assets, so those are (re)built here.
    wanted = [p for n in figures if document["figures"][n]["composer"] == "rows_1to1"
              for p in panel_ids(document["figures"][n])]
    if args.verify:
        wanted = [p for n in figures for p in panel_ids(document["figures"][n])]
    if wanted:
        cmd = [sys.executable, CODE / "make_panels.py", *wanted, "--force"]
        if args.verify:
            cmd.append("--verify")
        if sh(cmd):
            return 1

    cmd = [sys.executable, CODE / "assemble_figures.py"] + \
        [x for n in figures for x in ("--figure", n)]
    if args.check:
        cmd.append("--check")
    return sh(cmd)


if __name__ == "__main__":
    sys.exit(main())
