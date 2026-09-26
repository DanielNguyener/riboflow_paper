#!/usr/bin/env python3
"""Regenerate the manuscript's panel assets from `config/panel_manifest.yaml`.

`fit:` panels re-render until their ink is `fit.width_pt` wide and record `<stem>.clip.json`.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "config" / "panel_manifest.yaml"
REFERENCES = REPO / "figures" / "panel_references"

sys.path.insert(0, str(REPO / "code" / "common"))
import inputs  # noqa: E402

PDF_CREATIONDATE = re.compile(rb"/CreationDate \(D:[0-9+'\-Z]+\)")

log = inputs.make_log("make_panels")

def load_manifest(path=MANIFEST):
    import yaml
    document = yaml.safe_load(path.read_text())
    expected = "riboflow_paper/panel-manifest/2"
    if document.get("schema_version") != expected:
        raise SystemExit("%s has schema_version %r, expected %r"
                         % (path, document.get("schema_version"), expected))
    # Figure keys are labels, not integers: main-text figures are digits, supporting
    # figures are `S<n>`. YAML reads a bare digit as int, so normalise once here.
    document["figures"] = {str(k): v for k, v in document["figures"].items()}
    for panel in document["panels"]:
        if panel.get("figure") is not None:
            panel["figure"] = str(panel["figure"])
    panels = [p for p in document["panels"] if p.get("generator")]
    outputs = [p["output"] for p in panels]
    duplicates = {o for o in outputs if outputs.count(o) > 1}
    if duplicates:
        raise SystemExit("the manifest declares the same output more than once: %s"
                         % ", ".join(sorted(duplicates)))
    return document, panels

def is_supporting(number):
    """`S1`, `S2`, ... are supporting-information figures; digits are main-text figures."""
    return str(number)[:1].upper() == "S"


def figure_order(number):
    """Sort key: main-text figures by number, then supporting figures by number."""
    number = str(number)
    return (1, int(number[1:])) if is_supporting(number) else (0, int(number))


def figure_stem(number):
    """The published file stem PLOS expects: `Fig2` for the main text, `S1_Fig` for SI."""
    number = str(number)
    return "%s_Fig" % number.upper() if is_supporting(number) else "Fig%s" % number


def figure_label(number):
    """How the figure is cited: `Figure 2`, `S1 Fig`."""
    number = str(number)
    return "%s Fig" % number.upper() if is_supporting(number) else "Figure %s" % number


def compare_to_reference(generated, reference):
    """Byte comparison, tolerating ONLY the PDF /CreationDate field."""
    if not reference.exists():
        return "NO_REFERENCE", "no stored reference at %s" % reference.relative_to(REPO)
    a, b = generated.read_bytes(), reference.read_bytes()
    if a == b:
        return "IDENTICAL", "bit-for-bit identical (%d bytes)" % len(a)
    if len(a) != len(b):
        return "DIFFERS", "size %d vs reference %d" % (len(a), len(b))
    offsets = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    span = PDF_CREATIONDATE.search(b)
    if span and all(span.start() <= i < span.end() for i in offsets):
        if PDF_CREATIONDATE.sub(b"", a) == PDF_CREATIONDATE.sub(b"", b):
            return ("IDENTICAL_MOD_TIMESTAMP",
                    "%d byte(s) differ, all inside /CreationDate" % len(offsets))
    return "DIFFERS", ("%d byte(s) differ, NOT confined to /CreationDate (first at %d)"
                       % (len(offsets), offsets[0]))

def _flag_pairs(args_block):
    """Manifest `args:` -> command-line flags. Lists become comma-separated values."""
    flags = []
    for key, value in (args_block or {}).items():
        flag = "--" + key.replace("_", "-")
        if isinstance(value, bool):
            if value:
                flags.append(flag)
        elif isinstance(value, (list, tuple)):
            if key == "ylim":
                flags += [flag] + [str(v) for v in value]
            else:
                flags += [flag, ",".join(str(v) for v in value)]
        else:
            flags += [flag, str(value)]
    return flags

def build_command(entry, formats, force, output=None, figsize=None):
    command = [sys.executable, str(REPO / entry["generator"])]
    for key, value in (entry.get("inputs") or {}).items():
        command += ["--" + key.replace("_", "-"), str(REPO / value)
                    if not str(value).startswith("/") else str(value)]
    command += _flag_pairs(entry.get("args"))
    if figsize is not None:
        command += ["--figsize", "%.6f" % figsize[0], "%.6f" % figsize[1]]
    command += ["--output", str(output if output is not None else REPO / entry["output"])]
    command += ["--format", ",".join(formats)]
    if force:
        command.append("--force")
    return command


def fit_and_run(entry, formats, force):
    """Render a `fit:` panel until its ink is the declared width; record the clip."""
    import json
    sys.path.insert(0, str(REPO / "code" / "panels"))
    import figure_io

    fit = entry["fit"]
    stem = REPO / entry["output"]
    pdf = stem.with_suffix(".pdf")
    if pdf.exists() and not force:
        return 1, "%s exists; pass --force (fitted panels are re-rendered several times)" % pdf
    def build(width_in, height_in):
        # --force on every iteration: the loop overwrites its own previous render.
        return build_command(entry, formats, True, stem, (width_in, height_in))

    clip = figure_io.fit_panel(entry["id"], build, str(pdf), fit["width_pt"],
                               fit["height_in"], start_w_in=fit["width_pt"] / 72.0)
    with open(str(stem) + ".clip.json", "w") as handle:
        json.dump({"x0": clip.x0, "y0": clip.y0, "x1": clip.x1, "y1": clip.y1,
                   "target_w_pt": fit["width_pt"]}, handle, indent=2)
    return 0, ""

def run(command):
    completed = subprocess.run(command, capture_output=True, text=True)
    return completed.returncode, (completed.stdout or "") + (completed.stderr or "")

def _failure_excerpt(output):
    """The first 6 lines of a failure AND the last 18, indented (summary head, remedy tail)."""
    lines = output.splitlines()
    if len(lines) <= 24:
        shown = lines
    else:
        shown = lines[:6] + ["    ... %d lines omitted ..." % (len(lines) - 24)] + lines[-18:]
    return "\n".join("    " + line for line in shown)

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("panels", nargs="*", help="panel ids, e.g. fig02A fig05B")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--format", dest="formats", default=None,
                        help="override the manifest's formats, e.g. pdf,svg,png")
    parser.add_argument("--verify", action="store_true",
                        help="compare each panel against its stored reference")
    parser.add_argument("--accept", action="store_true",
                        help="store the current renders as the references (explicit act)")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    document, panels = load_manifest(args.manifest)
    by_id = {p["id"]: p for p in panels}

    if args.all:
        wanted = [p["id"] for p in panels]
    elif args.panels:
        unknown = [p for p in args.panels if p not in by_id]
        if unknown:
            raise SystemExit("unknown panel(s): %s (have %s)"
                             % (", ".join(unknown), ", ".join(by_id)))
        wanted = args.panels
    else:
        raise SystemExit("nothing selected. Pass panel ids, or --all.")

    if sys.version_info[:2] != (3, 9):
        log("WARNING: running Python %d.%d; the published artifacts were produced on 3.9"
            % sys.version_info[:2])
    os.environ.setdefault("MPLBACKEND", "Agg")

    formats = tuple(args.formats.split(",")) if args.formats else \
        tuple(document["defaults"]["formats"])

    results = []
    for panel_id in wanted:
        entry = by_id[panel_id]
        log("%s (%s %s)" % (panel_id, figure_label(entry["figure"]), entry.get("panel", "")))
        if entry.get("fit"):
            code, output = fit_and_run(entry, formats, args.force)
        else:
            code, output = run(build_command(entry, formats, args.force))
        if code:
            log("  FAILED\n%s" % _failure_excerpt(output))
            results.append((panel_id, False, None, None))
            continue
        primary = REPO / (entry["output"] + "." + formats[0])
        verdict = detail = None
        if args.accept:
            REFERENCES.mkdir(parents=True, exist_ok=True)
            destination = REFERENCES / primary.name
            destination.write_bytes(primary.read_bytes())
            log("  accepted as reference: %s" % destination.relative_to(REPO))
        if args.verify:
            verdict, detail = compare_to_reference(primary, REFERENCES / primary.name)
            log("  vs reference: %s -- %s" % (verdict, detail))
        results.append((panel_id, True, verdict, detail))

    print()
    failed = [p for p, ok, _v, _d in results if not ok]
    if args.verify:
        failed += [p for p, ok, v, _d in results
                   if ok and v not in (None, "IDENTICAL", "IDENTICAL_MOD_TIMESTAMP")]
    for panel_id, ok, verdict, detail in results:
        line = "%-8s %s" % (panel_id, "OK  " if ok else "FAIL")
        if verdict:
            line += "  %s" % verdict
        print(line)
    print("\n%d/%d panels produced." % (len(results) - len(failed), len(results)))
    print("These are PANEL ASSETS. `python code/assemble_figures.py --all --check` composes")
    print("them into figures/published/<Fig<N>|S<N>_Fig>.tif.")
    return 1 if failed else 0

if __name__ == "__main__":
    sys.exit(main())
