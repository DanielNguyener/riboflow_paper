#!/usr/bin/env python3
"""Compose the published figures from panel assets, per the manifest's `figures:` block.

Held to the PLOS spec; `--check` verifies outputs. Figure 1 is the author's schematic.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"
STAGING = REPO / "results" / "panels" / "_plos_fitted"
OUT_DIR = REPO / "figures" / "published"
SAMPLES_CSV = REPO / "supporting_information" / "S1_Table" / "samples.csv"

PAGE_W_IN, PAGE_H_IN = 7.5, 8.75          # PLOS maxima
MIN_W_IN = 2.63
FONT_MIN_PT, FONT_MAX_PT = 8.0, 12.0
MAX_BYTES = 10 * 1024 ** 2
ALLOWED_FONTS = ("Arial", "Times", "Symbol")

sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / "panels"))


def log(message):
    print("[assemble] %s" % message, flush=True)


def run(cmd):
    completed = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError("panel generation failed:\n%s\n%s\n%s"
                           % (" ".join(str(c) for c in cmd), completed.stdout,
                              completed.stderr))


def load_manifest(path):
    import make_panels
    document, panels = make_panels.load_manifest(path)
    return document, {p["id"]: p for p in panels}


def generator_command(entry, defaults, output, extras=(), formats=("pdf",), figsize=None):
    """The manifest's command for `entry`, to `output`, plus composer-specific flags."""
    import make_panels
    return make_panels.build_command(entry, defaults, formats, True, output, figsize) \
        + [str(x) for x in extras]


def bold_font():
    """(fontname, fontfile) for the panel letters: the bold cut of the panels' own family."""
    import panel_style as ps
    from matplotlib import font_manager
    name, _ = ps.resolve_font()
    path = font_manager.findfont(font_manager.FontProperties(family=name, weight="bold"),
                                 fallback_to_default=False)
    return name.replace(" ", "") + "-Bold", path


def measure(pdf_path):
    import fitz
    with fitz.open(pdf_path) as doc:
        rect = doc[0].rect
    return rect.width, rect.height


def resolve_gsm(cell_line):
    """The GSM the manuscript's tables use in place of a cell-line name."""
    import pandas as pd
    samples = pd.read_csv(SAMPLES_CSV)
    row = samples[samples["cell_line"].astype(str).str.replace(" ", "_") == cell_line]
    if row.empty:
        raise SystemExit("cell line %r not found in %s" % (cell_line, SAMPLES_CSV))
    return str(row["ribo_GSM"].iloc[0])


# ── Figure 2: two 24 x 9 grids side by side ─────────────────────────────────────────────
#: Inches. Both grids share AXES so cells match; B drops the GSM labels and sits flush
#: right of A, so B's left margin is the gutter. B's right margin: colourbar; A's bottom: legend.
FIG2_AXES = (2.72, 6.85)
FIG2_MARGINS_A = (0.88, 1.50, 0.05, 0.25)     # left, bottom, right, top
FIG2_MARGINS_B = (0.20, 1.50, 0.87, 0.25)


def compose_fig02_stack(spec, by_id, defaults, pdf_out):
    import fitz
    STAGING.mkdir(parents=True, exist_ok=True)
    border = float(spec["raster"]["border_pt"])
    label_pt = float(spec["label_pt"])
    a_id, b_id = spec["panels"]
    a, b = STAGING / (a_id + "_plos"), STAGING / (b_id + "_plos")
    run(generator_command(by_id[a_id], defaults, a,
                          ["--type-scale", "base", "--legend-ncol", "1",
                           "--axes-size", *FIG2_AXES, "--margins", *FIG2_MARGINS_A]))
    run(generator_command(by_id[b_id], defaults, b,
                          ["--type-scale", "base", "--hide-ylabels",
                           "--axes-size", *FIG2_AXES, "--margins", *FIG2_MARGINS_B]))
    panels = [("A", a.with_suffix(".pdf"), FIG2_MARGINS_A[0]),
              ("B", b.with_suffix(".pdf"), FIG2_MARGINS_B[0])]
    fontname, fontfile = bold_font()
    docs = [fitz.open(p) for _, p, _ in panels]
    widths = [d[0].rect.width for d in docs]
    heights = [d[0].rect.height for d in docs]
    page_w = sum(widths) + 2 * border
    page_h = max(heights) + 2 * border
    if page_w > PAGE_W_IN * 72 + 0.01 or page_h > PAGE_H_IN * 72 + 0.01:
        raise SystemExit("Figure 2 composes to %.2f x %.2f in, over the PLOS maximum"
                         % (page_w / 72, page_h / 72))
    out = fitz.open()
    page = out.new_page(width=page_w, height=page_h)
    x = border
    for (letter, _, left_margin_in), doc, w, h in zip(panels, docs, widths, heights):
        page.show_pdf_page(fitz.Rect(x, border, x + w, border + h), doc, 0)
        # Letter in the panel's top margin: over A's left margin, over B's grid itself.
        lx = x + (2.0 if letter == "A" else left_margin_in * 72)
        page.insert_text((lx, border + label_pt + 2), letter, fontsize=label_pt,
                         fontname=fontname, fontfile=fontfile, color=(0, 0, 0))
        x += w
    out.save(pdf_out)
    for d in docs:
        d.close()


# ── Figure 3: A / B / C|D rows, each panel fitted to its slot ────────────────────────────
FIG3_LABEL_GUTTER = 14.0             # a strip down the left edge that holds the letters
FIG3_GAP_V = 5.0
FIG3_GAP_H = 6.0
FIT_TOL = 0.4
FIT_MAX_ITER = 6
FIG3_ROWS = [["A"], ["B"], ["C", "D"]]
#: C renders STACKED so its 24 GSM labels appear once and fit beside D.
#: Row heights in inches for C and D; A and B share what is left.
FIG3_HEIGHT_CD_IN = 3.30
FIG3_WIDTH_D_IN = 2.45
HIGHLIGHT_CELL_LINE = "HeLa"
GAPDH_NAME, GAPDH_TRANSCRIPT = "GAPDH", "ENST00000396861.5"
COMT_NAME, COMT_TRANSCRIPT = "COMT", "ENST00000361682.11"


def fig03_generators(spec, by_id, defaults, gsm):
    """letter -> (entry, extra flags, starting figsize, output stem)."""
    ids = dict(zip("ABCD", spec["panels"]))
    return {
        "A": (by_id[ids["A"]],
              ["--title", "%s (%s) - %s" % (COMT_NAME, COMT_TRANSCRIPT, gsm),
               "--labels", "minimal", "--title-correlations",
               "--record-json", STAGING / "fig03A_record.json"],
              (10.0, 4.8), STAGING / "fig03A_fitted"),
        "B": (by_id[ids["B"]],
              ["--title", "%s (%s) - %s" % (GAPDH_NAME, GAPDH_TRANSCRIPT, gsm),
               "--labels", "minimal", "--title-correlations",
               "--record-json", STAGING / "fig03B_record.json"],
              (10.0, 4.8), STAGING / "fig03B_fitted"),
        "C": (by_id[ids["C"]], ["--layout", "stacked"], (5.0, 3.3), STAGING / "fig03C_fitted"),
        # No --highlight-sample: panel D marks no cell line.
        "D": (by_id[ids["D"]], ["--short-labels", "--no-points"], (2.45, 3.3),
              STAGING / "fig03D_fitted"),
    }


def fig03_targets(page_w, page_h, margin):
    """Panel boxes (in pt) that tile exactly into the target page."""
    content_width = page_w - 2 * margin - FIG3_LABEL_GUTTER
    content_height = page_h - 2 * margin
    available_height = content_height - (len(FIG3_ROWS) - 1) * FIG3_GAP_V
    h_cd = FIG3_HEIGHT_CD_IN * 72.0
    h_ab = (available_height - h_cd) / 2.0
    w_d = FIG3_WIDTH_D_IN * 72.0
    return {"A": (content_width, h_ab), "B": (content_width, h_ab),
            "C": (content_width - FIG3_GAP_H - w_d, h_cd), "D": (w_d, h_cd)}


def fig03_fit(generators, defaults, name, target_w, target_h):
    entry, extras, initial_figsize, out_stem = generators[name]
    nominal_w_pt, nominal_h_pt = initial_figsize[0] * 72.0, initial_figsize[1] * 72.0
    pdf_path = actual_w = actual_h = None
    for iteration in range(1, FIT_MAX_ITER + 1):
        figsize = (nominal_w_pt / 72.0, nominal_h_pt / 72.0)
        run(generator_command(entry, defaults, out_stem, extras, figsize=figsize))
        pdf_path = out_stem.with_suffix(".pdf")
        actual_w, actual_h = measure(pdf_path)
        err_w, err_h = target_w - actual_w, target_h - actual_h
        print("  %s iter %d: figsize=(%.3f, %.3f)in -> %.2f x %.2f pt "
              "(target %.2f x %.2f pt, err %.2f, %.2f)"
              % (name, iteration, figsize[0], figsize[1], actual_w, actual_h,
                 target_w, target_h, err_w, err_h))
        if abs(err_w) <= FIT_TOL and abs(err_h) <= FIT_TOL:
            break
        nominal_w_pt += err_w
        nominal_h_pt += err_h
    return pdf_path, actual_w, actual_h


def write_fig03_annotations(path):
    """The numbers panels A and B no longer print on themselves, from the render records."""
    lines = ["# Figure 3, panels A and B: annotations removed from the image", "",
             "Written by `code/assemble_figures.py` from the panel generator's render record. "
             "Each panel is drawn with `--labels minimal`: the route names, the correlation "
             "box and the trim-boundary captions are left off the axes.", "",
             "- Upper track: P-site coverage, mirrored -- **genome** route drawn upward "
             "(green), **transcriptome** route drawn downward (red).",
             "- Lower track: footprint coverage, both routes overlaid (same colours).",
             "- Dashed vertical lines mark the CDS trim boundaries.", ""]
    for letter, stem in (("A", "fig03A"), ("B", "fig03B")):
        record = json.loads((STAGING / ("%s_record.json" % stem)).read_text())
        resolved, corr = record["resolved"], record["correlations"]
        x0, x1 = record["axis_window"]
        lines += ["## Panel %s -- %s (%s), %s" % (letter, resolved["gene_name"],
                                                   resolved["transcript_id"],
                                                   record["sample"]), "",
                  "- CDS window: start +%d nt to stop −%d nt (CDS positions %d–%d)"
                  % (record["trim"], record["trim"], x0, x1), "",
                  "| track | Spearman ρ | Pearson r |", "|---|---|---|"]
        for key, name in (("psite", "P-site coverage"), ("footprint", "footprint coverage")):
            lines.append("| %s | %.3f | %.3f |" % (name, corr[key]["spearman"],
                                                   corr[key]["pearson"]))
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def compose_fig03_fit(spec, by_id, defaults, pdf_out):
    import fitz
    STAGING.mkdir(parents=True, exist_ok=True)
    margin = float(spec["raster"]["border_pt"])
    label_pt = float(spec["label_pt"])
    page_w, page_h = (float(v) * 72.0 for v in spec["page_in"])
    targets = fig03_targets(page_w, page_h, margin)
    gsm = resolve_gsm(HIGHLIGHT_CELL_LINE)
    generators = fig03_generators(spec, by_id, defaults, gsm)
    fitted = {}
    for name in "ABCD":
        print("fitting panel %s to %.2f x %.2f pt" % (name, *targets[name]))
        pdf_path, w, h = fig03_fit(generators, defaults, name, *targets[name])
        fitted[name] = (fitz.open(pdf_path), w, h)

    content_width = page_w - 2 * margin - FIG3_LABEL_GUTTER
    out = fitz.open()
    page = out.new_page(width=page_w, height=page_h)
    fontname, fontfile = bold_font()
    y = margin
    for row_index, row in enumerate(FIG3_ROWS):
        widths = [fitted[name][1] for name in row]
        heights = [fitted[name][2] for name in row]
        row_width = sum(widths) + FIG3_GAP_H * (len(row) - 1)
        row_height = max(heights)
        x = margin + FIG3_LABEL_GUTTER + (content_width - row_width) / 2.0
        for name in row:
            doc, w, h = fitted[name]
            y_off = y + (row_height - h) / 2.0
            rect = fitz.Rect(x, y_off, x + w, y_off + h)
            page.show_pdf_page(rect, doc, 0)
            # Panel letter level with the panel's top edge: in the page's left gutter for
            # the first panel of a row, in the panel's own top-left corner otherwise.
            letter_x = margin + 2.0 if name == row[0] else rect.x0 + 2.0
            page.insert_text((letter_x, y + label_pt), name, fontsize=label_pt,
                             fontname=fontname, fontfile=fontfile, color=(0, 0, 0))
            x += w + FIG3_GAP_H
        y += row_height
        if row_index < len(FIG3_ROWS) - 1:
            y += FIG3_GAP_V
    slack = page_h - margin - y
    if abs(slack) > 2.0:
        print("warning: %.2f pt of unused vertical space before the bottom margin" % slack)
    out.save(pdf_out)
    out.close()
    for doc, _, _ in fitted.values():
        doc.close()
    if spec.get("annotations"):
        log("wrote %s" % write_fig03_annotations(REPO / spec["annotations"]))


# ── Figure 4: one generator page, the generator writes the TIFF ──────────────────────────
def compose_single_panel(spec, by_id, defaults, pdf_out, tif_out):
    STAGING.mkdir(parents=True, exist_ok=True)
    (panel_id,) = spec["panels"]
    stem = STAGING / panel_id
    run(generator_command(by_id[panel_id], defaults, stem, formats=("pdf", "tif")))
    shutil.copyfile(stem.with_suffix(".pdf"), pdf_out)
    shutil.copyfile(stem.with_suffix(".tif"), tif_out)


# ── Figures 5, 6: panels placed 1:1 by ink box ───────────────────────────────────────────
def panel_entry(entry):
    import figure_io
    pdf = REPO / (entry["output"] + ".pdf")
    if not pdf.exists():
        raise SystemExit("%s is missing; run `python code/make_panels.py %s` first"
                         % (pdf.relative_to(REPO), entry["id"]))
    clip = figure_io.ink_box_trim_bottom(str(pdf))
    clip_file = REPO / (entry["output"] + ".clip.json")
    if clip_file.exists():
        # keep the fitted horizontal clip (the trim helper re-measures; widths must agree)
        fitted = json.loads(clip_file.read_text())
        clip.x0, clip.x1 = fitted["x0"], fitted["x1"]
    return {"pdf": str(pdf), "clip": clip, "stem": entry["id"]}


def compose_rows_1to1(spec, by_id, defaults, pdf_out):
    import figure_io
    page = spec["page"]
    rows = [[panel_entry(by_id[stem]) for stem in row] for row in spec["rows"]]
    for row in rows:
        row_w = sum(p["clip"].width for p in row) + page["gutter_pt"] * (len(row) - 1)
        if row_w + 2 * page["margin_pt"] > page["max_width_pt"]:
            raise SystemExit("a row is %.1f pt of ink; reduce panel widths" % row_w)
    page_w, page_h = figure_io.compose(
        rows, spec["letters"], str(pdf_out), margin_pt=page["margin_pt"],
        gutter_pt=page["gutter_pt"], row_gap_pt=page["row_gap_pt"],
        letter_pt=page["letter_pt"])
    figure_io.check_page(pdf_out.name, page_w, page_h)


# ── PDF -> TIFF and the spec check ──────────────────────────────────────────────────────
def rasterise_matplotlib_pdf(pdf_path, tiff_path, dpi):
    """One RGB raster of the whole page (border included), written as the TIFF."""
    import fitz
    from PIL import Image
    with fitz.open(pdf_path) as doc:
        pix = doc[0].get_pixmap(dpi=dpi, alpha=False, colorspace=fitz.csRGB)
    image = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    image.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))


def write_png(tiff_path, png_path):
    """A PNG of the TIFF's pixels, for previews and co-authors without a TIFF viewer."""
    from PIL import Image
    with Image.open(tiff_path) as image:
        dpi = image.info.get("dpi", (300, 300))
        image.save(png_path, format="PNG", dpi=dpi, optimize=True)


def check(pdf_path, tiff_path, raster):
    """Re-open the outputs and verify every measurable clause of the spec."""
    import fitz
    from PIL import Image
    problems = []
    min_pt = float(raster.get("min_font_pt", FONT_MIN_PT))
    with fitz.open(pdf_path) as doc:
        page = doc[0]
        w_in, h_in = page.rect.width / 72, page.rect.height / 72
        fonts = set()
        for block in page.get_text("dict")["blocks"]:
            if block["type"] != 0:
                continue
            for line in block["lines"]:
                for span in line["spans"]:
                    fonts.add((span["font"], round(span["size"], 2)))
    if not MIN_W_IN <= w_in <= PAGE_W_IN + 1e-3:
        problems.append("width %.2f in outside [%.2f, %.2f]" % (w_in, MIN_W_IN, PAGE_W_IN))
    if h_in > PAGE_H_IN + 1e-3:
        problems.append("height %.2f in over %.2f" % (h_in, PAGE_H_IN))
    for font, size in sorted(fonts):
        if not any(font.startswith(a) for a in ALLOWED_FONTS):
            problems.append("font %s is not Arial/Times/Symbol" % font)
        if not min_pt - 0.05 <= size <= FONT_MAX_PT + 0.05:
            problems.append("%s at %.1f pt outside %g-12 pt" % (font, size, min_pt))
    with Image.open(tiff_path) as image:
        if image.mode != "RGB":
            problems.append("TIFF mode %s, not RGB" % image.mode)
        if image.info.get("compression") != "tiff_lzw":
            problems.append("TIFF compression %s, not LZW" % image.info.get("compression"))
        px_w, px_h = image.size
        x_dpi = round(image.info.get("dpi", (0, 0))[0])
        if not 300 <= x_dpi <= 600:
            problems.append("TIFF dpi %s outside 300-600" % x_dpi)
        if px_w > round(PAGE_W_IN * x_dpi) or px_h > round(PAGE_H_IN * x_dpi):
            problems.append("TIFF %dx%d px over %dx%d at %d dpi"
                            % (px_w, px_h, round(PAGE_W_IN * x_dpi),
                               round(PAGE_H_IN * x_dpi), x_dpi))
        if px_w < round(MIN_W_IN * x_dpi):
            problems.append("TIFF %d px narrower than %.2f in" % (px_w, MIN_W_IN))
        try:
            image.seek(1)
            problems.append("TIFF has more than one page")
        except EOFError:
            pass
    size = Path(tiff_path).stat().st_size
    if size > MAX_BYTES:
        problems.append("TIFF is %.1f MB, over 10 MB" % (size / 1024 ** 2))
    log("%s: %.2f x %.2f in, %dx%d px @ %d dpi, %.2f MB, fonts %s"
        % (Path(tiff_path).name, w_in, h_in, px_w, px_h, x_dpi, size / 1024 ** 2,
           ", ".join("%s %.0f" % f for f in sorted(fonts))))
    return problems


COMPOSERS = ("fig02_stack", "fig03_fit", "single_panel", "rows_1to1")
RASTERS = ("matplotlib_pdf", "generator_tiff", "fitz")


def assemble(number, spec, by_id, defaults, output_dir, png=False):
    import figure_io
    composer, raster = spec["composer"], spec["raster"]
    if composer not in COMPOSERS or raster["kind"] not in RASTERS:
        raise SystemExit("figure %s: unknown composer/raster %s/%s"
                         % (number, composer, raster["kind"]))
    stem = Path(spec["output"]).name
    pdf = output_dir / ("%s_plos.pdf" % stem)
    tif = output_dir / ("%s.tif" % stem)
    log("Figure %s (%s) -> %s" % (number, composer, tif))
    if composer == "fig02_stack":
        compose_fig02_stack(spec, by_id, defaults, pdf)
    elif composer == "fig03_fit":
        compose_fig03_fit(spec, by_id, defaults, pdf)
    elif composer == "single_panel":
        compose_single_panel(spec, by_id, defaults, pdf, tif)
    else:
        compose_rows_1to1(spec, by_id, defaults, pdf)
    if raster["kind"] == "matplotlib_pdf":
        rasterise_matplotlib_pdf(pdf, tif, int(raster["dpi"]))
    elif raster["kind"] == "fitz":
        figure_io.export_tiff(str(pdf), str(tif), float(raster["dpi"]))
    if png:
        write_png(tif, output_dir / ("%s.png" % stem))
    return pdf, tif


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--figure", type=int, action="append", help="2, 3, 4, 5 or 6")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--manifest", type=Path, default=REPO / "config" / "panel_manifest.yaml")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--check", action="store_true", help="verify the outputs against the spec")
    parser.add_argument("--png", action="store_true", help="also write Fig<N>.png, the same pixels")
    args = parser.parse_args(argv)
    document, by_id = load_manifest(args.manifest)
    figures_block = document["figures"]
    if not args.all and not args.figure:
        parser.error("name a figure (--figure 5) or pass --all")
    unknown = [n for n in (args.figure or []) if n not in figures_block]
    if unknown:
        parser.error("no `figures:` entry for %s (have %s)"
                     % (unknown, sorted(figures_block)))
    figures = sorted(figures_block) if args.all else sorted(set(args.figure))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    failed = False
    for number in figures:
        spec = figures_block[number]
        pdf, tif = assemble(number, spec, by_id, document["defaults"], args.output_dir,
                            args.png)
        if args.check:
            for problem in check(pdf, tif, spec["raster"]):
                failed = True
                log("  FAIL Fig%s: %s" % (number, problem))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
