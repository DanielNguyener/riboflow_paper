"""Ink cropping, the fit loop, 1:1 PDF composition and the fitz TIFF export (Figures 5, 6).

The thresholds, dpi, tolerances and geometry are load-bearing: they reproduce the shipped
figures pixel for pixel. `assemble_figures.py` and `make_panels.py` import from here.
"""
from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def die(message):
    raise SystemExit("error: %s" % message)

#: PLOS Computational Biology's box: 789-2250 px wide, <= 2625 px tall at 300 dpi.
PAGE_MAX_W_PT = 540.0
PAGE_MAX_H_PT = 630.0
PAGE_MIN_W_PT = 789 / 300.0 * 72.0

MARGIN_PT = 4.0
GUTTER_PT = 6.0
ROW_GAP_PT = 10.0
LETTER_PT = 10.0
LETTER_BAND_PT = LETTER_PT * 1.25

FIT_TOL_PT = 0.4
FIT_MAX_ITER = 6
#: A hit within this is accepted after the loop: ink is measured at 144 dpi, so widths
#: quantise to 0.5 pt and an exact landing is not always reachable.
FIT_SETTLE_PT = 1.0


def ink_box(pdf, dpi=144, pad_pt=1.0):
    """The panel's drawn extent HORIZONTALLY, over its full page height.

    Vertical bounds are deliberately NOT trimmed -- cropping would break the cohort
    panels' declared-from-page-top row alignment.
    """
    import fitz
    import numpy as np
    from PIL import Image

    with fitz.open(pdf) as document:
        page = document[0]
        pixmap = page.get_pixmap(dpi=dpi, alpha=False)
        page_rect = fitz.Rect(page.rect)
    grey = np.array(Image.open(io.BytesIO(pixmap.tobytes("png"))).convert("L"))
    drawn = np.argwhere(grey < 250)
    if not len(drawn):
        return page_rect
    scale = dpi / 72.0
    left = drawn[:, 1].min() / scale
    right = (drawn[:, 1].max() + 1) / scale
    return fitz.Rect(max(page_rect.x0, left - pad_pt), page_rect.y0,
                     min(page_rect.x1, right + pad_pt), page_rect.y1)


def ink_box_trim_bottom(pdf, dpi=144, pad_pt=1.0):
    """`ink_box`, and also trim the BOTTOM -- but never the top (see `ink_box`)."""
    import fitz
    import numpy as np
    from PIL import Image

    box = ink_box(pdf, dpi=dpi, pad_pt=pad_pt)
    with fitz.open(pdf) as document:
        page = document[0]
        pixmap = page.get_pixmap(dpi=dpi, alpha=False)
        page_rect = fitz.Rect(page.rect)
    grey = np.array(Image.open(io.BytesIO(pixmap.tobytes("png"))).convert("L"))
    drawn = np.argwhere(grey < 250)
    if not len(drawn):
        return box
    bottom = (drawn[:, 0].max() + 1) / (dpi / 72.0)
    return fitz.Rect(box.x0, box.y0, box.x1, min(page_rect.y1, bottom + pad_pt))


def render(command, label, cwd=None):
    result = subprocess.run(command, cwd=str(cwd or REPO), capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout + result.stderr)
        die("%s failed (exit %d)" % (label, result.returncode))
    return result


def fit_panel(label, build_command, out_pdf, target_w_pt, height_in, start_w_in):
    """Render until the panel's INK is `target_w_pt` wide.

    `bbox_inches="tight"` makes the page differ from the requested figsize, so iterate.
    """
    nominal_w_pt = start_w_in * 72.0
    best = None
    for iteration in range(1, FIT_MAX_ITER + 1):
        render(build_command(nominal_w_pt / 72.0, height_in), label)
        box = ink_box(out_pdf)
        error = target_w_pt - box.width
        print("[fit] %-7s fit %d: ink %.2f pt, %+0.2f from target"
              % (label, iteration, box.width, error))
        if abs(error) <= FIT_TOL_PT:
            return box
        if best is None or abs(error) < abs(best[1]):
            best = (box, error, nominal_w_pt)
        nominal_w_pt += error
    if best is not None and abs(best[1]) <= FIT_SETTLE_PT:
        print("[fit] %-7s settling at %+0.2f pt (below the %.1f pt device-pixel floor)"
              % (label, best[1], FIT_SETTLE_PT))
        render(build_command(best[2] / 72.0, height_in), label)
        return ink_box(out_pdf)
    die("%s did not reach %.1f pt of ink in %d attempts (best error %+.2f pt). Adjust its "
        "height, or raise FIT_MAX_ITER."
        % (label, target_w_pt, FIT_MAX_ITER, best[1] if best else float("nan")))


def arial_bold():
    """The bold face for panel letters, resolved rather than hard-coded."""
    from matplotlib import font_manager
    from matplotlib.font_manager import FontProperties
    return font_manager.findfont(FontProperties(family="Arial", weight="bold"),
                                 fallback_to_default=False)


def compose(rows, letters, output_pdf, margin_pt=MARGIN_PT, gutter_pt=GUTTER_PT,
            row_gap_pt=ROW_GAP_PT, letter_pt=LETTER_PT):
    """Place the panels 1:1 with bold letters in a reserved band above each row.

    `rows` is a list of rows, each a list of {"pdf", "clip", "stem"}.
    """
    import fitz

    letter_band_pt = letter_pt * 1.25
    row_widths = [sum(p["clip"].width for p in row) + gutter_pt * (len(row) - 1)
                  for row in rows]
    row_heights = [max(p["clip"].height for p in row) for row in rows]
    page_w = max(row_widths) + 2 * margin_pt
    page_h = (sum(row_heights) + letter_band_pt * len(rows)
              + row_gap_pt * (len(rows) - 1) + 2 * margin_pt)

    fontfile = arial_bold()
    out = fitz.open()
    page = out.new_page(width=page_w, height=page_h)
    letter_iter = iter(letters)
    y = margin_pt
    for row, row_w, row_h in zip(rows, row_widths, row_heights):
        x = margin_pt + (max(row_widths) - row_w) / 2.0
        for panel in row:
            clip = panel["clip"]
            target = fitz.Rect(x, y + letter_band_pt,
                               x + clip.width, y + letter_band_pt + clip.height)
            with fitz.open(panel["pdf"]) as source:
                page.show_pdf_page(target, source, 0, clip=clip)
            letter = next(letter_iter, None)
            if letter:
                page.insert_text(fitz.Point(x, y + letter_pt), letter,
                                 fontsize=letter_pt, fontname="Arial-Bold",
                                 fontfile=fontfile)
            print("[compose] %-28s %-8s %7.1f x %7.1f pt at (%.1f, %.1f)"
                  % (os.path.basename(output_pdf), panel["stem"],
                     clip.width, clip.height, x, y + letter_band_pt))
            x += clip.width + gutter_pt
        y += letter_band_pt + row_h + row_gap_pt

    out.save(output_pdf)
    out.close()
    return page_w, page_h


def export_tiff(pdf_path, tiff_path, dpi=600.0):
    """PLOS submission raster: flattened RGB, LZW, no alpha, < 10 MB or retried at 300."""
    import fitz
    from PIL import Image

    for attempt_dpi in (dpi, 300.0):
        with fitz.open(pdf_path) as document:
            pixmap = document[0].get_pixmap(dpi=int(attempt_dpi), alpha=False)
        image = Image.frombytes("RGB", (pixmap.width, pixmap.height), pixmap.samples)
        image.save(tiff_path, format="TIFF", compression="tiff_lzw",
                   dpi=(attempt_dpi, attempt_dpi))
        size_mb = os.path.getsize(tiff_path) / 1e6
        if size_mb < 10.0:
            print("[tiff] wrote %s  (%d x %d px, %.0f dpi, RGB, LZW, %.2f MB)"
                  % (tiff_path, image.width, image.height, attempt_dpi, size_mb))
            return
        print("[tiff] %s is %.1f MB at %.0f dpi; retrying at 300"
              % (tiff_path, size_mb, attempt_dpi))
    die("%s exceeds 10 MB even at 300 dpi" % tiff_path)


def check_page(name, page_w, page_h):
    px_w, px_h = page_w / 72.0 * 300, page_h / 72.0 * 300
    ok = 789 <= px_w <= 2250 and px_h <= 2625
    print("[plos] %s page %.2f x %.2f pt = %.0f x %.0f px at 300 dpi  %s"
          % (name, page_w, page_h, px_w, px_h, "OK" if ok else "OUT OF BOX"))
    if not ok:
        die("%s is outside PLOS's figure box" % name)
