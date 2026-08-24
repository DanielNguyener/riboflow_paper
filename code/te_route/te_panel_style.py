"""Figure style: the PLOS print type scale, the page limits, and saving.

Mirrors panels/panel_style.py's conventions, reimplemented so this folder imports no
outside code. Numbers follow PLOS Computational Biology's figure requirements.
"""
from __future__ import annotations

import io
from pathlib import Path

#: PLOS permits Arial, Times or Symbol only.
FONT_FAMILY = "Arial"

#: The print scale, all within PLOS's 8-12 point window.
FONT_LABEL = 11
FONT_TICK = 10
FONT_ANNOTATION = 10
FONT_PANEL_LETTER = 12

#: Linear scale of strokes against the 11-pt base the line widths were drawn at. Marker `s` is
#: an AREA, so a dot keeps its apparent size only when scaled by the SQUARE of this.
TYPE_SCALE = FONT_LABEL / 11.0
MARKER_AREA = TYPE_SCALE ** 2

#: PLOS never wants a stroke that vanishes at print; every scaled line width is clamped here.
MIN_LINEWIDTH = 0.5

#: Page limits in inches, and the pixel limits they become at SAVE_DPI.
PAGE_WIDTH_MAX = 7.5
PAGE_HEIGHT_MAX = 8.75
SAVE_DPI = 300
MAX_PIXELS = (2250, 2625)
MAX_BYTES = 10 * 1024 * 1024
BORDER_PT = 2.0

SPEARMAN_FILL, SPEARMAN_LINE = "#8fb4d6", "#2c5f8a"
PEARSON_FILL, PEARSON_LINE = "#e3ab74", "#b3651a"

LEGEND_PAD_PT = 6.0

#: A framed key, shared by B and C, so a sample marker is not read as data.
KEY_FRAME = dict(frameon=True, fancybox=False, framealpha=1.0, facecolor="white",
                 edgecolor="#000000", borderpad=0.6, handletextpad=0.5, borderaxespad=0.5)


FORMATS = ("pdf", "svg", "png", "tif")


def lw(points):
    """A line width scaled to the print type size, never thinner than MIN_LINEWIDTH."""
    return max(points * TYPE_SCALE, MIN_LINEWIDTH)


def apply_rcparams():
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT_FAMILY],
        "font.size": FONT_TICK,
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT_FAMILY,
        "mathtext.it": FONT_FAMILY + ":italic",
        "mathtext.bf": FONT_FAMILY + ":bold",
        "pdf.fonttype": 42,      # editable text in the PDF
        "svg.fonttype": "none",
        "legend.handlelength": 1.4,
        "legend.labelspacing": 0.3,
        "legend.borderaxespad": 0.3,
    })


def resolve_formats(spec):
    if not spec:
        return ("pdf",)
    formats = tuple(f.strip().lower() for f in spec.split(",") if f.strip())
    formats = tuple("tif" if f == "tiff" else f for f in formats)
    unknown = [f for f in formats if f not in FORMATS]
    if unknown:
        raise SystemExit("unknown format(s): %s" % ", ".join(unknown))
    return formats


def below_axis_anchor(axis, pad_pt=LEGEND_PAD_PT):
    """Axes-fraction y clearing everything already drawn under `axis`.

    Measured, not chosen: a constant that clears the labels at one box size collides at another.
    """
    figure = axis.figure
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    candidates = list(axis.get_xticklabels(which="both"))
    candidates += [axis.xaxis.get_offset_text(), axis.xaxis.label]
    lowest = axis.get_window_extent(renderer).y0
    for artist in candidates:
        if artist.get_visible() and artist.get_text():
            lowest = min(lowest, artist.get_window_extent(renderer).y0)
    pad_px = pad_pt * figure.dpi / 72.0
    return float(axis.transAxes.inverted().transform((0.0, lowest - pad_px))[1])


def _save_tiff(figure, destination, tight):
    """A flattened RGB, LZW-compressed TIFF with a 2-pt white border, checked against PLOS.

    Not `figure.savefig(.tif)`: that would write an alpha channel, which PLOS rejects.
    """
    from PIL import Image, ImageOps

    buffer = io.BytesIO()
    kwargs = {"bbox_inches": "tight"} if tight else {}
    figure.savefig(buffer, format="png", dpi=SAVE_DPI, facecolor="white", **kwargs)
    buffer.seek(0)
    rgba = Image.open(buffer).convert("RGBA")
    white = Image.new("RGBA", rgba.size, (255, 255, 255, 255))
    rgb = Image.alpha_composite(white, rgba).convert("RGB")
    border = int(round(BORDER_PT / 72.0 * SAVE_DPI))
    rgb = ImageOps.expand(rgb, border=border, fill=(255, 255, 255))
    rgb.save(destination, format="TIFF", compression="tiff_lzw", dpi=(SAVE_DPI, SAVE_DPI))

    w, h = rgb.size
    size = destination.stat().st_size
    problems = []
    if w > MAX_PIXELS[0] or h > MAX_PIXELS[1]:
        problems.append("%d x %d px exceeds PLOS %d x %d px (%.2f x %.2f in at %d dpi)"
                        % (w, h, MAX_PIXELS[0], MAX_PIXELS[1], w / SAVE_DPI, h / SAVE_DPI,
                           SAVE_DPI))
    if size > MAX_BYTES:
        problems.append("%.1f MB exceeds PLOS 10 MB" % (size / 1024.0 ** 2))
    if problems:
        raise SystemExit("%s: %s" % (destination, "; ".join(problems)))


def save(figure, output, formats=("pdf",), force=False, tight=True):
    """Write one file per format. Refuses to overwrite without `force`."""
    stem = Path(output)
    if stem.suffix.lstrip(".") in FORMATS:
        stem = stem.with_suffix("")
    stem.parent.mkdir(parents=True, exist_ok=True)
    destinations = [stem.with_suffix("." + fmt) for fmt in formats]
    existing = [p for p in destinations if p.exists()]
    if existing and not force:
        raise SystemExit("refusing to overwrite:\n%s\nPass --force."
                         % "\n".join("    %s" % p for p in existing))
    for fmt, destination in zip(formats, destinations):
        if fmt == "tif":
            _save_tiff(figure, destination, tight)
            continue
        kwargs = {"bbox_inches": "tight"} if tight else {}
        figure.savefig(destination, dpi=SAVE_DPI, **kwargs)
    return destinations
