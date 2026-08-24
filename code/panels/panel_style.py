#!/usr/bin/env python3
"""Shared styling and output handling for the independent panel generators."""
from __future__ import annotations

from pathlib import Path

#: Arial per journal requirements; fallbacks are metric-compatible so a machine without
#: Arial renders at the same widths (matplotlib's DejaVu default is wider and shifts layout).
FONT_FAMILY = "Arial"
FONT_FALLBACKS = ("Helvetica", "Liberation Sans", "Arimo", "TeX Gyre Heros", "Nimbus Sans",
                  "DejaVu Sans")
FONT_TITLE = 11
FONT_LABEL = 11
FONT_TICK = 9
FONT_ANNOTATION = 9
#: tick size on purpose -- an inset is a legend, not a second axis.
FONT_INSET = 8

#: The single ENLARGED type scale, for panels reproduced large (heatmaps, ribo-vs-RNA).
FONT_LABEL_LARGE = 14
FONT_TICK_LARGE = 12
FONT_ANNOTATION_LARGE = 11

#: Linear enlargement of the *_LARGE scale. Marker `s` is an AREA -- scale it by the square.
LARGE_SCALE = FONT_LABEL_LARGE / FONT_LABEL          # ~1.27
LARGE_MARKER_AREA = LARGE_SCALE ** 2                 # ~1.62

#: Plot-box edge, in inches, for enlarged-scale panels. It is the AXES box, not the page,
#: so plot areas line up at assembly regardless of label size.
LARGE_AXES_BOX = 5.0

GENOME = "#3a923a"
TXOME = "#cc3d3d"
SPEARMAN_FILL = "#8fb4d6"
SPEARMAN_LINE = "#2c5f8a"
PEARSON_FILL = "#e3ab74"
PEARSON_LINE = "#b3651a"
MISSING_HATCH = "#bbbbbb"

DEFAULT_FORMATS = ("pdf",)
#: Applied to EVERY format: in a PDF it sets the resolution of rasterised artists
#: (the scatter clouds are rasterised on purpose); vector content is unaffected.
SAVE_DPI = 300
PNG_DPI = SAVE_DPI

LEGEND_PAD_PT = 8.0

def apply_rcparams():
    """One font family and one set of sizes, with editable text in the PDF."""
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT_FAMILY] + list(FONT_FALLBACKS),
        "font.size": FONT_TICK,
        "axes.titlesize": FONT_TITLE,
        "axes.labelsize": FONT_LABEL,
        "xtick.labelsize": FONT_TICK,
        "ytick.labelsize": FONT_TICK,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        # Mathtext ignores font.sans-serif -- point mathtext.* at the same family.
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT_FAMILY,
        "mathtext.it": "%s:italic" % FONT_FAMILY,
        "mathtext.bf": "%s:bold" % FONT_FAMILY,
    })

def resolve_font(family=FONT_FAMILY):
    """`(name, path)` of the first available family, walking FONT_FALLBACKS.

    matplotlib substitutes a missing family silently; this errors instead.
    """
    from matplotlib import font_manager

    for name in [family] + [f for f in FONT_FALLBACKS if f != family]:
        try:
            path = font_manager.findfont(font_manager.FontProperties(family=name),
                                         fallback_to_default=False)
        except Exception:
            continue
        return name, path
    raise SystemExit(
        "none of the panel fonts is installed: %s.\nInstall one (Arial ships with macOS "
        "and Office; Liberation Sans and Arimo are open metric-compatible substitutes) or "
        "change FONT_FAMILY in code/panels/panel_style.py."
        % ", ".join([family] + list(FONT_FALLBACKS)))

def resolve_formats(spec):
    """'pdf,svg,png' -> ('pdf', 'svg', 'png'), validated."""
    if not spec:
        return DEFAULT_FORMATS
    formats = tuple(f.strip().lower() for f in spec.split(",") if f.strip())
    unknown = [f for f in formats if f not in ("pdf", "svg", "png")]
    if unknown:
        raise SystemExit("unknown output format(s): %s. Choose from pdf, svg, png."
                         % ", ".join(unknown))
    return formats

def below_axis_anchor(axis, pad_pt=LEGEND_PAD_PT):
    """Axes-fraction y that clears everything already drawn under `axis`.

    Call after `tight_layout` (the layout pass moves the axes); returns a negative y.
    """
    figure = axis.figure
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()

    candidates = list(axis.get_xticklabels(which="both"))
    candidates.append(axis.xaxis.get_offset_text())
    candidates.append(axis.xaxis.label)
    lowest = axis.get_window_extent(renderer).y0
    for artist in candidates:
        if not artist.get_visible() or not artist.get_text():
            continue
        lowest = min(lowest, artist.get_window_extent(renderer).y0)

    pad_px = pad_pt * figure.dpi / 72.0
    return float(axis.transAxes.inverted().transform((0.0, lowest - pad_px))[1])

def legend_below(axis, handles=None, labels=None, pad_pt=LEGEND_PAD_PT, **kwargs):
    """A legend centred under the x axis, positioned clear of the x label.

    The legend lives outside the axes -- hand it to `save(..., extra_artists=[legend])`.
    """
    kwargs.setdefault("loc", "upper center")
    kwargs.setdefault("frameon", False)
    kwargs["bbox_to_anchor"] = (kwargs.pop("x", 0.5), below_axis_anchor(axis, pad_pt))
    kwargs.setdefault("bbox_transform", axis.transAxes)
    if handles is None:
        return axis.legend(**kwargs)
    if labels is None:
        return axis.legend(handles=handles, **kwargs)
    return axis.legend(handles, labels, **kwargs)

def save(figure, output, formats=DEFAULT_FORMATS, force=False, extra_artists=None,
         tight=True):
    """Write one figure in each requested format; refuses to overwrite without `force`.

    `tight=False` exports at exactly `figsize` (no crop) so co-placed panels scale equally.
    """
    output = Path(output)
    stem = output.with_suffix("") if output.suffix.lstrip(".") in ("pdf", "svg", "png") \
        else output
    stem.parent.mkdir(parents=True, exist_ok=True)

    destinations = [stem.with_suffix("." + fmt) for fmt in formats]
    existing = [p for p in destinations if p.exists()]
    if existing and not force:
        raise SystemExit(
            "refusing to overwrite %d existing file(s):\n%s\nPass --force if that is "
            "what you want." % (len(existing), "\n".join("    %s" % p for p in existing)))

    written = []
    for fmt, destination in zip(formats, destinations):
        kwargs = {"bbox_inches": "tight"} if tight else {}
        if extra_artists and tight:
            kwargs["bbox_extra_artists"] = list(extra_artists)
        kwargs["dpi"] = PNG_DPI if fmt == "png" else SAVE_DPI
        figure.savefig(destination, **kwargs)
        written.append(destination)
    return written

def require_columns(frame, columns, source):
    """Validate an input table's schema, naming every missing column at once."""
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise SystemExit(
            "%s is missing %d required column(s): %s\nPresent: %s"
            % (source, len(missing), ", ".join(missing), ", ".join(map(str, frame.columns))))
    return frame
