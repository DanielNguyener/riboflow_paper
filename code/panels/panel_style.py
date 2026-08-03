#!/usr/bin/env python3
"""Shared styling and output handling for the independent panel generators."""
from __future__ import annotations

from pathlib import Path

FONT_FAMILY = "DejaVu Sans"
FONT_TITLE = 11
FONT_LABEL = 11
FONT_TICK = 9
FONT_ANNOTATION = 9
#: tick size on purpose -- an inset is a legend, not a second axis.
FONT_INSET = 8

GENOME = "#3a923a"
TXOME = "#cc3d3d"
SPEARMAN_FILL = "#8fb4d6"
SPEARMAN_LINE = "#2c5f8a"
PEARSON_FILL = "#e3ab74"
PEARSON_LINE = "#b3651a"
MISSING_HATCH = "#bbbbbb"

DEFAULT_FORMATS = ("pdf",)
PNG_DPI = 200

LEGEND_PAD_PT = 8.0

def apply_rcparams():
    """One font family and one set of sizes, with editable text in the PDF."""
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT_FAMILY],
        "font.size": FONT_TICK,
        "axes.titlesize": FONT_TITLE,
        "axes.labelsize": FONT_LABEL,
        "xtick.labelsize": FONT_TICK,
        "ytick.labelsize": FONT_TICK,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })

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

    A hand-picked `bbox_to_anchor=(0.5, -0.16)` is an *axes fraction*, so the gap it opens
    scales with the axes height -- while the tick labels and the x-axis label below it
    occupy a fixed number of points. On a tall axes the constant is a generous gap; on a
    short one (two horizontal bars, say) it is smaller than a two-line x label, and the
    legend lands on top of the label. This measures the tick labels, the axis offset text
    and the x label as rendered, then converts a pad in points into the fraction that
    clears the lowest of them.

    Call it after the figure's layout pass (`tight_layout`), because that pass moves the
    axes. Returns a negative number: the y for `loc="upper center"` in `transAxes`.
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
    """A legend centred under the x axis, positioned clear of the x label -- not guessed.

    Thin wrapper over `Axes.legend` that supplies `loc`/`bbox_to_anchor` from
    `below_axis_anchor`; any other legend keyword passes straight through. The returned
    legend lives outside the axes, so hand it to `save(..., extra_artists=[legend])` or the
    export will crop it.
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

def save(figure, output, formats=DEFAULT_FORMATS, force=False, extra_artists=None):
    """Write one figure in each requested format. REFUSES to overwrite without `force`.

    `output` may carry an extension or not; the stem is what matters. Returns the paths
    written, so a caller can record exactly what it produced.
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
        kwargs = {"bbox_inches": "tight"}
        if extra_artists:
            kwargs["bbox_extra_artists"] = list(extra_artists)
        if fmt == "png":
            kwargs["dpi"] = PNG_DPI
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
