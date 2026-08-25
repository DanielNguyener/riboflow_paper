"""Panel layout: a legend placed under an axis must not land on the x-axis label.

A legend anchored at a constant *axes fraction* (`bbox_to_anchor=(0.5, -0.16)`) collides
with a two-line x label on a short axes, because tick labels and the label occupy a
constant number of *points*. `panel_style.below_axis_anchor` measures instead. These tests
pin that at every axes height, and pin the one type scale and the shared cohort plot box
the Figure 5 panels depend on. Pure matplotlib; no BAM, no annotation, no reference render.

Run with `python` (3.9).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PANELS = REPO / "code" / "panels"

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib.patches import Patch                              # noqa: E402

def _load(name):
    """Import a `code/panels/` module by path, with `code/panels` importable.

    The generators do `sys.path.insert(0, HERE); import panel_style`, so the directory has
    to be on the path for the module to execute at all.
    """
    if name in sys.modules:
        return sys.modules[name]
    if str(PANELS) not in sys.path:
        sys.path.insert(0, str(PANELS))
    spec = importlib.util.spec_from_file_location(name, str(PANELS / ("%s.py" % name)))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def ps():
    return _load("panel_style")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ── helpers ───────────────────────────────────────────────────────────────────
def _below_axis_extents(axis):
    """Display-space bboxes of everything the panel draws under the x axis."""
    renderer = axis.figure.canvas.get_renderer()
    artists = list(axis.get_xticklabels(which="both"))
    artists.append(axis.xaxis.get_offset_text())
    artists.append(axis.xaxis.label)
    return [a.get_window_extent(renderer) for a in artists
            if a.get_visible() and a.get_text()]


def _assert_legend_clears_labels(axis, legend, pad_pt):
    """The legend's top must sit at or below the lowest label, by at least `pad_pt`."""
    figure = axis.figure
    figure.canvas.draw()
    legend_box = legend.get_window_extent(figure.canvas.get_renderer())
    lowest = min(box.y0 for box in _below_axis_extents(axis))
    gap_pt = (lowest - legend_box.y1) * 72.0 / figure.dpi
    assert gap_pt >= pad_pt - 0.5, (
        "legend top is %.2f pt from the lowest x-axis label (wanted >= %.2f); "
        "negative means they overlap" % (gap_pt, pad_pt))
    for box in _below_axis_extents(axis):
        assert not legend_box.overlaps(box)


def _two_line_axes(height_in, xlabel="% of transcriptome-\nassigned reads"):
    figure, axis = plt.subplots(figsize=(3.4, height_in))
    axis.barh([0, 1], [40, 90])
    axis.set_xlim(0, 100)
    axis.set_xlabel(xlabel)
    figure.tight_layout()
    return figure, axis


# ── below_axis_anchor ─────────────────────────────────────────────────────────
@pytest.mark.parametrize("height_in", [1.8, 2.4, 3.8, 6.0, 10.0])
def test_anchor_clears_a_two_line_label_at_every_height(ps, height_in):
    """The measured anchor clears the label whether the axes is short or tall."""
    figure, axis = _two_line_axes(height_in)
    anchor = ps.below_axis_anchor(axis)
    assert anchor < 0.0, "the anchor must be below the axes"

    y_display = axis.transAxes.transform((0.0, anchor))[1]
    lowest = min(box.y0 for box in _below_axis_extents(axis))
    gap_pt = (lowest - y_display) * 72.0 / figure.dpi
    assert gap_pt == pytest.approx(ps.LEGEND_PAD_PT, abs=0.5)


def test_the_constant_fraction_collides_and_the_measurement_does_not(ps):
    """The regression itself: on a 3.8 in axes with a two-line x label, an anchor at
    -0.16 lands ON the label; the measured anchor clears it. If -0.16 ever stops colliding
    here, the measurement is no longer buying anything and this test should be revisited
    rather than deleted.
    """
    figure, axis = _two_line_axes(3.8)
    figure.canvas.draw()
    lowest = min(box.y0 for box in _below_axis_extents(axis))
    old = axis.transAxes.transform((0.0, -0.16))[1]
    assert old > lowest, "the old constant offset was supposed to collide, but did not"
    new = axis.transAxes.transform((0.0, ps.below_axis_anchor(axis)))[1]
    assert new < lowest


def test_anchor_ignores_empty_labels(ps):
    """A blank x label must not drag the anchor down to a stray zero-size extent."""
    figure, axis = _two_line_axes(3.0, xlabel="")
    with_ticks = ps.below_axis_anchor(axis)
    axis.set_xlabel("one line")
    figure.canvas.draw()
    with_label = ps.below_axis_anchor(axis)
    assert with_label < with_ticks, "adding a label must push the anchor further down"


def test_pad_is_in_points_not_fractions(ps):
    """Doubling the pad moves the anchor by exactly that many points, at any height."""
    for height_in in (2.0, 6.0):
        figure, axis = _two_line_axes(height_in)
        near = axis.transAxes.transform((0.0, ps.below_axis_anchor(axis, 4.0)))[1]
        far = axis.transAxes.transform((0.0, ps.below_axis_anchor(axis, 16.0)))[1]
        moved_pt = (near - far) * 72.0 / figure.dpi
        assert moved_pt == pytest.approx(12.0, abs=0.25)
        plt.close(figure)


# ── legend_below ──────────────────────────────────────────────────────────────
def test_legend_below_places_and_does_not_overlap(ps):
    figure, axis = _two_line_axes(3.8)
    legend = ps.legend_below(axis, handles=[Patch(color="#a6d96a", label="kept"),
                                            Patch(color="#dddddd", label="other")])
    _assert_legend_clears_labels(axis, legend, ps.LEGEND_PAD_PT)


def test_legend_below_passes_keywords_through(ps):
    figure, axis = _two_line_axes(3.8)
    legend = ps.legend_below(axis, handles=[Patch(color="#a6d96a", label="kept")],
                             fontsize=7, ncol=1, title="fate")
    assert legend.get_title().get_text() == "fate"
    assert legend.get_frame_on() is False
    assert legend.get_texts()[0].get_fontsize() == 7


def test_legend_below_accepts_handles_and_labels_positionally(ps):
    figure, axis = _two_line_axes(3.8)
    legend = ps.legend_below(axis, [Patch(color="#a6d96a")], ["kept"])
    assert [t.get_text() for t in legend.get_texts()] == ["kept"]


# ── the export keeps the legend ───────────────────────────────────────────────
def test_export_contains_a_legend_placed_below(ps, tmp_path):
    """Moving the legend down must not crop it away: `bbox_inches='tight'` grows the saved
    bounds to hold it, so the export is strictly taller with the legend than without."""
    figure, axis = _two_line_axes(3.8)
    legend = ps.legend_below(axis, handles=[Patch(color="#a6d96a", label="kept"),
                                            Patch(color="#dddddd", label="other")])
    written = ps.save(figure, tmp_path / "panel", ("png",), extra_artists=[legend])
    height_with = plt.imread(str(written[0])).shape[0]
    legend.remove()
    cropped = ps.save(figure, tmp_path / "cropped", ("png",))
    height_without = plt.imread(str(cropped[0])).shape[0]
    assert height_with > height_without


# ── one type scale, enforced ─────────────────────────────────────────────────

def test_no_panel_hard_codes_a_font_size():
    """Every size comes from `panel_style`, so the panels share one type scale.

    Figure 2 carried six literals -- 7 pt ticks and 9 pt axis labels against the shared
    9 / 11, plus 6.5 and 5.3 pt in-cell text -- so its type rendered visibly smaller than
    Figures 3-5 once the panels were assembled at one scale. Nothing caught it: each panel
    rendered fine on its own, and the mismatch only existed *between* panels.

    A literal is the failure mode, not a particular wrong number: `fontsize=9` and
    `fontsize=ps.FONT_TICK` look identical until the constant moves and one panel silently
    keeps the old value. So the rule is the absence of literals, which is checkable.
    """
    import re

    panels = REPO / "code" / "panels"
    literal = re.compile(r"\b(?:fontsize|labelsize|titlesize)\s*=\s*[0-9]")
    offenders = []
    for path in sorted(panels.glob("*.py")):
        if path.name == "panel_style.py":
            continue        # where the sizes are DEFINED, and described in prose
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if literal.search(line) and not line.lstrip().startswith("#"):
                offenders.append("%s:%d: %s" % (path.name, number, line.strip()))
    assert not offenders, (
        "hard-coded font size(s) -- use a panel_style constant so every panel moves "
        "together:\n  " + "\n  ".join(offenders))


def test_every_generator_applies_the_shared_rcparams():
    """A generator that skips `apply_rcparams` gets matplotlib's default font, not ours.

    Checked by import rather than by grep: a helper module that only draws into an axes its
    caller created does not need to call it, so the test asks which modules actually build
    a figure.
    """
    panels = REPO / "code" / "panels"
    offenders = []
    for path in sorted(panels.glob("*.py")):
        if path.name in ("panel_style.py",):
            continue
        text = path.read_text()
        builds_figure = "plt.subplots(" in text or "plt.figure(" in text
        if builds_figure and "apply_rcparams" not in text:
            offenders.append(path.name)
    assert not offenders, (
        "these modules build a figure without applying the shared rcParams: %s"
        % ", ".join(offenders))


# ── the stacked pair: figure 5 A and B ────────────────────────────────────────
TAXONOMY = REPO / "data" / "read_taxonomy" / "taxonomy" / "taxonomy_all.tsv"
SAMPLES = REPO / "supporting_information" / "S1_Table" / "samples.csv"


def _axes_box_pt(figure):
    position = figure.axes[0].get_position()
    width_in, height_in = figure.get_size_inches()
    return (position.width * width_in * 72.0, position.height * height_in * 72.0,
            position.y0 * height_in * 72.0, position.y1 * height_in * 72.0)


TIE_MASTER = REPO / "data" / "read_taxonomy" / "multimap_biotype" / "multimap_tie_biotype_all.tsv"
REACH_MASTER = REPO / "data" / "read_taxonomy" / "reach" / "genome_anchored_reach_all.tsv"


def _cohort_figures():
    """The four cohort panels, drawn at their shipped defaults: {letter: figure}."""
    samples = SAMPLES if SAMPLES.exists() else None
    panel_a = _load("plot_route_read_counts")
    panel_b = _load("plot_read_id_union")
    panel_c = _load("plot_multimap_biotype")
    panel_d = _load("plot_nonselected_isoform_reach")
    side = _load("_fig05_side_panel")
    return {
        "A": panel_a.draw(panel_a.prepare(TAXONOMY, samples))[0],
        "B": panel_b.draw(panel_b.prepare(TAXONOMY, samples), show_labels=False)[0],
        # C and D share one helper and pass their own labels in from `main`; the geometry
        # under test is the helper's, so the test calls it the way they do.
        "C": side.draw_side_panel(
            panel_c.prepare(TIE_MASTER, TAXONOMY, samples)["values"], [""] * 24, "#1a7d1a",
            "protein-coding\npseudogene tie", "% of multimapper\nreads (primary)")[0],
        "D": side.draw_side_panel(
            panel_d.prepare(REACH_MASTER, TAXONOMY, samples)["values"], [""] * 24, "#7fb9da",
            "alternative\nisoform exon", "% of genome-only\nunique reads")[0],
    }


@pytest.mark.skipif(not (TAXONOMY.exists() and TIE_MASTER.exists() and REACH_MASTER.exists()),
                    reason="needs the shipped read-taxonomy masters")
def test_fig05_cohort_panels_share_one_plot_box():
    """A, B, C and D draw the same 24 cell lines, and only A carries the tick labels.

    So every row has to sit at the same height once the panels are stacked side by side, and
    the assembler stacks them by their page TOP. Two things must therefore match: the box
    height (or the rows sit at different pitches) and the box's offset below the page top (or
    equal pitches still start at different heights -- C and D have a title band above their
    axes and A does not). Both come from `fig05_common`. A quarter point is a rounding
    allowance, not room to drift.
    """
    boxes = {}
    for letter, figure in _cohort_figures().items():
        position = figure.axes[0].get_position()
        width_in, height_in = figure.get_size_inches()
        page_pt = height_in * 72.0
        boxes[letter] = (position.height * height_in * 72.0,          # box height
                         page_pt - position.y1 * page_pt)             # page top -> box top
        plt.close(figure)

    heights = {k: v[0] for k, v in boxes.items()}
    offsets = {k: v[1] for k, v in boxes.items()}
    assert max(heights.values()) - min(heights.values()) < 0.25, (
        "plot boxes differ in height, so the 24 rows sit at different pitches: %s"
        % {k: "%.2f pt (%.2f/row)" % (v, v / 24) for k, v in heights.items()})
    assert max(offsets.values()) - min(offsets.values()) < 0.25, (
        "plot boxes are equal but start at different heights below the page top: %s"
        % {k: "%.2f pt" % v for k, v in offsets.items()})
