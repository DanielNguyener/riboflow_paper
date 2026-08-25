"""The P-site offset detector, tested directly.

`get_offset_periodicity` decides every `psite_offset` this repository ships, so its frame
rule, its threshold scan and both of its fallbacks are asserted here on hand-built metagenes
rather than inferred from a downstream table.
"""
from __future__ import annotations

import pytest

from psite_offset import get_offset_periodicity, ribotish_get_offset


def metagene(offset=12, upstream_peak=None, downstream_reads=600, noise=5,
             window=30, frame_purity=1.0):
    """A 5'-end metagene with its P-site `offset` nt from the read's 5' end.

    5'-ends of in-frame reads sit at positions congruent to `-offset` mod 3, so the
    downstream pileup is phased to that residue class. `frame_purity` dilutes it.
    """
    counts = {p: noise for p in range(-40, window)}
    counts[-offset] = upstream_peak if upstream_peak is not None else 400
    frame = (-offset) % 3
    in_frame = int(round(downstream_reads * frame_purity))
    spread = downstream_reads - in_frame
    for p in range(0, window):
        if p % 3 == frame:
            counts[p] = counts.get(p, 0) + in_frame // (window // 3)
        else:
            counts[p] = counts.get(p, 0) + spread // (2 * (window // 3) or 1)
    return counts


# ── the canonical call ───────────────────────────────────────────────────────

@pytest.mark.parametrize("offset", [11, 12, 13])
def test_a_clean_metagene_returns_its_own_offset(offset):
    assert get_offset_periodicity(metagene(offset=offset)) == offset


def test_the_offset_is_positive_and_measured_from_the_read_5_prime_end():
    """The convention is a distance ALONG THE READ, so the return is positive even though
    the metagene positions it is read off are negative."""
    assert get_offset_periodicity(metagene(offset=12)) == 12


# ── the frame comes from downstream periodicity, not one start peak ──────────

def test_a_bimodal_start_peak_does_not_flip_the_frame():
    """Two near-tied start positions in DIFFERENT frames is exactly the case the plain
    argmax rule gets wrong: whichever happens to be one read taller decides the frame."""
    counts = metagene(offset=12)
    counts[-11] = counts[-12] + 1         # taller, and -11 % 3 != -12 % 3
    assert get_offset_periodicity(counts) == 12
    assert ribotish_get_offset(counts) == 11, "the argmax rule is what this guards against"


# ── the two documented fallbacks ─────────────────────────────────────────────

def test_a_thin_downstream_signal_falls_back_to_the_argmax_rule():
    counts = metagene(offset=12, downstream_reads=30)     # below min_down = 200
    assert get_offset_periodicity(counts) == ribotish_get_offset(counts)


def test_an_unphased_downstream_signal_falls_back_to_the_argmax_rule():
    counts = metagene(offset=12, frame_purity=0.34)       # below dom_frac = 0.40
    assert get_offset_periodicity(counts) == ribotish_get_offset(counts)


def test_the_fallback_thresholds_are_the_documented_ones():
    """200 reads and a 0.40 dominant fraction, just above the 1/3 no-periodicity floor."""
    import inspect
    signature = inspect.signature(get_offset_periodicity)
    assert signature.parameters["min_down"].default == 200
    assert signature.parameters["dom_frac"].default == 0.40
    assert signature.parameters["win_codons"].default == 10


# ── degenerate input is refused, never guessed ───────────────────────────────

def test_no_upstream_signal_returns_the_default():
    assert get_offset_periodicity({p: 0 for p in range(-40, 30)}) == 12
    assert get_offset_periodicity({}) == 12


def test_a_custom_default_is_honoured():
    assert get_offset_periodicity({}, default=15) == 15


def test_the_scan_stays_inside_the_flank():
    """A peak outside +/- flank of defOffset must not be picked up: the search window is
    what keeps a spurious upstream pileup from becoming the offset."""
    counts = metagene(offset=12)
    counts[-30] = 10_000
    assert get_offset_periodicity(counts) == 12


# ── the two routes must run the identical algorithm ──────────────────────────

def test_both_route_steps_import_the_same_detector():
    """A drift between the genome and transcriptome offset calls would turn the route
    comparison into a comparison of offset tables."""
    from pathlib import Path
    qc = Path(__file__).resolve().parents[1] / "code" / "ribo_seq_qc"
    for name in ("01_readlen_psite_qc.py", "01t_readlen_psite_qc_transcriptome.py"):
        text = (qc / name).read_text()
        assert "from psite_offset import" in text or "import psite_offset" in text, name
        assert "def ribotish_get_offset" not in text, "%s re-implements the detector" % name
        assert "def get_offset_periodicity" not in text, "%s re-implements the detector" % name
