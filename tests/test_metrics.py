"""The two correlation definitions, and why they are written out rather than imported.

`compute_coverage_concordance` defines both metrics used everywhere in this repository:

    _spear(g, t)     Spearman on RAW counts       -- rank-invariant, so normalization
                                                     choices cannot move it
    _pe_log2(g, t)   Pearson on log2(count + 1)   -- centred by hand in float64

The log2 transform is not cosmetic: per-base coverage spans four orders of magnitude, and
a Pearson on raw counts would be dominated by a handful of peak positions. The pseudocount
keeps zeros finite.
"""
from __future__ import annotations

import numpy as np
import pytest

import compute_coverage_concordance as ccc
from conftest import TX_MINUS, TX_PLUS, cds_interior


# ── Spearman on raw counts ───────────────────────────────────────────────────

def test_spearman_of_a_vector_with_itself_is_one():
    x = np.array([0, 1, 5, 2, 9, 3], dtype=np.int32)
    assert ccc._spear(x, x) == pytest.approx(1.0)


def test_spearman_is_invariant_to_a_monotonic_rescaling():
    """This is the property the metric is chosen for: any per-sample normalization is a
    positive monotone transform and therefore cannot change it."""
    x = np.array([0, 1, 5, 2, 9, 3], dtype=np.int32)
    y = np.array([1, 3, 8, 4, 20, 6], dtype=np.int32)
    assert ccc._spear(x, y) == pytest.approx(ccc._spear(x, 7 * y))
    assert ccc._spear(x, y) == pytest.approx(ccc._spear(x, y.astype(float) / 1000.0))


def test_spearman_of_an_exact_reversal_is_minus_one():
    x = np.arange(10, dtype=np.int32)
    assert ccc._spear(x, x[::-1].copy()) == pytest.approx(-1.0)


def test_spearman_matches_scipy_on_raw_counts():
    scipy_stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(0)
    x = rng.integers(0, 40, 300).astype(np.int32)
    y = rng.integers(0, 40, 300).astype(np.int32)
    assert ccc._spear(x, y) == pytest.approx(
        scipy_stats.spearmanr(x, y).correlation, abs=1e-12)


# ── Pearson on log2(count + 1) ───────────────────────────────────────────────

def test_pearson_log2_of_a_vector_with_itself_is_one():
    x = np.array([0, 1, 5, 2, 9, 3], dtype=np.int32)
    assert ccc._pe_log2(x, x) == pytest.approx(1.0)


def test_pearson_log2_is_computed_on_the_log_scale_not_the_raw_one():
    """A single large peak dominates a raw Pearson and barely moves a log2 one. If these
    two agreed, the transform would not be being applied."""
    base = np.array([1, 2, 3, 4, 5] * 20, dtype=np.int32)
    other = base.copy()
    other[0] = 10_000
    spiked = base.copy()
    spiked[0] = 10_000
    raw = float(np.corrcoef(base, other)[0, 1])
    assert ccc._pe_log2(base, other) != pytest.approx(raw, abs=1e-6)
    assert ccc._pe_log2(spiked, other) == pytest.approx(1.0)


def test_the_pseudocount_keeps_all_zero_positions_finite():
    x = np.zeros(50, dtype=np.int32)
    y = np.zeros(50, dtype=np.int32)
    y[:5] = 3
    value = ccc._pe_log2(x, y)
    assert value == 0.0 or np.isnan(value), "a constant vector has no correlation"
    assert not np.isinf(value)


def test_the_pseudocount_is_a_parameter_not_a_constant():
    x = np.array([0, 1, 2, 8, 0, 4], dtype=np.int32)
    y = np.array([0, 2, 1, 9, 1, 3], dtype=np.int32)
    assert ccc._pe_log2(x, y, pseudocount=1.0) != \
        pytest.approx(ccc._pe_log2(x, y, pseudocount=0.1))


def test_pearson_log2_matches_a_hand_written_reference():
    """The implementation centres in float64 by hand. This reproduces that arithmetic
    independently, so the test would catch a change of formulation rather than agreeing
    with it by construction."""
    x = np.array([0, 1, 2, 8, 0, 4], dtype=np.int32)
    y = np.array([0, 2, 1, 9, 1, 3], dtype=np.int32)
    a = np.log2(x.astype(np.float64) + 1.0)
    b = np.log2(y.astype(np.float64) + 1.0)
    a -= a.mean()
    b -= b.mean()
    expected = float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))
    assert ccc._pe_log2(x, y) == pytest.approx(expected, abs=1e-15)


# ── degenerate input ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("metric", [ccc._spear, ccc._pe_log2])
def test_an_empty_vector_does_not_raise(metric):
    empty = np.zeros(0, dtype=np.int32)
    value = metric(empty, empty)
    assert np.isnan(value) or value == 0.0


@pytest.mark.parametrize("metric", [ccc._spear, ccc._pe_log2])
def test_a_constant_vector_has_no_correlation(metric):
    x = np.ones(20, dtype=np.int32)
    y = np.arange(20, dtype=np.int32)
    value = metric(x, y)
    assert np.isnan(value) or value == 0.0, "a zero-variance vector cannot correlate"


# ── against the real pipeline output ─────────────────────────────────────────

def test_identical_routes_give_a_perfect_correlation(coverage):
    """TX_MINUS receives the same coverage from both routes, so both metrics must be
    exactly 1. This is the end-to-end check that the two routes are on one coordinate."""
    genome = cds_interior(coverage, TX_MINUS, "genome_psite")
    txome = cds_interior(coverage, TX_MINUS, "txome_psite")
    assert ccc._spear(genome, txome) == pytest.approx(1.0)
    assert ccc._pe_log2(genome, txome) == pytest.approx(1.0)


def test_the_metrics_are_finite_on_a_real_transcript(coverage):
    genome = cds_interior(coverage, TX_PLUS, "genome_footprint")
    txome = cds_interior(coverage, TX_PLUS, "txome_footprint")
    assert np.isfinite(ccc._spear(genome, txome))
    assert np.isfinite(ccc._pe_log2(genome, txome))


def test_the_codon_view_groups_bases_into_triplets():
    """`_codon_view` reshapes a per-base vector into one row per codon, so column j is
    frame j. It does NOT aggregate: the three frames stay separable, which is the whole
    reason a codon view exists."""
    x = np.arange(9, dtype=np.int32)
    codons = ccc._codon_view(x)
    assert codons.shape == (3, 3)
    assert codons.tolist() == [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    assert list(codons[:, 0]) == [0, 3, 6], "column 0 is frame 0"


def test_the_codon_view_discards_a_trailing_partial_codon():
    """A partial codon at the 3' end has no frame-2 base, so including it would bias the
    frame comparison toward frames 0 and 1."""
    x = np.ones(10, dtype=np.int32)
    codons = ccc._codon_view(x)
    assert codons.shape == (3, 3)
    assert codons.size == 9, "the tenth base is dropped, not padded"
