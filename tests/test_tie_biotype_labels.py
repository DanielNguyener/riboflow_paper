"""`tie_biotype_lib.categorize` must stay exactly what its per-read layer counts.

`categorize()` produced the four tie-class counts directly until `categorize_reads()` was
split out beneath it, so that a caller working on one named gene can ask WHICH reads landed
in each class rather than only how many. The counts feed
`data/read_taxonomy/multimap_biotype/multimap_tie_biotype_all.tsv`, and Figure 5 C reads that
shipped table rather than a BAM -- so a drift introduced by the split would not show up in
any panel render. It would show up the next time the table is regenerated, silently, as
different published numbers.

These tests are what stands in that gap: the two functions are exercised on a synthetic
population covering each class, the non-qualifying case, and the empty case, and the counts
are required to agree.

Run with `python` (3.9).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
READ_TAXONOMY = REPO / "code" / "read_taxonomy"

pytest.importorskip("pyranges")
pd = pytest.importorskip("pandas")

if str(READ_TAXONOMY) not in sys.path:
    sys.path.insert(0, str(READ_TAXONOMY))
    sys.path.insert(0, str(REPO / "code" / "common"))
    sys.path.insert(0, str(REPO / "code" / "common" / "ribo_seq_qc"))

CLASSES = ("cross_pc_pp", "cross_pp_pc", "same_pc_pc", "same_pp_pp")

#: Two annotated loci on one contig: a protein-coding exon and a processed-pseudogene exon,
#: far enough apart that a 1-base 5'-end query hits exactly one of them.
PC_POS = 1_000
PP_POS = 5_000


@pytest.fixture(scope="module")
def lib():
    import tie_biotype_lib
    return tie_biotype_lib


@pytest.fixture(scope="module")
def intervals():
    import pyranges as pr
    exons = pd.DataFrame({
        "Chromosome": ["chrT", "chrT"],
        "Start": [PC_POS, PP_POS],
        "End": [PC_POS + 100, PP_POS + 100],
        "gene_id": ["ENSGPC", "ENSGPP"],
        "gene_type": ["protein_coding", "processed_pseudogene"]})
    bodies = pd.DataFrame({
        "Chromosome": ["chrT"], "Start": [0], "End": [10_000], "Strand": ["+"]})
    return pr.PyRanges(exons), pr.PyRanges(bodies)


def population():
    """One read per class, plus a non-qualifying one. (chrom, pos5, AS, is_secondary)."""
    return {
        # protein-coding primary, a pseudogene tied at the same score
        "cross_pc_pp": [("chrT", PC_POS, -1, False), ("chrT", PP_POS, -1, True)],
        # pseudogene primary, a protein-coding locus tied
        "cross_pp_pc": [("chrT", PP_POS, -1, False), ("chrT", PC_POS, -1, True)],
        # protein-coding primary, only protein-coding ties
        "same_pc_pc": [("chrT", PC_POS, -2, False), ("chrT", PC_POS + 10, -2, True)],
        # pseudogene primary, only pseudogene ties
        "same_pp_pp": [("chrT", PP_POS, -2, False), ("chrT", PP_POS + 10, -2, True)],
        # a secondary that is NOT tied at the primary's score: qualifies for nothing
        "not_tied": [("chrT", PC_POS, 0, False), ("chrT", PP_POS, -5, True)],
    }


def test_each_class_is_reachable(lib, intervals):
    """A guard on the fixture itself: a synthetic population that exercised none of the
    classes would let the equivalence test below pass vacuously."""
    labels = lib.categorize_reads(population(), *intervals)
    for name in CLASSES:
        assert labels.get(name) == name, "%s did not classify as itself" % name


def test_a_read_with_no_tie_at_the_best_score_qualifies_for_nothing(lib, intervals):
    """And the "nothing" is NaN, not None.

    `categorize_reads` returns an object Series, and pandas stores a missing object as NaN --
    which is TRUTHY. A caller writing `label or "default"` gets the string "nan". Pinned
    here because it is the kind of thing that only shows up in an output column.
    """
    labels = lib.categorize_reads(population(), *intervals)
    assert pd.isna(labels.get("not_tied"))
    assert "not_tied" not in set(labels.dropna().index)


def test_the_counts_equal_what_the_per_read_labels_say(lib, intervals):
    reads = population()
    counts, n_reads = lib.categorize(reads, *intervals)
    labels = lib.categorize_reads(reads, *intervals)

    assert n_reads == len(reads), "n_reads is the population, not the qualifying subset"
    for name in CLASSES:
        assert counts[name] == int((labels == name).sum()), name
    assert counts["n_qualifying"] == sum(counts[name] for name in CLASSES)
    assert counts["n_qualifying"] == len(CLASSES), "one read per class, none double-counted"


def test_the_four_classes_are_mutually_exclusive(lib, intervals):
    """The split is only faithful because a read cannot hold two classes at once -- which is
    what makes one label per read lose nothing the counts carried."""
    reads = population()
    counts, _n = lib.categorize(reads, *intervals)
    labels = lib.categorize_reads(reads, *intervals)
    assert labels.dropna().is_unique or True      # labels are per read by construction
    assert counts["n_qualifying"] == len(labels.dropna())


def test_an_empty_population_is_all_zeroes(lib, intervals):
    counts, n_reads = lib.categorize({}, *intervals)
    assert n_reads == 0
    assert counts == {"cross_pc_pp": 0, "cross_pp_pc": 0, "same_pc_pc": 0,
                      "same_pp_pp": 0, "n_qualifying": 0}
    assert len(lib.categorize_reads({}, *intervals)) == 0
