"""Gene-level read membership and the two Figure 5A mechanisms, on synthetic data.

A genome alignment places its read in a gene's union only if it is the read's primary or a
secondary whose AS equals the primary's, and it overlaps the gene's full span. Three places
apply that rule and must agree: the indexed fetch Figure 5A uses
(`gene_read_partition_lib.fetch_gene_candidates` + `resolve_genome_side`), the read-state
slice Figure 6 uses (`read_state.ReadState.gene_side`) and the locus status Figure 5B labels
its transcriptome track with (`build_locus_data.locus_ribo_reads`).

The mechanisms are tested at the gene, not at whichever placement STAR called primary: the
protein-coding / processed-pseudogene tie (`gene_pseudogene_tie`) and the alternative exon
(`alt_exon_overlap`).

Run with `python` (3.9).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import _write_bam

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"

CHROM = "chrG"
GENE = (1_000, 3_000)                 # the gene's full span, half-open
RUN = "SRRT"                          # names `<run>.<serial>` pack exactly in the store
TOP, LOW = 50, 40                     # alignment scores


def _name(serial):
    return "%s.%d" % (RUN, serial)


#: serial -> (role, on the genome side, primary NH or None when there is no primary)
READS = {
    1: ("primary elsewhere, TIED secondary at the gene (read B)", True, 2),
    2: ("primary elsewhere, LOWER-scoring secondary at the gene", False, 2),
    3: ("secondary at the gene, no primary record", False, None),
    4: ("unique, at the gene", True, 1),
    5: ("unique, spliced from before the gene into it", True, 1),
    6: ("unique, elsewhere", False, 1),
    7: ("primary at the gene, tied secondary elsewhere (read A)", True, 2),
    8: ("primary elsewhere, secondary at the gene with no AS", False, 2),
}


def _genome_records():
    def r(serial, pos, cigar="30M", secondary=False, nh=1, score=TOP):
        return {"ref": CHROM, "pos": pos, "cigar": cigar, "secondary": secondary, "nh": nh,
                "as_": score, "name": _name(serial)}
    return [
        r(1, 8_000, nh=2), r(1, 1_500, secondary=True, nh=2),
        r(2, 8_200, nh=2), r(2, 1_550, secondary=True, nh=2, score=LOW),
        r(3, 1_600, secondary=True, nh=2),
        r(4, 1_800),
        r(5, 900, cigar="10M500N20M"),        # blocks [900, 910) and [1410, 1430)
        r(6, 9_000),
        r(7, 1_200, nh=2), r(7, 8_500, secondary=True, nh=2),
        r(8, 8_700, nh=2), r(8, 1_700, secondary=True, nh=2, score=None),
    ]


@pytest.fixture(scope="module")
def bams(tmp_path_factory):
    root = tmp_path_factory.mktemp("membership")
    genome = _write_bam(root / "genome.bam", [(CHROM, 20_000)], _genome_records())
    txome = _write_bam(root / "txome.bam", [("ENSTT0001.1", 2_000)],
                       [{"ref": "ENSTT0001.1", "pos": 100, "cigar": "30M",
                         "name": _name(serial)} for serial in (1, 3)])
    return genome, txome


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def lib():
    return _load(CODE / "alignment_fate" / "gene_read_partition_lib.py",
                 "gene_read_partition_lib")


@pytest.fixture(scope="module")
def tie_lib():
    for entry in (CODE / "read_taxonomy", CODE / "common", CODE / "common" / "ribo_seq_qc"):
        if str(entry) not in sys.path:
            sys.path.insert(0, str(entry))
    import tie_biotype_lib
    return tie_biotype_lib


@pytest.fixture(scope="module")
def expected():
    return {_name(s) for s, (_role, on_side, _nh) in READS.items() if on_side}


def _primaries(bam):
    import pysam
    out = {}
    with pysam.AlignmentFile(str(bam)) as handle:
        for read in handle.fetch(until_eof=True):
            if not read.is_secondary:
                out[read.query_name] = (read.reference_name, "+", read.get_blocks(),
                                        read.get_tag("NH"), read.get_tag("AS"))
    return out


# ── membership: one rule, three paths ─────────────────────────────────────────

def test_fetch_keeps_primary_and_tied_secondaries(lib, bams, expected):
    candidates = lib.fetch_gene_candidates(bams[0], CHROM, *GENE)
    side, joined_by, audit = lib.resolve_genome_side(candidates, _primaries(bams[0]))
    assert side == expected
    assert joined_by[_name(1)] == "tied_secondary"
    assert joined_by[_name(7)] == "primary"
    assert audit == {"n_joined_by_tied_secondary_only": 1, "n_secondary_only_no_primary": 1,
                     "n_secondary_only_lower_score": 1, "n_missing_as": 1}


def test_read_state_gene_side_matches_the_fetch(bams, expected, tmp_path):
    sys.path.insert(0, str(CODE / "clustering"))
    import read_state
    path = read_state.build("T", bams[0], bams[1], {"ENSTT0001": "ENSTT0001.1"},
                            tmp_path / "state.h5", log=lambda _m: None)
    with read_state.ReadState(path) as state:
        index = {int(key) & ((1 << 40) - 1): i for i, key in enumerate(state.read_key)}
        side = set(state.gene_side(CHROM, *GENE).tolist())
        assert {_name(serial) for serial, i in index.items() if i in side} == expected
        for serial, (_role, _side, nh) in READS.items():
            i = index[serial]
            assert (i in state.genome_present) == (nh is not None), serial
            assert (i in state.genome_unique) == (nh == 1), serial
        assert len(state.records_dict([index[1]])[index[1]]) == 2


def test_locus_status_uses_the_same_rule(bams, expected):
    locus = _load(CODE / "alignment_fate" / "build_locus_data.py", "build_locus_data")
    _members, status, audit = locus.locus_ribo_reads(
        str(bams[0]), CHROM, *GENE, lengths={30}, offsets={30: 12})
    assert set(status) == expected
    assert status == {_name(s): nh for s, (_r, side, nh) in READS.items() if side}
    assert audit == {"joined_by_tied_secondary_only": 1, "secondary_only_no_primary": 1,
                     "secondary_only_lower_score": 1, "missing_as": 1}
    # Reads A and B tie the same way; which one STAR made primary does not matter.
    assert (locus.txome_read_layer(_name(1), status)
            == locus.txome_read_layer(_name(7), status) == "shared_multi")
    assert locus.txome_read_layer(_name(2), status) == "txome_only"


def test_qualifies(lib):
    assert lib.qualifies(False, lib.MISSING_AS, None)          # the primary always
    assert lib.qualifies(True, TOP, TOP)
    assert not lib.qualifies(True, LOW, TOP)
    assert not lib.qualifies(True, TOP, None)                  # no primary record
    assert not lib.qualifies(True, lib.MISSING_AS, lib.MISSING_AS)


# ── the protein-coding / processed-pseudogene tie, at THIS gene ───────────────

GENE_EXON, PP_EXON, OTHER_PC_EXON = (1_000, 1_100), (5_000, 5_100), (7_000, 7_100)


@pytest.fixture(scope="module")
def annotation():
    import pyranges as pr
    exons = pd.DataFrame({
        "Chromosome": [CHROM] * 4,
        "Start": [GENE_EXON[0], 1_500, PP_EXON[0], OTHER_PC_EXON[0]],
        "End": [GENE_EXON[1], 1_600, PP_EXON[1], OTHER_PC_EXON[1]],
        "gene_id": ["ENSGG.1", "ENSGG.1", "ENSGP.1", "ENSGO.1"],
        "gene_type": ["protein_coding", "protein_coding", "processed_pseudogene",
                      "protein_coding"]})
    return {
        "exon_gene_df": exons,
        "exon_pr": pr.PyRanges(exons),
        "gene_body_pr": pr.PyRanges(pd.DataFrame(
            {"Chromosome": [CHROM], "Start": [0], "End": [20_000], "Strand": ["+"]})),
        # Selected transcript: [1000, 1100) and [2000, 2100); [1500, 1600) is omitted.
        "table": {"ENSTG.1": {"gene_id": "ENSGG.1", "g_start": np.array([1_000, 2_000]),
                              "g_end": np.array([1_100, 2_100])}},
    }


def _rec(pos5, score=TOP, secondary=True):
    return (CHROM, pos5, score, secondary)


def test_pseudogene_tie_is_tested_at_the_gene(lib, tie_lib, annotation):
    records = {
        "A": [_rec(1_050, secondary=False), _rec(5_050)],        # primary on the gene
        "B": [_rec(5_050, secondary=False), _rec(1_050)],        # primary on the pseudogene
        "other_pc": [_rec(7_050, secondary=False), _rec(5_050), _rec(1_050, score=LOW)],
        "three": [_rec(5_060, secondary=False), _rec(1_050), _rec(5_050)],
        "lower": [_rec(1_050, secondary=False), _rec(5_050, score=LOW)],
        "gene_only": [_rec(1_050, secondary=False), _rec(1_060)],
    }
    ties = lib.gene_pseudogene_tie(tie_lib, annotation, "ENSGG.1", records)
    got = {q: (tie, n_top) for q, (tie, n_top, _detail) in ties.items()}
    assert got == {"A": (True, 2), "B": (True, 2), "other_pc": (False, 2),
                   "three": (True, 3), "lower": (False, 1), "gene_only": (False, 2)}


# ── the alternative exon: omitted exonic sequence, nothing else ──────────────

@pytest.mark.parametrize("strand, blocks, want", [
    ("+", [(1_520, 1_550)], True),                      # on the omitted exon
    ("-", [(1_520, 1_550)], True),                      # strand does not enter
    ("+", [(1_590, 1_620)], True),                      # partly on it
    ("+", [(1_000, 1_050), (2_000, 2_050)], False),     # an unannotated junction, selected exons
    ("+", [(1_080, 1_150)], False),                     # runs into the intron, not the omitted exon
    ("-", [(1_000, 1_030)], False),                     # strand mismatch alone is not enough
])
def test_alt_exon_overlap(lib, annotation, strand, blocks, want):
    primary = {"r": (CHROM, strand, blocks, 1, TOP)}
    assert lib.alt_exon_overlap(annotation, "ENSTG.1", primary, ["r"]) == {"r": want}


# ── the chain and the fold refuse what can no longer happen ──────────────────

def test_chain_rejects_a_genome_side_read_without_a_primary(lib):
    with pytest.raises(lib.PartitionError, match="no primary genome alignment"):
        lib.classify_union((None, None, None), {}, "T", (CHROM, *GENE), {"orphan"}, set(),
                           {}, set(), set(), {}, {})


@pytest.mark.parametrize("category", (
    "genome_multi_primary_in_gene_pseudogene_tie", "genome_multi_primary_pseudogene",
    "genome_multi_primary_elsewhere_other", "genome_multi_primary_lost_in_dedup",
    "genome_unique_absent_nonselected_isoform_exon"))
def test_retired_categories_fail_validation(lib, category):
    assert category not in lib.PARTITION_CATEGORIES
    fold = _load(CODE / "panels" / "plot_gene_read_partition.py", "plot_gene_read_partition")
    with pytest.raises(SystemExit, match="unknown chain category"):
        fold._route7_segment(category, True)
