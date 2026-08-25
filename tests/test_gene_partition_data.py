"""Figure 6A's compact table: the shipped seven-segment partition carries the validated
counts, sums to each gene's union, and the fold program refuses anything else.

No per-read dump is needed: the builder's `check_expected` is exercised on structures
built here. Run with `python` (3.9).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "data" / "alignment_fate" / "gene_partition_route7.tsv"
META = REPO / "data" / "alignment_fate" / "gene_partition_route7.json"


def _load_builder():
    path = REPO / "code" / "alignment_fate" / "build_gene_partition_data.py"
    spec = importlib.util.spec_from_file_location("build_gene_partition_data", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def builder():
    return _load_builder()


@pytest.fixture(scope="module")
def rows():
    with open(TABLE) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def meta():
    return json.loads(META.read_text())


def test_the_table_holds_the_validated_counts(rows, builder):
    counts = {}
    for row in rows:
        counts.setdefault(row["gene_name"], {})[row["segment_key"]] = int(row["n_reads"])
    assert counts == builder.EXPECTED_COUNTS


def test_segments_sum_to_the_union_and_percentages_to_one_hundred(rows, builder):
    by_gene = {}
    for row in rows:
        entry = by_gene.setdefault(row["gene_name"], {"union": int(row["n_union"]),
                                                      "n": 0, "pct": 0.0})
        entry["n"] += int(row["n_reads"])
        entry["pct"] += float(row["pct_of_union"])
    assert {g: e["union"] for g, e in by_gene.items()} == builder.EXPECTED_UNION
    for gene, entry in by_gene.items():
        assert entry["n"] == entry["union"], gene
        assert abs(entry["pct"] - 100.0) < 1e-6, gene


def test_gene_order_and_identity(rows, meta, builder):
    order = []
    for row in rows:
        if row["gene_name"] not in order:
            order.append(row["gene_name"])
    assert tuple(order) == builder.GENE_ORDER
    assert meta["gene_order"] == list(builder.GENE_ORDER)
    assert meta["gsm"] == "GSM2100602"
    assert meta["transcripts"]["LRRFIP1"] == "ENST00000308482.14"
    assert [s["key"] for s in meta["segments"]] == \
        sorted({row["segment_key"] for row in rows},
               key=[s["key"] for s in meta["segments"]].index)


def test_the_builder_refuses_a_drifted_count(builder, meta):
    entries = []
    for gene in builder.GENE_ORDER:
        counts = dict(builder.EXPECTED_COUNTS[gene])
        entries.append({"gene_name": gene, "n_union": builder.EXPECTED_UNION[gene],
                        "counts": counts})
    builder.check_expected({"entries": entries})          # the shipped numbers pass
    entries[0]["counts"]["r7_txonly"] += 1
    with pytest.raises(SystemExit):
        builder.check_expected({"entries": entries})
    entries[0]["counts"]["r7_txonly"] -= 1
    entries.reverse()
    with pytest.raises(SystemExit):
        builder.check_expected({"entries": entries})


def test_the_locus_artifact_describes_the_caption_facts():
    locus = json.loads((REPO / "data" / "alignment_fate" / "locus_LRRFIP1.json").read_text())
    assert locus["gene"] == "LRRFIP1" and locus["strand"] == "+"
    assert locus["selected_transcript"] == "ENST00000308482.14"
    assert locus["alternative_transcript"] == "ENST00000244815.9"
    assert locus["n_absent_nt"] == 3619
    assert locus["txome_min_mapq"] == 42
    assert locus["signal"] == "psite"
    import numpy as np
    with np.load(REPO / "data" / "alignment_fate" / "locus_LRRFIP1.npz") as data:
        assert set(data.files) >= {"genomic_position", "genome_cov", "txome_cov",
                                   "sel_exons", "alt_exons", "absent_blocks"}
        n = len(data["genomic_position"])
        assert len(data["genome_cov"]) == n == len(data["txome_cov"])
        assert float(data["genome_cov"].sum()) == locus["counts"]["genome_cov_total"]
