#!/usr/bin/env python3
"""Genome-MULTIMAPPER vs transcriptome concordance."""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import concordance_lib as cl
fc = cl.fc

OUTDIR = fc.output_root() / "read_taxonomy" / "multimap_concordance"

_MISSING_AS = -(10 ** 9)

def read_genome_multi_records(bam_path, target_qnames):
    """qname -> list of (chrom, strand, pos5, n_blocks, blk_min, blk_max, AS) for EVERY
    reported genome alignment (primary + secondary) of each multimapper (NH>1) read whose
    qname is in `target_qnames`. Full BAM pass (no early exit): a read's NH loci sit at
    different coordinates in the sorted file, so all records are only guaranteed after EOF.

    pos5 = 5'-most genomic base (leftmost fwd / rightmost-1 rev), matching
    concordance_lib.read_genome_unique. blocks from get_blocks() (splits on CIGAR N-ops).
    """
    out = defaultdict(list)
    primary_multi = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_supplementary:
            continue
        q = r.query_name
        if q not in target_qnames:
            continue
        try:
            nh = r.get_tag("NH")
        except KeyError:
            nh = None
        if not r.is_secondary and nh is not None and nh > 1:
            primary_multi.add(q)
        if nh is None or nh <= 1:
            continue
        blocks = r.get_blocks()
        if not blocks:
            continue
        strand = "-" if r.is_reverse else "+"
        pos5 = blocks[0][0] if strand == "+" else blocks[-1][1] - 1
        AS = r.get_tag("AS") if r.has_tag("AS") else _MISSING_AS
        out[q].append((r.reference_name, strand, pos5, len(blocks),
                       blocks[0][0], blocks[-1][1], int(AS)))
    bam.close()
    return {q: recs for q, recs in out.items() if q in primary_multi}

def classify_records(records_by_qname, txome_dict, transcript_payload):
    """Per-read match summary DataFrame (index = qname) for reads whose txome primary
    resolves to an APPRIS transcript in `txome_dict`.

    Columns: n_records, max_as, n_at_max, ambiguous, any_match, match_at_best,
    match_only_lower. `match` per genome record == concordance_lib's `concordant` rule
    (within the expected transcript's gene body + strand + splice-junction count + exact
    5' coordinate); records outside the expected gene body are simply non-matching (no
    biotype sub-classification is needed here).
    """
    table = transcript_payload["table"]

    rows = []
    for q, recs in records_by_qname.items():
        tv = txome_dict.get(q)
        if tv is None:
            continue
        tid, tx_pos, tx_len = tv
        for (chrom, strand, pos5, nblk, gmin, gmax, AS) in recs:
            rows.append((q, tid, tx_pos, tx_len, chrom, strand, pos5, nblk, gmin, gmax, AS))

    cols = ["qname", "tid", "tx_pos", "tx_len", "g_chrom", "g_strand",
            "g_pos5", "g_nblocks", "g_min", "g_max", "AS"]
    if not rows:
        return pd.DataFrame(columns=["n_records", "max_as", "n_at_max", "ambiguous",
                                     "any_match", "match_at_best", "match_only_lower"])
    df = pd.DataFrame(rows, columns=cols)

    match = np.zeros(len(df), dtype=bool)
    for tid, idx in df.groupby("tid", sort=False).groups.items():
        t = table.get(tid)
        if t is None:
            continue
        pos = df.index.get_indexer(idx)
        sub = df.loc[idx]

        cum_start = t["cum_start"]; g_start = t["g_start"]; g_end = t["g_end"]
        n_exons = len(cum_start)
        tx_p = sub["tx_pos"].to_numpy(); tx_l = sub["tx_len"].to_numpy()
        idx_start = np.clip(np.searchsorted(cum_start, tx_p, side="right") - 1, 0, n_exons - 1)
        idx_end = np.clip(np.searchsorted(cum_start, tx_p + tx_l - 1, side="right") - 1, 0, n_exons - 1)
        expected_njunc = idx_end - idx_start
        offset = tx_p - cum_start[idx_start]
        if t["strand"] == "+":
            expected_pos5 = g_start[idx_start] + offset
        else:
            expected_pos5 = g_end[idx_start] - 1 - offset

        obs_chrom = sub["g_chrom"].to_numpy()
        obs_strand = sub["g_strand"].to_numpy()
        obs_pos5 = sub["g_pos5"].to_numpy(dtype=np.int64)
        obs_min = sub["g_min"].to_numpy(); obs_max = sub["g_max"].to_numpy()
        obs_njunc = sub["g_nblocks"].to_numpy() - 1

        same_chrom = obs_chrom == t["chrom"]
        within_body = same_chrom & (obs_min < t["body_end"]) & (obs_max > t["body_start"])
        strand_ok = obs_strand == t["strand"]
        splice_ok = obs_njunc == expected_njunc
        coord_ok = np.abs(obs_pos5 - expected_pos5.astype(np.int64)) <= cl.COORD_TOL
        match[pos] = within_body & strand_ok & splice_ok & coord_ok

    df["match"] = match
    df["max_as"] = df.groupby("qname", sort=False)["AS"].transform("max")
    df["at_max"] = df["AS"] == df["max_as"]
    df["match_best"] = df["match"] & df["at_max"]

    summ = df.groupby("qname", sort=False).agg(
        n_records=("AS", "size"),
        max_as=("AS", "max"),
        n_at_max=("at_max", "sum"),
        any_match=("match", "any"),
        match_at_best=("match_best", "any"),
    )
    summ["n_at_max"] = summ["n_at_max"].astype(int)
    summ["ambiguous"] = summ["n_at_max"] >= 2
    summ["match_only_lower"] = summ["any_match"] & ~summ["match_at_best"]
    return summ
