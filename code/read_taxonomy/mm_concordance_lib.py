#!/usr/bin/env python3
"""Genome-MULTIMAPPER vs transcriptome concordance."""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

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
    """qname -> [(chrom, strand, pos5, n_blocks, blk_min, blk_max, AS)] for every reported
    genome alignment of each NH>1 read in `target_qnames`.

    Full BAM pass, no early exit: a read's NH loci are scattered through the sorted file.
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
