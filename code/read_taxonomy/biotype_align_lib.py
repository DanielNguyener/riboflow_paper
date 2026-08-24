#!/usr/bin/env python3
"""GENCODE gene_type of genome-MULTIMAPPER ALIGNMENTS."""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import biotype_lib as bl
import mm_concordance_lib as mm
import taxonomy_lib as tl
cl, fc = bl.cl, bl.fc

OUTDIR = fc.output_root() / "read_taxonomy" / "multimap_biotype"

def population_records(sample, log=print):
    """qname -> [genome-locus records] for the dark-green population (gM_tU + gM_tM).

    Values are the (chrom, strand, pos5, n_blocks, blk_min, blk_max, AS) tuples from
    `mm_concordance_lib.read_genome_multi_records`.
    """
    log(f"[{sample}] reading txome BAM (present qname set)...")
    t_all, _ = tl.status_sets(fc.txome_bam(sample), "txome")
    log(f"[{sample}] n_txome_present={len(t_all):,}; enumerating genome multimapper loci...")
    records_by_qname = mm.read_genome_multi_records(fc.genome_bam(sample), t_all)
    return records_by_qname
