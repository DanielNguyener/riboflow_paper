#!/usr/bin/env python3
"""The multi-mapping read taxonomy."""
from __future__ import annotations

import sys
from pathlib import Path


_HERE = Path(__file__).resolve().parent
_COMMON = _HERE.parent / "common"
for _entry in (str(_HERE), str(_COMMON), str(_COMMON / "ribo_seq_qc")):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
import bam_inputs as fc

GENOME_STATES = ("unique", "multi", "absent")
TXOME_STATES = ("present", "absent")
#: five cells of the read ID taxonomy  genome status x transcriptome presence minus
#: the empty absent absent cell  transcriptome presence is a primary alignment in the
#: post dedup BAM  RiboFlow_v2 already filtered that BAM at MAPQ >= 10 so no further
#: threshold here  the MAPQ >= 42 rule belongs to the coverage and TE analyses
CELLS = tuple((g, t) for g in GENOME_STATES for t in TXOME_STATES
              if not (g == "absent" and t == "absent"))

#: one taxonomy cell -> its manuscript category key  categories.KEYS is the one definition
CELL_KEY = {("unique", "present"): "sh_u", ("multi", "present"): "sh_m",
            ("unique", "absent"): "go_u", ("multi", "absent"): "go_m",
            ("absent", "present"): "to"}

def cell_key(genome_status, txome_status):
    """`sh_u`, `go_m`, ...: the column stem of one taxonomy cell."""
    return CELL_KEY[(genome_status, txome_status)]
