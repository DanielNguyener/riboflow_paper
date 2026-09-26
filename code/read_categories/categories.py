#!/usr/bin/env python3
"""The five read categories of Figures 4-6, defined once.

Each library's read-ID union is partitioned into shared, genome-only and
transcriptome-only reads; shared and genome-only reads are further split by whether they
map uniquely to the genome (`NH == 1`). The manuscript's names, always in this order:

    SH-U  shared genome-unique          SH-M  shared genome-multimapped
    GO-U  genome-only unique            GO-M  genome-only multimapped
    TO    transcriptome-only

Two transcriptome-PRESENCE rules exist, and each figure names the one it uses. Both take
a read from the post-dedup transcriptome BAM (RiboFlow_v2 filtered it at MAPQ >= 10),
never the MAPQ >= 42 rule of S1 Fig / Figures 2-3:

  * ANY PRIMARY (Figure 4's library-wide categories, the 5B locus): the read has any
    primary transcriptome alignment.
  * SELECTED TRANSCRIPT (Figure 5A, Figure 6): the read's primary transcriptome
    alignment resolves to a selected reference transcript (APPRIS principal).

The two agree on HeLa's total (7,686,204 reads) because the transcriptome reference
contains only the selected transcripts; they are still different rules, kept separate.
"""
from __future__ import annotations

#: Canonical keys, in the manuscript's order.
KEYS = ("sh_u", "sh_m", "go_u", "go_m", "to")

#: The abbreviations the figures print.
ABBR = {"sh_u": "SH-U", "sh_m": "SH-M", "go_u": "GO-U", "go_m": "GO-M", "to": "TO"}

#: The captions' long names.
LONG = {"sh_u": "shared genome-unique",
        "sh_m": "shared genome-multimapped",
        "go_u": "genome-only unique",
        "go_m": "genome-only multimapped",
        "to": "transcriptome-only"}

#: One colour per category, shared by every panel that draws them.
COLOR = {"sh_u": "#a6d96a", "sh_m": "#1a7d1a",
         "go_u": "#7fb9da", "go_m": "#0d57a1", "to": "#cc3d3d"}

#: (abbreviation, colour) in order -- the key most panels draw.
KEY = tuple((ABBR[k], COLOR[k]) for k in KEYS)


def n_col(key):
    """Count-column name for a category, e.g. 'n_sh_u'."""
    return "n_%s" % key


def pct_col(key):
    """Percentage-column name for a category, e.g. 'pct_sh_u'."""
    return "pct_%s" % key


#: Sentinel for an alignment with no AS tag: far below any real score, so an unscored
#: alignment can never tie with a scored one.
MISSING_AS = -(10 ** 9)


def qualifies(is_secondary, score, primary_score):
    """Can this genome alignment place its read at a gene? (tied best-scoring rule)

    The primary always; a secondary only when it and the primary both carry an AS and the
    two are equal. A read with no primary record, or an alignment with no AS, is never
    tied: that is reported, not assumed. `read_state.ReadState` applies the same rule
    vectorised over the whole store.
    """
    if not is_secondary:
        return True
    return (primary_score is not None and primary_score != MISSING_AS
            and score != MISSING_AS and score == primary_score)
