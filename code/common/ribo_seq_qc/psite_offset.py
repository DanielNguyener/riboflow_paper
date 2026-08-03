#!/usr/bin/env python3
"""The P-site offset detector behind every `psite_offset` value this repository ships."""

def ribotish_get_offset(counts_by_pos, defOffset=12, flank=6, default=12):
    """
    Port of ribotish get_offset() (ribo.py:1111).

    1. Locate TIS peak in upstream region → derive frame (peak_pos % 3).
    2. Compute noise threshold from far-upstream vs CDS-body in-frame counts.
    3. Scan ±flank nt around defOffset, restricted to the TIS frame; return
       the first in-frame position that exceeds the threshold.

    Returns a positive offset (nt from 5' end) or `default` on failure.
    """
    CODON = 3

    upstream = {p: c for p, c in counts_by_pos.items() if p < 0}
    if not upstream or max(upstream.values()) == 0:
        return default

    tis_pos = max(upstream, key=upstream.get)
    frame   = tis_pos % 3

    threshold_pos = -defOffset - flank

    a0 = [counts_by_pos.get(p, 0) for p in sorted(counts_by_pos)
          if p % CODON == frame and p <= threshold_pos]
    a1 = [counts_by_pos.get(p, 0) for p in sorted(counts_by_pos)
          if p % CODON == frame and p >= threshold_pos]
    ai = [counts_by_pos.get(p, 0) for p in sorted(counts_by_pos)
          if p % CODON == frame and threshold_pos < p < -defOffset + flank]

    if not a0 or not a1 or not ai:
        return default

    a0m = max(a0)
    a1m = max(a1)
    th  = a0m + int((a1m - a0m) / 6.0)

    for p in range(-defOffset - flank + 1, -defOffset + flank):
        if p % CODON != frame:
            continue
        if counts_by_pos.get(p, 0) > th:
            return -p

    return default

def get_offset_periodicity(counts_by_pos, defOffset=12, flank=6, default=12,
                           win_codons=10, min_down=200, dom_frac=0.40):
    """Frame-robust P-site offset: like ribotish, but the reading FRAME comes from
    DOWNSTREAM 3-nt periodicity instead of the single tallest start-peak position.

    Rationale: `ribotish_get_offset` sets `frame = argmax_upstream % 3`, which flips on
    a near-tie when the start metagene is bimodal (e.g. HEK293T 26 nt: positions -9 and
    -8 within ~1%, different frames -> +9 vs +11). The downstream 5'-end pileup is phased
    to a single residue class (P-sites are frame-0, so 5'-ends sit at p ≡ -offset mod 3),
    and the dominant residue over the first `win_codons` codons of 5'-end positions is a
    far more stable frame signal than one start position. We take that dominant residue as
    the frame, then run the IDENTICAL ribotish magnitude/threshold scan constrained to it
    — so the return convention (canonical ~+11/+12) and the shared genome/transcriptome
    implementation are preserved. Uses only `counts_by_pos` — no extra I/O.

    The frame is scored on the 5'-END window (p in [0,win)), NOT the P-site window
    (p+offset): the latter (tried and rejected) pulls start-proximal mixed-frame reads
    into the window for long read lengths and destabilises the call on 31/32 nt reads.

    Falls back to plain `ribotish_get_offset` when the downstream signal is too thin
    (< `min_down` reads) or not clearly phased (dominant-frame fraction < `dom_frac`,
    just above the 1/3 no-periodicity baseline), so this is a strict superset of the
    old behavior on clean/ambiguous-free lengths.
    """
    CODON = 3
    WIN = win_codons * CODON

    upstream = {p: c for p, c in counts_by_pos.items() if p < 0}
    if not upstream or max(upstream.values()) == 0:
        return default

    # downstream 5'-end residue histogram (P-sites land in CDS body for these)
    mass = [0, 0, 0]
    for p, c in counts_by_pos.items():
        if 0 <= p < WIN:
            mass[p % CODON] += c
    total = sum(mass)
    if total < min_down or max(mass) / total < dom_frac:
        return ribotish_get_offset(counts_by_pos, defOffset, flank, default)
    frame = max(range(CODON), key=lambda r: mass[r])

    threshold_pos = -defOffset - flank
    a0 = [counts_by_pos.get(p, 0) for p in sorted(counts_by_pos)
          if p % CODON == frame and p <= threshold_pos]
    a1 = [counts_by_pos.get(p, 0) for p in sorted(counts_by_pos)
          if p % CODON == frame and p >= threshold_pos]
    ai = [counts_by_pos.get(p, 0) for p in sorted(counts_by_pos)
          if p % CODON == frame and threshold_pos < p < -defOffset + flank]
    if not a0 or not a1 or not ai:
        return default
    th = max(a0) + int((max(a1) - max(a0)) / 6.0)
    for p in range(-defOffset - flank + 1, -defOffset + flank):
        if p % CODON != frame:
            continue
        if counts_by_pos.get(p, 0) > th:
            return -p
    return default
