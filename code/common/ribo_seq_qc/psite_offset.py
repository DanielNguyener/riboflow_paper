#!/usr/bin/env python3
"""The P-site offset detector behind every `psite_offset` value this repository ships."""

def ribotish_get_offset(counts_by_pos, defOffset=12, flank=6, default=12):
    """Port of ribotish get_offset() (ribo.py:1111).

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
    """Frame-robust P-site offset: ribotish scan, but the frame comes from downstream
    3-nt periodicity (scored on the 5'-END window, NOT the P-site window — the latter
    destabilises 31/32 nt calls); falls back to `ribotish_get_offset` on thin/unphased signal.
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
