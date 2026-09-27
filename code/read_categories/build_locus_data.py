#!/usr/bin/env python3
"""BAMs + GTF + APPRIS + QC tables -> the locus artifact (`locus_<GENE>.{npz,json}`).

Per-base P-site coverage of both routes over the merged exons of the selected isoform and
the best-supported alternative isoform (the one carrying the most genome-only unique reads
on non-selected sequence). Writes <output>.npz (coverage vectors + exon blocks) and
<output>.json (metadata). Introns are drawn as a constant 90-unit gap. The BAM template
and cigar-aware P-site are re-implemented here, identically to `common/` and `coverage/`.
The transcriptome half counts every primary alignment on the selected transcript: the
post-dedup BAM is already RiboFlow_v2's MAPQ >= 10 set, and no stricter cut is applied
(the same rule as Figures 4-6).

Each route's coverage is also split into the read populations of the Figure 5A partition
(`panels/plot_gene_categories.ROUTE7_KEY`), so the two panels share one colour
vocabulary. "Shared" is read-level presence in both BAMs, exactly as in 5A:
  genome half   = shared_unique      (genome NH==1, read present in the transcriptome BAM)
                + shared_multi       (genome NH>1, present in the transcriptome BAM)
                + genome_only        (genome NH==1, absent from the transcriptome BAM)
                + genome_only_multi  (genome NH>1, absent from the transcriptome BAM)
  txome half    = shared_unique (txome primary on the selected transcript, a top-score
                                 genome placement at the locus, NH==1)
                + shared_multi  (..., a top-score genome placement at the locus, NH>1)
                + txome_only    (..., no top-score genome placement at the locus)
A top-score placement is the read's genome primary or a secondary whose AS equals the
primary's (any read length): 5A's gene-membership rule. NH is the primary's, as in Figures
3-4. A lower-scoring secondary at the locus places nothing; neither does a secondary whose
read has no primary record or no AS. Genome-track coverage places each member read AT MOST
ONCE, from its qualifying (top-score, at-locus) placements only: identical P-sites collapse
to one count, disagreeing P-sites omit the read from position-resolved coverage (counted
per layer in the JSON), and no placement is ever picked arbitrarily.

Published locus: LRRFIP1, chr2:237,627,586-237,781,643 (+), selected ENST00000308482.14,
alternative ENST00000244815.9 (3,619 nt absent from the selected reference).
Run with `python` (3.9).
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import sys

import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))
import inputs as paths  # noqa: E402
from inputs import die  # noqa: E402
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "coverage"))
import psite_placement  # noqa: E402
from intervals import merge, subtract  # noqa: E402

#: txome route uniqueness  bowtie2 emits no NH tag  project wide rule
#: fixed width in plotted units of dashed intron connector  recorded in artifact
INTRON_GAP = 90.0


# ── QC tables ────────────────────────────────────────────────────────────────

def read_window_and_offsets(qc_csv, sample):
    """The sample's selected read lengths and their P-site offsets, from the QC table.

    `psite_placement.load_offsets` is the one reader of the QC master tables."""
    offsets = psite_placement.load_offsets(qc_csv, sample)
    return set(offsets), offsets


def psite_reference_position(read, offset):
    """The read's P-site as a reference coordinate. The one P-site rule: `cigar_aware`.

    Offset walked along the READ; returns None when it lands in an insertion or off the
    alignment — such reads are counted and dropped, never approximated.
    """
    target = (read.query_length - 1 - offset) if read.is_reverse else offset
    if target < 0 or target >= (read.query_length or 0):
        return None
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if query_pos == target:
            return ref_pos
    return None


# ── annotation ───────────────────────────────────────────────────────────────

def selected_transcript(appris_path, gene_name):
    """(transcript id, full header, length) of the gene's transcriptome-reference entry."""
    hits = []
    with open(appris_path) as handle:
        for line in handle:
            header = line.split("\t")[0]
            fields = header.split("|")
            if len(fields) > 6 and fields[5] == gene_name:
                hits.append((fields[0], header, int(fields[6])))
    if not hits:
        die("%s has no transcript in the transcriptome reference" % gene_name)
    if len(hits) > 1:
        die("%s maps to %d reference transcripts: %s"
            % (gene_name, len(hits), [h[0] for h in hits]))
    return hits[0]


def gene_transcripts(gtf_path, gene_name):
    """{transcript_id: [(start, end), ...]} (0-based half-open) for one gene, + chrom, strand."""
    exons = collections.defaultdict(list)
    chrom = strand = None
    needle = 'gene_name "%s"' % gene_name
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line[0] == "#" or needle not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "exon":
                continue
            attributes = fields[8]
            if needle not in attributes:
                continue
            tid = attributes.split('transcript_id "', 1)[1].split('"', 1)[0]
            exons[tid].append((int(fields[3]) - 1, int(fields[4])))
            chrom, strand = fields[0], fields[6]
    if not exons:
        die("%s has no exon in %s" % (gene_name, gtf_path))
    return {t: sorted(v) for t, v in exons.items()}, chrom, strand


def exonic_bases(blocks, strand):
    """Every exonic base of the merged blocks, in 5'->3' plot order (the SplicedAxis order)."""
    ordered = list(blocks) if strand == "+" else list(reversed(blocks))
    out = []
    for start, end in ordered:
        out.append(np.arange(start, end) if strand == "+"
                   else np.arange(end - 1, start - 1, -1))
    return np.concatenate(out) if out else np.array([], dtype=int)


# ── reads ────────────────────────────────────────────────────────────────────

def locus_ribo_reads(bam_path, chrom, start, end, lengths, offsets):
    """One region fetch (+ one genome pass) -> (members, status, audit).

    `status` holds the primary NH of every read with a TOP-SCORE placement overlapping the
    locus -- its primary, or a secondary tied with the primary's AS -- with no length
    filter: Figure 5A's gene-membership rule
    (`gene_read_partition_lib.resolve_genome_side`), reused to label the transcriptome
    track's reads. A secondary-only read needs its primary's AS, which may lie anywhere, so
    those reads get one full genome pass.

    `members` restricts `status` to reads in the read-length window and carries, per read,
    every QUALIFYING placement at the locus -- the primary when it lies here, plus every
    secondary whose AS ties the primary's. A lower-scoring secondary places nothing, and no
    placement is ever picked arbitrarily: {qname: {"nh": int,
    "alignments": [(blocks, psite or None), ...]}}. The P-site is cigar_aware per
    placement; the caller collapses a read's placements to at most one count.
    """
    primary_here = {}                                  # qname to score blocks psite
    status = {}
    secondaries = collections.defaultdict(list)        # qname to list of score qlen blocks psite
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        if not bam.has_index():
            die("%s has no index; the locus view is a region fetch" % bam_path)
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_supplementary:
                continue
            score = int(read.get_tag("AS")) if read.has_tag("AS") else None
            qlen = read.query_length or read.infer_query_length() or 0
            blocks = read.get_blocks()
            psite = (psite_reference_position(read, offsets[qlen])
                     if qlen in lengths and blocks else None)
            if read.is_secondary:
                secondaries[read.query_name].append((score, qlen, blocks, psite))
                continue
            try:
                nh = int(read.get_tag("NH"))
            except KeyError:
                die("%s has an alignment with no NH tag" % bam_path)
            if read.query_name in status:
                die("%s: read %s has two primary alignments" % (bam_path, read.query_name))
            status[read.query_name] = nh
            if blocks:
                primary_here[read.query_name] = (score, qlen, blocks, psite)
    finally:
        bam.close()

    audit = {"joined_by_tied_secondary_only": 0, "secondary_only_no_primary": 0,
             "secondary_only_lower_score": 0, "missing_as": 0}
    primary_score_of = {q: rec[0] for q, rec in primary_here.items()}
    pending = set(secondaries) - set(status)
    primaries = genome_primaries(bam_path, pending) if pending else {}
    for qname in sorted(pending):
        scores = [score for score, _l, _b, _p in secondaries[qname]]
        found = primaries.get(qname)
        if found is None:
            audit["secondary_only_no_primary"] += 1
            continue
        nh, primary_score = found
        missing = primary_score is None or any(score is None for score in scores)
        if missing:
            audit["missing_as"] += 1
        if primary_score is not None and primary_score in scores:
            status[qname] = nh
            primary_score_of[qname] = primary_score
            audit["joined_by_tied_secondary_only"] += 1
        elif not missing:
            audit["secondary_only_lower_score"] += 1

    members = {}
    for qname, nh in status.items():
        qualifying = []
        if qname in primary_here:
            _score, qlen, blocks, psite = primary_here[qname]
            if qlen in lengths:
                qualifying.append((blocks, psite))
        top = primary_score_of.get(qname)
        for score, qlen, blocks, psite in secondaries.get(qname, ()):
            if qlen in lengths and blocks and score is not None and score == top:
                qualifying.append((blocks, psite))
        if qualifying:
            members[qname] = {"nh": nh, "alignments": qualifying}
    return members, status, audit


def genome_primaries(bam_path, qnames):
    """{qname: (NH, AS or None)} of the primary genome record of each read in `qnames`."""
    found = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_name in qnames:
                found[read.query_name] = (int(read.get_tag("NH")),
                                          int(read.get_tag("AS")) if read.has_tag("AS")
                                          else None)
    finally:
        bam.close()
    return found


def txome_present(bam_path, wanted):
    """Which of `wanted` appear anywhere in the transcriptome BAM (a full pass, on purpose)."""
    found = set()
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_name in wanted:
                found.add(read.query_name)
    finally:
        bam.close()
    return found


def txome_coverage(bam_path, reference, lengths, offsets, tx2genome, signal):
    """Per-read genomic positions from the transcriptome route's primary alignments on
    `reference` (every one: the post-dedup BAM is already RiboFlow_v2's MAPQ >= 10 set).

    Returns ({qname: [genomic position, ...]}, n_reads, dropped): one position per read for
    the P-site signal, every covered base for the footprint signal. Kept per read so the
    track can be split by each read's genome status.
    """
    hits = collections.defaultdict(list)
    n_reads = 0
    dropped = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        if not bam.has_index():
            die("%s has no index; the transcriptome track is a region fetch" % bam_path)
        if reference not in set(bam.references):
            die("%s is not a reference in %s" % (reference, bam_path))
        for read in bam.fetch(reference):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_length not in lengths:
                continue
            n_reads += 1
            if signal == "psite":
                psite = psite_reference_position(read, offsets[read.query_length])
                if psite is None or not (0 <= psite < len(tx2genome)):
                    dropped += 1
                    continue
                hits[read.query_name].append(int(tx2genome[psite]))
            else:
                lo = max(read.reference_start, 0)
                hi = min(read.reference_end, len(tx2genome))
                hits[read.query_name].extend(int(pos) for pos in tx2genome[lo:hi])
    finally:
        bam.close()
    return dict(hits), n_reads, dropped


def depth_over(hits, positions):
    """Depth at `positions` from {qname: [genomic position, ...]}."""
    depth = collections.Counter()
    for placed in hits.values():
        for pos in placed:
            depth[pos] += 1
    return np.array([depth.get(int(p), 0) for p in positions], dtype=float)


#: figure 5A populations each track splits into  stacking order  nearest baseline first
#: keys match panels/plot_gene_categories.ROUTE7_KEY wording
#: SH-U SH-M GO-U GO-M on genome track  SH-U SH-M TO on transcriptome track
GENOME_LAYERS = ("shared_unique", "shared_multi", "genome_only", "genome_only_multi")
TXOME_LAYERS = ("shared_unique", "shared_multi", "txome_only")


def txome_read_layer(qname, genome_status):
    """Which 5A population a transcriptome-track read belongs to, from its primary NH when it
    has a top-score genome placement at the locus (none: transcriptome-only at the gene)."""
    nh = genome_status.get(qname)
    if nh is None:
        return "txome_only"
    return "shared_unique" if nh == 1 else "shared_multi"


def member_positions(entry, signal):
    """The one set of genomic positions a member read contributes, or its exclusion.

    Every qualifying (top-score, at-locus) placement of the read votes; identical votes
    collapse to one. Returns None when the placements disagree (POSITION-AMBIGUOUS: the
    read is omitted from position-resolved coverage, never assigned arbitrarily), [] when
    no placement resolves a P-site, else the positions counted exactly once.
    """
    if signal == "psite":
        votes = {int(psite) for _blocks, psite in entry["alignments"] if psite is not None}
        if len(votes) > 1:
            return None
        return [votes.pop()] if votes else []
    votes = {tuple(sorted({pos for start, end in blocks for pos in range(start, end)}))
             for blocks, _psite in entry["alignments"]}
    if len(votes) > 1:
        return None
    return list(votes.pop()) if votes else []


def layer_coverage(members, layer_qnames, positions, signal):
    """(depth vector, n placed, n position-ambiguous, n P-site-unresolved) for one layer.

    One read contributes at most once: `member_positions` collapses tied placements and
    excludes disagreeing ones.
    """
    depth = collections.Counter()
    placed = ambiguous = unresolved = 0
    for qname in layer_qnames:
        positions_of_read = member_positions(members[qname], signal)
        if positions_of_read is None:
            ambiguous += 1
            continue
        if not positions_of_read:
            unresolved += 1
            continue
        placed += 1
        for pos in positions_of_read:
            depth[pos] += 1
    vector = np.array([depth.get(int(p), 0) for p in positions], dtype=float)
    return vector, placed, ambiguous, unresolved


# ── driver ───────────────────────────────────────────────────────────────────

def build(gene, sample, inputs, qc_genome, qc_txome, signal):
    lengths, offsets = read_window_and_offsets(qc_genome, sample)
    t_lengths, t_offsets = read_window_and_offsets(qc_txome, sample)
    print("[locus] %s signal=%s   genome lengths %s offsets %s"
          % (sample, signal, sorted(lengths), sorted(set(offsets.values()))))
    print("[locus]          transcriptome lengths %s offsets %s"
          % (sorted(t_lengths), sorted(set(t_offsets.values()))))

    sel_tid, sel_header, sel_len = selected_transcript(inputs["appris"], gene)
    transcripts, chrom, strand = gene_transcripts(inputs["gtf"], gene)
    if sel_tid not in transcripts:
        die("selected transcript %s is not in the GTF for %s" % (sel_tid, gene))
    sel_exons = merge(transcripts[sel_tid])
    span_start = min(s for exons in transcripts.values() for s, _e in exons)
    span_end = max(e for exons in transcripts.values() for _s, e in exons)
    print("[locus] %s %s:%d-%d (%s), %d annotated transcripts; selected %s (%d exons)"
          % (gene, chrom, span_start, span_end, strand, len(transcripts), sel_tid,
             len(sel_exons)))

    tx2genome = exonic_bases(sel_exons, strand)
    if len(tx2genome) != sel_len:
        die("%s: GTF exons give a spliced length of %d, but the transcriptome reference "
            "declares %d -- the projection would be wrong" % (sel_tid, len(tx2genome), sel_len))

    members, genome_status, membership_audit = locus_ribo_reads(
        inputs["ribo_genome"], chrom, span_start, span_end, lengths, offsets)
    genome_unique = {q for q, entry in members.items() if entry["nh"] == 1}
    genome_multi = set(members) - genome_unique
    print("[locus] %d member reads at the locus in the window (%d unique, %d multimapping;"
          " %d read ids with a top-score placement at the locus at any length)"
          % (len(members), len(genome_unique), len(genome_multi), len(genome_status)))
    print("[locus]   membership: %s" % ", ".join("%s %d" % kv for kv in membership_audit.items()))

    txome_hits, n_txome, t_dropped = txome_coverage(
        inputs["ribo_txome"], sel_header, t_lengths, t_offsets, tx2genome, signal)
    print("[locus] %d transcriptome-route reads on %s, projected onto the genome axis"
          % (n_txome, sel_tid))
    if t_dropped:
        print("[locus]   %d transcriptome reads dropped for the same reason" % t_dropped)

    print("[locus] one transcriptome pass to find which genome reads are genome-only ...")
    present = txome_present(inputs["ribo_txome"], set(members))
    genome_only = {q for q in genome_unique if q not in present}
    print("[locus] genome-only uniquely mapping: %d of %d unique"
          % (len(genome_only), len(genome_unique)))

    # alternative isoform  most genome only unique reads on non selected sequence
    scores = {}
    for tid, exons in transcripts.items():
        if tid == sel_tid:
            continue
        extra = subtract(merge(exons), sel_exons)
        if not extra:
            continue
        n = 0
        for qname in genome_only:
            blocks, _psite = members[qname]["alignments"][0]   # NH==1 so the primary
            if any(b_s < e_e and e_s < b_e
                   for b_s, b_e in blocks for e_s, e_e in extra):
                n += 1
        scores[tid] = n
    if not scores or max(scores.values()) == 0:
        die("no annotated transcript of %s carries genome-only reads outside the "
            "selected isoform; this gene is not an example of the effect" % gene)
    alt_tid = max(scores, key=lambda t: (scores[t], len(transcripts[t])))
    print("[locus] alternative isoform support (genome-only reads on non-selected "
          "sequence):")
    for tid, n in sorted(scores.items(), key=lambda kv: -kv[1])[:5]:
        print("          %-22s %6d%s" % (tid, n, "   <- chosen" if tid == alt_tid else ""))

    alt_exons = merge(transcripts[alt_tid])
    absent_blocks = subtract(alt_exons, sel_exons)
    print("[locus] %s adds %d nt of exonic sequence absent from %s"
          % (alt_tid, sum(e - s for s, e in absent_blocks), sel_tid))

    # per base vectors exactly as panel draws them  union exonic bases in 5' to 3' order
    union = merge(list(sel_exons) + list(alt_exons))
    gs = exonic_bases(union, strand)
    txome_cov = depth_over(txome_hits, gs)

    # tracks split into figure 5A populations  genome track unique vs multi is primary NH
    # shared vs genome only is transcriptome presence  multimapper positions come only
    # from top score placements at locus  one count per read  disagreeing placements
    # omitted see member_positions  transcriptome track uses NH of each read with a top
    # score genome placement at locus  from the status fetch
    genome_members = {"shared_unique": genome_unique - genome_only,
                      "shared_multi": genome_multi & present,
                      "genome_only": genome_only,
                      "genome_only_multi": genome_multi - present}
    genome_layers, genome_placed, genome_ambiguous, genome_unresolved = {}, {}, {}, {}
    for name in GENOME_LAYERS:
        (genome_layers[name], genome_placed[name],
         genome_ambiguous[name], genome_unresolved[name]) = layer_coverage(
            members, genome_members[name], gs, signal)
    genome_cov = sum(genome_layers.values())
    g_dropped = sum(genome_unresolved.values())
    txome_members = {name: set() for name in TXOME_LAYERS}
    for qname in txome_hits:
        txome_members[txome_read_layer(qname, genome_status)].add(qname)
    txome_layers = {name: depth_over({q: txome_hits[q] for q in members_}, gs)
                    for name, members_ in txome_members.items()}
    if not np.array_equal(sum(genome_layers.values()), genome_cov):
        die("the genome layers do not sum to the genome track")
    if not np.array_equal(sum(txome_layers.values()), txome_cov):
        die("the transcriptome layers do not sum to the transcriptome track")
    if signal == "psite":
        # no read counted twice within a track  one P-site count per placed read
        if float(genome_cov.sum()) > sum(genome_placed.values()):
            die("the genome track carries more P-site counts than placed reads")
        if float(txome_cov.sum()) != float(sum(len(v) for v in txome_hits.values())):
            die("the transcriptome track carries more P-site counts than reads")
        if any(len(v) != 1 for v in txome_hits.values()):
            die("a transcriptome read carries more than one P-site")
    for name in GENOME_LAYERS:
        print("[locus] genome track %-18s %6d reads, %d placed, %d position-ambiguous "
              "omitted, %d without a resolvable P-site"
              % (name, len(genome_members[name]), genome_placed[name],
                 genome_ambiguous[name], genome_unresolved[name]))
    print("[locus] transcriptome track by 5A population: %s"
          % ", ".join("%s %d reads" % (n, len(txome_members[n])) for n in TXOME_LAYERS))

    arrays = {
        "genomic_position": gs.astype(np.int64),
        "genome_cov": genome_cov.astype(np.float64),
        "txome_cov": txome_cov.astype(np.float64),
        "sel_exons": np.array(sel_exons, dtype=np.int64).reshape(-1, 2),
        "alt_exons": np.array(alt_exons, dtype=np.int64).reshape(-1, 2),
        "absent_blocks": np.array(absent_blocks, dtype=np.int64).reshape(-1, 2),
    }
    for name in GENOME_LAYERS:
        arrays["genome_cov_%s" % name] = genome_layers[name].astype(np.float64)
    for name in TXOME_LAYERS:
        arrays["txome_cov_%s" % name] = txome_layers[name].astype(np.float64)
    meta = {
        "gene": gene, "chrom": chrom, "strand": strand,
        "locus": {"start": int(union[0][0]), "end": int(union[-1][1]),
                  "label": "%s:%s-%s (%s)" % (chrom, format(union[0][0], ","),
                                              format(union[-1][1], ","), strand)},
        "selected_transcript": sel_tid, "selected_length": sel_len,
        "alternative_transcript": alt_tid,
        "alternative_gene": gene,
        "alternative_chosen_by": "forced" if alt_transcript else
                                 "most genome-only unique reads on non-selected sequence",
        "alternative_support_top5": dict(sorted(scores.items(), key=lambda kv: -kv[1])[:5]),
        "n_absent_nt": int(sum(e - s for s, e in absent_blocks)),
        "sample": sample,
        "signal": signal,
        "intron_gap_plot_units": INTRON_GAP,
        "txome_reads": "every primary alignment in the post-dedup BAM (RiboFlow_v2 MAPQ>=10)",
        "genome_window": {"lengths": sorted(lengths),
                          "offsets": {str(k): v for k, v in sorted(offsets.items())}},
        "txome_window": {"lengths": sorted(t_lengths),
                         "offsets": {str(k): v for k, v in sorted(t_offsets.items())}},
        "counts": {"genome_member_reads_in_window": len(members),
                   "status_membership_audit": membership_audit,
                   "genome_unique": len(genome_unique),
                   "genome_multi": len(genome_multi),
                   "genome_only_unique": len(genome_only),
                   "genome_psite_dropped": g_dropped,
                   "txome_on_selected": n_txome,
                   "txome_psite_dropped": t_dropped,
                   "genome_cov_total": float(genome_cov.sum()),
                   "txome_cov_total": float(txome_cov.sum()),
                   "genome_track_reads_by_layer":
                       {n: len(genome_members[n]) for n in GENOME_LAYERS},
                   "genome_track_placed_by_layer":
                       {n: genome_placed[n] for n in GENOME_LAYERS},
                   "genome_track_position_ambiguous_by_layer":
                       {n: genome_ambiguous[n] for n in GENOME_LAYERS},
                   "genome_track_psite_unresolved_by_layer":
                       {n: genome_unresolved[n] for n in GENOME_LAYERS},
                   "txome_track_reads_by_layer":
                       {n: len(txome_members[n]) for n in TXOME_LAYERS},
                   "genome_track_cov_by_layer":
                       {n: float(genome_layers[n].sum()) for n in GENOME_LAYERS},
                   "txome_track_cov_by_layer":
                       {n: float(txome_layers[n].sum()) for n in TXOME_LAYERS}},
        "coverage_semantics": {
            "genome_cov": "genome route, reads with a top-score placement at the locus in "
                          "the read-length window, at most one count per read at its "
                          "P-site (cigar_aware). Only placements tied for the read's "
                          "highest genomic score qualify; identical P-sites from several "
                          "qualifying placements collapse to one count; disagreeing "
                          "P-sites omit the read from position-resolved coverage "
                          "(counted in genome_track_position_ambiguous_by_layer), never "
                          "an arbitrary pick",
            "txome_cov": "transcriptome route, every primary on the selected transcript "
                         "(post-dedup BAM, RiboFlow_v2 MAPQ>=10), P-site projected through "
                         "the transcript->genome map",
            "order": "5'->3' over the merged exons of both models (SplicedAxis order)",
            "layers": {
                "genome_cov_shared_unique": "genome_cov, NH==1, reads present in the "
                                            "transcriptome BAM (partition: shared, unique)",
                "genome_cov_shared_multi": "genome_cov, NH>1, reads present in the "
                                           "transcriptome BAM (partition: shared, "
                                           "multimapping)",
                "genome_cov_genome_only": "genome_cov, NH==1, reads absent from the "
                                          "transcriptome BAM (partition: genome-only, "
                                          "unique)",
                "genome_cov_genome_only_multi": "genome_cov, NH>1, reads absent from the "
                                                "transcriptome BAM (partition: "
                                                "genome-only, multimapping)",
                "txome_cov_shared_unique": "txome_cov, reads with a top-score genome "
                                           "placement at the locus and NH==1 (partition: "
                                           "shared, unique)",
                "txome_cov_shared_multi": "txome_cov, reads with a top-score genome "
                                          "placement at the locus and NH>1 (partition: "
                                          "shared, multimapping)",
                "txome_cov_txome_only": "txome_cov, reads with no top-score genome "
                                        "placement at the locus (partition: "
                                        "transcriptome-only at gene)",
                "stacking": "GENOME_LAYERS / TXOME_LAYERS order, baseline first; each "
                            "route's layers sum exactly to its track"}},
    }
    return arrays, meta


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bams")
    parser.add_argument("--gtf")
    parser.add_argument("--appris")
    parser.add_argument("--qc-genome", required=True)
    parser.add_argument("--qc-txome", required=True)
    parser.add_argument("--output")
    args = parser.parse_args(argv)

    # published locus  lrrfip1 in the hela example library
    gene, sample, gsm = "LRRFIP1", "HeLa", "GSM2100602"
    inputs = paths.resolve_external_inputs(args.bams, args.gtf, args.appris, sample)
    qc_genome = str(paths.repo_path(args.qc_genome))
    qc_txome = str(paths.repo_path(args.qc_txome))
    for qc in (qc_genome, qc_txome):
        if not os.path.exists(qc):
            die("QC table missing: %s" % qc)
    stem = args.output or os.path.join(paths.REPO, "results", "read_categories",
                                       "locus_%s" % gene)

    arrays, meta = build(gene, sample, inputs, qc_genome, qc_txome, "psite")
    meta["gsm"] = gsm
    meta["inputs"] = {key: {"file": os.path.basename(inputs[key]),
                            "sha256": paths.sha256_of(inputs[key])}
                      for key in ("ribo_genome", "ribo_txome", "gtf", "appris")}
    meta["inputs"]["qc_genome"] = {"file": os.path.relpath(qc_genome, str(paths.REPO)),
                                   "sha256": paths.sha256_of(qc_genome)}
    meta["inputs"]["qc_txome"] = {"file": os.path.relpath(qc_txome, str(paths.REPO)),
                                  "sha256": paths.sha256_of(qc_txome)}
    meta["builder"] = "code/read_categories/build_locus_data.py"

    os.makedirs(os.path.dirname(stem), exist_ok=True)
    np.savez(stem + ".npz", **arrays)          # uncompressed for deterministic bytes
    with open(stem + ".json", "w") as handle:
        handle.write(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print("[locus] wrote %s.npz / .json" % stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
