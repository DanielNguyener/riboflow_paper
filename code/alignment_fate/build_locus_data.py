#!/usr/bin/env python3
"""BAMs + GTF + APPRIS + QC tables -> the locus artifact (`locus_<GENE>.{npz,json}`).

Per-base P-site coverage of both routes over the merged exons of the selected isoform and
the best-supported alternative isoform (the one carrying the most genome-only unique reads
on non-selected sequence). Writes <output>.npz (coverage vectors + exon blocks) and
<output>.json (metadata). Introns are drawn as a constant 90-unit gap. The BAM template,
MAPQ >= 42 rule and cigar-aware P-site are re-implemented here, identically to `common/`
and `coverage/`.

Published locus: LRRFIP1, chr2:237,627,586-237,781,643 (+), selected ENST00000308482.14,
alternative ENST00000244815.9 (3,619 nt absent from the selected reference).
Run with `python` (3.9).
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))
import inputs as paths  # noqa: E402
from inputs import die  # noqa: E402

#: The shipped QC tables: the one input that legitimately defaults into the repository.
QC_GENOME_DEFAULT = "data/ribo_seq_qc/genome/tables/readlen_window_qc.csv"
QC_TXOME_DEFAULT = "data/ribo_seq_qc/transcriptome/tables/readlen_window_qc.csv"

#: Transcriptome-route uniqueness (bowtie2 emits no NH tag): the project-wide rule.
TXOME_MIN_MAPQ = 42
#: Fixed width, in plotted units, of the dashed intron connector; recorded in the artifact.
INTRON_GAP = 90.0


# ── QC tables ────────────────────────────────────────────────────────────────

def read_window_and_offsets(qc_csv, sample):
    """The sample's selected read lengths and their P-site offsets, from the QC table."""
    offsets = {}
    with open(qc_csv) as handle:
        for row in csv.DictReader(handle):
            if row["sample"] != sample or row["in_phase1"] != "True":
                continue
            length = int(row["read_length"])
            raw = row.get("psite_offset", "")
            if raw in ("", "NA", "None"):
                die("%s: read length %d is selected but has no P-site offset in %s"
                    % (sample, length, qc_csv))
            offsets[length] = int(float(raw))
    if not offsets:
        die("no selected read lengths for %s in %s" % (sample, qc_csv))
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


def transcript_to_genome(exons, strand):
    ordered = sorted(exons) if strand == "+" else sorted(exons, reverse=True)
    out = []
    for start, end in ordered:
        out.append(np.arange(start, end) if strand == "+"
                   else np.arange(end - 1, start - 1, -1))
    return np.concatenate(out) if out else np.array([], dtype=int)


def gene_transcripts(gtf_path, gene_name):
    """{transcript_id: [(start, end), ...]} (0-based half-open) for one gene, + chrom, strand."""
    exons = collections.defaultdict(list)
    chrom = strand = None
    needle = 'gene_name "%s"' % gene_name
    with open(gtf_path) as handle:
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


def merge(intervals):
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(s, e) for s, e in out]


def subtract(intervals, holes):
    out = []
    for start, end in intervals:
        pieces = [(start, end)]
        for h_start, h_end in holes:
            nxt = []
            for p_start, p_end in pieces:
                if h_end <= p_start or h_start >= p_end:
                    nxt.append((p_start, p_end))
                    continue
                if p_start < h_start:
                    nxt.append((p_start, h_start))
                if h_end < p_end:
                    nxt.append((h_end, p_end))
            pieces = nxt
        out.extend(pieces)
    return [(s, e) for s, e in out if e > s]


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
    """{qname: (blocks, nh, psite)} for every primary ribo alignment at the locus, in window."""
    out = {}
    dropped = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        if not bam.has_index():
            die("%s has no index; the locus view is a region fetch" % bam_path)
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_length not in lengths:
                continue
            blocks = read.get_blocks()
            if not blocks:
                continue
            try:
                nh = int(read.get_tag("NH"))
            except KeyError:
                die("%s has an alignment with no NH tag" % bam_path)
            psite = psite_reference_position(read, offsets[read.query_length])
            if psite is None:
                dropped += 1
            out[read.query_name] = (blocks, nh, psite)
    finally:
        bam.close()
    return out, dropped


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


def txome_unique_coverage(bam_path, reference, lengths, offsets, tx2genome, signal,
                          min_mapq=TXOME_MIN_MAPQ):
    """Genomic per-base depth from the transcriptome route's uniquely mapping reads."""
    depth = collections.Counter()
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
            if read.mapping_quality < min_mapq:
                continue
            if read.query_length not in lengths:
                continue
            n_reads += 1
            if signal == "psite":
                psite = psite_reference_position(read, offsets[read.query_length])
                if psite is None or not (0 <= psite < len(tx2genome)):
                    dropped += 1
                    continue
                depth[int(tx2genome[psite])] += 1
            else:
                lo = max(read.reference_start, 0)
                hi = min(read.reference_end, len(tx2genome))
                for pos in tx2genome[lo:hi]:
                    depth[int(pos)] += 1
    finally:
        bam.close()
    return depth, n_reads, dropped


def coverage_over(reads, positions, signal):
    """Depth at `positions` from {qname: (blocks, nh, psite)}; psite = one count per read."""
    depth = collections.Counter()
    for blocks, _nh, psite in reads.values():
        if signal == "psite":
            if psite is not None:
                depth[int(psite)] += 1
        else:
            for start, end in blocks:
                for pos in range(start, end):
                    depth[pos] += 1
    return np.array([depth.get(int(p), 0) for p in positions], dtype=float)


# ── driver ───────────────────────────────────────────────────────────────────

def build(gene, sample, inputs, qc_genome, qc_txome, signal, alt_transcript=None):
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

    tx2genome = transcript_to_genome(sel_exons, strand)
    if len(tx2genome) != sel_len:
        die("%s: GTF exons give a spliced length of %d, but the transcriptome reference "
            "declares %d -- the projection would be wrong" % (sel_tid, len(tx2genome), sel_len))

    reads, g_dropped = locus_ribo_reads(inputs["ribo_genome"], chrom, span_start, span_end,
                                       lengths, offsets)
    genome_unique = {q for q, (_b, nh, _p) in reads.items() if nh == 1}
    print("[locus] %d ribo reads at the locus in the window, %d uniquely mapping"
          % (len(reads), len(genome_unique)))
    if g_dropped:
        print("[locus]   %d genome reads have no resolvable P-site and are dropped" % g_dropped)

    txome_depth, n_txome, t_dropped = txome_unique_coverage(
        inputs["ribo_txome"], sel_header, t_lengths, t_offsets, tx2genome, signal)
    print("[locus] %d transcriptome-route reads on %s (MAPQ >= %d), projected onto the "
          "genome axis" % (n_txome, sel_tid, TXOME_MIN_MAPQ))
    if t_dropped:
        print("[locus]   %d transcriptome reads dropped for the same reason" % t_dropped)

    print("[locus] one transcriptome pass to find which genome reads are genome-only ...")
    present = txome_present(inputs["ribo_txome"], set(reads))
    genome_only = {q for q in genome_unique if q not in present}
    print("[locus] genome-only uniquely mapping: %d of %d unique"
          % (len(genome_only), len(genome_unique)))

    # Alternative isoform: most genome-only unique reads on non-selected sequence.
    scores = {}
    if alt_transcript:
        alt_tid = alt_transcript
        if alt_tid not in transcripts:
            die("%s is not an annotated transcript of %s" % (alt_tid, gene))
    else:
        for tid, exons in transcripts.items():
            if tid == sel_tid:
                continue
            extra = subtract(merge(exons), sel_exons)
            if not extra:
                continue
            n = 0
            for qname in genome_only:
                blocks, _nh, _psite = reads[qname]
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

    # Per-base vectors exactly as the panel draws them: union exonic bases, 5'->3' order.
    union = merge(list(sel_exons) + list(alt_exons))
    gs = exonic_bases(union, strand)
    genome_cov = coverage_over({q: reads[q] for q in genome_unique}, gs, signal)
    txome_cov = np.array([txome_depth.get(int(p), 0) for p in gs], dtype=float)

    arrays = {
        "genomic_position": gs.astype(np.int64),
        "genome_cov": genome_cov.astype(np.float64),
        "txome_cov": txome_cov.astype(np.float64),
        "sel_exons": np.array(sel_exons, dtype=np.int64).reshape(-1, 2),
        "alt_exons": np.array(alt_exons, dtype=np.int64).reshape(-1, 2),
        "absent_blocks": np.array(absent_blocks, dtype=np.int64).reshape(-1, 2),
    }
    meta = {
        "gene": gene, "chrom": chrom, "strand": strand,
        "locus": {"start": int(union[0][0]), "end": int(union[-1][1]),
                  "label": "%s:%s-%s (%s)" % (chrom, format(union[0][0], ","),
                                              format(union[-1][1], ","), strand)},
        "selected_transcript": sel_tid, "selected_length": sel_len,
        "alternative_transcript": alt_tid,
        "alternative_chosen_by": "forced" if alt_transcript else
                                 "most genome-only unique reads on non-selected sequence",
        "alternative_support_top5": dict(sorted(scores.items(), key=lambda kv: -kv[1])[:5]),
        "n_absent_nt": int(sum(e - s for s, e in absent_blocks)),
        "sample": sample,
        "signal": signal,
        "intron_gap_plot_units": INTRON_GAP,
        "txome_min_mapq": TXOME_MIN_MAPQ,
        "genome_window": {"lengths": sorted(lengths),
                          "offsets": {str(k): v for k, v in sorted(offsets.items())}},
        "txome_window": {"lengths": sorted(t_lengths),
                         "offsets": {str(k): v for k, v in sorted(t_offsets.items())}},
        "counts": {"genome_reads_at_locus_in_window": len(reads),
                   "genome_unique": len(genome_unique),
                   "genome_only_unique": len(genome_only),
                   "genome_psite_dropped": g_dropped,
                   "txome_unique_on_selected": n_txome,
                   "txome_psite_dropped": t_dropped,
                   "genome_cov_total": float(genome_cov.sum()),
                   "txome_cov_total": float(txome_cov.sum())},
        "coverage_semantics": {
            "genome_cov": "genome route, NH==1 primaries in the read-length window, one "
                          "count per read at its P-site (cigar_aware)",
            "txome_cov": "transcriptome route, MAPQ>=42 primaries on the selected "
                         "transcript, P-site projected through the transcript->genome map",
            "order": "5'->3' over the merged exons of both models (SplicedAxis order)"},
    }
    return arrays, meta


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gene", required=True, help="gene NAME, e.g. LRRFIP1")
    parser.add_argument("--sample", default="HeLa")
    parser.add_argument("--gsm", default="GSM2100602")
    parser.add_argument("--bams", help="RiboFlow output root (else RIBOFLOW_PAPER_BAMS)")
    parser.add_argument("--gtf", help="GENCODE GTF (else RIBOFLOW_PAPER_GTF)")
    parser.add_argument("--appris", help="reference lengths TSV (else RIBOFLOW_PAPER_APPRIS)")
    parser.add_argument("--qc-genome", default=QC_GENOME_DEFAULT,
                        help="repo-relative or absolute; the genome route's QC table")
    parser.add_argument("--qc-txome", default=QC_TXOME_DEFAULT,
                        help="the transcriptome route's QC table")
    parser.add_argument("--signal", choices=("psite", "footprint"), default="psite")
    parser.add_argument("--alt-transcript", help="force the alternative isoform")
    parser.add_argument("--output", help="output stem; default results/alignment_fate/locus_<GENE>")
    parser.add_argument("--record-input-paths", action="store_true",
                        help="also record the absolute input paths in the JSON")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    inputs = paths.resolve_external_inputs(args.bams, args.gtf, args.appris, args.sample)
    qc_genome = str(paths.repo_path(args.qc_genome))
    qc_txome = str(paths.repo_path(args.qc_txome))
    for qc in (qc_genome, qc_txome):
        if not os.path.exists(qc):
            die("QC table missing: %s" % qc)
    stem = args.output or os.path.join(paths.REPO, "results", "alignment_fate",
                                       "locus_%s" % args.gene)
    if os.path.exists(stem + ".npz") and not args.force:
        die("%s.npz exists; pass --force" % stem)

    arrays, meta = build(args.gene, args.sample, inputs, qc_genome, qc_txome, args.signal,
                         args.alt_transcript)
    meta["gsm"] = args.gsm
    meta["inputs"] = {key: {"file": os.path.basename(inputs[key]),
                            "sha256": paths.sha256_of(inputs[key])}
                      for key in ("ribo_genome", "ribo_txome", "gtf", "appris")}
    meta["inputs"]["qc_genome"] = {"file": os.path.relpath(qc_genome, str(paths.REPO)),
                                   "sha256": paths.sha256_of(qc_genome)}
    meta["inputs"]["qc_txome"] = {"file": os.path.relpath(qc_txome, str(paths.REPO)),
                                  "sha256": paths.sha256_of(qc_txome)}
    if args.record_input_paths:
        for key in ("ribo_genome", "ribo_txome", "gtf", "appris"):
            meta["inputs"][key]["path"] = os.path.abspath(inputs[key])
    meta["builder"] = "code/alignment_fate/build_locus_data.py"

    os.makedirs(os.path.dirname(stem), exist_ok=True)
    np.savez(stem + ".npz", **arrays)          # uncompressed: deterministic bytes
    with open(stem + ".json", "w") as handle:
        handle.write(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print("[locus] wrote %s.npz / .json" % stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
