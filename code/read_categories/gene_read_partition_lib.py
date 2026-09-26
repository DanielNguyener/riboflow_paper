#!/usr/bin/env python3
"""Gene-anchored read-ID partition: every read at a gene, on either route, in one chain.

Denominator = the UNION of read IDs at the gene on either route. A read enters through the
genome route when a QUALIFYING genome alignment -- its primary, or a secondary whose AS
equals the primary's (a tied-best placement) -- overlaps the gene's full multi-isoform span;
it enters through the transcriptome route when its transcriptome primary is the gene's
selected transcript. A read counts once per union; a read with tied-best placements at two
genes is in both unions, so unions are never summed as a partition of the library.

Genome state is the primary's, as in Figures 3-4: unique = NH == 1, multi = NH > 1,
absent = no primary. The two mechanisms are tested at the gene being classified, never at
whichever placement STAR happened to call primary:
  pseudogene tie   a top-score placement whose 5' base is in an exon of THIS gene, and a
                   different top-score placement whose 5' base is a processed pseudogene
  alternative exon a genome-only unique read with an aligned block on exonic sequence of
                   this gene that its selected transcript omits
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

if str(REPO / "code" / "common") not in sys.path:
    sys.path.insert(0, str(REPO / "code" / "common"))
from intervals import merge as _merge, subtract as _subtract  # noqa: E402

if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from categories import MISSING_AS, qualifies  # noqa: E402

#: The chain, in evaluation order; one read gets exactly one of these. The order nests by
#: definedness — `classify_gU_tA` must only see transcriptome-ABSENT reads.
PARTITION_CATEGORIES = (
    "txome_only_genome_absent",
    "txome_only_genome_elsewhere",
    "genome_multi_top_at_gene_pseudogene_tie",
    "genome_multi_top_at_gene_no_pseudogene_tie",
    "genome_unique_txome_present",
    "genome_unique_absent_omitted_exon",
    "genome_unique_absent_splice_junction",
    "genome_unique_absent_representable",
    "genome_unique_absent_pseudogene",
    "genome_unique_absent_other",
)


#: `reach_lib` category -> this chain's label, for a genome-only unique read that does NOT
#: overlap the gene's omitted exonic sequence. `reach_lib` (Figure 4D) is not a gene-level
#: test -- it calls a strand mismatch or an unrepresentable block "nonselected isoform
#: exon" -- so here it only names the rest; its `nonselected_isoform_exon` falls to
#: `genome_unique_absent_other`. The raw label survives in the `--dump-reads` table.
REACH_TO_CATEGORY = {
    "splice_junction_absent": "genome_unique_absent_splice_junction",
    "representable_not_present_in_dedup_bam": "genome_unique_absent_representable",
    "pseudogene": "genome_unique_absent_pseudogene",
}

TIDY_COLUMNS = ["sample", "gene_id", "gene_name", "transcript_id", "category",
                "n_reads", "pct_of_union"]
WIDE_COLUMNS = ["sample", "gene_id", "gene_name", "transcript_id",
                "gene_chromosome", "gene_start", "gene_end",
                "n_union", "n_genome_side", "n_txome_side", "n_shared",
                "n_genome_only", "n_txome_only", "n_genome_unique", "n_genome_multi",
                "n_joined_by_tied_secondary_only", "n_secondary_only_no_primary",
                "n_secondary_only_lower_score", "n_missing_as"]




class PartitionError(RuntimeError):
    pass


def top_placements(recs):
    """The read's qualifying records, from [(chrom, pos5, AS, is_secondary), ...]."""
    primary_score = next((r[2] for r in recs if not r[3]), None)
    return [r for r in recs if qualifies(r[3], r[2], primary_score)]


def load_libraries():
    """The read-taxonomy libraries, imported by path without disturbing sys.path."""
    sys.path.insert(0, str(REPO / "code"))
    from common.inputs import import_from
    dirs = (REPO / "code" / "read_categories", REPO / "code" / "common",
            REPO / "code" / "common" / "ribo_seq_qc")
    return tuple(import_from(dirs[0], name, dirs[1:])
                 for name in ("reference_lib", "reach_lib", "tie_biotype_lib"))


def build_span_table(exon_gene_df):
    """Version-stripped gene id -> (chromosome, start, end, n_chromosomes), once for
    every gene; the per-call aggregation rescanned a 1.4-million-row frame each time."""
    base = exon_gene_df["gene_id"].astype(str).str.split(".").str[0]
    return exon_gene_df.assign(_base=base).groupby("_base").agg(
        chrom=("Chromosome", "first"), start=("Start", "min"),
        end=("End", "max"), n_chrom=("Chromosome", "nunique"))


def gene_locus(annotation, gene_id):
    """The gene's full genomic span, over EVERY annotated isoform.

    Not the selected transcript's own span: nonselected-isoform exons lie outside it.
    """
    base = str(gene_id).split(".", 1)[0]
    try:
        row = annotation["spans"].loc[base]
    except KeyError:
        raise PartitionError("gene %r has no exon in the annotation" % gene_id)
    if int(row["n_chrom"]) != 1:
        exon_gene_df = annotation["exon_gene_df"]
        members = exon_gene_df.loc[
            exon_gene_df["gene_id"].astype(str).str.split(".").str[0] == base,
            "Chromosome"].unique()
        raise PartitionError("gene %r spans %d chromosomes (%s); this is not handled"
                             % (gene_id, len(members), ", ".join(map(str, members))))
    return str(row["chrom"]), int(row["start"]), int(row["end"])


def fetch_gene_candidates(genome_bam, chrom, start, end):
    """{read ID: [(is_secondary, AS), ...]} for every alignment overlapping [start, end).

    Indexed fetch, strand-agnostic. Candidates only: whether a secondary qualifies needs the
    primary's AS, which may sit anywhere in the genome (`resolve_genome_side`).
    """
    import pysam

    candidates = defaultdict(list)
    bam = pysam.AlignmentFile(str(genome_bam), "rb")
    try:
        if not bam.has_index():
            raise PartitionError(
                "%s has no index. The gene-anchored denominator is an indexed region "
                "fetch; index it with `samtools index`." % genome_bam)
        for read in bam.fetch(str(chrom), int(start), int(end)):
            if read.is_unmapped or read.is_supplementary:
                continue
            score = int(read.get_tag("AS")) if read.has_tag("AS") else MISSING_AS
            candidates[read.query_name].append((bool(read.is_secondary), score))
    finally:
        bam.close()
    return dict(candidates)


def resolve_genome_side(candidates, primary):
    """Candidates -> (genome side, {read: "primary" | "tied_secondary"}, audit counts).

    `primary` maps read -> (chrom, strand, blocks, NH, AS) of its primary record.
    """
    side, joined_by = set(), {}
    audit = dict.fromkeys(("n_joined_by_tied_secondary_only", "n_secondary_only_no_primary",
                           "n_secondary_only_lower_score", "n_missing_as"), 0)
    for qname, alignments in candidates.items():
        record = primary.get(qname)
        primary_score = record[4] if record else None
        missing = (primary_score == MISSING_AS
                   or any(score == MISSING_AS for _sec, score in alignments))
        if missing:
            audit["n_missing_as"] += 1
        if any(not sec for sec, _score in alignments):
            side.add(qname)
            joined_by[qname] = "primary"
        elif any(qualifies(sec, score, primary_score) for sec, score in alignments):
            side.add(qname)
            joined_by[qname] = "tied_secondary"
            audit["n_joined_by_tied_secondary_only"] += 1
        elif record is None:
            audit["n_secondary_only_no_primary"] += 1
        elif not missing:
            audit["n_secondary_only_lower_score"] += 1
    return side, joined_by, audit


def collect_genome_state(genome_bam, target_qnames):
    """One full genome pass -> (present, unique, primary, records) for `target_qnames`.

    A full pass, not a fetch: a multimapper's other loci are anywhere in the genome and the
    tie test needs all of them.
    """
    import pysam

    sys.path.insert(0, str(REPO / "code" / "common"))
    import bam_inputs

    present, unique = set(), set()
    primary = {}
    records = defaultdict(list)

    bam = pysam.AlignmentFile(str(genome_bam), "rb")
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                continue
            qname = read.query_name
            if qname not in target_qnames:
                continue
            try:
                nh = int(read.get_tag("NH"))
            except KeyError:
                raise PartitionError(
                    "%s has an alignment with no NH tag. This analysis is about "
                    "multimapping, so a BAM that cannot report NH cannot answer it."
                    % genome_bam)
            blocks = read.get_blocks()
            # A missing AS is a sentinel far below any real score, never None and never 0:
            # `tie_biotype_lib` compares scores for equality to find ties, so None crashes it
            # and 0 would let an unscored alignment tie with a real one. Same convention as
            # `tie_biotype_lib._MISSING_AS`. STAR always emits AS, so this is unreachable on
            # the published cohort and cannot move its numbers.
            score = int(read.get_tag("AS")) if read.has_tag("AS") else MISSING_AS
            if not read.is_secondary:
                present.add(qname)
                if bam_inputs.is_unique_genome_read(read):
                    unique.add(qname)
                if blocks:
                    primary[qname] = (read.reference_name,
                                      "-" if read.is_reverse else "+", blocks, nh, score)
            if not blocks:
                continue
            pos5 = blocks[-1][1] - 1 if read.is_reverse else blocks[0][0]
            records[qname].append(
                (read.reference_name, pos5, score, bool(read.is_secondary)))
    finally:
        bam.close()
    return present, unique, primary, dict(records)


def _in_locus(record, chrom, start, end):
    """Does this primary record overlap the gene locus? Strand-agnostic."""
    if record is None:
        return False
    r_chrom, _strand, blocks, _nh, _as = record
    if str(r_chrom) != str(chrom):
        return False
    return min(b[0] for b in blocks) < end and max(b[1] for b in blocks) > start




def _overlaps(blocks, intervals):
    return any(b_start < i_end and i_start < b_end
               for b_start, b_end in blocks for i_start, i_end in intervals)


def gene_exons(annotation, gene_id):
    """(chrom, gene_type, merged exons) of one gene over EVERY annotated transcript.

    Indexed once per annotation (the per-gene filter would otherwise scan every exon row).
    """
    index = annotation.get("_gene_exons")
    if index is None:
        frame = annotation["exon_gene_df"]
        base = frame["gene_id"].astype(str).str.split(".").str[0]
        index = {}
        for gene, rows in frame.groupby(base, sort=False):
            index[gene] = (str(rows["Chromosome"].iat[0]), str(rows["gene_type"].iat[0]),
                           _merge(zip(rows["Start"].tolist(), rows["End"].tolist())))
        annotation["_gene_exons"] = index
    return index.get(str(gene_id).split(".", 1)[0])


def omitted_exonic_sequence(annotation, tid):
    """(chrom, intervals): the gene's exonic sequence its selected transcript omits."""
    info = annotation["table"][tid]
    chrom, _type, exons = gene_exons(annotation, info["gene_id"])
    selected = list(zip(info["g_start"].tolist(), info["g_end"].tolist()))
    return chrom, _subtract(exons, selected)


def alt_exon_overlap(annotation, tid, primary, qnames):
    """{read: bool}: an aligned block of the read's primary lies on omitted exonic sequence.

    Strand-agnostic, and indifferent to junctions and to whether the rest of the alignment
    fits the selected transcript: the only question is the omitted sequence.
    """
    chrom, omitted = omitted_exonic_sequence(annotation, tid)
    out = {}
    for qname in qnames:
        record = primary.get(qname)
        out[qname] = bool(record is not None and str(record[0]) == chrom
                          and _overlaps(record[2], omitted))
    return out


def gene_pseudogene_tie(tie_biotype_lib, annotation, gene_id, records_by_qname):
    """{read: (tie, n_top, detail)}: a protein-coding / processed-pseudogene tie AT THIS GENE.

    Among the read's top-score placements (`top_placements`), one must have its 5' base in
    an exon of `gene_id` and a DIFFERENT one its 5' base on a processed pseudogene (the Fig 4
    5'-base biotype join, `tie_biotype_lib.classify_loci_frame`). Which of them STAR made
    primary does not enter; a tie with another protein-coding gene is not a tie here.
    """
    import bisect

    chrom, gene_type, exons = gene_exons(annotation, gene_id)
    starts = [e[0] for e in exons]

    def in_gene(record):
        if str(record[0]) != chrom:
            return False
        i = bisect.bisect_right(starts, record[1]) - 1
        return i >= 0 and exons[i][0] <= record[1] < exons[i][1]

    tops = {q: top_placements(recs) for q, recs in records_by_qname.items()}
    rows = [(i, q, r[0], int(r[1]), int(r[2]), bool(r[3]))
            for i, (q, r) in enumerate((q, r) for q, recs in tops.items() for r in recs)]
    frame = pd.DataFrame(rows, columns=list(tie_biotype_lib.LOCUS_FRAME_COLUMNS))
    biotypes = (tie_biotype_lib.classify_loci_frame(
        frame, annotation["exon_pr"], annotation["gene_body_pr"])["biotype"].tolist()
        if rows else [])

    out, cursor = {}, 0
    for qname, top in tops.items():
        types = biotypes[cursor:cursor + len(top)]
        cursor += len(top)
        at_gene = {i for i, r in enumerate(top) if in_gene(r)}
        on_pp = {i for i, t in enumerate(types) if t == tie_biotype_lib.PP}
        tie = gene_type == tie_biotype_lib.PC and any(
            i != j for i in at_gene for j in on_pp)
        detail = ";".join("%s:%d%s%s%s" % (r[0], r[1], "" if r[3] else "*",
                                           "@gene" if i in at_gene else "",
                                           "@pp" if i in on_pp else "")
                          for i, r in enumerate(top))
        out[qname] = (tie, len(top), detail)
    return out


def classify_union(libs, annotation, tid, locus, genome_side, txome_side,
                   txome_present, genome_present, genome_unique, primary, records):
    """The chain. Returns (labels: Series qname -> category, detail: per-read mechanism data).

    `genome_side` must already hold only reads with a qualifying placement at the gene
    (`resolve_genome_side`, `ReadState.gene_side`). Transcriptome presence is membership in
    `txome_present` (a primary alignment in the post-dedup BAM, which RiboFlow_v2 filtered at
    MAPQ >= 10); no MAPQ test is made here.
    """
    _reference_lib, reach_lib, tie_biotype_lib = libs

    genome_side = set(genome_side)
    txome_side = set(txome_side)
    union = genome_side | txome_side
    labels = pd.Series(index=sorted(union), dtype=object)
    tie_detail, reach_label, alt_exon = {}, {}, {}

    # 1-2. Transcriptome route here, no top-score genome placement HERE: absent vs elsewhere.
    for qname in txome_side - genome_side:
        labels[qname] = ("txome_only_genome_elsewhere" if qname in genome_present
                         else "txome_only_genome_absent")

    at_gene = sorted(genome_side)
    # A qualifying placement needs a primary (a secondary ties only against the primary's AS).
    orphan = [q for q in at_gene if q not in genome_present]
    if orphan:
        raise PartitionError(
            "%d genome-side read(s) have no primary genome alignment, e.g. %s; a secondary "
            "is never tied without one" % (len(orphan), orphan[:5]))
    uniq = [q for q in at_gene if q in genome_unique]
    multi = [q for q in at_gene if q not in genome_unique]

    # 3-4. Genome multimappers with a top-score placement here: the tie is tested at THIS
    # gene, over every top-score placement, whichever of them is primary.
    if multi:
        ties = gene_pseudogene_tie(tie_biotype_lib, annotation,
                                   annotation["table"][tid]["gene_id"],
                                   {q: records.get(q, []) for q in multi})
        for qname in multi:
            tie, n_top, detail = ties[qname]
            tie_detail[qname] = (tie, n_top, detail)
            labels[qname] = ("genome_multi_top_at_gene_pseudogene_tie" if tie
                             else "genome_multi_top_at_gene_no_pseudogene_tie")

    # 5. Genome-unique and present on the transcriptome route (on this transcript or any).
    for qname in uniq:
        if qname in txome_present:
            labels[qname] = "genome_unique_txome_present"

    # 6-10. Genome-unique, absent from every transcriptome primary. The alternative exon is
    # the gene-level overlap test; `reach_lib` only names the rest.
    absent = [q for q in uniq if q not in txome_present]
    if absent:
        alt_exon = alt_exon_overlap(annotation, tid, primary, absent)
        blocks = {q: (primary[q][0], primary[q][1], primary[q][2])
                  for q in absent if q in primary}
        reach = reach_lib.classify_gU_tA(
            absent, blocks, annotation["exon_pr"], annotation["exon_gene_df"],
            annotation["gene_body_pr"], annotation["table"], annotation["gene2tid"],
            annotation["omitted_genes"])
        for qname in absent:
            reach_label[qname] = reach.get(qname, "")
            labels[qname] = ("genome_unique_absent_omitted_exon" if alt_exon[qname]
                             else REACH_TO_CATEGORY.get(reach_label[qname],
                                                        "genome_unique_absent_other"))

    unlabelled = labels[labels.isna()]
    if len(unlabelled):
        raise PartitionError(
            "%d read(s) fell through the chain, e.g. %s"
            % (len(unlabelled), list(unlabelled.index[:5])))

    detail = {"tie": tie_detail, "reach_label": reach_label, "alt_exon": alt_exon}
    return labels, detail


def load_annotation(libs):
    """Everything the two category systems need, built once for all genes."""
    import pyranges as pr

    reference_lib, reach_lib, _tie = libs
    payload = reference_lib.build_transcript_table()
    table = payload["table"]
    exon_gene_df = reference_lib.build_exon_gene_table()
    return {
        "table": table,
        "base2ver": payload["base2ver"],
        "exon_pr": reference_lib.load_exon_gene_pr(),
        "exon_gene_df": exon_gene_df,
        "spans": build_span_table(exon_gene_df),
        "gene_body_pr": pr.PyRanges(
            reach_lib.fc.config.load_all_gene_bodies().reset_index(drop=True)),
        "gene2tid": reach_lib.gene_to_transcript_map(table),
        "omitted_genes": reach_lib.omitted_pc_genes(
            exon_gene_df, set(v["gene_id"] for v in table.values())),
    }


class FateError(RuntimeError):
    pass


def resolve_transcripts(table, gene_ids=(), transcript_ids=(), coverage=None):
    """Gene and/or transcript IDs -> versioned APPRIS transcript IDs present in `table`.

    Refuses to guess when a gene maps to more than one candidate, listing them.
    """
    resolved = []
    for tid in transcript_ids:
        if tid in table:
            resolved.append(tid)
            continue
        base = tid.split(".", 1)[0]
        hits = [t for t in table if t.split(".", 1)[0] == base]
        if len(hits) == 1:
            resolved.append(hits[0])
        elif not hits:
            raise FateError("transcript %r is not in the APPRIS transcript table" % tid)
        else:
            raise FateError("transcript %r is ambiguous: %s" % (tid, ", ".join(hits)))

    for gene_id in gene_ids:
        if coverage is not None:
            index = coverage.resolve_gene(gene_id)
            resolved.append(coverage.transcript_info(index)["transcript_id"])
            continue
        base = gene_id.split(".", 1)[0]
        hits = [tid for tid, entry in table.items()
                if str(entry.get("gene_id", "")).split(".", 1)[0] == base]
        if len(hits) == 1:
            resolved.append(hits[0])
        elif not hits:
            raise FateError("gene %r has no transcript in the APPRIS table" % gene_id)
        else:
            raise FateError(
                "gene %r maps to %d transcripts and no unique one resolves it; pass the "
                "versioned transcript ID instead.\n%s"
                % (gene_id, len(hits), "\n".join("    %s" % h for h in sorted(hits))))
    return list(dict.fromkeys(resolved))


def compute_partition(sample, genome_bam, txome_bam, gene_ids=(), transcript_ids=(),
                      coverage=None, log=lambda _m: None):
    """The whole computation. Returns (wide, tidy, dump)."""
    sys.path.insert(0, str(HERE))

    libs = load_libraries()
    reference_lib = libs[0]
    annotation = load_annotation(libs)
    table = annotation["table"]

    tids = resolve_transcripts(
        table, gene_ids, transcript_ids, coverage)
    if not tids:
        raise PartitionError("no transcripts requested")

    log("reading the transcriptome BAM")
    txome_present, _all = reference_lib.read_txome_primary(txome_bam, annotation["base2ver"])
    txome_side = {tid: set() for tid in tids}
    for qname, hit in txome_present.items():
        if hit in txome_side:
            txome_side[hit].add(qname)

    names, genes = _display(table, tids, coverage)
    log("fetching the genome-side candidates of each gene")
    loci, candidates = {}, {}
    for tid in tids:
        loci[tid] = gene_locus(annotation, table[tid]["gene_id"])
        candidates[tid] = fetch_gene_candidates(genome_bam, *loci[tid])

    targets = set()
    for tid in tids:
        targets |= set(candidates[tid]) | txome_side[tid]
    log("one genome pass over %d read ids" % len(targets))
    g_present, g_unique, primary, records = collect_genome_state(genome_bam, targets)

    genome_side, joined_by, audit = {}, {}, {}
    for tid in tids:
        genome_side[tid], joined_by[tid], audit[tid] = resolve_genome_side(
            candidates[tid], primary)
        log("  %-8s %s:%d-%d  %d genome-side (%d by a tied secondary only), %d "
            "transcriptome-side read ids; not placed: %d secondary-only without a primary, "
            "%d secondary-only below the primary's AS; %d with a missing AS"
            % (names[tid] or tid, loci[tid][0], loci[tid][1], loci[tid][2],
               len(genome_side[tid]), audit[tid]["n_joined_by_tied_secondary_only"],
               len(txome_side[tid]), audit[tid]["n_secondary_only_no_primary"],
               audit[tid]["n_secondary_only_lower_score"], audit[tid]["n_missing_as"]))

    wide_rows, tidy_rows, dump_rows = [], [], []
    for tid in tids:
        labels, detail = classify_union(
            libs, annotation, tid, loci[tid], genome_side[tid], txome_side[tid],
            txome_present, g_present, g_unique, primary, records)
        counts = labels.value_counts()
        n_union = int(len(labels))
        if int(counts.sum()) != n_union:
            raise PartitionError("%s %s: labels (%d) do not cover the union (%d)"
                                 % (sample, tid, int(counts.sum()), n_union))
        unknown = set(counts.index) - set(PARTITION_CATEGORIES)
        if unknown:
            raise PartitionError("undeclared category/ies: %s" % sorted(unknown))

        # Every category is emitted, zero included.
        for category in PARTITION_CATEGORIES:
            n_reads = int(counts.get(category, 0))
            tidy_rows.append({
                "sample": sample, "gene_id": genes[tid], "gene_name": names[tid],
                "transcript_id": tid, "category": category, "n_reads": n_reads,
                "pct_of_union": (100.0 * n_reads / n_union) if n_union else 0.0})

        shared = genome_side[tid] & txome_side[tid]
        wide_rows.append({
            "sample": sample, "gene_id": genes[tid], "gene_name": names[tid],
            "transcript_id": tid, "gene_chromosome": loci[tid][0],
            "gene_start": loci[tid][1], "gene_end": loci[tid][2],
            "n_union": n_union, "n_genome_side": len(genome_side[tid]),
            "n_txome_side": len(txome_side[tid]), "n_shared": len(shared),
            "n_genome_only": len(genome_side[tid] - txome_side[tid]),
            "n_txome_only": len(txome_side[tid] - genome_side[tid]),
            "n_genome_unique": len(genome_side[tid] & g_unique),
            "n_genome_multi": len(genome_side[tid] & (g_present - g_unique)),
            **audit[tid]})

        for qname, category in labels.items():
            record = primary.get(qname)
            blocks = record[2] if record else None
            tie, n_top, placements = detail["tie"].get(qname, (False, -1, ""))
            dump_rows.append({
                "sample": sample, "gene_name": names[tid], "transcript_id": tid,
                "read_id": qname, "category": category,
                "genome_placement_at_gene": joined_by[tid].get(qname, ""),
                "pseudogene_tie_at_gene": bool(tie),
                "n_top_placements": n_top,
                "top_placements": placements,
                "alt_exon_overlap": bool(detail["alt_exon"].get(qname, False)),
                "reach_category": detail["reach_label"].get(qname, ""),
                "genome_primary_nh": record[3] if record else -1,
                "genome_primary_chromosome": record[0] if record else "",
                "genome_primary_strand": record[1] if record else "",
                "genome_primary_start": min(b[0] for b in blocks) if blocks else -1,
                "genome_primary_end": max(b[1] for b in blocks) if blocks else -1,
                "genome_primary_in_gene": _in_locus(record, *loci[tid]),
                "n_genome_loci": len(records.get(qname, ())),
                "on_genome_side": qname in genome_side[tid],
                "txome_primary_transcript": txome_present.get(qname, "")})

    tidy = pd.DataFrame(tidy_rows, columns=TIDY_COLUMNS)
    wide = pd.DataFrame(wide_rows, columns=WIDE_COLUMNS)
    return wide, tidy, pd.DataFrame(dump_rows)


def _display(table, tids, coverage):
    """Gene names and ids, preferring the coverage file when one is supplied."""
    names, genes = {}, {}
    for tid in tids:
        name = table[tid].get("gene_name", "")
        names[tid] = "" if name is None or isinstance(name, float) else str(name)
        genes[tid] = str(table[tid].get("gene_id", "") or "")
    if coverage is not None:
        for tid in tids:
            try:
                index = coverage.index_of_transcript(tid)
            except Exception:
                continue
            info = coverage.transcript_info(index)
            names[tid] = info["gene_name"]
            genes[tid] = info["gene_id"]
    return names, genes
