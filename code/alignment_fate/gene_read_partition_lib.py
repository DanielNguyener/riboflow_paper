#!/usr/bin/env python3
"""Gene-anchored read-ID partition: every read at a gene, on either route, in one chain.

Denominator = the UNION of read IDs at the gene's full multi-isoform locus on either route,
partitioned by one priority chain reusing `tie_biotype_lib` + `reach_lib` categories.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

#: The chain, in evaluation order; one read gets exactly one of these. The order nests by
#: definedness — `classify_gU_tA` must only see transcriptome-ABSENT reads.
PARTITION_CATEGORIES = (
    "txome_only_genome_absent",
    "txome_only_genome_elsewhere",
    "genome_multi_primary_in_gene_pseudogene_tie",
    "genome_multi_primary_in_gene_no_pseudogene_tie",
    "genome_multi_primary_pseudogene",
    "genome_multi_primary_elsewhere_other",
    "genome_multi_primary_lost_in_dedup",
    "genome_unique_shared_concordant",
    "genome_unique_shared_discordant",
    "genome_unique_txome_other_transcript",
    "genome_unique_txome_multimapped",
    "genome_unique_absent_nonselected_isoform_exon",
    "genome_unique_absent_splice_junction",
    "genome_unique_absent_representable",
    "genome_unique_absent_pseudogene",
    "genome_unique_absent_other",
)

CATEGORY_LABEL = {
    "txome_only_genome_absent":
        "transcriptome-only: no genome alignment at all",
    "txome_only_genome_elsewhere":
        "transcriptome-only: genome placed it elsewhere",
    "genome_multi_primary_in_gene_pseudogene_tie":
        "genome-multi: primary here, tied with a processed pseudogene",
    "genome_multi_primary_in_gene_no_pseudogene_tie":
        "genome-multi: primary here, no pseudogene tie",
    "genome_multi_primary_pseudogene":
        "genome-multi: primary on a processed pseudogene",
    "genome_multi_primary_elsewhere_other":
        "genome-multi: primary elsewhere, other biotype",
    "genome_multi_primary_lost_in_dedup":
        "genome-multi: primary record removed by deduplication",
    "genome_unique_shared_concordant":
        "genome-unique: same locus on both routes",
    "genome_unique_shared_discordant":
        "genome-unique: different locus on the two routes",
    "genome_unique_txome_other_transcript":
        "genome-unique: transcriptome assigned another transcript",
    "genome_unique_txome_multimapped":
        "genome-unique: transcriptome multimapped",
    "genome_unique_absent_nonselected_isoform_exon":
        "genome-only: exon of a nonselected isoform",
    "genome_unique_absent_splice_junction":
        "genome-only: junction absent from the selected isoform",
    "genome_unique_absent_representable":
        "genome-only: representable, absent from the dedup'd transcriptome BAM",
    "genome_unique_absent_pseudogene":
        "genome-only: pseudogene",
    "genome_unique_absent_other":
        "genome-only: intronic, intergenic or other biotype",
}

#: `reach_lib` category -> this chain's label; the uncollapsed label survives in the
#: `--dump-reads` table, so the fold is presentational and reversible.
REACH_TO_CATEGORY = {
    "nonselected_isoform_exon": "genome_unique_absent_nonselected_isoform_exon",
    "splice_junction_absent": "genome_unique_absent_splice_junction",
    "representable_not_present_in_dedup_bam": "genome_unique_absent_representable",
    "pseudogene": "genome_unique_absent_pseudogene",
}

#: Reach labels impossible on a transcriptome-ABSENT population; checked, not trusted.
REACH_IMPOSSIBLE = ("shared_unique_concordant", "shared_unique_discordant",
                    "genome_unique_transcriptome_multimapped")

#: `tie_biotype_lib` class -> category, for a read whose primary is NOT at the anchor gene;
#: both `*_pp_*`-primary classes map to the pseudogene bucket.
TIE_PRIMARY_ELSEWHERE = {
    "cross_pp_pc": "genome_multi_primary_pseudogene",
    "same_pp_pp": "genome_multi_primary_pseudogene",
}

TIDY_COLUMNS = ["sample", "gene_id", "gene_name", "transcript_id", "category",
                "n_reads", "pct_of_union"]
WIDE_COLUMNS = ["sample", "gene_id", "gene_name", "transcript_id",
                "gene_chromosome", "gene_start", "gene_end",
                "n_union", "n_genome_side", "n_txome_side", "n_shared",
                "n_genome_only", "n_txome_only", "n_genome_unique", "n_genome_multi"]


class PartitionError(RuntimeError):
    pass


def load_libraries():
    """Import the read-taxonomy libraries by path, without disturbing sys.path."""
    saved = list(sys.path)
    try:
        for directory in (str(REPO / "code" / "read_taxonomy"),
                          str(REPO / "code" / "common"),
                          str(REPO / "code" / "common" / "ribo_seq_qc")):
            if directory not in sys.path:
                sys.path.insert(0, directory)
        import concordance_lib
        import reach_lib
        import tie_biotype_lib
        return concordance_lib, reach_lib, tie_biotype_lib
    finally:
        sys.path[:] = saved + [p for p in sys.path if p not in saved]


def gene_locus(exon_gene_df, gene_id):
    """The gene's full genomic span, over EVERY annotated isoform.

    Not the selected transcript's own span: nonselected-isoform exons lie outside it.
    """
    base = str(gene_id).split(".", 1)[0]
    rows = exon_gene_df[exon_gene_df["gene_id"].astype(str).str.split(".").str[0] == base]
    if rows.empty:
        raise PartitionError("gene %r has no exon in the annotation" % gene_id)
    chroms = rows["Chromosome"].unique()
    if len(chroms) != 1:
        raise PartitionError("gene %r spans %d chromosomes (%s); this is not handled"
                             % (gene_id, len(chroms), ", ".join(map(str, chroms))))
    return str(chroms[0]), int(rows["Start"].min()), int(rows["End"].max())


def fetch_gene_qnames(genome_bam, chrom, start, end):
    """Every read ID with a reported genome alignment overlapping [start, end) on `chrom`.

    Indexed fetch, secondaries count, strand-agnostic — a multimapper whose primary sits on
    a pseudogene elsewhere must still belong to the gene.
    """
    import pysam

    qnames = set()
    bam = pysam.AlignmentFile(str(genome_bam), "rb")
    try:
        if not bam.has_index():
            raise PartitionError(
                "%s has no index. The gene-anchored denominator is an indexed region "
                "fetch; index it with `samtools index`." % genome_bam)
        for read in bam.fetch(str(chrom), int(start), int(end)):
            if read.is_unmapped or read.is_supplementary:
                continue
            qnames.add(read.query_name)
    finally:
        bam.close()
    return qnames


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
            score = int(read.get_tag("AS")) if read.has_tag("AS") else None
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


def _projection_frame(qnames, txome_population, primary):
    """The frame `transcript_fate_lib._project_match` consumes, for genome-unique reads."""
    rows = []
    for qname in qnames:
        record = primary.get(qname)
        if record is None:
            continue
        chrom, strand, blocks, _nh, _as = record
        _tid, tx_pos, tx_len = txome_population[qname]
        pos5 = blocks[-1][1] - 1 if strand == "-" else blocks[0][0]
        rows.append((qname, tx_pos, tx_len, chrom, strand, pos5, len(blocks),
                     min(b[0] for b in blocks), max(b[1] for b in blocks)))
    return pd.DataFrame(rows, columns=[
        "qname", "tx_pos", "tx_len", "g_chrom", "g_strand", "g_pos5", "g_nblocks",
        "g_min", "g_max"])


def classify_union(libs, annotation, tid, locus, genome_side, txome_side, txome_unique,
                   txome_present, genome_present, genome_unique, primary, records):
    """The chain. Returns (labels: Series qname -> category, detail: raw tie/reach labels)."""
    concordance_lib, reach_lib, tie_biotype_lib = libs
    import transcript_fate_lib

    chrom, start, end = locus
    genome_side = set(genome_side)
    txome_side = set(txome_side)
    union = genome_side | txome_side
    labels = pd.Series(index=sorted(union), dtype=object)
    tie_class, reach_label = {}, {}

    # 1-2. Transcriptome route here, no genome alignment HERE: split absent vs elsewhere.
    for qname in txome_side - genome_side:
        labels[qname] = ("txome_only_genome_elsewhere" if qname in genome_present
                         else "txome_only_genome_absent")

    at_gene = sorted(genome_side)
    uniq = [q for q in at_gene if q in genome_unique]
    multi = [q for q in at_gene if q in genome_present and q not in genome_unique]
    # A post-dedup BAM can keep a multimapper's secondaries after collapsing its primary:
    # no primary to classify, so named, not dropped.
    orphan = [q for q in at_gene if q not in genome_present]
    for qname in orphan:
        labels[qname] = "genome_multi_primary_lost_in_dedup"

    # 3-6. Genome multimappers: anchor-gene orientation first, tie-biotype class second —
    # the raw class names describe the PRIMARY locus, which is often not the anchor gene.
    if multi:
        classes = tie_biotype_lib.categorize_reads(
            {q: records[q] for q in multi if q in records},
            annotation["exon_pr"], annotation["gene_body_pr"])
        for qname in multi:
            # Missing class is NaN, not None; NaN is truthy, so `klass or ""` would write "nan".
            klass = classes.get(qname) if len(classes) else None
            if klass is None or pd.isna(klass):
                klass = None
            tie_class[qname] = klass or ""
            if _in_locus(primary.get(qname), chrom, start, end):
                labels[qname] = (
                    "genome_multi_primary_in_gene_pseudogene_tie" if klass == "cross_pc_pp"
                    else "genome_multi_primary_in_gene_no_pseudogene_tie")
            else:
                labels[qname] = TIE_PRIMARY_ELSEWHERE.get(
                    klass, "genome_multi_primary_elsewhere_other")

    # 7-8. Genome-unique, transcriptome-assigned to THIS transcript: project against it.
    shared_here = [q for q in uniq if q in txome_unique and q in txome_side]
    if shared_here:
        frame = _projection_frame(shared_here, txome_present, primary)
        if not frame.empty:
            match = transcript_fate_lib._project_match(
                concordance_lib, frame, annotation["table"][tid])
            for qname, ok in zip(frame["qname"], match):
                labels[qname] = ("genome_unique_shared_concordant" if ok
                                 else "genome_unique_shared_discordant")

    # 9-10. Genome-unique, transcriptome-present but not confidently on this transcript.
    for qname in uniq:
        if qname in txome_unique and qname not in txome_side:
            labels[qname] = "genome_unique_txome_other_transcript"
        elif qname in txome_present and qname not in txome_unique:
            labels[qname] = "genome_unique_txome_multimapped"

    # 11-16. Genome-unique, absent from every transcriptome primary: the reach classifier.
    absent = [q for q in uniq if q not in txome_present]
    if absent:
        blocks = {q: (primary[q][0], primary[q][1], primary[q][2])
                  for q in absent if q in primary}
        reach = reach_lib.classify_gU_tA(
            absent, blocks, annotation["exon_pr"], annotation["exon_gene_df"],
            annotation["gene_body_pr"], annotation["table"], annotation["gene2tid"],
            annotation["omitted_genes"])
        impossible = [c for c in REACH_IMPOSSIBLE if (reach == c).any()]
        if impossible:
            raise PartitionError(
                "the reach classifier returned %s on a transcriptome-absent population; "
                "the population was built wrong" % ", ".join(impossible))
        for qname in absent:
            reach_label[qname] = reach.get(qname, "")
            labels[qname] = REACH_TO_CATEGORY.get(
                reach_label[qname], "genome_unique_absent_other")

    unlabelled = labels[labels.isna()]
    if len(unlabelled):
        raise PartitionError(
            "%d read(s) fell through the chain, e.g. %s"
            % (len(unlabelled), list(unlabelled.index[:5])))

    detail = {"tie_class": tie_class, "reach_label": reach_label}
    return labels, detail


def load_annotation(libs):
    """Everything the two category systems need, built once for all genes."""
    import pyranges as pr

    concordance_lib, reach_lib, _tie = libs
    payload = concordance_lib.build_transcript_table()
    table = payload["table"]
    exon_gene_df = concordance_lib.build_exon_gene_table()
    return {
        "table": table,
        "base2ver": payload["base2ver"],
        "exon_pr": concordance_lib.load_exon_gene_pr(),
        "exon_gene_df": exon_gene_df,
        "gene_body_pr": pr.PyRanges(
            reach_lib.fc.config.load_all_gene_bodies().reset_index(drop=True)),
        "gene2tid": reach_lib.gene_to_transcript_map(table),
        "omitted_genes": reach_lib.omitted_pc_genes(
            exon_gene_df, set(v["gene_id"] for v in table.values())),
    }


def compute_partition(sample, genome_bam, txome_bam, gene_ids=(), transcript_ids=(),
                      coverage=None, log=lambda _m: None):
    """The whole computation. Returns (wide, tidy, dump)."""
    sys.path.insert(0, str(HERE))
    import transcript_fate_lib

    libs = load_libraries()
    concordance_lib = libs[0]
    annotation = load_annotation(libs)
    table = annotation["table"]

    tids = transcript_fate_lib.resolve_transcripts(
        table, gene_ids, transcript_ids, coverage)
    if not tids:
        raise PartitionError("no transcripts requested")

    log("reading the transcriptome BAM")
    txome_present, txome_unique, _all = concordance_lib.read_txome_primary(
        txome_bam, annotation["base2ver"])
    txome_side = {tid: set() for tid in tids}
    for qname, value in txome_present.items():
        if value[0] in txome_side:
            txome_side[value[0]].add(qname)

    names, genes = _display(table, tids, coverage)
    log("fetching the genome side of each gene")
    loci, genome_side = {}, {}
    for tid in tids:
        loci[tid] = gene_locus(annotation["exon_gene_df"], table[tid]["gene_id"])
        genome_side[tid] = fetch_gene_qnames(genome_bam, *loci[tid])
        log("  %-8s %s:%d-%d  %d genome-side, %d transcriptome-side read ids"
            % (names[tid] or tid, loci[tid][0], loci[tid][1], loci[tid][2],
               len(genome_side[tid]), len(txome_side[tid])))

    targets = set()
    for tid in tids:
        targets |= genome_side[tid] | txome_side[tid]
    log("one genome pass over %d read ids" % len(targets))
    g_present, g_unique, primary, records = collect_genome_state(genome_bam, targets)

    wide_rows, tidy_rows, dump_rows = [], [], []
    for tid in tids:
        labels, detail = classify_union(
            libs, annotation, tid, loci[tid], genome_side[tid], txome_side[tid],
            txome_unique, txome_present, g_present, g_unique, primary, records)
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
            "n_genome_multi": len(genome_side[tid] & (g_present - g_unique))})

        for qname, category in labels.items():
            record = primary.get(qname)
            blocks = record[2] if record else None
            txome_hit = txome_present.get(qname)
            dump_rows.append({
                "sample": sample, "gene_name": names[tid], "transcript_id": tid,
                "read_id": qname, "category": category,
                "tie_biotype_class": detail["tie_class"].get(qname, ""),
                "reach_category": detail["reach_label"].get(qname, ""),
                "genome_primary_nh": record[3] if record else -1,
                "genome_primary_chromosome": record[0] if record else "",
                "genome_primary_strand": record[1] if record else "",
                "genome_primary_start": min(b[0] for b in blocks) if blocks else -1,
                "genome_primary_end": max(b[1] for b in blocks) if blocks else -1,
                "genome_primary_in_gene": _in_locus(record, *loci[tid]),
                "n_genome_loci": len(records.get(qname, ())),
                "on_genome_side": qname in genome_side[tid],
                "txome_primary_transcript": txome_hit[0] if txome_hit else "",
                "txome_unique": qname in txome_unique})

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
