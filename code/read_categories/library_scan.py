#!/usr/bin/env python3
"""The per-library category tables of Figure 4 (definitions and row builders).

The per-read collection runs inside the shared ribo pass
(`code/coverage/build_shared_coverage.py --only ...,categories`): its BAM loops feed
`collect_txome_record` / `collect_genome_record` (transcriptome presence is any primary
alignment; see categories.py), and `stage_rows` computes the three Figure 4 tables from
that single scan:

    taxonomy_all.tsv                 the five category counts per library (4A, 4B)
    multimap_tie_biotype_all.tsv     protein-coding-pseudogene ties (4C)
    genome_anchored_reach_all.tsv    omitted alternative-exon overlap (4D)

Each library stages one row per table; the cohort driver concatenates the staged rows
into the three master tables with `aggregate`.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reach_lib as rl  # noqa: E402
import reference_lib as cl  # noqa: E402
import taxonomy_lib as tl  # noqa: E402
import tie_biotype_lib as tie  # noqa: E402
from categories import MISSING_AS  # noqa: E402
fc = tl.fc

def _out(*parts):
    return fc.output_root().joinpath("read_categories", *parts)

TABLES = {
    "taxonomy": {
        "staging": _out("_staging_taxonomy"),
        "master": _out("taxonomy_all.tsv"),
    },
    "reach": {
        "staging": _out("_staging_reach"),
        "master": _out("genome_anchored_reach_all.tsv"),
    },
    "tie_biotype": {
        "staging": _out("_staging_tie_biotype"),
        "master": _out("multimap_tie_biotype_all.tsv"),
    },
}

CATS = ["cross_pc_pp", "cross_pp_pc", "same_pc_pc", "same_pp_pp"]


# ── the one scan: fed record by record from the shared ribo pass ──────────────
#
# The transcriptome loop collects the presence set; it must be COMPLETE before the genome
# loop starts (the genome collector tests membership). The genome loop then collects, in
# one traversal: the primary and unique read-id sets (taxonomy), every genome locus of
# the multimapping presence-set reads (the tie test), and the primary blocks of the
# unique reads outside the presence set (the reach test). Secondaries are needed: the
# caller must not filter them out before calling the genome collector.

def new_state():
    return {"t_all": set(), "g_all": set(), "g_uniq": set(),
            "tie_records": defaultdict(list), "primary_multi": set(),
            "gUtA_blocks": {}}

def collect_txome_record(state, r):
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        return
    state["t_all"].add(r.query_name)

def collect_genome_record(state, r):
    if r.is_unmapped or r.is_supplementary:
        return
    q = r.query_name
    t_all = state["t_all"]

    if not r.is_secondary:
        state["g_all"].add(q)
        if fc.is_unique_genome_read(r):
            state["g_uniq"].add(q)
            if q not in t_all:
                blocks = r.get_blocks()
                if blocks:
                    strand = "-" if r.is_reverse else "+"
                    state["gUtA_blocks"][q] = (r.reference_name, strand, blocks)

    if q in t_all:
        try:
            nh = r.get_tag("NH")
        except KeyError:
            nh = None
        if not r.is_secondary and nh is not None and nh > 1:
            state["primary_multi"].add(q)
        if nh is None or nh <= 1:
            return
        blocks = r.get_blocks()
        if not blocks:
            return
        pos5 = blocks[0][0] if not r.is_reverse else blocks[-1][1] - 1
        score = int(r.get_tag("AS")) if r.has_tag("AS") else MISSING_AS
        state["tie_records"][q].append((r.reference_name, pos5, score,
                                        bool(r.is_secondary)))

def finish_state(state):
    """Restrict the tie loci to reads whose PRIMARY is multimapping; drop the helpers."""
    tie_records = state.pop("tie_records")
    primary_multi = state.pop("primary_multi")
    state["tie_records"] = {q: recs for q, recs in tie_records.items()
                            if q in primary_multi}
    return state


# ── taxonomy: one row of the five categories per library (4A, 4B) ─────────────

def build_row(sample, counts, n_universe):
    row = {"sample": sample, "n_universe": n_universe}
    for g, t in tl.CELLS:
        n = counts[(g, t)]
        key = tl.cell_key(g, t)
        row[f"n_{key}"] = n
        row[f"pct_{key}"] = 100.0 * n / n_universe if n_universe else float("nan")
    for g in tl.GENOME_STATES:
        row[f"n_genome_{g}"] = sum(counts[(g, t)] for t in tl.TXOME_STATES if (g, t) in counts)
    for t in tl.TXOME_STATES:
        row[f"n_txome_{t}"] = sum(counts[(g, t)] for g in tl.GENOME_STATES if (g, t) in counts)
    return row


def taxonomy_row(sample, state, log=print):
    g_all, g_uniq, t_all = state["g_all"], state["g_uniq"], state["t_all"]
    log(f"  [{sample}] genome mapped={len(g_all):,} unique={len(g_uniq):,} | "
        f"txome mapped={len(t_all):,}")

    inter = len(g_all & t_all)
    smaller = min(len(g_all), len(t_all))
    frac = inter / smaller if smaller else 0.0
    assert frac > 0.05, (
        f"[{sample}] QNAME-namespace mismatch: genome n txome intersection {inter:,} "
        f"= {frac:.1%} of min({len(g_all):,},{len(t_all):,}) — suffix/dedup bug suspected")
    if frac < 0.60:
        log(f"  [{sample}] WARNING: low qname overlap ({frac:.1%} of smaller route) — "
            f"asymmetric read recovery; expect a large `absent` fraction")

    counts = {cell: 0 for cell in tl.CELLS}
    for q in g_all:
        gs = "unique" if q in g_uniq else "multi"
        counts[(gs, "present" if q in t_all else "absent")] += 1
    counts[("absent", "present")] = len(t_all - g_all)

    n_universe = len(g_all | t_all)
    assert sum(counts.values()) == n_universe, \
        f"[{sample}] cell sum {sum(counts.values()):,} != n_universe {n_universe:,}"

    row = build_row(sample, counts, n_universe)
    log(f"[{sample}] n_universe={n_universe:,}  " + "  ".join(
        "%s=%.2f%%" % (tl.cell_key(g, t), row["pct_" + tl.cell_key(g, t)])
        for g, t in tl.CELLS))
    return row, counts


# ── tie_biotype: protein-coding-pseudogene ties per library (4C) ──────────────

def tie_row(sample, state, tax_counts, log=print):
    recs = state["tie_records"]
    expected = tax_counts[("multi", "present")]
    assert len(recs) == expected, (
        f"[{sample}] tie population {len(recs)} disagrees with the taxonomy gM_tP count "
        f"{expected}")

    exon_pr = cl.load_exon_gene_pr()
    gene_pr = cl.gene_body_pr()
    counts, n_reads = tie.categorize(recs, exon_pr, gene_pr)
    row = {"sample": sample, "n_reads": n_reads, **{c: counts[c] for c in CATS},
           "n_qualifying": counts["n_qualifying"]}
    for c in CATS + ["n_qualifying"]:
        row[f"pct_{c}"] = 100.0 * counts.get(c, row.get(c)) / n_reads if n_reads else float("nan")
    row["pct_n_qualifying"] = 100.0 * counts["n_qualifying"] / n_reads if n_reads else float("nan")

    log(f"[{sample}] n_reads={n_reads:,} qualifying={counts['n_qualifying']:,} "
        f"({row['pct_n_qualifying']:.1f}%): "
        f"cross_pc_pp {counts['cross_pc_pp']:,}, cross_pp_pc {counts['cross_pp_pc']:,}, "
        f"same_pc_pc {counts['same_pc_pc']:,}, same_pp_pp {counts['same_pp_pp']:,}")
    return row


# ── reach: omitted alternative-exon overlap per library (4D) ──────────────────

def reach_row(sample, state, tax_counts, log=print):
    n_gUtP = tax_counts[("unique", "present")]
    n_gUtA = tax_counts[("unique", "absent")]
    n_genome_unique = n_gUtP + n_gUtA
    log(f"[{sample}] N_G={n_genome_unique} gU_tP={n_gUtP} gU_tA={n_gUtA}")

    gUtA_qnames = state["g_uniq"] - state["t_all"]
    genome_blocks = state["gUtA_blocks"]
    assert len(gUtA_qnames) == n_gUtA, (
        f"[{sample}] gU_tA set ({len(gUtA_qnames)}) disagrees with the taxonomy count "
        f"({n_gUtA})")

    transcript_payload = cl.build_transcript_table()
    table = transcript_payload["table"]
    selected_genes = set(v["gene_id"] for v in table.values())
    gene2tid = rl.gene_to_transcript_map(table)
    exon_gene_df = cl.build_exon_gene_table()
    exon_gene_pr = cl.load_exon_gene_pr()
    all_gene_body_pr = _all_gene_body_pr()
    omitted_genes = rl.omitted_pc_genes(exon_gene_df, selected_genes)

    log(f"[{sample}] classifying gU_tA reads...")
    labels = rl.classify_gU_tA(gUtA_qnames, genome_blocks, exon_gene_pr, exon_gene_df,
                               all_gene_body_pr, table, gene2tid, omitted_genes)
    counts = labels.value_counts()

    log(f"[{sample}] testing direct overlap with omitted exonic sequence...")
    omitted_index = rl.build_omitted_exon_index(exon_gene_df, table)
    overlap_qnames = rl.omitted_exon_overlap_qnames(gUtA_qnames, genome_blocks, omitted_index)

    row = {"sample": sample, "n_genome_unique": n_genome_unique,
           "n_gU_tP": n_gUtP, "n_gU_tA": n_gUtA}
    for cat in rl.REACH_CATEGORIES:
        row[f"n_{cat}"] = int(counts.get(cat, 0))
    check_sum = n_gUtP + sum(row[f"n_{c}"] for c in rl.REACH_CATEGORIES)
    assert check_sum == n_genome_unique, f"partition does not sum to N_G ({check_sum} != {n_genome_unique})"
    for cat in rl.REACH_CATEGORIES:
        row[f"pct_{cat}"] = 100.0 * row[f"n_{cat}"] / n_gUtA if n_gUtA else float("nan")
    # Figure 4D: NOT a partition slice. The gene-level alternative-exon test of Figure 5A
    # applied to the whole gU_tA population; each read ID counts at most once.
    row["n_omitted_exon_overlap"] = len(overlap_qnames)
    row["pct_omitted_exon_overlap"] = (
        100.0 * len(overlap_qnames) / n_gUtA if n_gUtA else float("nan"))
    log(f"[{sample}] omitted-exon overlap {len(overlap_qnames)} / {n_gUtA} "
        f"({row['pct_omitted_exon_overlap']:.2f}%)")
    return row


def _all_gene_body_pr():
    import pyranges as pr
    df = fc.config.load_all_gene_bodies()
    return pr.PyRanges(df.reset_index(drop=True))


# ── the per-library worker ─────────────────────────────────────────────────────

def _stage(analysis, sample, row):
    staging = TABLES[analysis]["staging"]
    staging.mkdir(parents=True, exist_ok=True)
    dest = staging / f"{sample}.tsv"
    pd.DataFrame([row]).to_csv(dest, sep="\t", index=False)
    print(f"[{sample}] wrote {dest}", flush=True)


def stage_rows(sample, state, log=print):
    """One library: the three table rows from a finished scan state, staged."""
    tax, tax_counts = taxonomy_row(sample, state, log)
    _stage("taxonomy", sample, tax)
    _stage("tie_biotype", sample, tie_row(sample, state, tax_counts, log))
    _stage("reach", sample, reach_row(sample, state, tax_counts, log))


def aggregate(analysis, samples):
    """Concatenate the per-sample staging TSVs into the analysis's master table."""
    spec = TABLES[analysis]
    staging, master = spec["staging"], spec["master"]
    frames = [pd.read_csv(staging / ("%s.tsv" % s), sep="\t")
              for s in samples if (staging / ("%s.tsv" % s)).exists()]
    if not frames:
        print("no staging rows in %s -- nothing aggregated" % staging, file=sys.stderr)
        return None
    table = pd.concat(frames, ignore_index=True).sort_values("sample")
    master.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(master, sep="\t", index=False, lineterminator="\n")
    print("\nwrote %s (%d samples)" % (master, table["sample"].nunique()))
    return master
