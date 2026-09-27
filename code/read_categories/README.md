# The five read categories (Figures 4–6)

Each library's read-ID union is partitioned into shared, genome-only and
transcriptome-only (TO) reads; shared and genome-only reads are further split by whether
they map uniquely to the genome (`NH == 1`): **SH-U, SH-M, GO-U, GO-M, TO**.
`categories.py` is the one definition — keys, printed abbreviations, colours, the
`MISSING_AS` sentinel, the tied best-scoring rule, and the two transcriptome-presence
rules (any primary for Figure 4 and the 5B locus; selected-transcript for Figure 5A and
Figure 6). Transcriptome presence here is always the post-dedup BAM (RiboFlow_v2's
MAPQ ≥ 10), never the MAPQ ≥ 42 rule of S1 Fig / Figures 2–3.

Everything in this directory is run by `code/make_tables.py`; the stages, in figure
order:

## Figure 4 — the cohort, one row per library (`taxonomy`, `multimap_biotype`, `reach`)

`library_scan.py <analysis>` scans the two post-dedup BAMs of every library (a
`--sample` worker per library) and writes one master table each under
`results/read_categories/` (shipped in `data/read_categories/`):

| analysis | master table | feeds |
|---|---|---|
| `taxonomy` | `taxonomy_all.tsv` — the five category counts and percentages per library | 4A, 4B, and the shared row order of all four panels |
| `tie_biotype` | `multimap_tie_biotype_all.tsv` — protein-coding–pseudogene ties among shared genome-multimapping reads | 4C |
| `reach` | `genome_anchored_reach_all.tsv` — omitted alternative-exon overlap of genome-only unique reads | 4D |

Libraries: `taxonomy_lib.py` (the per-read genome/transcriptome state),
`tie_biotype_lib.py` (the tie test, anchored on the primary), `reach_lib.py` (where
genome-only unique reads fall relative to the selected transcript), `reference_lib.py`
(the APPRIS/GTF annotation tables and gene bodies with biotype).

## Figure 5 (`gene_partition`, `locus`)

* `build_gene_categories.py` — 5A. Runs the ten-category chain
  (`gene_read_partition_lib.py`) over COMT, GAPDH and LRRFIP1, writes the per-read dump,
  and folds it through the route-explicit seven segments (`ROUTE7_SEGMENTS`) into
  `gene_partition_route7.tsv/.json`, asserting the published counts. A gene's union is
  every read with a tied best-scoring genomic placement at the gene or a transcriptome
  primary on its selected transcript.
* `build_locus_data.py` — 5B. The LRRFIP1 locus P-site coverage, split by category, with
  the selected and alternative isoform models (`locus_LRRFIP1.npz/.json`). Reads the
  BAMs directly: the locus needs read lengths, offsets and tied secondaries' positions.

## Figure 6 (`clustering`)

Chain, all under `results/clustering/` with the stem `HeLa.post_dedup`:

1. `read_state.py` — the two BAM passes, once, into `read_state.h5` (~20 min, 250 MB).
   A durable product like the coverage HDF5: reused when present.
2. `build_gene_counts.py` — the five category counts for every APPRIS-selected gene from
   the store (`gene_counts.tsv`), then the clustering input (`gene_counts_filtered.tsv`:
   status ok and `n_union >= MIN_UNION = 100`). `--verify N` cross-checks N genes
   against the untouched `compute_partition` on the BAMs.
3. `ward_cluster.R` — Euclidean distance on the five proportions, Ward's linkage, the
   k = 4 cut (`clusters_k4.tsv`, `cluster_centroids.tsv`, `tree_merge.tsv`); base R.
   Cluster 1 is the most concordant centroid by construction.
4. `cluster_validation.py {pseudogene_counts|omitted_sequence|reference_duplication}` —
   the per-gene annotation tables behind panels 6E–G, from the GTF (and APPRIS) alone.

`panels/plot_read_category_clusters.py` draws Figure 6 from those tables;
`sensitivity_min_union.py` (gitignored extra) sweeps the `MIN_UNION` cutoff.
