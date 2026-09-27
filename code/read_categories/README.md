# read_categories — Figures 4, 5 and 6

Every read ID that aligns on either route falls into one of five categories:
**SH-U, SH-M, GO-U, GO-M, TO** (shared / genome-only / transcriptome-only, unique /
multimapping in the genome). All scripts here are run by `code/make_tables.py`.

## Scripts

| script | what it does |
|---|---|
| `categories.py` | The one definition of the five categories: keys, names, colours, and the shared alignment rules. Everything else imports it. |
| `library_scan.py` | Figure 4. Scans each library's two BAMs and writes one table per analysis: `taxonomy` (the five category counts, panels 4A/4B), `tie_biotype` (pseudogene ties, 4C), `reach` (omitted alternative-exon overlap, 4D). |
| `taxonomy_lib.py` | Which reads are in which BAM, and whether they map uniquely. |
| `tie_biotype_lib.py` | The pseudogene-tie test: a multimapper whose best placements sit on a protein-coding gene and a processed pseudogene. |
| `reach_lib.py` | Where genome-only unique reads fall relative to the selected transcript. |
| `reference_lib.py` | The APPRIS/GTF annotation tables the other scripts share. |
| `build_gene_categories.py` | Figure 5A. Classifies every read at COMT, GAPDH and LRRFIP1 and writes `gene_partition_route7.tsv/.json`. |
| `build_locus_data.py` | Figure 5B. LRRFIP1 P-site coverage split by category, with the selected and alternative isoform models (`locus_LRRFIP1.npz/.json`). |
| `gene_read_partition_lib.py` | The per-gene classification chain both Figure 5A and Figure 6 use. |
| `read_state.py` | Reads both HeLa BAMs once into `read_state.h5` (~20 min, 250 MB) so Figure 6 can ask about every gene without re-reading them. Reused when present. |
| `build_gene_counts.py` | Figure 6, step 1. The five category counts for every gene (`gene_counts.tsv`), then the clustering input (`gene_counts_filtered.tsv`: at least 100 union reads). |
| `ward_cluster.R` | Figure 6, step 2. Ward clustering of the five proportions, cut at k = 4 (`clusters_k4.tsv`, `cluster_centroids.tsv`, `tree_merge.tsv`). Base R. |
| `cluster_validation.py` | Figure 6, step 3. Three per-gene annotation tables for panels 6E–G: pseudogene counts, omitted sequence, reference duplication. |

## Outputs

Figure 4 and 5 tables go to `results/read_categories/` (shipped in
`data/read_categories/`); Figure 6 tables go to `results/clustering/` (shipped in
`data/clustering/`). `code/panels/plot_cohort_panels.py`, `plot_gene_categories.py`,
`plot_locus_coverage.py` and `plot_read_category_clusters.py` draw the figures from them.

For these figures a read counts as "on the transcriptome route" when it has a primary
alignment in the post-dedup BAM (RiboFlow_v2's MAPQ ≥ 10 filter), not the MAPQ ≥ 42 rule
used by S1 Fig and Figures 2–3.
