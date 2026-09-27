# read_categories

Figures 4, 5 and 6. Every read ID that aligns on either route falls into one of five
categories: SH-U, SH-M, GO-U, GO-M and TO (shared, genome-only or transcriptome-only,
and unique or multimapping in the genome). All scripts here are run by
`code/make_tables.py`.

| script | what it does |
|---|---|
| `categories.py` | Defines the five categories once: keys, names, colors and the shared alignment rules. Everything else imports it. |
| `library_scan.py` | Figure 4. Reads each library's two BAMs once and writes all three tables from that one scan: taxonomy (the five counts, panels 4A and 4B), tie_biotype (pseudogene ties, 4C) and reach (omitted alternative-exon overlap, 4D). |
| `taxonomy_lib.py` | Which reads are in which BAM, and whether they map uniquely. |
| `tie_biotype_lib.py` | The pseudogene-tie test: a multimapper whose best placements sit on a protein-coding gene and a processed pseudogene. |
| `reach_lib.py` | Where genome-only unique reads fall relative to the selected transcript. |
| `reference_lib.py` | The APPRIS and GTF annotation tables the other scripts share. |
| `build_gene_categories.py` | Figure 5A. Classifies every read at COMT, GAPDH and LRRFIP1 and writes gene_partition_route7.tsv and .json. |
| `build_locus_data.py` | Figure 5B. LRRFIP1 P-site coverage split by category, with the selected and alternative isoform models. |
| `gene_read_partition_lib.py` | The per-gene classification chain that Figure 5A and Figure 6 share. |
| `read_state.py` | Reads both HeLa BAMs once into read_state.h5 (about 20 minutes, 250 MB) so Figure 6 can ask about every gene without re-reading them. Reused when present. |
| `build_gene_counts.py` | Figure 6, step 1. The five counts for every gene (gene_counts.tsv), then the clustering input (gene_counts_filtered.tsv: at least 100 union reads). |
| `ward_cluster.R` | Figure 6, step 2. Ward clustering of the five proportions, cut at k = 4. Base R. |
| `cluster_validation.py` | Figure 6, step 3. Three per-gene annotation tables for panels 6E to 6G: pseudogene counts, omitted sequence, reference duplication. |

Figure 4 and 5 tables go to `results/read_categories/` (shipped in
`data/read_categories/`). Figure 6 tables go to `results/clustering/` (shipped in
`data/clustering/`). The panels are drawn by `code/panels/plot_cohort_panels.py`,
`plot_gene_categories.py`, `plot_locus_coverage.py` and `plot_read_category_clusters.py`.

For these figures a read counts as present on the transcriptome route when it has a
primary alignment in the post-dedup BAM (the MAPQ >= 10 filter RiboFlow_v2 applies).
S1 Fig and Figures 2 and 3 use MAPQ >= 42 instead.
