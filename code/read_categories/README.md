# read_categories

Figures 4, 5 and 6. Every read ID that aligns on either route falls into one of five
categories: SH-U, SH-M, GO-U, GO-M and TO (shared, genome-only or transcriptome-only,
and unique or multimapping in the genome). All scripts run via `code/make_tables.py`.

| script | what it does |
|---|---|
| `categories.py` | The five categories defined once: keys, names, colors, shared alignment rules. Everything else imports it. |
| `library_scan.py` | Figure 4. Per-record collectors and row builders for the three tables: taxonomy (the five counts, 4A/4B), tie_biotype (pseudogene ties, 4C), reach (omitted alternative-exon overlap, 4D). Collection runs inside the shared ribo pass (`code/coverage/build_shared_coverage.py --only ...,categories`, the `ribo_pass` stage); the cohort driver aggregates the staged rows. |
| `taxonomy_lib.py` | Which reads are in which BAM, and whether they map uniquely. |
| `tie_biotype_lib.py` | Pseudogene-tie test: a multimapper whose best placements sit on a protein-coding gene and a processed pseudogene. |
| `reach_lib.py` | Where genome-only unique reads fall relative to the selected transcript. |
| `reference_lib.py` | The shared APPRIS and GTF annotation tables. |
| `build_gene_categories.py` | Figure 5A. Classifies every read at COMT, GAPDH and LRRFIP1; writes gene_partition_route7.tsv/.json. |
| `build_locus_data.py` | Figure 5B. LRRFIP1 P-site coverage split by category, with selected and alternative isoform models. |
| `gene_read_partition_lib.py` | Per-gene classification chain shared by Figures 5A and 6. |
| `read_state.py` | Reads both HeLa BAMs once into read_state.h5 (~20 min, 250 MB); reused when present. |
| `build_gene_counts.py` | Figure 6, step 1. Five counts per gene (gene_counts.tsv), then the clustering input (gene_counts_filtered.tsv: >= 100 union reads). |
| `ward_cluster.R` | Figure 6, step 2. Ward clustering of the five proportions, cut at k = 4. Base R. |
| `cluster_validation.py` | Figure 6, step 3. Per-gene tables for panels 6E-6G: pseudogene counts, omitted sequence, reference duplication. |

Figure 4 and 5 tables: `results/read_categories/` (shipped in `data/read_categories/`).
Figure 6 tables: `results/clustering/` (shipped in `data/clustering/`). Panels:
`code/panels/plot_cohort_panels.py`, `plot_gene_categories.py`, `plot_locus_coverage.py`,
`plot_read_category_clusters.py`.

Transcriptome presence here means a primary alignment in the post-dedup BAM (the
MAPQ >= 10 filter RiboFlow_v2 applies). S1 Fig and Figures 2-3 use MAPQ >= 42 instead.
