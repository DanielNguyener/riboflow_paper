# data

The shipped analysis tables, one directory per `code/` function, at the same relative path a
regenerated table takes under `results/` (`make_tables.py --into-data` is the prefix swap).

| | built by | read by |
|---|---|---|
| `annotation/orf_catalog.tsv` | `orf_catalog` | `te_stats` (gene names) |
| `ribo_seq_qc/{genome,transcriptome}/tables/*.csv`, `offsets/` | `qc`, `offsets` | Figure 2; every stage that needs the read-length window |
| `coverage/concordance/*` | `concordance` | Figure 3C–D |
| `ribo_rna/counts/*.csv` | `te_counts` | `te_normalize` |
| `te_route/normalized/size_factors.csv` | `te_normalize` | — (one factor per assay and library) |
| `te_route/tables/*.tsv` | `te_stats` | Figure 4 |
| `te_route/housekeeping/*.csv` | HRT Atlas v1.0 (external) | Figure 4C labels |
| `read_taxonomy/**/*.tsv` | `taxonomy`, `reach`, `multimap_biotype` | Figure 5 |
| `alignment_fate/gene_partition_route7.*`, `locus_LRRFIP1.*` | `gene_partition`, `locus` | Figure 6 |

Not shipped: `results/coverage/<sample>.shared_coverage.h5` (~25 MB each; Figure 3A/B need
HeLa's) and the per-read partition dump (`results/alignment_fate/`).
