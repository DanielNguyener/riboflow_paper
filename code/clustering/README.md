# Clustering genes by how their reads split between the two routes (Figure 6)

**Question.** For every APPRIS-selected gene of one sample, the reads in its union fall
into five fates: present on both routes and genome-unique, present on both and
genome-multimapped, genome-only unique, genome-only multimapped, transcriptome-only. Do
genes fall into a small number of composition types, and what are they?

This directory is the `clustering` stage of `code/make_tables.py` (HeLa, post-dedup). Its
tables ship under `data/clustering/`, and `code/panels/plot_read_fate_clusters.py` draws
Figure 6 from them (`config/panel_manifest.yaml`, panel `fig06`). Everything is written
flat under `results/clustering/` with the stem `<sample>.<label>` (`HeLa.post_dedup`).

## The stage

```
python code/make_tables.py --bams DIR --gtf G --appris A --stages clustering [--into-data]
```

runs, in order (`R=results/clustering`, `S=HeLa.post_dedup`):

| step | program | writes |
|---|---|---|
| 0 | `read_state.py --sample HeLa` (~20 min, once; reused when present) | `$R/$S.read_state.h5`: the two BAM passes `compute_partition` would otherwise repeat for every gene |
| 1 | `build_gene_counts.py --sample HeLa` (~55 min) | `$R/$S.gene_counts.tsv`: one row per gene, `n_union` and the five counts; `$S.gene_mechanisms.tsv`: the two mechanism segments the seven-way split attributes (side table, not used by the figure) |
| 2 | `filter_genes.py` | `$R/$S.gene_counts_filtered.tsv`: `status == ok` and `n_union >= 100` (`filter_genes.MIN_UNION`, the one cutoff) |
| 3 | `ward_cluster.R --k 4` (base R) | `$S.tree.rds`, `$S.tree_merge.tsv` (R's merge matrix and heights), `$S.tree.tsv` (cophenetic correlation), `$S.hclust_sweep.tsv` (k = 2..8: within-SS, between/total, average silhouette, merge height), `$S.clusters_k4.tsv` (per gene: cluster and the five shares in %), `$S.cluster_centroids.tsv` (per cluster: n and the member-mean composition in %) |
| 4 | `validate_cluster_pseudogene_counts.py` | `$S.pseudogene_counts_{genes,by_cluster,tests}.tsv`, `_sources.json`: GENCODE v34 2-way consensus pseudogenes naming each gene as parent (processed / unprocessed / other / unclassified) |
| 5 | `validate_cluster_omitted_sequence.py` | `$S.omitted_sequence_{genes,by_cluster,tests}.tsv`, `_sources.json`: `1 - merged exonic nt of the selected transcript / merged exonic nt of every transcript of the gene` (and the CDS twin) |
| 6 | `validate_cluster_reference_duplication.py` | `$S.reference_duplication_{entries,by_cluster,tests}.tsv`, `_sources.json`: per reference entry, the exonic bases shared with the selected transcript of another gene on the same strand |

Shipped (`--into-data`): the count table, the cut (`clusters_k4`, `cluster_centroids`),
the tree (`tree_merge`) and the three per-gene validation tables, plus the consensus
pseudogene GTF (`data/clustering/gencode.v34.2wayconspseudos.gtf.gz`, GENCODE release 34,
an external input like the housekeeping lists).

## The five metrics

Each gene's union of read IDs is classified by the Figure 5A machinery, unchanged:
`alignment_fate/gene_read_partition_lib.classify_union` (the ten-category chain) and
`panels/plot_gene_read_partition._route7_segment` (the fold to the seven bar segments).
The seven segments are folded once more, keeping where a read went and dropping the
mechanism the seven-way split attributes to it:

| metric | route-7 segments |
|---|---|
| `n_shared_genome_unique` | `r7_shared_unique` |
| `n_shared_genome_multimapped` | `r7_shared_multi_pp` + `r7_shared_multi_other` |
| `n_genome_only_unique` | `r7_gonly_unique_omit` + `r7_gonly_unique_other` |
| `n_genome_only_multimapped` | `r7_gonly_multi` |
| `n_transcriptome_only` | `r7_txonly` |

`build_gene_counts.py` asserts that every segment is folded exactly once and that the
five sum to `n_union` on every row, compares the seven for COMT/GAPDH/LRRFIP1 to the
shipped Figure 5A counts, and writes only the five to the count table. `--verify N` runs
the untouched `compute_partition` on N genes and asserts every segment count matches.

The gene-assignment rule is Figure 5A's own: a read joins a gene's union if a top-score
genome placement -- its primary, or a secondary whose AS equals the primary's -- overlaps
the gene's full multi-isoform genomic span, or if its transcriptome primary is the gene's
selected APPRIS transcript. A lower-scoring secondary never places a read, and neither does
a secondary whose read has no primary record or no AS. A read with tied-best placements at
two genes is in both unions: the question is which genes show a read-assignment pattern,
not which gene owns a multimapper. Genome-unique / multimapped is the primary's NH, as in
Figures 3 and 4. "Shared" means present on the transcriptome route on ANY APPRIS transcript,
not necessarily the gene's selected one (the Figure 5A rule, unchanged). Unions
overlap between genes, so `n_union` does not sum to the library size. Genes are grouped on
`transcript_id`, never on symbol; PAR genes are `status = excluded`; genes with an empty
union are `status = zero_union` (a composition of nothing is undefined, not zero).

## The clustering

Each gene is its five counts divided by `n_union`; the distance is Euclidean between those
shares and the tree is Ward's (`stats::hclust`, `ward.D2`). Rows are sorted by `n_union`
descending then symbol, the order every table keeps. The tree is cut at k = 4: the
average silhouette (on a seeded 4,000-gene subsample) is flat over k = 2..4 and k = 4 is
the partition by dominant fate. Clusters are relabelled so that cluster 1 has the most
concordant centroid (highest shared/unique share). HeLa: 11,864 genes, one cluster of
9,699 concordant genes (83 % shared/unique at the centroid) and three of 811 / 690 / 664
genes each dominated by one discordant fate (shared/multimapped 68 %, genome-only/unique
68 %, genome-only/multimapped 73 %); average silhouette 0.711 at k = 4, the maximum over
k = 2..8.

## Reading Figure 6

Fate labels are abbreviated as in Figures 4-5: SH shared (both routes), GO genome-only, TO
transcriptome-only; U genome-unique (`NH == 1`), M genome-multimapping.

`plot_read_fate_clusters.py` draws, top to bottom: **A** the dendrogram, every subtree
below the cut in the colour of the fate that dominates its cluster's centroid (the Fig 4B
colours); **B** the five shares of every gene, genes in leaf order, with the cluster strip and the three Figure 5A genes (COMT, GAPDH, LRRFIP1) pointed out;
**C** per fate, the % of each gene's union it holds, one box per cluster; **D** the
centroids as stacked bars; **E** the % of each cluster's genes with >= 1 consensus
pseudogene; **F** the % of each gene's exonic sequence its reference transcript omits;
**G** the % of each gene's exonic sequence duplicated by another reference entry on the
same strand. The letters are drawn by the generator: the page is one panel to
`assemble_figures.py`.

## Extra: cutoff sensitivity (not part of the figure)

`sensitivity_min_union.py` (gitignored) reruns the filter and the Ward tree at
`n_union >= 50 / 100 / 200` into `results/clustering/sensitivity/min<N>/` and compares
the k = 4 cuts on the genes every threshold keeps (sizes, dominant fates, adjusted Rand
index and label agreement against the 100 run).
