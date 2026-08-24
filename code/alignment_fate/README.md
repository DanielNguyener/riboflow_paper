# alignment_fate — every read at a gene, on either route (Figure 6)

| | |
|---|---|
| `build_gene_read_partition.py` (+ `gene_read_partition_lib.py`, `transcript_fate_lib.py`) | the union of read IDs at a gene on either route, partitioned by one priority chain; `--dump-reads` writes one row per read |
| `build_gene_partition_data.py` | folds the per-read dump through `panels/plot_gene_read_partition.prepare_route_explicit` into the seven-segment table `gene_partition_route7.{tsv,json}`; refuses counts other than the validated ones |
| `build_locus_data.py` | BAMs + GTF + APPRIS + QC tables → `locus_<GENE>.{npz,json}`: per-base P-site coverage of both routes over the merged exons of the selected and the best-supported alternative isoform |

```bash
python code/make_tables.py --bams DIR --gtf G --appris A --stages gene_partition,locus
```

Figure 6A segments (unit: read IDs; denominator: the union at the gene). "Shared" and
"genome-only" are global BAM-presence terms; "transcriptome-only" is gene-local.

| key | meaning |
|---|---|
| `r7_shared_unique` | shared, genome-unique |
| `r7_shared_multi_pp` | shared, genome-multimapping, protein-coding–pseudogene tie |
| `r7_shared_multi_other` | shared, other genome-multimapping |
| `r7_gonly_unique_omit` | genome-only, genome-unique, on an exon the selected isoform omits |
| `r7_gonly_unique_other` | genome-only, other genome-unique |
| `r7_gonly_multi` | genome-only, genome-multimapping |
| `r7_txonly` | transcriptome-only at the gene |

Validated counts (GSM2100602): COMT 1084/0/34/105/37/26/9 (union 1,295); GAPDH
1057/2207/805/63/49/64/115 (4,360); LRRFIP1 281/40/18/755/27/388/16 (1,525).

Figure 6B locus: chr2:237,627,586–237,781,643 (+), selected `ENST00000308482.14`,
alternative `ENST00000244815.9` (3,619 nt absent from the selected reference; chosen as
the isoform carrying the most genome-only unique reads on non-selected sequence). Introns
are drawn as a constant 90-unit gap. `build_locus_data.py` re-implements the BAM template,
MAPQ ≥ 42 rule and cigar-aware P-site locally, identically to `common/` and `coverage/`.
