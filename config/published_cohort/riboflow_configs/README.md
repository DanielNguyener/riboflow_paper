# RiboFlow_v2 configurations for the 24 published samples

One YAML per cell line, plus [`../run_riboflow_cohort.sh`](../run_riboflow_cohort.sh):
the configuration behind the Figures 2-6 alignments. RiboFlow_v2 is not in this
repository. The performance benchmark was a different run, recorded in
[`benchmark/runs/`](../../../benchmark/runs/).

## Placeholders

Absolute paths were replaced with placeholders. 

| Placeholder                                                                                    | Meaning                                                                                                                          |
| ---------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `${REFERENCES_FOR_RIBOFLOW}`                                                                 | a clone of`ribosomeprofiling/references_for_riboflow` (rRNA filter, APPRIS transcriptome, regions BED, transcript-lengths TSV) |
| `${STAR_INDEX_DIR}`                                                                          | the STAR GRCh38 genome index; some YAMLs append`/GRCh38`                                                                       |
| `${FASTQ_DIR}` | input FASTQ root, laid out `$FASTQ_DIR/{ribo,rna}/<GSM>/<SRR>_1.fastq.gz` |                                                                                                                                  |
| `${RIBOFLOW_REPO}`                                                                           | the RiboFlow_v2 checkout containing`main.nf`                                                                                   |
| `${NXF_PROFILE}`                                                                             | Nextflow profile; the script defaults to`lonestar6`                                                                            |

## Notes

RiboFlow's write-time filters are `mapping_quality_cutoff: 0` with
`samtools_filter_arguments: "-F 2052"` on both genome routes, and
`mapping_quality_cutoff: 10` on both transcriptome routes. `-F 2052` keeps secondary
alignments, so the genome BAMs retain multimappers; the transcriptome routes default to
`-F 2308`, which drops them.

The filters this repository applies when reading the BAMs are different, and depend on the
BAM class rather than on ribo versus RNA-seq: a genome read must be primary with `NH == 1`,
a transcriptome read primary with MAPQ >= 42 (`code/common/bam_inputs.py`). `read_taxonomy` and `alignment_fate` apply the same
`NH == 1` test to classify reads

## Per-sample FASTQ counts

| Sample | Ribo FASTQ | RNA FASTQ |
|---|---|---|
| A2780 | 1 | 1 |
| A549 | 1 | 1 |
| BJ | 1 | 1 |
| BRx-142 | 4 | 4 |
| cardiac_fibroblasts | 1 | 1 |
| Cybrid_Cells | 1 | 1 |
| early_neurons | 2 | 2 |
| fibroblast | 1 | 1 |
| H1-hESC | 1 | 1 |
| H9-hESC | 2 | 2 |
| HEK293 | 1 | 1 |
| HEK293T | 1 | 1 |
| HeLa | 1 | 1 |
| hESC | 1 | 1 |
| Huh7 | 1 | 1 |
| LuCaP-PDX | 1 | 1 |
| MCF10A | 1 | 1 |
| megakaryocytes | 6 | 3 |
| neurons | 1 | 1 |
| normal_prostate | 1 | 1 |
| RD-CCL-136 | 3 | 3 |
| U-251 | 1 | 1 |
| U-343 | 1 | 1 |
| U2OS | 1 | 1 |

## Running

```bash
export RIBOFLOW_REPO=/path/to/RiboFlow_v2
export FASTQ_DIR=/path/to/fastqs
export REFERENCES_FOR_RIBOFLOW=/path/to/references_for_riboflow
export STAR_INDEX_DIR=/path/to/star_index
export NXF_PROFILE=              # unset for a local run

cp config/published_cohort/riboflow_configs/*.yaml "$RIBOFLOW_REPO/RiboFlow_YAMLs_cohort/"
bash config/published_cohort/run_riboflow_cohort.sh
```

Finished samples are skipped and `-resume` continues partial ones; `FORCE=1` re-submits
every sample using the Nextflow work cache.
