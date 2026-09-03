# Cohort sample manifest - schema `riboflow_paper/cohort-manifest/1`

`config/cohort_manifest.tsv` lists every sample in the cohort, its identifiers, and where its
four alignment classes and their indexes live. One row per sample.

Which stage reads which BAM

| Stage    | ribo genome | ribo txome | RNA genome | RNA txome |
| -------- | :---------: | :--------: | :--------: | :-------: |
| Figure 3 |     ●     |     ●     |            |          |
| Figure 4 |     ●     |     ●     |     ●     |    ●    |
| Figure 5 |     ●     |     ●     |            |          |
| Figure6  |     ●     |     ●     |            |          |

## Default RiboFlow-genome layout

```
{sample_id}/genome/alignment_ribo/merged/{sample_id}.post_dedup.bam
{sample_id}/transcriptome/alignment_ribo/merged/{sample_id}.transcriptome.post_dedup.bam
{sample_id}/rnaseq/genome/alignment_ribo/merged/{sample_id}.rnaseq.post_dedup.bam
{sample_id}/rnaseq/transcriptome/alignment_ribo/merged/{sample_id}.rnaseq.transcriptome.post_dedup.bam
```

Validation

```bash
python code/coverage/build_cohort_coverage.py --manifest config/cohort_manifest.tsv \
    --bams /path/to/riboflow/output --validate
```

## Selecting samples

```bash
python code/coverage/build_cohort_coverage.py --samples HeLa,A549 ...
python code/coverage/build_cohort_coverage.py --all ...
```
