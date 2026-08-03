# Cohort sample manifest — schema `riboflow_paper/cohort-manifest/1`

`config/cohort_manifest.tsv` declares, for every sample in the cohort, its identifiers and the
location of **all four alignment classes** and their indexes. It is the single place the pipeline
learns where a sample's data lives.

## Columns

| Column | Meaning |
|---|---|
| `schema_version` | `riboflow_paper/cohort-manifest/1`. Present on every row so a single row is self-describing; a reader that does not recognise the value must refuse rather than guess. |
| `sample_id` | The key used in every output table's `sample` column and as the BAM directory name. Cell-line names containing spaces use underscores (`Cybrid Cells` → `Cybrid_Cells`). |
| `cell_line` | Display name, spaces intact. **Display metadata only** — never a join key. |
| `ribo_gsm`, `rna_gsm` | GEO accessions for the ribosome-profiling and matched RNA-seq libraries. `ribo_gsm` is what the figures use for axis labels. |
| `ribo_gse` | GEO series for the ribo library. |
| `ribo_genome_bam`, `ribo_genome_bai` | Ribo-seq aligned to the genome (STAR), and its index. |
| `ribo_txome_bam`, `ribo_txome_bai` | Ribo-seq aligned to the transcriptome (Bowtie2), and its index. |
| `rna_genome_bam`, `rna_genome_bai` | RNA-seq aligned to the genome, and its index. |
| `rna_txome_bam`, `rna_txome_bai` | RNA-seq aligned to the transcriptome, and its index. |

Paths may be **absolute**, or **relative to a `--bams` root** given on the command line. Relative is
what ships, so the manifest is portable across machines; only the root changes.

An index column is required because several stages call `pysam`'s region `fetch()`, which cannot
work without one. Naming it explicitly means a missing or stale index is reported by `--validate`
instead of surfacing as an opaque failure mid-run.

## Which stage reads which BAM

| Stage | ribo genome | ribo txome | RNA genome | RNA txome |
|---|:--:|:--:|:--:|:--:|
| shared coverage (Figure 3) | ● | ● | | |
| ribo-vs-RNA route (Figure 4) | ● | ● | ● | ● |
| read taxonomy, concordance, reach, ties (Figure 5) | ● | ● | | |
| alignment fates (Figure 5 E) | ● | ● | | |

## The default RiboFlow-genome layout

The shipped manifest was generated from `docs/accessions.tsv` and these four templates, then every
path was checked to exist:

```
{sample_id}/genome/alignment_ribo/merged/{sample_id}.post_dedup.bam
{sample_id}/transcriptome/alignment_ribo/merged/{sample_id}.transcriptome.post_dedup.bam
{sample_id}/rnaseq/genome/alignment_ribo/merged/{sample_id}.rnaseq.post_dedup.bam
{sample_id}/rnaseq/transcriptome/alignment_ribo/merged/{sample_id}.rnaseq.transcriptome.post_dedup.bam
```

with `.bai` appended for each index — 24 samples, 192 paths, every one checked to exist.

A tree shaped differently needs no code change and no environment variable — edit the manifest.

## Validation

```bash
python code/coverage/build_cohort_coverage.py --manifest config/cohort_manifest.tsv \
    --bams /path/to/riboflow/output --validate
```

Checks, before any compute starts: the `schema_version` on every row is recognised; `sample_id` is
unique and non-empty; every requested BAM and index exists and is non-empty; every BAM is readable
by `pysam` and carries the `@SQ` dictionary the stage needs. It reports **all** problems at once
rather than stopping at the first.

## No implicit whole-cohort run

`build_shared_coverage.py` processes **exactly one sample**. The cohort driver requires an explicit
selection:

```bash
python code/coverage/build_cohort_coverage.py --samples HeLa,A549 ...   # named samples
python code/coverage/build_cohort_coverage.py --all ...                 # every row, opt-in
```


