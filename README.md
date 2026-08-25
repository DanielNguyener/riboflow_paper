# riboflow_paper

Analysis and figure code for the RiboFlow_v2 manuscript: ribosome profiling and matched
RNA-seq from 24 human cell lines, aligned to the genome and to the transcriptome, with the
two alignment routes compared. From RiboFlow_v2 alignments, the code produces the analysis
tables, Figures 2–6, Table 1, and S1 Table.

Read processing is done by the separate
[RiboFlow_v2](https://github.com/ribosomeprofiling/riboflow) pipeline. The alignments used
here were produced at commit
[`e5e041c6`](https://github.com/ribosomeprofiling/riboflow/commit/e5e041c6fa842c27fabe46d2ca87d8aff3696874)
with the configurations in
[`config/published_cohort/`](config/published_cohort/riboflow_configs/README.md).

## Layout

| | |
|---|---|
| `code/` | `make_tables.py` (BAMs → tables), `make_panels.py` (tables → panels), `assemble_figures.py` (panels → figures), `make_figures.py` (all three); one subdirectory per analysis: `ribo_seq_qc/` (Fig 2), `coverage/` (Fig 3), `ribo_rna/` + `te_route/` (Fig 4, R), `read_taxonomy/` (Fig 5), `alignment_fate/` (Fig 6), `panels/`, `common/` |
| `config/` | `panel_manifest.yaml` (panels, figures, composition), `cohort_manifest.tsv` (+ `.schema.md`), `inputs.example.yaml`, `published_cohort/` |
| `data/` | shipped analysis tables, one directory per `code/` subdirectory |
| `results/` | regenerated output (not versioned, except `coverage/coverage_checksums.tsv`) |
| `figures/` | `panel_references/*.pdf` and `published/Fig{2..6}.{tif,_plos.pdf}` |
| `docs/` | `methods_te_route.md` (Figure 4 statistics), `hdf5_schema.md` (coverage file format), `numeric_claims.tsv` (every published number and its source), `accessions.tsv` |
| `benchmark/` | Nextflow traces and scenario definitions behind Table 1 |
| `supporting_information/S1_Table/` | `samples.csv` and its generator |
| `tests/` | test suite |

## Installation

Python 3.9, R ≥ 4 (base only, for `code/te_route/*.R`), and the Arial font.

```bash
pip install -r requirements.txt -r requirements-dev.txt
```

## Reproducing the tables

```bash
python code/make_tables.py --bams DIR --gtf GTF --appris APPRIS --all --into-data
```

Writes to `data/` (`results/` without `--into-data`). `--validate` lists what would run;
`--stages` selects stages. Inputs may also be given by `RIBOFLOW_PAPER_{BAMS,GTF,APPRIS}` or
`config/local.yaml` (see `config/inputs.example.yaml`).

## Reproducing the figures

```bash
python code/make_figures.py --all --check
```

Renders the panels from `data/` and writes `figures/published/Fig{2..6}.{tif,_plos.pdf}`.
Figure 3A/3B also need `results/coverage/HeLa.shared_coverage.h5`, built from the GSM2100602 BAMs
by `code/coverage/build_shared_coverage.py` ([`docs/hdf5_schema.md`](docs/hdf5_schema.md)).
`python code/make_panels.py --all --verify` compares panels with `figures/panel_references/`;
`python benchmark/summarize_benchmarks.py --check` recomputes Table 1.

## External inputs

| input | identifier | use |
|---|---|---|
| RiboFlow_v2 alignments, 24 cell lines | Zenodo [10.5281/zenodo.22083992](https://doi.org/10.5281/zenodo.22083992) | `--bams DIR`; layout in `config/cohort_manifest.tsv` |
| Sequencing data | GEO accessions in [`docs/accessions.tsv`](docs/accessions.tsv) | input to RiboFlow_v2 |
| GENCODE annotation | release 34 (GRCh38) | `--gtf` |
| APPRIS principal-isoform transcriptome | [`references_for_riboflow`](https://github.com/ribosomeprofiling/references_for_riboflow), `transcriptome/human/v2` | RiboFlow_v2 reference; its transcript-lengths table is `--appris` (not redistributed) |
| Housekeeping gene lists | HRT Atlas v1.0 (`data/te_route/housekeeping/`) | Figure 4C labels |
| Sample QC table | [ribobaser](https://github.com/CenikLab/ribobaser) | S1 Table generator (not redistributed) |

## Tests

```bash
python -m pytest tests -q
```

Tests needing the coverage HDF5 (`RIBOFLOW_PAPER_COVERAGE_H5`), BAMs
(`RIBOFLOW_PAPER_BAMS`), `Rscript`, or the S1 Table inputs (`RIBOBASER_RDA`,
`RIBOFLOW_PAPER_S1_XLSX`) are skipped when those are absent.

## Code and data availability

Code: this repository (MIT). Tables: `data/`. Alignments: Zenodo 10.5281/zenodo.22083992.
Raw reads: GEO ([`docs/accessions.tsv`](docs/accessions.tsv)).

## Citation and license

Citation information will be added at release. License: [`LICENSE`](LICENSE).
