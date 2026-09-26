# riboflow_paper

[![DOI](https://zenodo.org/badge/1321237103.svg)](https://doi.org/10.5281/zenodo.22102431)

Analysis and figure code for the RiboFlow_v2 manuscript: ribosome profiling and matched
RNA-seq from 24 human cell lines, aligned to the genome and to the transcriptome, with the
two alignment routes compared. From RiboFlow_v2 alignments, the code produces the analysis
tables, Figures 2–6 and S1 Fig, the performance benchmark, and S1 Table.

Read processing is done by the separate
[RiboFlow_v2](https://github.com/ribosomeprofiling/riboflow) pipeline, with the
configurations in
[`config/published_cohort/`](config/published_cohort/riboflow_configs/README.md).

## Layout

| | |
|---|---|
| `code/` | `make_tables.py` (BAMs → tables), `make_panels.py` (tables → panels), `assemble_figures.py` (panels → figures), `make_figures.py` (all three); one subdirectory per analysis: `ribo_seq_qc/` (S1 Fig), `coverage/` (Fig 2), `ribo_rna/` + `te_route/` (Fig 3, R), `read_taxonomy/` (Fig 4), `alignment_fate/` (Fig 5), `clustering/` (Fig 6, R), `panels/`, `common/` |
| `config/` | `panel_manifest.yaml` (panels, figures, composition), `cohort_manifest.tsv` (+ `.schema.md`), `inputs.example.yaml`, `published_cohort/` |
| `data/` | shipped analysis tables, one directory per `code/` subdirectory |
| `results/` | regenerated output |
| `figures/` | `panel_references/*.pdf` and `published/{Fig2,Fig3,Fig4,Fig5,Fig6,S1_Fig}.{tif,_plos.pdf}` |
| `docs/` | `methods_te_route.md` (Figure 3 statistics), `hdf5_schema.md` (coverage file format), `numeric_claims.tsv` (every published number and its source), `accessions.tsv` |
| `benchmark/` | performance benchmark: Nextflow traces, scenario definitions and results ([`benchmark/README.md`](benchmark/README.md)) |
| `supporting_information/S1_Table/` | `samples.csv` and its generator |
| `tests/` | test suite |

## Installation

Python 3.9, R ≥ 4 (base only, for `code/te_route/*.R` and `code/clustering/ward_cluster.R`), and the Arial font.

```bash
pip install -r requirements.txt -r requirements-dev.txt
```

## Reproducing the tables

```bash
python code/make_tables.py --bams DIR --gtf GTF --appris APPRIS --all --into-data
```

Writes to `data/` (`results/` without `--into-data`).

## Reproducing the figures

```bash
python code/make_figures.py --all --check
```

Renders the panels from `data/` and writes `figures/published/{Fig2,Fig3,Fig4,Fig5,Fig6,S1_Fig}.{tif,_plos.pdf}`.

Figure 2A/2B also need `results/coverage/HeLa.shared_coverage.h5`, built from the GSM2100602 BAMs
by `code/coverage/build_shared_coverage.py` ([`docs/hdf5_schema.md`](docs/hdf5_schema.md)).

`python code/make_panels.py --all --verify` compares panels with `figures/panel_references/`;

`python benchmark/summarize_benchmarks.py --check` recomputes the performance benchmark
results ([`benchmark/README.md`](benchmark/README.md)).

## External inputs

| input | identifier | use |
|---|---|---|
| RiboFlow_v2 alignments, 24 cell lines | Zenodo [10.5281/zenodo.22083992](https://doi.org/10.5281/zenodo.22083992) | `--bams DIR`; layout in `config/cohort_manifest.tsv` |
| Sequencing data | GEO accessions in [`docs/accessions.tsv`](docs/accessions.tsv) | input to RiboFlow_v2 |
| GENCODE annotation | release 34 (GRCh38) | `--gtf` |
| APPRIS principal-isoform transcriptome | [`references_for_riboflow`](https://github.com/ribosomeprofiling/references_for_riboflow), `transcriptome/human/v2` | RiboFlow_v2 reference; its transcript-lengths table is `--appris` |
| Housekeeping gene lists | HRT Atlas v1.0 (`data/te_route/housekeeping/`) | Figure 3C labels |
| Consensus pseudogenes | GENCODE release 34, Yale-UCSC 2-way consensus set (`data/clustering/gencode.v34.2wayconspseudos.gtf.gz`) | Figure 6E |
| Sample QC table | [ribobaser](https://github.com/CenikLab/ribobaser) | S1 Table generator (not redistributed) |



## Code and data availability

Code: this repository (MIT), archived at Zenodo
[10.5281/zenodo.22102431](https://doi.org/10.5281/zenodo.22102431). Tables: `data/`.
Alignments: Zenodo [10.5281/zenodo.22083992](https://doi.org/10.5281/zenodo.22083992).
Raw reads: GEO ([`docs/accessions.tsv`](docs/accessions.tsv)).

## Citation and license

Citation information will be added at release. License: [`LICENSE`](LICENSE).
