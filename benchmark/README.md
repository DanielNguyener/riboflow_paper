# RiboFlow_v2 performance benchmark

End-to-end computational performance of the
[RiboFlow_v2](https://github.com/ribosomeprofiling/riboflow) pipeline on public data

## Results

Mean of three replicate runs per workload, under a 16-CPU / 64-GB scheduling profile.

| Workload | Input reads | Wall time (s) | CPU time (h) | CPU utilization (%) | Max task RSS (GiB) |
|---|---:|---:|---:|---:|---:|
| 1 matched Ribo-seq/RNA-seq pair  | 14,187,683 |   516 | 0.942 | 41.10 | 30.02 |
| 3 matched Ribo-seq/RNA-seq pairs | 62,536,429 | 1,646 | 3.878 | 53.00 | 30.08 |
| Change                           | ×4.41 | ×3.19 | ×4.12 | +11.90 pp | +0.22 % |


## Details

- **Data.** [GSE269734](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269734) The one-pair run is GSM8325903 (Ribo-seq) with GSM8325891 (RNA-seq). The three-pair run
  adds GSM8325907/GSM8325895 and GSM8325911/GSM8325899. The SRR-level file lists are in
  [`scenarios/`](scenarios).
- **Processing Stages.** Adapter and quality trimming, rRNA/tRNA filtering, genome and transcriptome
  alignment, Ribo-seq UMI extraction and UMIcollapse deduplication, and `.ribo` file generation.
  FastQC and bigWig generation are excluded (see
  [`scenarios/base_gse269734.yaml`](scenarios/base_gse269734.yaml)).
- **Host and scheduling profile** (recorded in [`host.json`](host.json)). AMD Ryzen Threadripper
  3990X, Linux 6.8, Nextflow 26.04.3, against a prebuilt GRCh38 STAR index and the APPRIS v2
  transcriptome. Execution used 16 CPUs (`taskset:0-15`) with a 64-GB Nextflow
  ceiling. Per-process CPU and memory reservations were held constant across both workloads by
  [`config/fixed_resources.config`](config/fixed_resources.config).
- **Measurement.** Timing and memory come from each run's Nextflow task trace by [`config/benchmark.config`](config/benchmark.config)

## What is in this directory

| Path | |
|---|---|
| `runs/<scenario>.rep<N>/trace.txt` | Nextflow task trace|
| `runs/<scenario>.rep<N>/wall_seconds` | wall clock |
| `runs/<scenario>.rep<N>/exit_code` | wrapper exit status |
| `runs/<scenario>.rep<N>/params.yaml` | pipeline configuration for that run |
| `runs/<scenario>.rep<N>/inputs.json` | the FASTQ files consumed, with sizes |
| `runs/<scenario>.rep<N>/manifest.json` | scenario, replicate, run order, description |
| `host.json` | machine, Nextflow version, CPU/memory pinning, reference set |
| `reads.csv` | per-run read counts at each pipeline stage, both routes |
| `benchmark_summary.csv` | the results table above |
| `individual_runs.csv` | one row per run, before averaging |
| `scenarios/` | the two workload definitions and the shared base configuration |
| `config/` | the measurement and fixed-resource Nextflow config layers |
| `summarize_benchmarks.py` | regenerates both CSVs from the traces |
| `run_benchmark.sh` | used to execute the benchmark |

## Producing the results table

```bash
python summarize_benchmarks.py --check
```

Recomputes `benchmark_summary.csv` and `individual_runs.csv` from `runs/*/trace.txt`,
`host.json` and `reads.csv`
