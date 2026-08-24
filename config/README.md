# config

| | |
|---|---|
| `panel_manifest.yaml` | every panel (generator, inputs, output) and every figure (panels, composition, raster) — the only place figure numbers appear |
| `cohort_manifest.tsv` (+ `.schema.md`) | the 24 samples and their four BAMs, relative to `--bams` |
| `inputs.example.yaml` | template for `local.yaml` (gitignored): machine-local BAM root, GTF, APPRIS |
| `published_cohort/` | the RiboFlow_v2 per-sample configs and run script that produced the BAMs |
