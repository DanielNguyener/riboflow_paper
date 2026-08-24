# The shared-coverage HDF5

`schema = "riboflow_paper/shared-coverage/3"`, one file per sample,
`results/coverage/<sample>.shared_coverage.h5`: genome-route and transcriptome-route coverage
on one transcript coordinate. ~25 MB per sample (19,736 transcripts, 70,500,740 positions).

```bash
python code/coverage/coverage_schema.py --validate results/coverage/HeLa.shared_coverage.h5
```

The file stores the four coverage arrays and each transcript's CDS bounds. Everything else a
reader needs (regions, event counts, CDS coverage keys) is a function of those and is
computed on read (`code/coverage/coverage_schema.py`, `CoverageFile`).

## Coordinate

Position *i* of a transcript's slice is position *i* (0-based) of its complete 5′→3′ GTF
exons, concatenated; spliced exon length equals the transcriptome reference length for every
stored transcript (asserted at build time). Genome alignments are projected into this space
through those exons; transcriptome alignments are placed directly. Transcripts are stored in
sorted `transcript_id` order and concatenated; `coverage_offset` is the running sum of
`transcript_len`.

## Root attributes

| attribute                          | value                                                                                                                                                                             |
| ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `schema`, `schema_version`     | `riboflow_paper/shared-coverage/3`, `3`                                                                                                                                       |
| `sample`, `assay`              | e.g.`HeLa`, `ribo` \| `rna`                                                                                                                                                 |
| `routes`                         | `['genome', 'transcriptome']`                                                                                                                                                   |
| `coordinate_system`              | `transcript_5p_to_3p`                                                                                                                                                           |
| `psite_placement`                | `cigar_aware` (the genome route's rule; transcriptome P-sites are `reference_start + offset`)                                                                                 |
| `paper_cds_trim`                 | `15` — nt excluded at each CDS end by every consumer                                                                                                                           |
| `n_transcripts`, `n_positions` | `19736`, `70500740`                                                                                                                                                           |
| `created_utc`                    | build stamp                                                                                                                                                                       |
| `provenance`                     | JSON: every input by name/size/SHA-256, the parameters, both assignment policies, the read-length offsets used, library versions, the code digest and the (path-redacted) command |

## `/transcripts` — one row per transcript

| column                                        | dtype             | meaning                                                                 |
| --------------------------------------------- | ----------------- | ----------------------------------------------------------------------- |
| `transcript_id`, `gene_id`, `gene_name` | fixed-width bytes | identity (GENCODE; the gene name is display metadata)                   |
| `transcript_len`                            | int32             | spliced length                                                          |
| `cds_start`, `cds_end`                    | int32             | CDS as [start, end), stop codon in UTR3;`-1, -1` when there is no CDS |
| `coverage_offset`                           | int64             | first index of the transcript in every`/coverage` array               |

Regions derive as `UTR5 = [0, cds_start)`, `CDS = [cds_start, cds_end)`,
`UTR3 = [cds_end, transcript_len)` (empty parts omitted).

## `/coverage` — four `int32` arrays, each `(n_positions,)`

| dataset              | route         | measure                              |
| -------------------- | ------------- | ------------------------------------ |
| `genome_psite`     | genome        | one count per read at its P-site     |
| `txome_psite`      | transcriptome | one count per read at its P-site     |
| `genome_footprint` | genome        | depth over each read's aligned bases |
| `txome_footprint`  | transcriptome | depth over each read's aligned bases |

Storage: chunk 65,536, gzip 9, shuffle. ~23 MB of the file.

## Derived on read

| quantity                       | definition                                                                                                 |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------- |
| event count                    | sum of a track over the whole transcript (`CoverageFile.event_counts`)                                   |
| coverage state                 | `no_reads_assigned` (count 0) / `reads_outside_requested_slice` (count > 0, slice sum 0) / `covered` |
| CDS coverage key (concordance) | P-site: any count in the untrimmed CDS; footprint: any count in the trimmed interior (`window_sums`)     |

## History

| version | change                                                                                                                                     |
| ------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| 1 → 2  | CIGAR-aware P-site placement became the only rule                                                                                          |
| 2 → 3  | dropped the exon map, region and junction-bin tables, offsets group, provenance sub-group and 25 per-transcript columns (~46 MB → ~25 MB) |
