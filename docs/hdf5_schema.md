# The shared coverage file

One HDF5 file per sample at `results/coverage/<sample>.shared_coverage.h5`, about 25 MB.
It holds coverage from both alignment routes on one shared transcript coordinate
(19,736 transcripts, 70,500,740 positions). The schema name is
`riboflow_paper/shared-coverage/3`.

To check a file:

```bash
python code/coverage/coverage_schema.py --validate results/coverage/HeLa.shared_coverage.h5
```

The file stores only the four coverage arrays and each transcript's CDS bounds.
Everything else a reader needs (regions, event counts, coverage keys) is computed when
the file is read, by `CoverageFile` in `code/coverage/coverage_schema.py`.

## The coordinate

Position i of a transcript is position i (0-based) of its complete 5' to 3' exons joined
together. This spliced length equals the transcriptome reference length for every stored
transcript, and the build asserts that. Genome alignments are projected into this space
through the exons; transcriptome alignments are placed directly. Transcripts are stored
in sorted transcript_id order, and `coverage_offset` says where each one starts in the
arrays.

## Root attributes

| attribute | value |
|---|---|
| `schema`, `schema_version` | `riboflow_paper/shared-coverage/3`, `3` |
| `sample`, `assay` | for example `HeLa`, `ribo` or `rna` |
| `routes` | `['genome', 'transcriptome']` |
| `coordinate_system` | `transcript_5p_to_3p` |
| `psite_placement` | `cigar_aware` for the genome route; transcriptome P-sites are `reference_start + offset` |
| `paper_cds_trim` | `15`, the nt excluded at each CDS end by every consumer |
| `n_transcripts`, `n_positions` | `19736`, `70500740` |
| `created_utc` | build time |
| `provenance` | JSON: every input by name, size and SHA-256, the parameters, both assignment policies, the read-length offsets, library versions, and the command |

## /transcripts, one row per transcript

| column | type | meaning |
|---|---|---|
| `transcript_id`, `gene_id`, `gene_name` | bytes | GENCODE identity |
| `transcript_len` | int32 | spliced length |
| `cds_start`, `cds_end` | int32 | CDS as [start, end); the stop codon is in the UTR3; `-1, -1` when there is no CDS |
| `coverage_offset` | int64 | first index of the transcript in every coverage array |

Regions follow from the CDS bounds: UTR5 is [0, cds_start), CDS is [cds_start, cds_end),
UTR3 is [cds_end, transcript_len).

## /coverage, four int32 arrays

| dataset | route | measure |
|---|---|---|
| `genome_psite` | genome | one count per read at its P-site |
| `txome_psite` | transcriptome | one count per read at its P-site |
| `genome_footprint` | genome | depth over each read's aligned bases |
| `txome_footprint` | transcriptome | depth over each read's aligned bases |

Each array has one entry per position. Storage is chunked (65,536), gzip 9, shuffle.

## Computed on read

| quantity | definition |
|---|---|
| event count | sum of a track over the whole transcript |
| coverage state | `no_reads_assigned` (count 0), `reads_outside_requested_slice` (count > 0 but the slice sums to 0), or `covered` |
| CDS coverage key | P-site: any count in the untrimmed CDS; footprint: any count in the trimmed interior |
