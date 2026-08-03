# The shared-coverage HDF5

`schema = "riboflow_paper/shared-coverage/2"`, one file per sample,
`results/coverage/<sample>.shared_coverage.h5`. Genome-route and transcriptome-route
coverage on one transcript coordinate, with the annotation and provenance needed to read it
without any other file.

Everything below is measured against the 
HeLa file on disk (19,736 transcripts, 46.2 MB).

```bash
python code/coverage/coverage_schema.py --validate results/coverage/HeLa.shared_coverage.h5
```
---

## Root attributes

| attribute | value |
|---|---|
| `schema` | `riboflow_paper/shared-coverage/2` |
| `sample` | e.g. `HeLa` |
| `assay` | `ribo` \| `rna` |
| `routes` | `['genome' 'transcriptome']` |
| `reference_name` | e.g. `appris_human_v2_selected` |
| `coordinate_system` | `transcript_5p_to_3p` |
| `exon_source` | `gencode_exon_features` |
| `stop_codon_assignment` | `utr3` |
| `psite_placement` | `cigar_aware` — the **genome** route's rule. Transcriptome P-sites are `reference_start + offset` (no introns on a transcript reference; Bowtie2 `--norc`) |
| `paper_cds_trim` | `15` |
| `n_transcripts` | `19736` |
| `n_positions` | `70500740` — Σ `transcript_len` |
| `created_utc`, `generator` | build stamp |

> Position *i* of a transcript's slice is position *i* (0-based) of its complete 5′→3′ GTF
> exons, concatenated. Spliced exon length equals the transcriptome reference length for
> every stored transcript (asserted at build time). Genome alignments are projected into this
> space through those exons; transcriptome alignments are placed directly. Region membership
> (UTR5 / CDS / UTR3) is an interval overlay in `/regions`.

## `/coverage` — four `int32` arrays, each `(n_positions,)`

| dataset | route | measure |
|---|---|---|
| `genome_psite` | genome | psite |
| `txome_psite` | transcriptome | psite |
| `genome_footprint` | genome | footprint |
| `txome_footprint` | transcriptome | footprint |

Each dataset carries its own `route`, `measure` and `assay` attributes — no name-parsing
needed to know which aligner produced a number.

A transcript's slice is `[coverage_offset, coverage_offset + transcript_len)`. Storage:
`chunk = 65536`, `gzip_level = 9`, `shuffle = True` (group attributes).

## `/transcripts` — 19,736 rows, sorted by `transcript_id`

| group | columns | dtype |
|---|---|---|
| identity | `transcript_id`, `gene_id`, `transcript_name`, `gene_name`, `chrom`, `strand` | vlen utf-8 |
| geometry | `transcript_len`, `n_exons`, `cds_len_gtf`, `n_cds_exons` | `int64` |
| offsets | `coverage_offset`, `exon_offset`, `region_offset`, `bin_offset` | `int64` |
| facts | `in_transcriptome_reference`, `length_filtered`, `has_annotated_stop`, `cds_divisible_by_3` | `bool` |
| coverage state | `hist_cds_{genome,txome}_{psite,footprint}_key` | `bool` |
| counts | `n_{genome,txome}_psite_events`, `n_{genome,txome}_footprint_bases` | `int64` |

On the published reference: `length_filtered` is false for **17,071** transcripts,
`has_annotated_stop` false for **5**, `cds_divisible_by_3` false for **25**.

`length_filtered` only asks whether a transcript's UTRs are long enough for a metagene
window — it is a QC convenience, not a universe. Figure 4's universe is every APPRIS
transcript with a canonical CDS that the transcriptome reference also names, not the
17,071 that pass `length_filtered`.

There is **no `appris_category`**: no formal APPRIS rank reaches this pipeline.

### Why four coverage states instead of one `has_coverage` flag

A single boolean would merge states that mean different things.
`coverage_schema.describe_coverage_state` names which one holds:

| state | meaning |
|---|---|
| `covered` | reads assigned, and the requested slice is non-zero |
| `no_reads_assigned` | no read was assigned to this transcript on this route |
| `reads_outside_requested_slice` | reads were assigned, but all fell outside the slice asked for |
| `slice_all_zero` | the slice is in range and genuinely zero |

`hist_cds_*_key` is asymmetric by design: the P-site accumulator keys a transcript with any
in-CDS read even if every one fell in the trimmed zone, while the footprint accumulator only
keys non-zero vectors. Reproducing published per-transcript row counts needs this flag rather
than re-deriving it.

## `/exons` — 207,370 rows

`transcript_index`, `exon_index`, `chrom`, `g_start`, `g_end`, `tx_start`, `tx_end`.

`exon_index` is 0-based in 5′→3′ (strand) order; `g_start`/`g_end` are genomic 0-based
half-open. A transcript position round-trips to the genome as:

```
g = e.g_start + off        if strand == "+"
g = e.g_end - 1 - off      if strand == "-"
```

Build-time assertions per transcript: exons are never region-split, `tx_start`/`tx_end` are
contiguous and gapless, `Σ(g_end − g_start) == transcript_len`.

`/genomic_segments` is reserved for a future region-split table, so `/exons` can keep meaning
*whole exons*.

## `/regions` — 58,512 rows

`transcript_index`, `label`, `source`, plus three coordinate pairs:
`raw_header_start_1based` / `raw_header_end_1based` (1-based inclusive, stop **inside** the
CDS), `raw_bed_start` / `raw_bed_end` (0-based half-open, stop in UTR3), and the normalized
`start` / `end` used downstream. Absent values are `-1`; absent regions are omitted rather
than stored as a zero-length interval at position 0.

Labels are only `UTR5`, `CDS`, `UTR3`, and tile `[0, transcript_len)` exactly.

## `/ribo_region_bins` — 95,256 rows

`transcript_index`, `label`, `ribopy_alias`, `start`, `end`. A **derived** five-way overlay,
not annotation. Group attributes record how it was made:

```
algorithm          ribopy_get_extended_boundaries
left_span          35
right_span         10
start_site_source  header_cds_start_minus_1
stop_site_source   header_cds_end_stop_inclusive
parameter_source   cli_default | cli | <path to .ribo> sha256=<...>
```

## `/offsets/{genome,transcriptome}`

`read_length` and `psite_offset`, both `int32` — the selected window and offsets actually
used for this sample and route.

## `/provenance`

A JSON blob under the `json` attribute. Parameters that change the numbers are also promoted
to attributes for inspection without parsing: `paper_cds_trim`, `genome_min_mapq` (4),
`txome_min_mapq` (42), `left_span`, `right_span`, `psite_placement`,
`appris_principal_ranks_consumed`, plus `command`, `paths_redacted`, `code_version`.

`psite_placement` describes the **genome** route only — it's the sole route where the rule
can matter, since transcriptome P-sites are `reference_start + offset`.

The blob also records: every consumed input **by content** — name, byte count, SHA-256 (plus
`.bai` digests and `@SQ` dictionaries for BAMs) — both assignment policies with their five
counts each, library versions, and the SHA-256 of every module that defines a number
(`build_shared_coverage.py`, `coverage_schema.py`, `psite_placement.py`,
`transcript_coords.py`, `transcript_regions.py`).

**Identity is by content, not by path.** Absolute paths are recorded only if
`--record-input-paths` is passed — a path names the machine, the project layout and often the
operator. Reuse is fail-closed: `--reuse-existing` refuses on any manifest difference; there
is no force-reuse flag.

## Size

| | |
|---|---|
| full coordinate | 70,500,740 positions, 1,128 MB uncompressed |
| CDS interior only | 33,232,821 positions, 532 MB |
| `/coverage` on disk | 23.2 MB |
| `/exons` · `/transcripts` · `/regions` · `/ribo_region_bins` | 2.8 · 1.3 · 1.0 · 1.0 MB |
| sum of datasets | 29.3 MB |
| **file on disk** | **46.2 MB** |

The 29.3 → 46.2 MB gap is the HDF5 variable-length string heap (~700,000 short strings across
the metadata tables). Encoding them as integer codes against a lookup table would recover
most of it (~30 MB), but hasn't been done — it would make the file less self-describing to a
plain `h5py` reader. Noted here so the trade-off is visible rather than accidental.

Storage defaults are gzip 9, chunk 65,536 `int32`, shuffle on — 23.2 MB above (2.06% of raw).
Shuffle helps because three of every four bytes of an `int32` count are zero.
