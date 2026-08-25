"""Shared fixtures for the coverage-pipeline test suite.

Everything here is SYNTHETIC and self-contained: a three-transcript annotation, a matching
pair of indexed genome / transcriptome BAMs, a GENCODE-format GTF, an APPRIS
transcript-lengths file, and a QC master supplying the read-length window and P-site
offsets. No real BAM, no GENCODE download and no APPRIS table is needed, so the whole
suite runs in seconds on a clean checkout -- and it runs the REAL chain

    GTF + APPRIS -> transcript coordinate -> regions -> coverage HDF5 -> concordance

rather than a mock of it.

THE GEOMETRY
------------
Chosen so every coordinate rule the pipeline depends on is exercised and so the expected
vectors can be written down by hand. Transcript coordinates run 5'->3'; the CDS is an
interval overlay inside a longer transcript, so CDS-relative and transcript-relative
positions are never accidentally the same number.

  TX_PLUS   ENSTPLUS0001.1  chr1 '+'  transcript 200 nt = UTR5 30 + CDS 96 + UTR3 74
                            two exons: chr1:970-1060 (90 nt), chr1:2000-2110 (110 nt)
                            CDS exons:  chr1:1000-1060 (60), chr1:2000-2036 (36)
                            so cds_rel 0..59 sit in exon 1 and 60..95 in exon 2
  TX_MINUS  ENSTMINUS002.1  chr2 '-'  transcript 180 nt = UTR5 20 + CDS 90 + UTR3 70
                            one exon: chr2:4930-5110; CDS chr2:5000-5090
                            cds_rel 0 is at genomic 5089 and runs DOWNWARD
  TX_SHORT  ENSTSHORT003.1  chr3 '+'  transcript 80 nt = UTR5 10 + CDS 24 + UTR3 46
                            CDS 24 nt is SHORTER than 2*trim, so its trimmed CDS
                            interior is empty on both routes

The two routes are wired so TX_PLUS and TX_MINUS receive the SAME cds-relative coverage
from both, which makes a perfect-correlation assertion meaningful, while individual reads
are added on one side only to exercise each filter.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
CODE = REPO / "code"
COVERAGE_CODE = CODE / "coverage"

# Import the pipeline modules the ordinary way: they are plain modules in one directory
# that import each other by bare name. No sys.modules pre-registration, no loader tricks.
for _entry in (COVERAGE_CODE, CODE / "panels", CODE / "common",
               CODE / "common" / "ribo_seq_qc"):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))


# ── the suite must not write into the repository ──────────────────────────────
# `results/` is the pipeline's output root, and it is also where a test that forgets
# `tmp_path` will land. That happened: with the documented `RIBOFLOW_PAPER_*` variables
# exported, a no-argument entry-point check ran a real stage and left nine files behind
# (a catalog, two dead pickles and six annotation caches). Nothing failed -- the suite
# passed and the working tree was quietly different afterwards.
#
# This guard makes that failure loud. It records `results/` before any test runs and
# compares afterwards, so the next test that writes there names itself.

def _results_inventory():
    results = REPO / "results"
    if not results.exists():
        return {}
    return {str(p.relative_to(REPO)): p.stat().st_size
            for p in results.rglob("*") if p.is_file()}


@pytest.fixture(scope="session", autouse=True)
def repository_results_are_not_written_to():
    """Fail the session if any test creates, deletes or resizes a file under `results/`."""
    before = _results_inventory()
    yield
    after = _results_inventory()
    created = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    resized = sorted(k for k in set(before) & set(after) if before[k] != after[k])
    assert not (created or removed or resized), (
        "the test suite modified the repository's results/ tree -- tests that need "
        "outputs must use tmp_path.\n  created: %s\n  removed: %s\n  resized: %s"
        % (created, removed, resized))

# ── the synthetic annotation ──────────────────────────────────────────────────
TX_PLUS = "ENSTPLUS0001.1"
TX_MINUS = "ENSTMINUS002.1"
TX_SHORT = "ENSTSHORT003.1"

TRIM = 15
READ_LEN = 30          # in the phase-1 window
READ_LEN2 = 28         # also selected, SAME P-site offset -- mirrors HeLa, whose four
                       # selected lengths all use offset 11. Two lengths sharing an offset
                       # is how a position-deduplicated BAM can legitimately put more than
                       # one P-site on a single base.
OTHER_LEN = 31         # present in the QC master but NOT selected
OFFSET = 12
SAMPLE = "SYN"

# transcript -> (chrom, strand, [CDS (start, end) …] 5'->3', cds_total, trimmed interior)
GEOMETRY = {
    TX_PLUS: ("chr1", "+", [(1000, 1060), (2000, 2036)], 96, 66),
    TX_MINUS: ("chr2", "-", [(5000, 5090)], 90, 60),
    TX_SHORT: ("chr3", "+", [(100, 124)], 24, 0),
}

# transcript -> full exon spans (genomic, 5'->3'), covering UTR5 + CDS + UTR3.
# GENCODE `exon` features span the whole mature transcript, not just the CDS; the
# coordinate builder asserts that the spliced exon length equals the reference length,
# so a CDS-only exon set here would (correctly) be rejected.
EXONS = {
    TX_PLUS: [(970, 1060), (2000, 2110)],       # 90 + 110 = 200
    TX_MINUS: [(4930, 5110)],                   # 180, read 5'->3' from the high end
    TX_SHORT: [(90, 170)],                      # 80
}

# transcript -> (reference name, 5'UTR length). The reference name carries the
# `|UTR5:..|CDS:..|UTR3:..|` header the region parser reads.
TXOME_REFS = {
    TX_PLUS: ("ENSTPLUS0001.1|ENSGPLUS.1|-|-|PLUSGENE-201|PLUSGENE|200|UTR5:1-30|"
              "CDS:31-126|UTR3:127-200|", 30),
    TX_MINUS: ("ENSTMINUS002.1|ENSGMINUS.1|-|-|MINUSGENE-201|MINUSGENE|180|UTR5:1-20|"
               "CDS:21-110|UTR3:111-180|", 20),
    TX_SHORT: ("ENSTSHORT003.1|ENSGSHORT.1|-|-|SHORTGENE-201|SHORTGENE|80|UTR5:1-10|"
               "CDS:11-34|UTR3:35-80|", 10),
}
TXOME_LENGTHS = {TX_PLUS: 200, TX_MINUS: 180, TX_SHORT: 80}
GENOME_CONTIGS = [("chr1", 10000), ("chr2", 10000), ("chr3", 10000)]

GENE_NAMES = {TX_PLUS: "PLUSGENE", TX_MINUS: "MINUSGENE", TX_SHORT: "SHORTGENE"}
# The gene id stored in the coverage file comes from the transcriptome reference HEADER
# (field 2), not from the GTF, so the fixture must use the same string the header carries.
GENE_IDS = {TX_PLUS: "ENSGPLUS.1", TX_MINUS: "ENSGMINUS.1", TX_SHORT: "ENSGSHORT.1"}

#: Transcript-coordinate CDS start for each transcript == its 5'UTR length.
CDS_START = {tid: TXOME_REFS[tid][1] for tid in TXOME_REFS}


def cds_rel_to_genomic(tid, rel):
    """Where cds-relative position `rel` sits in the genome, 5'->3' along the CDS."""
    _chrom, strand, exons, total, _interior = GEOMETRY[tid]
    assert 0 <= rel < total
    walked = 0
    for start, end in exons:
        n = end - start
        if rel < walked + n:
            within = rel - walked
            return start + within if strand == "+" else end - 1 - within
        walked += n
    raise AssertionError("unreachable")


def txome_pos(tid, rel):
    """Where cds-relative position `rel` sits on the transcript coordinate."""
    return CDS_START[tid] + rel


# ── BAM builders ──────────────────────────────────────────────────────────────
def _write_bam(path: Path, contigs, records, index=True):
    """Write a coordinate-sorted, optionally indexed BAM.

    Each record is a dict: ref, pos (0-based), cigar, mapq, reverse, secondary,
    supplementary, unmapped, name, nh.

    `nh` is written as an `NH` tag when present and omitted otherwise, because the two
    references differ: STAR tags every genome alignment with its multiplicity, and the
    genome uniqueness rule (`NH == 1`) reads it, while bowtie2's transcriptome BAMs carry
    no `NH` at all and are judged on MAPQ. A fixture that tagged both would hide that.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": sn, "LN": ln} for sn, ln in contigs]}
    order = {sn: i for i, (sn, _ln) in enumerate(contigs)}
    records = sorted(records, key=lambda r: (order[r["ref"]], r["pos"]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, r in enumerate(records):
            segment = pysam.AlignedSegment(out.header)
            segment.query_name = r.get("name", "read%d" % i)
            segment.reference_id = order[r["ref"]]
            segment.reference_start = r["pos"]
            segment.cigarstring = r["cigar"]
            n = sum(int(x) for x in _cigar_query_lengths(r["cigar"]))
            segment.query_sequence = "A" * n
            segment.query_qualities = pysam.qualitystring_to_array("I" * n)
            segment.mapping_quality = r.get("mapq", 255)
            if r.get("nh") is not None:
                segment.set_tag("NH", int(r["nh"]))
            flag = 0
            if r.get("reverse"):
                flag |= 16
            if r.get("secondary"):
                flag |= 256
            if r.get("supplementary"):
                flag |= 2048
            if r.get("unmapped"):
                flag |= 4
            segment.flag = flag
            out.write(segment)
    if index:
        pysam.index(str(path))
    return path


def _cigar_query_lengths(cigar):
    """Query-consuming lengths in a CIGAR string (M/I/S/=/X, not N/D)."""
    number, out = "", []
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            if char in "MIS=X":
                out.append(number)
            number = ""
    return out


def genome_read(tid, rel, *, mapq=255, nh=1, length=READ_LEN, offset=OFFSET, **kw):
    """A pure-match genome read whose P-site lands on cds-relative `rel` of `tid`.

    For an alignment that is a single run of matches the P-site is `reference_start +
    offset` on '+' and `reference_end - 1 - offset` on '-', so the read is placed
    backwards from the wanted P-site. (Spliced and clipped reads are built explicitly
    where they are needed; for those the two are NOT the same, which is the point.)
    """
    chrom, strand, _exons, _total, _interior = GEOMETRY[tid]
    genomic = cds_rel_to_genomic(tid, rel)
    if strand == "+":
        return dict(ref=chrom, pos=genomic - offset, cigar="%dM" % length, mapq=mapq,
                    nh=nh, **kw)
    return dict(ref=chrom, pos=genomic + offset + 1 - length, cigar="%dM" % length,
                mapq=mapq, nh=nh, reverse=True, **kw)


def txome_read(tid, rel, *, mapq=42, length=READ_LEN, offset=OFFSET, **kw):
    """A transcriptome read whose P-site lands on cds-relative `rel` of `tid`."""
    return dict(ref=TXOME_REFS[tid][0], pos=txome_pos(tid, rel) - offset,
                cigar="%dM" % length, mapq=mapq, **kw)


# ── the fixture cohort ────────────────────────────────────────────────────────
# (transcript, cds_rel of the P-site, read length). The two rel-20 entries are DIFFERENT
# lengths on purpose: in a position-deduplicated BAM two reads of the same length cannot
# share a start position, but two lengths sharing a P-site offset can stack on one base --
# which is why a real sample's vectors top out at the number of selected lengths.
#   rel 20  TX_PLUS   interior hit, two read lengths -> count 2
#   rel 40  TX_PLUS   single interior hit
#   rel  5  TX_PLUS   inside the 5' trim zone -> P-site outside the interior
#   rel 30  TX_MINUS  minus-strand interior hit
#   rel 50  TX_MINUS  minus-strand interior hit
PSITE_PLAN = [(TX_PLUS, 20, READ_LEN), (TX_PLUS, 20, READ_LEN2),
              (TX_PLUS, 40, READ_LEN), (TX_PLUS, 5, READ_LEN),
              (TX_MINUS, 30, READ_LEN), (TX_MINUS, 50, READ_LEN)]

#: The junction-spanning genome read: 10 nt in TX_PLUS exon 1 (cds_rel 50..59) then 20 nt
#: in exon 2 (cds_rel 60..79). Walking 12 positions along the REFERENCE from 1050 gives
#: 1062, which is inside the intron -- a base the read does not cover. Walking 12 ALIGNED
#: READ BASES gives genomic 2002, i.e. cds_rel 62. The pipeline does the latter, so this
#: read IS counted, at cds_rel 62.
JUNCTION_READ = dict(ref="chr1", pos=1050, cigar="10M940N20M", mapq=255, nh=1,
                     name="junction_spanner")
JUNCTION_PSITE_CDS_REL = 62


def build_genome_bam(path: Path):
    records = [genome_read(tid, rel, length=n) for tid, rel, n in PSITE_PLAN]
    # --- reads that must be EXCLUDED, one per filter -------------------------
    # NH:i:2 is what excludes it. The MAPQ is STAR's own value for a 2-locus read, kept
    # so the record stays realistic, but nothing reads it any more.
    records.append(genome_read(TX_PLUS, 25, nh=2, mapq=3, name="excl_multimapper"))
    records.append(genome_read(TX_PLUS, 26, length=OTHER_LEN, name="excl_readlen"))
    records.append(genome_read(TX_PLUS, 27, secondary=True, name="excl_secondary"))
    records.append(genome_read(TX_PLUS, 28, supplementary=True, name="excl_supplementary"))
    records.append(dict(JUNCTION_READ))
    # Entirely intergenic: neither route may see it.
    records.append(dict(ref="chr1", pos=8000, cigar="30M", mapq=255, nh=1,
                        name="intergenic"))
    # TX_SHORT gets coverage on both routes; its trimmed CDS interior is still empty,
    # for being shorter than 2*trim rather than for lack of reads.
    records.append(genome_read(TX_SHORT, 12, name="short_tx"))
    return _write_bam(path, GENOME_CONTIGS, records)


def build_txome_bam(path: Path):
    records = [txome_read(tid, rel, length=n) for tid, rel, n in PSITE_PLAN]
    records.append(txome_read(TX_PLUS, 25, mapq=41, name="excl_low_mapq"))
    records.append(txome_read(TX_PLUS, 26, length=OTHER_LEN, name="excl_readlen"))
    records.append(txome_read(TX_PLUS, 27, secondary=True, name="excl_secondary"))
    # Wholly inside the 5'UTR: a real position on the transcript, outside the CDS.
    records.append(dict(ref=TXOME_REFS[TX_PLUS][0], pos=0, cigar="30M", mapq=42,
                        name="utr5_only"))
    # Wholly inside the 3'UTR.
    records.append(dict(ref=TXOME_REFS[TX_PLUS][0], pos=130, cigar="30M", mapq=42,
                        name="utr3_only"))
    # Straddles the CDS start: transcript 20..50, so 10 nt of UTR5 and 20 nt of CDS.
    records.append(dict(ref=TXOME_REFS[TX_PLUS][0], pos=20, cigar="30M", mapq=42,
                        name="clip5"))
    # Straddles the CDS end: transcript 110..140, so 16 nt of CDS and 14 nt of UTR3.
    records.append(dict(ref=TXOME_REFS[TX_PLUS][0], pos=110, cigar="30M", mapq=42,
                        name="clip3"))
    records.append(txome_read(TX_SHORT, 12, name="short_tx"))
    contigs = [(TXOME_REFS[t][0], TXOME_LENGTHS[t])
               for t in (TX_PLUS, TX_MINUS, TX_SHORT)]
    return _write_bam(path, contigs, records)


# ── annotation builders ───────────────────────────────────────────────────────
def build_synthetic_gtf(path: Path, exons=None, geometry=None, gene_names=None) -> Path:
    """A GENCODE-format GTF for the synthetic cohort.

    Emits `gene`, `exon`, `CDS` and `UTR` features. The `exon` features span the WHOLE
    mature transcript (UTR5 + CDS + UTR3): that is what GENCODE does, and the coordinate
    builder asserts the spliced exon length equals the transcriptome reference length, so
    a CDS-only exon set would be rejected -- correctly.

    GTF coordinates are 1-based inclusive; the tables here hold 0-based half-open pairs,
    so `Start` is written as `start + 1` and `End` as `end`.
    """
    exons = exons or EXONS
    geometry = geometry or GEOMETRY
    gene_names = gene_names or GENE_NAMES
    lines = ["##description: synthetic test annotation, not GENCODE",
             "##provider: riboflow_paper/tests"]

    def attrs(tid):
        return ('gene_id "%s"; transcript_id "%s"; gene_type "protein_coding"; '
                'gene_name "%s"; transcript_type "protein_coding";'
                % (GENE_IDS[tid], tid, gene_names[tid]))

    for tid, spans in exons.items():
        chrom, strand, cds_spans, _total, _interior = geometry[tid]
        lo = min(s for s, _e in spans)
        hi = max(e for _s, e in spans)
        lines.append("\t".join([chrom, "test", "gene", str(lo + 1), str(hi),
                                ".", strand, ".", attrs(tid)]))
        for start, end in spans:
            lines.append("\t".join([chrom, "test", "exon", str(start + 1), str(end),
                                    ".", strand, ".", attrs(tid)]))
        for start, end in cds_spans:
            lines.append("\t".join([chrom, "test", "CDS", str(start + 1), str(end),
                                    ".", strand, "0", attrs(tid)]))
        # UTR features: everything in the exons that is not CDS.
        cds_lo = min(s for s, _e in cds_spans)
        cds_hi = max(e for _s, e in cds_spans)
        for start, end in spans:
            for u_start, u_end in ((start, min(end, cds_lo)), (max(start, cds_hi), end)):
                if u_end > u_start:
                    lines.append("\t".join([chrom, "test", "UTR", str(u_start + 1),
                                            str(u_end), ".", strand, ".", attrs(tid)]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def build_synthetic_appris(path: Path, transcripts=None) -> Path:
    """The APPRIS transcript-lengths file: `<reference name>TAB<length>` per isoform.

    The reference name is exactly the string the synthetic transcriptome BAM uses, so the
    two sides of the cohort cannot drift apart, and it carries the
    `|UTR5:..|CDS:..|UTR3:..|` header the region parser reads.
    """
    transcripts = transcripts or (TX_PLUS, TX_MINUS, TX_SHORT)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join("%s\t%d" % (TXOME_REFS[t][0], TXOME_LENGTHS[t])
                              for t in transcripts) + "\n")
    return path


def build_regions_bed(path: Path, transcripts=None) -> Path:
    """The `actual_regions.bed` cross-check: 0-based half-open, in TRANSCRIPT coordinates,
    with the stop codon assigned to UTR3.

    These transcripts have no annotated stop codon, so the BED CDS is the header CDS
    shifted to 0-based -- i.e. it agrees with the header on both ends.
    """
    transcripts = transcripts or (TX_PLUS, TX_MINUS, TX_SHORT)
    rows = []
    for tid in transcripts:
        name, utr5 = TXOME_REFS[tid][0], CDS_START[tid]
        total_cds = GEOMETRY[tid][3]
        length = TXOME_LENGTHS[tid]
        for label, start, end in (("UTR5", 0, utr5),
                                  ("CDS", utr5, utr5 + total_cds),
                                  ("UTR3", utr5 + total_cds, length)):
            if end > start:
                rows.append("\t".join([name, str(start), str(end), label, "0", "+"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n")
    return path


def build_qc_master(path: Path, sample=SAMPLE, offset=OFFSET):
    """A minimal `readlen_window_qc.csv`: two phase-1 lengths sharing an offset, one
    excluded length, and a second sample that must never be picked up."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"sample": sample, "read_length": READ_LEN, "n_reads": 1000,
         "in_phase1": True, "psite_offset": offset, "frame0_pct": 70.0},
        {"sample": sample, "read_length": READ_LEN2, "n_reads": 800,
         "in_phase1": True, "psite_offset": offset, "frame0_pct": 68.0},
        {"sample": sample, "read_length": OTHER_LEN, "n_reads": 10,
         "in_phase1": False, "psite_offset": offset, "frame0_pct": 40.0},
        {"sample": "OTHER_SAMPLE", "read_length": READ_LEN, "n_reads": 500,
         "in_phase1": True, "psite_offset": 11, "frame0_pct": 65.0},
    ]).to_csv(path, index=False)
    return path


# A second geometry, identical in shape but with TX_PLUS's second exon moved. Used to
# write DIFFERENT content to the SAME GTF path: the annotation change a path comparison
# cannot see and a digest comparison must.
VARIANT_EXONS = dict(EXONS)
VARIANT_EXONS[TX_PLUS] = [(970, 1060), (3000, 3110)]
VARIANT_GEOMETRY = dict(GEOMETRY)
VARIANT_GEOMETRY[TX_PLUS] = ("chr1", "+", [(1000, 1060), (3000, 3036)], 96, 66)


# ── the build fixtures ────────────────────────────────────────────────────────
class Inputs:
    """Every path a coverage build needs, in one place."""

    def __init__(self, root: Path):
        self.root = root
        self.genome_bam = build_genome_bam(root / "syn.genome.bam")
        self.txome_bam = build_txome_bam(root / "syn.txome.bam")
        self.gtf = build_synthetic_gtf(root / "annotation" / "synthetic.gtf")
        self.appris = build_synthetic_appris(root / "annotation" / "appris.tsv")
        self.regions = build_regions_bed(root / "annotation" / "regions.bed")
        self.qc_genome = build_qc_master(root / "qc_genome.csv")
        self.qc_txome = build_qc_master(root / "qc_txome.csv")
        self.output = root / "out"


def build_config(inputs: Inputs, **overrides):
    """An argparse-shaped config for `build_shared_coverage.build()`."""
    import argparse

    values = dict(
        sample=SAMPLE,
        genome_bam=inputs.genome_bam, txome_bam=inputs.txome_bam,
        gtf=inputs.gtf, appris=inputs.appris, regions=inputs.regions,
        qc_genome=inputs.qc_genome, qc_txome=inputs.qc_txome,
        output=inputs.output, trim=TRIM,
        assay="ribo", reference_name="synthetic_v1",
        chunk=1 << 12, gzip_level=1, shuffle=True,
        hash_bams=False, record_input_paths=False,
    )
    values.update(overrides)
    return argparse.Namespace(**values)


@pytest.fixture
def inputs(tmp_path):
    return Inputs(tmp_path)


@pytest.fixture
def coverage_path(inputs):
    """A real coverage HDF5, built from the synthetic BAMs by the real builder."""
    import build_shared_coverage

    path, _report = build_shared_coverage.build(build_config(inputs))
    return path


@pytest.fixture
def coverage(coverage_path):
    """An open, validated `CoverageFile` over the synthetic build."""
    import coverage_schema

    handle = coverage_schema.open_coverage(coverage_path)
    yield handle
    handle.close()


# ── expected vectors, written down by hand ───────────────────────────────────
def expected_psite(tid, trim=TRIM):
    """Hand-computed CDS-interior P-site vector, on the GENOME route.

    The junction spanner contributes at cds_rel 62 because placement walks the READ: see
    JUNCTION_READ. Under a CIGAR-unaware rule it would land in the intron and be lost.
    """
    total = GEOMETRY[tid][3]
    full = np.zeros(total, dtype=np.int32)
    for other, rel, _n in PSITE_PLAN:
        if other == tid:
            full[rel] += 1
    if tid == TX_PLUS:
        full[JUNCTION_PSITE_CDS_REL] += 1
    return full[trim:total - trim] if trim else full


def expected_psite_txome(tid, trim=TRIM):
    """The transcriptome route sees the same plan minus the junction spanner, which has
    no transcriptome counterpart -- a spliced reference contains no introns to span."""
    total = GEOMETRY[tid][3]
    full = np.zeros(total, dtype=np.int32)
    for other, rel, _n in PSITE_PLAN:
        if other == tid:
            full[rel] += 1
    return full[trim:total - trim] if trim else full


def expected_footprint_genome(tid, trim=TRIM):
    """Hand-computed CDS-interior footprint vector on the genome route.

    A footprint contributes to every CDS base it covers; bases outside the interior are
    trimmed away afterwards, so a boundary-spanning read is PARTIALLY counted.
    """
    total = GEOMETRY[tid][3]
    full = np.zeros(total, dtype=np.int64)
    spans = {
        # (cds_rel start, length) of every kept read's footprint, per transcript
        TX_PLUS: [(8, READ_LEN), (8, READ_LEN2),   # the two rel-20 reads: 20 - 12 = 8
                  (28, READ_LEN),                  # rel 40 -> 40 - 12 = 28
                  (0, 23),                         # rel 5 -> starts at -7, clipped to 0
                  (50, 10), (60, 20)],             # the junction spanner, two blocks
        TX_MINUS: [(18, 30),                       # rel 30 on '-' -> mirrored to 18
                   (38, 30)],                      # rel 50 -> 38
    }[tid]
    for start, length in spans:
        full[start:start + length] += 1
    return full[trim:total - trim].astype(np.int32) if trim else full.astype(np.int32)


def expected_footprint_txome(tid, trim=TRIM):
    """Hand-computed CDS-interior footprint vector on the transcriptome route.

    The same plan as the genome route, minus the junction spanner, plus the two
    CDS-boundary-straddling reads clipped to the CDS. The two wholly-UTR reads cover no
    CDS base and so do not appear here (they DO appear in the UTR regions of the stored
    full-transcript vector, which is what makes the CDS slice a slice and not a filter).
    """
    total = GEOMETRY[tid][3]
    full = np.zeros(total, dtype=np.int64)
    spans = {
        TX_PLUS: [(8, READ_LEN), (8, READ_LEN2),
                  (28, READ_LEN),
                  (0, 23),
                  (0, 20),                 # clip5: transcript 20..50 -> cds_rel 0..19
                  (80, 16)],               # clip3: transcript 110..140 -> cds_rel 80..95
        TX_MINUS: [(18, 30), (38, 30)],
    }[tid]
    for start, length in spans:
        full[start:min(start + length, total)] += 1
    return full[trim:total - trim].astype(np.int32) if trim else full.astype(np.int32)


def cds_interior(coverage, tid, signal, trim=TRIM):
    """The stored full-transcript vector, sliced to the trimmed CDS interior."""
    index = coverage.index_of_transcript(tid)
    start, end = CDS_START[tid], CDS_START[tid] + GEOMETRY[tid][3]
    return coverage.get_track(index, signal)[start + trim:end - trim]


# ── the coverage builder, loaded by path ─────────────────────────────────────

@pytest.fixture(scope="session")
def bsc():
    """`build_shared_coverage`, for the accumulator and summation unit tests."""
    import importlib.util
    directory = Path(__file__).resolve().parents[1] / "code" / "coverage"
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
    spec = importlib.util.spec_from_file_location(
        "build_shared_coverage", directory / "build_shared_coverage.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

