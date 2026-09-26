#!/usr/bin/env python3
"""Transcript resolution and genome-locus projection shared by the gene-partition stage.

`resolve_transcripts` turns gene/transcript IDs into versioned APPRIS transcript IDs and
`_project_match` tests genome loci for concordance with a transcript; both are consumed by
`gene_read_partition_lib`.
"""
from __future__ import annotations


class FateError(RuntimeError):
    pass

def resolve_transcripts(table, gene_ids=(), transcript_ids=(), coverage=None):
    """Gene and/or transcript IDs -> versioned transcript IDs present in `table`.

    Refuses to guess when a gene maps to more than one candidate, listing them.
    """
    resolved = []
    for tid in transcript_ids:
        if tid in table:
            resolved.append(tid)
            continue
        base = tid.split(".", 1)[0]
        hits = [t for t in table if t.split(".", 1)[0] == base]
        if len(hits) == 1:
            resolved.append(hits[0])
        elif not hits:
            raise FateError("transcript %r is not in the APPRIS transcript table" % tid)
        else:
            raise FateError("transcript %r is ambiguous: %s" % (tid, ", ".join(hits)))

    for gene_id in gene_ids:
        if coverage is not None:
            index = coverage.resolve_gene(gene_id)
            resolved.append(coverage.transcript_info(index)["transcript_id"])
            continue
        base = gene_id.split(".", 1)[0]
        hits = [tid for tid, entry in table.items()
                if str(entry.get("gene_id", "")).split(".", 1)[0] == base]
        if len(hits) == 1:
            resolved.append(hits[0])
        elif not hits:
            raise FateError("gene %r has no transcript in the APPRIS table" % gene_id)
        else:
            raise FateError(
                "gene %r maps to %d transcripts and no unique one resolves it. Pass "
                "--transcript-id to choose; this is not guessed.\n%s"
                % (gene_id, len(hits), "\n".join("    %s" % h for h in sorted(hits))))

    seen, ordered = set(), []
    for tid in resolved:
        if tid not in seen:
            seen.add(tid)
            ordered.append(tid)
    return ordered
