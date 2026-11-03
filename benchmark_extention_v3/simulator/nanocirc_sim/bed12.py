"""Writes BED12 fields straight from a Candidate's Block list, so the real
exon/intron structure always makes it into the output."""
from __future__ import annotations

_RGB_BY_TYPE = {
    "ecirna": "27,158,119",
    "eicirna": "217,95,2",
    "cirna": "117,112,179",
    "antisense": "230,171,2",
    "intergenic": "231,41,138",
}


def bed12_fields(candidate) -> list:
    """Strict BED12 fields (12 columns) for one candidate."""
    blocks = candidate.blocks
    circ_start = blocks[0].start
    circ_end = blocks[-1].end
    block_sizes = ",".join(str(b.length) for b in blocks) + ","
    block_starts = ",".join(str(b.start - circ_start) for b in blocks) + ","
    rgb = _RGB_BY_TYPE.get(candidate.circ_type, "0,0,0")
    score = int(min(1000, max(0, round(candidate.tpm * 10))))
    return [
        candidate.chrom,
        str(circ_start),
        str(circ_end),
        candidate.id,
        str(score),
        candidate.strand,
        str(circ_start),
        str(circ_end),
        rgb,
        str(len(blocks)),
        block_sizes,
        block_starts,
    ]


def bed12_line(candidate) -> str:
    return "\t".join(bed12_fields(candidate))
