"""The exon/intron block structure of a circRNA, as an ordered list of Blocks."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Block:
    start: int  # 0-based genomic coordinate
    end: int  # half-open
    kind: str  # "exon" | "intron_retained" | "region"

    @property
    def length(self) -> int:
        return self.end - self.start


def total_length(blocks: list[Block]) -> int:
    return sum(b.length for b in blocks)


def validate_blocks(blocks: list[Block]) -> None:
    if not blocks:
        raise ValueError("a circRNA must have at least one block")
    for a, b in zip(blocks, blocks[1:]):
        if b.start < a.end:
            raise ValueError(f"blocks overlap or are unsorted: {a} then {b}")


def block_signature(blocks: list[Block]) -> tuple:
    """Hashable fingerprint for deduplicating isoforms from the same locus."""
    return tuple((b.start, b.end, b.kind) for b in blocks)


def span(blocks: list[Block]) -> tuple[int, int]:
    return blocks[0].start, blocks[-1].end


def with_shifted_boundary(blocks: list[Block], new_start: int, new_end: int):
    """Replace the first block's start and last block's end. Returns None if
    that would invert or overlap the blocks."""
    if not blocks:
        return None
    new_blocks = list(blocks)
    if len(new_blocks) == 1:
        if new_start >= new_end:
            return None
        new_blocks[0] = Block(new_start, new_end, new_blocks[0].kind)
    else:
        first, last = new_blocks[0], new_blocks[-1]
        if new_start >= first.end or new_end <= last.start:
            return None
        new_blocks[0] = Block(new_start, first.end, first.kind)
        new_blocks[-1] = Block(last.start, new_end, last.kind)
    try:
        validate_blocks(new_blocks)
    except ValueError:
        return None
    return new_blocks


def length_in_bounds(blocks: list[Block], min_length: int, max_length: int) -> bool:
    """True if the circRNA's total length is within [min_length, max_length]."""
    length = total_length(blocks)
    return min_length <= length <= max_length


def assemble_sequence(genome, chrom: str, blocks: list[Block], strand: str) -> str:
    """Concatenate block sequences in transcription order."""
    ordered = blocks if strand == "+" else list(reversed(blocks))
    return "".join(genome.fetch(chrom, b.start, b.end, strand) for b in ordered)
