"""Splice-site (donor/acceptor) sampling and boundary search.

Donor/acceptor side depends on strand: on `+`, end-side is donor and
start-side is acceptor; on `-` it's swapped, read via reverse complement."""
from __future__ import annotations

from .genome import reverse_complement

CANONICAL_MOTIFS = {
    "GT-AG": ("GT", "AG"),
    "GC-AG": ("GC", "AG"),
    "AT-AC": ("AT", "AC"),
}

SPLICE_TYPES = ["GT-AG", "GC-AG", "AT-AC", "other"]


def classify_splice_type(donor: str, acceptor: str) -> str:
    for name, (d, a) in CANONICAL_MOTIFS.items():
        if donor == d and acceptor == a:
            return name
    return "other"


def is_canonical(splice_type: str) -> bool:
    return splice_type != "other"


def sample_desired_splice_type(rng, distribution: dict) -> str:
    types = SPLICE_TYPES
    probs = [distribution.get(t, 0.0) for t in types]
    total = sum(probs)
    if total <= 0:
        return "GT-AG"
    probs = [p / total for p in probs]
    return rng.choice(types, p=probs)


def _end_side_motif(genome, chrom: str, pos: int, strand: str) -> str:
    """The 2bp flanking the circRNA just after `pos`, transcript-oriented."""
    motif = genome.fetch(chrom, pos, pos + 2, "+")
    return reverse_complement(motif) if strand == "-" else motif


def _start_side_motif(genome, chrom: str, pos: int, strand: str) -> str:
    """The 2bp flanking the circRNA just before `pos`, transcript-oriented."""
    motif = genome.fetch(chrom, pos - 2, pos, "+")
    return reverse_complement(motif) if strand == "-" else motif


def donor_acceptor(genome, chrom: str, start: int, end: int, strand: str):
    """Return (donor, acceptor) dinucleotides for a circRNA at [start, end) on strand."""
    end_motif = _end_side_motif(genome, chrom, end, strand)
    start_motif = _start_side_motif(genome, chrom, start, strand)
    if strand == "-":
        return start_motif, end_motif  # donor=start-side, acceptor=end-side
    return end_motif, start_motif  # donor=end-side, acceptor=start-side


def _shifts(center: int, max_adjustment: int):
    """Deterministic shift order: 0, -1, +1, -2, +2, ... (nearest first, stable tie-break)."""
    yield center
    for d in range(1, max_adjustment + 1):
        yield center - d
        yield center + d


def _search_end_side(genome, chrom: str, end: int, strand: str, target_motif: str, max_adjustment: int) -> int:
    window_start = end - max_adjustment
    window = genome.fetch(chrom, window_start, end + max_adjustment + 2, "+")
    for e in _shifts(end, max_adjustment):
        i = e - window_start
        if 0 <= i and i + 2 <= len(window):
            motif = window[i:i + 2]
            if strand == "-":
                motif = reverse_complement(motif)
            if motif == target_motif:
                return e
    return end


def _search_start_side(genome, chrom: str, start: int, strand: str, target_motif: str, max_adjustment: int) -> int:
    window_start = start - max_adjustment - 2
    window = genome.fetch(chrom, window_start, start + max_adjustment, "+")
    for s in _shifts(start, max_adjustment):
        i = (s - 2) - window_start
        if 0 <= i and i + 2 <= len(window):
            motif = window[i:i + 2]
            if strand == "-":
                motif = reverse_complement(motif)
            if motif == target_motif:
                return s
    return start


def find_boundary_for_desired_type(genome, chrom: str, start: int, end: int, strand: str,
                                    desired_type: str, max_adjustment: int = 10):
    """Independently shift `start` and `end` within +/-max_adjustment to
    realize `desired_type`. Returns (new_start, new_end, achieved_splice_type)."""
    target = CANONICAL_MOTIFS.get(desired_type)
    if target is None:
        donor, acceptor = donor_acceptor(genome, chrom, start, end, strand)
        return start, end, classify_splice_type(donor, acceptor)

    target_donor, target_acceptor = target
    end_target = target_donor if strand == "+" else target_acceptor
    start_target = target_acceptor if strand == "+" else target_donor

    new_end = _search_end_side(genome, chrom, end, strand, end_target, max_adjustment)
    new_start = _search_start_side(genome, chrom, start, strand, start_target, max_adjustment)

    donor, acceptor = donor_acceptor(genome, chrom, new_start, new_end, strand)
    achieved = classify_splice_type(donor, acceptor)
    return new_start, new_end, achieved


def shift_all_junctions_independently(genome, chrom: str, strand: str, exon_starts: list, exon_ends: list,
                                       rng, dist_splice: dict, max_adjustment: int = 10):
    """Canonicalize every junction of a multi-block circRNA independently on
    `strand`: the outer BSJ plus every internal exon-exon junction.
    Returns (new_starts, new_ends, outer_achieved_splice_type), or None if
    any shift would invert or overlap the blocks it touches."""
    n = len(exon_starts)
    starts = list(exon_starts)
    ends = list(exon_ends)

    desired = sample_desired_splice_type(rng, dist_splice)
    new_start, new_end, outer_achieved = find_boundary_for_desired_type(
        genome, chrom, starts[0], ends[-1], strand, desired, max_adjustment)
    starts[0], ends[-1] = new_start, new_end

    for j in range(n - 1):
        desired_j = sample_desired_splice_type(rng, dist_splice)
        new_s, new_e, _ = find_boundary_for_desired_type(
            genome, chrom, starts[j + 1], ends[j], strand, desired_j, max_adjustment)
        starts[j + 1], ends[j] = new_s, new_e

    for k in range(n):
        if starts[k] >= ends[k]:
            return None
    for k in range(n - 1):
        if ends[k] > starts[k + 1]:
            return None
    return starts, ends, outer_achieved
