"""Intergenic circRNAs: no gene overlap on either strand."""
from __future__ import annotations

from ..blocks import Block
from ..candidate import Candidate
from ..sampling import sample_length
from ..splice_sites import find_boundary_for_desired_type, sample_desired_splice_type

CIRC_TYPE = "intergenic"


def make_attempt_fn(rng, genome, annotation, ann_index, cfg):
    dist_length = cfg.length_distributions[CIRC_TYPE]
    dist_splice = cfg.splice_site_distributions[CIRC_TYPE]
    max_adjustment = cfg.splice_search_max_adjustment.get(CIRC_TYPE, 10)
    chroms = [c for c, length in genome.chrom_lengths.items() if length >= 10_000]

    def attempt():
        if not chroms:
            return None
        chrom = chroms[int(rng.integers(len(chroms)))]
        chrom_len = genome.chrom_lengths[chrom]

        length = sample_length(rng, dist_length)
        if length >= chrom_len - 2:
            return None
        start = int(rng.integers(0, chrom_len - length))
        end = start + length
        # no biological basis for a strand; fetch_strand is just for sequence
        # extraction, reported strand stays "." (unknown)
        fetch_strand = "+" if rng.random() < 0.5 else "-"

        if not ann_index.is_intergenic(chrom, start, end):
            return None

        desired = sample_desired_splice_type(rng, dist_splice)
        new_start, new_end, achieved = find_boundary_for_desired_type(
            genome, chrom, start, end, fetch_strand, desired, max_adjustment
        )
        if not (cfg.min_circle_length <= new_end - new_start <= cfg.max_circle_length):
            return None
        if not ann_index.is_intergenic(chrom, new_start, new_end):
            return None

        mature_seq = genome.fetch(chrom, new_start, new_end, fetch_strand)
        if not mature_seq:
            return None

        return Candidate(
            circ_type=CIRC_TYPE,
            chrom=chrom, start=new_start, end=new_end, strand=".",
            blocks=[Block(new_start, new_end, "region")], mature_seq=mature_seq,
            gene_id="NA", tx_id="NA",
            splice_type=achieved, exon_count=1,
        )

    return attempt
