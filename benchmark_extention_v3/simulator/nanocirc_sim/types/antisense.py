"""Antisense circRNAs: transcribed from the strand OPPOSITE an overlapping
gene's exons."""
from __future__ import annotations

from ..blocks import Block, assemble_sequence, length_in_bounds
from ..candidate import Candidate
from ..sampling import sample_exon_count, sample_start_idx_biased
from ..splice_sites import shift_all_junctions_independently

CIRC_TYPE = "antisense"


def _flip(strand: str) -> str:
    return "-" if strand == "+" else "+"


def make_attempt_fn(rng, genome, annotation, ann_index, cfg):
    dist_exon = cfg.exon_count_distributions[CIRC_TYPE]
    dist_splice = cfg.splice_site_distributions[CIRC_TYPE]
    skip_prob = cfg.exon_skip_prob.get(CIRC_TYPE, 0.0)
    max_adjustment = cfg.splice_search_max_adjustment.get(CIRC_TYPE, 10)

    def build(gene_id: str, tx_id: str, start_idx: int, end_idx: int):
        """Build from a fixed exon window on the host gene's coordinates,
        placed on the opposite strand."""
        tx = annotation.tx_structure[gene_id][tx_id]
        exons = tx.exons
        antisense_strand = _flip(tx.strand)

        # keep boundary exons; internal exons are independently skippable
        included_idx = [start_idx]
        n_skipped = 0
        for i in range(start_idx + 1, end_idx):
            if rng.random() < skip_prob:
                n_skipped += 1
                continue
            included_idx.append(i)
        if end_idx != start_idx:
            included_idx.append(end_idx)

        # search every junction independently for a canonical motif on antisense_strand
        exon_starts = [exons[i].start for i in included_idx]
        exon_ends = [exons[i].end for i in included_idx]
        shifted = shift_all_junctions_independently(
            genome, tx.chrom, antisense_strand, exon_starts, exon_ends, rng, dist_splice, max_adjustment
        )
        if shifted is None:
            return None
        starts, ends, achieved = shifted
        blocks = [Block(starts[k], ends[k], "exon") for k in range(len(starts))]
        if not length_in_bounds(blocks, cfg.min_circle_length, cfg.max_circle_length):
            return None

        # reject if a gene already exists on the circRNA's own (flipped) strand
        if ann_index.genes(tx.chrom, antisense_strand).overlaps(blocks[0].start, blocks[-1].end):
            return None

        mature_seq = assemble_sequence(genome, tx.chrom, blocks, antisense_strand)
        if not mature_seq:
            return None

        return Candidate(
            circ_type=CIRC_TYPE,
            chrom=tx.chrom, start=blocks[0].start, end=blocks[-1].end, strand=antisense_strand,
            blocks=blocks, mature_seq=mature_seq,
            gene_id=gene_id, tx_id=tx_id,
            splice_type=achieved, skipped_exons=n_skipped, retained_introns=0,
            exon_count=len(blocks),
            boundary_key=(start_idx, end_idx),
        )

    def attempt(gene_id: str, tx_id: str):
        tx = annotation.tx_structure[gene_id][tx_id]
        exons = tx.exons
        exon_count = sample_exon_count(rng, dist_exon)
        if exon_count > len(exons):
            return None
        start_idx = sample_start_idx_biased(rng, len(exons), exon_count, cfg.terminal_exon_prob)
        end_idx = start_idx + exon_count - 1
        return build(gene_id, tx_id, start_idx, end_idx)

    attempt.build_with_boundary = lambda gene_id, tx_id, boundary_key: build(gene_id, tx_id, *boundary_key)
    return attempt
