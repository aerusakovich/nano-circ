"""Exonic circRNAs (ecirna): fully-exonic mature sequence, with exon skipping."""
from __future__ import annotations

from ..blocks import Block, assemble_sequence, length_in_bounds, with_shifted_boundary
from ..candidate import Candidate
from ..sampling import sample_exon_count, sample_start_idx_biased
from ..splice_sites import find_boundary_for_desired_type, sample_desired_splice_type

CIRC_TYPE = "ecirna"


def make_attempt_fn(rng, genome, annotation, ann_index, cfg, circ_type: str = CIRC_TYPE):
    dist_exon = cfg.exon_count_distributions[circ_type]
    dist_splice = cfg.splice_site_distributions[circ_type]
    skip_prob = cfg.exon_skip_prob.get(circ_type, 0.0)
    max_adjustment = cfg.splice_search_max_adjustment.get(circ_type, 10)

    def build(gene_id: str, tx_id: str, start_idx: int, end_idx: int):
        """Build a candidate from a fixed exon window, independently
        re-rolling which internal exons are skipped."""
        tx = annotation.tx_structure[gene_id][tx_id]
        exons = tx.exons

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

        blocks = [Block(exons[i].start, exons[i].end, "exon") for i in included_idx]

        desired = sample_desired_splice_type(rng, dist_splice)
        new_start, new_end, achieved = find_boundary_for_desired_type(
            genome, tx.chrom, blocks[0].start, blocks[-1].end, tx.strand, desired, max_adjustment
        )
        blocks = with_shifted_boundary(blocks, new_start, new_end)
        if blocks is None:
            return None
        if not length_in_bounds(blocks, cfg.min_circle_length, cfg.max_circle_length):
            return None

        if not ann_index.genes(tx.chrom, tx.strand).overlaps(blocks[0].start, blocks[-1].end):
            return None

        # reject if final coordinates would actually classify as eicirna
        if ann_index.classify_final_type(tx.chrom, tx.strand, [(b.start, b.end) for b in blocks]) != circ_type:
            return None

        mature_seq = assemble_sequence(genome, tx.chrom, blocks, tx.strand)
        if not mature_seq:
            return None

        return Candidate(
            circ_type=circ_type,
            chrom=tx.chrom, start=blocks[0].start, end=blocks[-1].end, strand=tx.strand,
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
