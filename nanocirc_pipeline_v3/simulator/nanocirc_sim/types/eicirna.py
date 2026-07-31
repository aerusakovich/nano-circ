"""Exon-intron circRNAs (eicirna): exons plus at least one retained intron.
Each retained intron is its own Block, so it shows up as a real BED12 block."""
from __future__ import annotations

from ..blocks import Block, assemble_sequence, length_in_bounds, with_shifted_boundary
from ..candidate import Candidate
from ..sampling import sample_exon_count, sample_start_idx_biased
from ..splice_sites import find_boundary_for_desired_type, sample_desired_splice_type

CIRC_TYPE = "eicirna"


def make_attempt_fn(rng, genome, annotation, ann_index, cfg):
    dist_exon = cfg.exon_count_distributions[CIRC_TYPE]
    dist_splice = cfg.splice_site_distributions[CIRC_TYPE]
    retention_prob = cfg.intron_retention_prob.get(CIRC_TYPE, 0.7)
    max_adjustment = cfg.splice_search_max_adjustment.get(CIRC_TYPE, 10)

    def build(gene_id: str, tx_id: str, start_idx: int, end_idx: int):
        """Build a candidate from a fixed exon window, independently re-rolling
        which introns are retained."""
        tx = annotation.tx_structure[gene_id][tx_id]
        exons = tx.exons

        gaps = []  # (exon_index_before_gap, intron_start, intron_end)
        for i in range(start_idx, end_idx):
            s, e = exons[i].end, exons[i + 1].start
            if e > s:
                gaps.append((i, s, e))
        if not gaps:
            return None

        retain = {i: (rng.random() < retention_prob) for i, _, _ in gaps}
        if not any(retain.values()):
            forced_i = gaps[int(rng.integers(len(gaps)))][0]
            retain[forced_i] = True

        blocks = []
        for i in range(start_idx, end_idx + 1):
            blocks.append(Block(exons[i].start, exons[i].end, "exon"))
            if retain.get(i):
                gap = next(g for g in gaps if g[0] == i)
                blocks.append(Block(gap[1], gap[2], "intron_retained"))

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

        # reject if final coordinates would actually classify as ecirna
        if ann_index.classify_final_type(tx.chrom, tx.strand, [(b.start, b.end) for b in blocks]) != CIRC_TYPE:
            return None

        mature_seq = assemble_sequence(genome, tx.chrom, blocks, tx.strand)
        if not mature_seq:
            return None

        n_retained = sum(1 for b in blocks if b.kind == "intron_retained")
        n_exons = sum(1 for b in blocks if b.kind == "exon")
        return Candidate(
            circ_type=CIRC_TYPE,
            chrom=tx.chrom, start=blocks[0].start, end=blocks[-1].end, strand=tx.strand,
            blocks=blocks, mature_seq=mature_seq,
            gene_id=gene_id, tx_id=tx_id,
            splice_type=achieved, skipped_exons=0, retained_introns=n_retained,
            exon_count=n_exons,
            boundary_key=(start_idx, end_idx),
        )

    def attempt(gene_id: str, tx_id: str):
        tx = annotation.tx_structure[gene_id][tx_id]
        exons = tx.exons
        exon_count = max(2, sample_exon_count(rng, dist_exon))  # need >=1 intron gap
        if exon_count > len(exons):
            return None
        start_idx = sample_start_idx_biased(rng, len(exons), exon_count, cfg.terminal_exon_prob)
        end_idx = start_idx + exon_count - 1
        return build(gene_id, tx_id, start_idx, end_idx)

    attempt.build_with_boundary = lambda gene_id, tx_id, boundary_key: build(gene_id, tx_id, *boundary_key)
    return attempt
