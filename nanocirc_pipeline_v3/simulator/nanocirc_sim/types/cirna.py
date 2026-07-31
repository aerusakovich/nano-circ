"""Intronic circRNAs (ciRNA): exactly one intron per circRNA."""
from __future__ import annotations

from ..blocks import Block
from ..candidate import Candidate
from ..splice_sites import find_boundary_for_desired_type, sample_desired_splice_type

CIRC_TYPE = "cirna"


def make_attempt_fn(rng, genome, annotation, ann_index, cfg, valid_transcripts):
    dist_splice = cfg.splice_site_distributions[CIRC_TYPE]
    min_intron_length = cfg.min_intron_length
    max_adjustment = cfg.splice_search_max_adjustment.get(CIRC_TYPE, 10)

    def attempt():
        if not valid_transcripts:
            return None
        gene_id, tx_id = valid_transcripts[int(rng.integers(len(valid_transcripts)))]
        tx = annotation.tx_structure[gene_id][tx_id]
        introns = [iv for iv in tx.introns if iv[1] - iv[0] >= min_intron_length]
        if not introns:
            return None
        start, end = introns[int(rng.integers(len(introns)))]

        desired = sample_desired_splice_type(rng, dist_splice)
        new_start, new_end, achieved = find_boundary_for_desired_type(
            genome, tx.chrom, start, end, tx.strand, desired, max_adjustment
        )
        if not (cfg.min_circle_length <= new_end - new_start <= cfg.max_circle_length):
            return None

        if ann_index.classify_region(tx.chrom, tx.strand, new_start, new_end) != "cirna":
            return None

        mature_seq = genome.fetch(tx.chrom, new_start, new_end, tx.strand)
        if not mature_seq:
            return None

        return Candidate(
            circ_type=CIRC_TYPE,
            chrom=tx.chrom, start=new_start, end=new_end, strand=tx.strand,
            blocks=[Block(new_start, new_end, "region")], mature_seq=mature_seq,
            gene_id=gene_id, tx_id=tx_id,
            splice_type=achieved, exon_count=1,
        )

    return attempt
