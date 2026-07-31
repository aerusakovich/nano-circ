import numpy as np

from nanocirc_sim.blocks import Block
from nanocirc_sim.candidate import Candidate
from nanocirc_sim.config import SimConfig
from nanocirc_sim.isoforms import generate_with_isoforms


def _make_cand(chrom, start, end, boundary_key=None, kind="exon"):
    return Candidate(circ_type="ecirna", chrom=chrom, start=start, end=end, strand="+",
                      blocks=[Block(start, end, kind)], mature_seq="A" * (end - start),
                      boundary_key=boundary_key)


def test_global_dedup_rejects_identical_structure_from_different_loci():
    # Simulates the real bug: two different transcripts of the same gene
    # (different locus_id) that happen to produce a byte-identical circRNA.
    # It must only ever appear once in the results, regardless of which
    # locus "produced" it.
    rng = np.random.default_rng(0)
    cfg = SimConfig.default()
    cfg.same_bsj_isoform_prob = 0.0
    valid_loci = [("geneA", "tx1"), ("geneA", "tx2")]

    def attempt_fn(gene_id, tx_id):
        return _make_cand("chr1", 1000, 1100)  # always identical, regardless of locus

    results, attempts, _ = generate_with_isoforms(5, rng, valid_loci, attempt_fn, cfg, show_progress=False)
    assert len(results) == 1
    assert attempts >= 5  # kept retrying (and rejecting the duplicate) trying to reach target=5


def test_same_bsj_isoform_prob_one_forces_shared_boundary():
    rng = np.random.default_rng(0)
    cfg = SimConfig.default()
    cfg.same_bsj_isoform_prob = 1.0
    cfg.isoform_mean = 3.0
    cfg.isoform_max = 3
    valid_loci = [("geneA", "tx1")]

    counter = {"n": 0}

    def attempt_fn(gene_id, tx_id):
        counter["n"] += 1
        key = (counter["n"], counter["n"])
        return _make_cand("chr1", 1000 + counter["n"], 1100 + counter["n"], boundary_key=key)

    def build_with_boundary(gene_id, tx_id, boundary_key):
        counter["n"] += 1
        s, _ = boundary_key
        # different internal content (kind), but the SAME genomic boundary
        return _make_cand("chr1", 1000 + s, 1100 + s, boundary_key=boundary_key, kind=f"variant{counter['n']}")

    attempt_fn.build_with_boundary = build_with_boundary

    results, attempts, _ = generate_with_isoforms(3, rng, valid_loci, attempt_fn, cfg, show_progress=False)
    assert len(results) == 3
    boundary_keys = {c.boundary_key for c in results}
    assert len(boundary_keys) == 1, "with same_bsj_isoform_prob=1.0, all isoforms should share one BSJ boundary"
    # but they're still structurally distinct (different block 'kind'), so all got kept
    assert len({(c.start, c.end, c.blocks[0].kind) for c in results}) == 3


def test_same_bsj_isoform_prob_zero_never_reuses_boundary():
    rng = np.random.default_rng(0)
    cfg = SimConfig.default()
    cfg.same_bsj_isoform_prob = 0.0
    cfg.isoform_mean = 3.0
    cfg.isoform_max = 3
    valid_loci = [("geneA", "tx1")]

    counter = {"n": 0}

    def attempt_fn(gene_id, tx_id):
        counter["n"] += 1
        key = (counter["n"], counter["n"])
        return _make_cand("chr1", 1000 + counter["n"], 1100 + counter["n"], boundary_key=key)

    def build_with_boundary(gene_id, tx_id, boundary_key):
        raise AssertionError("build_with_boundary should never be called when same_bsj_isoform_prob=0.0")

    attempt_fn.build_with_boundary = build_with_boundary

    results, attempts, _ = generate_with_isoforms(3, rng, valid_loci, attempt_fn, cfg, show_progress=False)
    assert len(results) == 3
    assert len({c.boundary_key for c in results}) == 3
