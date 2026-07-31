import numpy as np

from nanocirc_sim.blocks import Block
from nanocirc_sim.candidate import Candidate
from nanocirc_sim.expression import assign_read_targets, assign_tpm


def _make_candidates(n, circ_type="ecirna"):
    return [
        Candidate(circ_type=circ_type, chrom="chr1", start=i * 10, end=i * 10 + 5, strand="+",
                  blocks=[Block(i * 10, i * 10 + 5, "region")], mature_seq="A" * 5)
        for i in range(n)
    ]


def test_read_targets_sum_exactly_to_budget():
    rng = np.random.default_rng(0)
    cands = _make_candidates(37)
    assign_tpm(rng, cands)
    assign_read_targets(rng, cands, total_reads=10_000)
    assert sum(c.n_reads_target for c in cands) == 10_000
    assert all(c.n_reads_target >= 0 for c in cands)


def test_zero_total_reads_gives_zero_targets():
    rng = np.random.default_rng(0)
    cands = _make_candidates(5)
    assign_tpm(rng, cands)
    assign_read_targets(rng, cands, total_reads=0)
    assert all(c.n_reads_target == 0 for c in cands)


def test_empty_candidate_list_is_safe():
    rng = np.random.default_rng(0)
    assign_tpm(rng, [])
    assign_read_targets(rng, [], total_reads=1000)  # must not raise


def test_min_reads_floor_applied_to_every_candidate():
    rng = np.random.default_rng(0)
    cands = _make_candidates(50)
    assign_tpm(rng, cands)
    # zero out TPM for most candidates so, without a floor, they'd get 0 reads
    for c in cands[1:]:
        c.tpm = 0.0
    assign_read_targets(rng, cands, total_reads=10_000, min_reads_per_circrna=30)
    assert all(c.n_reads_target >= 30 for c in cands)
    assert sum(c.n_reads_target for c in cands) == 10_000


def test_min_reads_floor_gracefully_degrades_when_unaffordable():
    rng = np.random.default_rng(0)
    cands = _make_candidates(10)
    assign_tpm(rng, cands)
    # requesting a floor of 100 each (1000 total) with only 50 reads to give out
    assign_read_targets(rng, cands, total_reads=50, min_reads_per_circrna=100)
    assert sum(c.n_reads_target for c in cands) == 50
    assert all(c.n_reads_target == 5 for c in cands)  # 50 // 10 == 5, applied evenly
