import numpy as np

from nanocirc_sim.sampling import sample_start_idx_biased


def test_terminal_exon_prob_zero_never_touches_first_or_last():
    rng = np.random.default_rng(0)
    n_exons, exon_count = 10, 3
    max_start_idx = n_exons - exon_count  # 7
    for _ in range(500):
        start_idx = sample_start_idx_biased(rng, n_exons, exon_count, terminal_exon_prob=0.0)
        end_idx = start_idx + exon_count - 1
        assert start_idx != 0, "should never use the first exon as start"
        assert end_idx != n_exons - 1, "should never use the last exon as end"
        assert 0 <= start_idx <= max_start_idx


def test_terminal_exon_prob_high_only_picks_terminal_windows():
    # terminal_exon_prob is clipped to 0.5 max (P(first) + P(last) <= 1),
    # so a very high request means: only the first/last windows are ever
    # drawn, never an internal one.
    rng = np.random.default_rng(0)
    n_exons, exon_count = 10, 3
    max_start_idx = n_exons - exon_count
    seen = set()
    for _ in range(500):
        seen.add(sample_start_idx_biased(rng, n_exons, exon_count, terminal_exon_prob=1.0))
    assert seen == {0, max_start_idx}


def test_falls_back_when_avoiding_terminal_exons_is_infeasible():
    rng = np.random.default_rng(0)
    # exon_count == n_exons: the only valid window necessarily includes both ends
    start_idx = sample_start_idx_biased(rng, n_exons=5, exon_count=5, terminal_exon_prob=0.0)
    assert start_idx == 0

    # exon_count == n_exons - 1: only two windows exist ([0..n-2]="first" or
    # [1..n-1]="last"), both touch a terminal exon: can't avoid either.
    seen = {sample_start_idx_biased(rng, n_exons=5, exon_count=4, terminal_exon_prob=0.0) for _ in range(200)}
    assert seen == {0, 1}


def test_empirical_rate_roughly_matches_requested_probability():
    rng = np.random.default_rng(1)
    n_exons, exon_count = 20, 3
    max_start_idx = n_exons - exon_count
    n = 20_000
    uses_first = uses_last = 0
    for _ in range(n):
        start_idx = sample_start_idx_biased(rng, n_exons, exon_count, terminal_exon_prob=0.075)
        end_idx = start_idx + exon_count - 1
        uses_first += start_idx == 0
        uses_last += end_idx == n_exons - 1
    assert abs(uses_first / n - 0.075) < 0.02
    assert abs(uses_last / n - 0.075) < 0.02
