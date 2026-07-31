"""Small, seeded-RNG-driven sampling helpers shared by the type generators."""
from __future__ import annotations

import numpy as np


def sample_exon_count(rng: np.random.Generator, dist: dict) -> int:
    keys = list(dist.keys())
    probs = np.array([dist[k] for k in keys], dtype=float)
    total = probs.sum()
    if total <= 0:
        return 1
    probs = probs / total
    choice = rng.choice(keys, p=probs)
    if choice == "5+":
        return int(rng.integers(5, 11))
    return int(choice)


def sample_start_idx_biased(rng: np.random.Generator, n_exons: int, exon_count: int,
                             terminal_exon_prob: float) -> int:
    """Choose the start exon index for an `exon_count`-long window out of
    `n_exons`, biased so P(touches first/last exon) == terminal_exon_prob."""
    max_start_idx = n_exons - exon_count
    if max_start_idx <= 0:
        return 0
    n_choices = max_start_idx + 1
    if n_choices <= 2:
        # every choice touches a terminal exon, nothing to bias
        return int(rng.integers(0, n_choices))

    p = min(max(terminal_exon_prob, 0.0), 0.5)
    weights = np.full(n_choices, (1.0 - 2 * p) / (n_choices - 2))
    weights[0] = p
    weights[-1] = p
    return int(rng.choice(n_choices, p=weights))


def sample_length(rng: np.random.Generator, dist: dict) -> int:
    """Skewed-quartile length sampler."""
    r = rng.random()
    if r < 0.25:
        val = rng.uniform(dist["min"], dist["q1"])
    elif r < 0.5:
        val = rng.uniform(dist["q1"], dist["median"])
    elif r < 0.75:
        val = rng.uniform(dist["median"], dist["q3"])
    else:
        if rng.random() < 0.9:
            val = rng.uniform(dist["q3"], dist["q3"] * 2)
        else:
            val = rng.uniform(dist["q3"] * 2, dist["max"])
    return max(1, int(val))


def sample_copy_number(rng: np.random.Generator, cfg) -> int:
    ranges = cfg.copy_number_ranges
    probs = cfg.copy_number_probabilities
    idx = int(rng.choice(len(ranges), p=probs))
    lo, hi = ranges[idx]
    if idx == len(ranges) - 1:
        base = rng.uniform(lo, lo + 40)
        if rng.random() < 0.1:
            return int(base + rng.exponential(scale=100))
        return int(base)
    return int(rng.uniform(lo, hi))


def rolling_circle_sequence(rng: np.random.Generator, mature_seq: str, copy_number: int) -> str:
    if len(mature_seq) > 1:
        start_pos = int(rng.integers(0, len(mature_seq)))
    else:
        start_pos = 0
    return mature_seq[start_pos:] + mature_seq * copy_number
