"""One seeded RNG threaded through the whole package for reproducibility."""
from __future__ import annotations

import numpy as np


class RNGFactory:
    def __init__(self, seed: int):
        self.seed = seed
        self._seed_sequence = np.random.SeedSequence(seed)
        self.generator = np.random.default_rng(self._seed_sequence)

    def spawn(self, n: int):
        """Return n independent child generators, deterministic given the master seed."""
        return [np.random.default_rng(s) for s in self._seed_sequence.spawn(n)]
