"""Assigns per-isoform TPM and a target read count (TPM-weighted multinomial
over total_reads, so the targets sum exactly to the budget)."""
from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)

_TPM_BETA_PARAMS = {
    "ecirna": (0.5, 2.0, 100.0),
    "eicirna": (0.5, 2.0, 80.0),
    "cirna": (0.3, 2.5, 60.0),
    "antisense": (0.4, 2.2, 70.0),
    "intergenic": (0.2, 3.0, 40.0),
}


def assign_tpm(rng: np.random.Generator, candidates: list) -> None:
    for c in candidates:
        a, b, scale = _TPM_BETA_PARAMS.get(c.circ_type, (0.5, 2.0, 100.0))
        c.tpm = float(rng.beta(a, b) * scale)


def assign_read_targets(rng: np.random.Generator, candidates: list, total_reads: int,
                         min_reads_per_circrna: int = 0) -> None:
    """Floors each candidate at `min_reads_per_circrna`, allocates the rest
    of `total_reads` by TPM-weighted multinomial."""
    if not candidates:
        return
    if total_reads <= 0:
        for c in candidates:
            c.n_reads_target = 0
        return

    n = len(candidates)
    floor = min(min_reads_per_circrna, total_reads // n)
    if floor < min_reads_per_circrna:
        logger.warning(
            "min_reads_per_circrna=%d not affordable for %d candidates within total_reads=%d; "
            "using floor=%d instead", min_reads_per_circrna, n, total_reads, floor,
        )
    remaining = total_reads - floor * n

    weights = np.array([max(c.tpm, 0.0) for c in candidates], dtype=float)
    total_w = weights.sum()
    if total_w <= 0:
        weights = np.ones(n)
        total_w = weights.sum()
    probs = weights / total_w
    extra = rng.multinomial(remaining, probs) if remaining > 0 else np.zeros(n, dtype=int)
    for c, e in zip(candidates, extra):
        c.n_reads_target = floor + int(e)
