"""Shared retry-until-satisfied loop used by every circRNA type generator."""
from __future__ import annotations

import logging

from tqdm import tqdm

logger = logging.getLogger(__name__)


def _default_dedup_key(cand):
    return (cand.chrom, cand.strand, tuple((b.start, b.end, b.kind) for b in cand.blocks))


def generate_until_satisfied(target_count: int, attempt_fn, max_attempts_multiplier: int = 50,
                              label: str = "circRNA", show_progress: bool = True, dedup: bool = True):
    """Calls `attempt_fn()` until `target_count` candidates are collected or
    the attempt budget runs out. Returns (results, attempts, rejection_rate)."""
    if target_count <= 0:
        return [], 0, 0.0

    max_attempts = max(target_count * max_attempts_multiplier, 200)
    results = []
    seen = set()
    attempts = 0
    pbar = tqdm(total=target_count, desc=label, disable=not show_progress)
    try:
        while len(results) < target_count and attempts < max_attempts:
            attempts += 1
            cand = attempt_fn()
            if cand is None:
                continue
            if dedup:
                key = _default_dedup_key(cand)
                if key in seen:
                    continue
                seen.add(key)
            results.append(cand)
            pbar.update(1)
    finally:
        pbar.close()

    rejection_rate = 1.0 - (len(results) / attempts) if attempts else 0.0
    if len(results) < target_count:
        logger.warning(
            "%s: requested %d, achieved %d after %d attempts (%.1f%% rejected, hit attempt budget)",
            label, target_count, len(results), attempts, 100 * rejection_rate,
        )
    else:
        logger.info(
            "%s: achieved requested %d in %d attempts (%.1f%% rejected)",
            label, target_count, attempts, 100 * rejection_rate,
        )
    return results, attempts, rejection_rate
