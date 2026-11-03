"""Generates multiple isoforms per locus, shared by ecirna/eicirna/antisense."""
from __future__ import annotations

import logging

from tqdm import tqdm

logger = logging.getLogger(__name__)


def sample_isoforms_per_locus(rng, isoform_mean: float, isoform_max: int) -> int:
    """Number of isoforms to draw from one locus, capped at isoform_max."""
    extra = rng.poisson(max(isoform_mean - 1.0, 0.0))
    return int(min(max(1, 1 + extra), isoform_max))


def _signature(cand):
    return (cand.chrom, cand.strand, tuple((b.start, b.end, b.kind) for b in cand.blocks))


def generate_with_isoforms(target_count: int, rng, valid_loci, attempt_fn, cfg, label: str = "circRNA",
                            show_progress: bool = True):
    """Returns (results, attempts, rejection_rate), same contract as
    types.base.generate_until_satisfied."""
    if target_count <= 0 or not valid_loci:
        return [], 0, 0.0

    max_attempts = max(target_count * cfg.max_attempts_multiplier, 200)
    results = []
    attempts = 0
    seen_signatures = set()  # global dedup across every locus
    locus_next_index = {}  # locus_id -> next isoform_index to assign
    locus_boundaries = {}  # locus_id -> list of boundary_key already used at that locus

    same_bsj_prob = getattr(cfg, "same_bsj_isoform_prob", 0.0)
    supports_boundary_reuse = hasattr(attempt_fn, "build_with_boundary")

    pbar = tqdm(total=target_count, desc=label, disable=not show_progress)
    try:
        while len(results) < target_count and attempts < max_attempts:
            gene_id, tx_id = valid_loci[int(rng.integers(len(valid_loci)))]
            locus_id = f"{gene_id}:{tx_id}"

            # caps this one visit only; the locus can be revisited later
            n_isoforms_wanted = sample_isoforms_per_locus(rng, cfg.isoform_mean, cfg.isoform_max)
            n_isoforms_wanted = min(n_isoforms_wanted, target_count - len(results))

            produced = 0
            local_attempts = 0
            local_budget = max(n_isoforms_wanted * 15, 15)
            while produced < n_isoforms_wanted and local_attempts < local_budget and attempts < max_attempts:
                local_attempts += 1
                attempts += 1

                prior_boundaries = locus_boundaries.get(locus_id)
                reuse_boundary = (
                    supports_boundary_reuse and prior_boundaries and rng.random() < same_bsj_prob
                )
                if reuse_boundary:
                    boundary_key = prior_boundaries[int(rng.integers(len(prior_boundaries)))]
                    cand = attempt_fn.build_with_boundary(gene_id, tx_id, boundary_key)
                else:
                    cand = attempt_fn(gene_id, tx_id)

                if cand is None:
                    continue
                sig = _signature(cand)
                if sig in seen_signatures:
                    continue
                seen_signatures.add(sig)

                if supports_boundary_reuse and cand.boundary_key is not None:
                    bucket = locus_boundaries.setdefault(locus_id, [])
                    if cand.boundary_key not in bucket:
                        bucket.append(cand.boundary_key)

                cand.locus_id = locus_id
                cand.isoform_index = locus_next_index.get(locus_id, 0)
                locus_next_index[locus_id] = cand.isoform_index + 1
                results.append(cand)
                produced += 1
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
