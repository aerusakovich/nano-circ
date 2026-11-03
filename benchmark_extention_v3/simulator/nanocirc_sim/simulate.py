"""Runs one full simulation: generate candidates, assign expression and
target reads, write outputs."""
from __future__ import annotations

import logging
import time

from . import expression as expr
from . import sampling
from . import writer
from .annotation import Annotation
from .genome import GenomeLoader
from .intervals import AnnotationIndex
from .isoforms import generate_with_isoforms
from .rng import RNGFactory
from .types import antisense, cirna, ecirna, eicirna, intergenic
from .types.base import generate_until_satisfied

logger = logging.getLogger(__name__)


def run_simulation(*, genome_path: str, gtf_path: str, config, outdir: str,
                    total: int = None, explicit_counts: dict = None, seed: int = None,
                    show_progress: bool = True) -> dict:
    t0 = time.time()
    config.validate()
    seed = config.seed if seed is None else seed
    rngf = RNGFactory(seed)
    rng = rngf.generator

    logger.info("Loading annotation from %s", gtf_path)
    annotation = Annotation.from_gtf(gtf_path)
    ann_index = AnnotationIndex(annotation)
    logger.info("Loading genome from %s", genome_path)
    genome = GenomeLoader(genome_path)

    counts = config.resolve_counts(total=total, explicit_counts=explicit_counts)
    logger.info("Target per-type counts: %s (total %d)", counts, sum(counts.values()))

    sense_loci = annotation.get_transcripts(min_exons=2, protein_coding_only=True)
    # antisense draws from the same locus pool; types/antisense.py flips the strand.
    antisense_loci = sense_loci
    if (counts.get("ecirna", 0) > 0 or counts.get("eicirna", 0) > 0 or counts.get("cirna", 0) > 0
            or counts.get("antisense", 0) > 0) and not sense_loci:
        logger.warning(
            "No protein-coding genes with >=2 exons found in GTF; "
            "ecirna/eicirna/cirna/antisense counts will be 0.")

    all_candidates = []
    type_stats = {}

    def _record(circ_type, results, attempts, rejection_rate, requested):
        type_stats[circ_type] = {
            "requested": requested, "achieved": len(results),
            "attempts": attempts, "rejection_rate": round(rejection_rate, 4),
        }
        if requested > 0 and len(results) < requested * config.shortfall_warn_threshold:
            logger.warning(
                "%s: only achieved %d/%d requested (below %.0f%% threshold)",
                circ_type, len(results), requested, config.shortfall_warn_threshold * 100,
            )
        all_candidates.extend(results)

    if counts.get("ecirna", 0) > 0:
        attempt = ecirna.make_attempt_fn(rng, genome, annotation, ann_index, config)
        res, att, rej = generate_with_isoforms(
            counts["ecirna"], rng, sense_loci, attempt, config, label="ecirna", show_progress=show_progress)
        _record("ecirna", res, att, rej, counts["ecirna"])

    if counts.get("eicirna", 0) > 0:
        attempt = eicirna.make_attempt_fn(rng, genome, annotation, ann_index, config)
        res, att, rej = generate_with_isoforms(
            counts["eicirna"], rng, sense_loci, attempt, config, label="eicirna", show_progress=show_progress)
        _record("eicirna", res, att, rej, counts["eicirna"])

    if counts.get("antisense", 0) > 0:
        attempt = antisense.make_attempt_fn(rng, genome, annotation, ann_index, config)
        res, att, rej = generate_with_isoforms(
            counts["antisense"], rng, antisense_loci, attempt, config, label="antisense", show_progress=show_progress)
        _record("antisense", res, att, rej, counts["antisense"])

    if counts.get("cirna", 0) > 0:
        attempt = cirna.make_attempt_fn(rng, genome, annotation, ann_index, config, sense_loci)
        res, att, rej = generate_until_satisfied(
            counts["cirna"], attempt, config.max_attempts_multiplier, label="cirna", show_progress=show_progress)
        _record("cirna", res, att, rej, counts["cirna"])

    if counts.get("intergenic", 0) > 0:
        attempt = intergenic.make_attempt_fn(rng, genome, annotation, ann_index, config)
        res, att, rej = generate_until_satisfied(
            counts["intergenic"], attempt, config.max_attempts_multiplier, label="intergenic",
            show_progress=show_progress)
        _record("intergenic", res, att, rej, counts["intergenic"])

    logger.info("Generated %d circRNAs total", len(all_candidates))

    expr.assign_tpm(rng, all_candidates)
    expr.assign_read_targets(rng, all_candidates, config.total_reads, config.min_reads_per_circrna)
    for c in all_candidates:
        c.copy_number = sampling.sample_copy_number(rng, config)
        c.pseudo_seq = sampling.rolling_circle_sequence(rng, c.mature_seq, c.copy_number)

    genome.close()

    elapsed = time.time() - t0
    paths = writer.write_all_outputs(
        all_candidates, outdir, config=config, seed=seed,
        inputs={"genome": genome_path, "gtf": gtf_path},
        type_stats=type_stats, elapsed_seconds=elapsed,
    )
    logger.info("Simulation complete in %.1fs. Outputs written to %s", elapsed, outdir)
    return {"candidates": all_candidates, "type_stats": type_stats, "paths": paths, "elapsed_seconds": elapsed}
