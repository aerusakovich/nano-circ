"""Fits config.py distributions from a real circRNA database (BED6/BED12)."""
from __future__ import annotations

import copy
import logging

import numpy as np

from .annotation import Annotation
from .config import CIRC_TYPES, SimConfig
from .genome import GenomeLoader
from .intervals import AnnotationIndex
from .splice_sites import SPLICE_TYPES, classify_splice_type, donor_acceptor

logger = logging.getLogger(__name__)

MIN_EXAMPLES_TO_FIT = 5
# bp tolerance for matching a real circRNA's boundary to an annotated exon boundary.
TERMINAL_EXON_MATCH_TOL = 15


def _parse_bed_record(line: str):
    """Parses one BED6/BED12 line. Returns (chrom, strand, blocks, block_count,
    has_block_info) or None."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6:
        return None
    try:
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    except ValueError:
        return None
    if end <= start:
        return None
    strand = fields[5] if fields[5] in ("+", "-") else "+"

    blocks = None
    block_count = None
    has_block_info = False
    if len(fields) >= 12:
        try:
            block_count = int(fields[9])
            sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            if block_count > 0 and len(sizes) == block_count and len(starts) == block_count:
                blocks = sorted((start + s, start + s + sz) for s, sz in zip(starts, sizes))
                has_block_info = True
        except (ValueError, IndexError):
            blocks = None
    if not blocks:
        blocks = [(start, end)]
        block_count = 1
    return chrom, strand, blocks, block_count, has_block_info


def _exon_count_bucket(n: int) -> str:
    return str(n) if 1 <= n <= 4 else "5+"


def _uses_terminal_exon(span_start: int, span_end: int, gene_ids: list, annotation) -> tuple:
    """Returns (matched, uses_first, uses_last) against real transcript exon boundaries."""
    tol = TERMINAL_EXON_MATCH_TOL
    matched = uses_first = uses_last = False
    for gene_id in gene_ids:
        for tx in annotation.tx_structure.get(gene_id, {}).values():
            exons = tx.exons
            if len(exons) < 2:
                continue
            near_any_start = any(abs(span_start - e.start) <= tol for e in exons)
            near_any_end = any(abs(span_end - e.end) <= tol for e in exons)
            if near_any_start and near_any_end:
                matched = True
                if abs(span_start - exons[0].start) <= tol:
                    uses_first = True
                if abs(span_end - exons[-1].end) <= tol:
                    uses_last = True
    return matched, uses_first, uses_last


def _quartiles(values) -> dict:
    arr = np.array(values, dtype=float)
    return {
        "min": float(arr.min()), "q1": float(np.percentile(arr, 25)),
        "median": float(np.percentile(arr, 50)), "q3": float(np.percentile(arr, 75)),
        "max": float(arr.max()),
    }


def fit_config_from_database(*, db_bed_path: str, genome_path: str, gtf_path: str,
                              base_config: SimConfig = None) -> SimConfig:
    """Classify every record in `db_bed_path` and return a new SimConfig
    refit from that real population."""
    cfg = copy.deepcopy(base_config) if base_config is not None else SimConfig.default()

    logger.info("dbfit: loading annotation %s", gtf_path)
    annotation = Annotation.from_gtf(gtf_path)
    ann_index = AnnotationIndex(annotation)
    genome = GenomeLoader(genome_path)

    lengths = {t: [] for t in CIRC_TYPES}
    exon_buckets = {t: [] for t in CIRC_TYPES}
    splice_counts = {t: {s: 0 for s in SPLICE_TYPES} for t in CIRC_TYPES}
    type_counts = {t: 0 for t in CIRC_TYPES}
    n_terminal_matched = n_uses_first = n_uses_last = 0
    n_lines = 0
    n_skipped = 0

    with open(db_bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            n_lines += 1
            parsed = _parse_bed_record(line)
            if parsed is None:
                n_skipped += 1
                continue
            chrom, strand, blocks, block_count, has_block_info = parsed
            span_start, span_end = blocks[0][0], blocks[-1][1]
            mature_len = sum(e - s for s, e in blocks)

            base_type = ann_index.classify_region_by_blocks(chrom, strand, blocks)
            if base_type in ("ecirna", "eicirna"):
                overlapping = ann_index.overlapping_genes(chrom, strand, span_start, span_end)
                matched, uses_first, uses_last = _uses_terminal_exon(span_start, span_end, overlapping, annotation)
                if matched:
                    n_terminal_matched += 1
                    n_uses_first += uses_first
                    n_uses_last += uses_last
            circ_type = ann_index.classify_final_type(chrom, strand, blocks, base_type=base_type)

            type_counts[circ_type] += 1
            lengths[circ_type].append(mature_len)
            if block_count and has_block_info:
                exon_buckets[circ_type].append(_exon_count_bucket(block_count))

            donor, acceptor = donor_acceptor(genome, chrom, span_start, span_end, strand)
            splice_counts[circ_type][classify_splice_type(donor, acceptor)] += 1

    genome.close()
    logger.info("dbfit: %d records read, %d skipped/unparseable, classified counts: %s",
                n_lines, n_skipped, type_counts)

    total_classified = sum(type_counts.values())
    if total_classified > 0:
        cfg.type_fractions = {t: type_counts[t] / total_classified for t in CIRC_TYPES}

    for t in CIRC_TYPES:
        if len(lengths[t]) >= MIN_EXAMPLES_TO_FIT:
            cfg.length_distributions[t] = _quartiles(lengths[t])
        else:
            logger.info("dbfit: only %d length examples for type=%s, keeping base config", len(lengths[t]), t)

        if len(exon_buckets[t]) >= MIN_EXAMPLES_TO_FIT:
            counts = {b: exon_buckets[t].count(b) for b in ["1", "2", "3", "4", "5+"]}
            total = sum(counts.values())
            cfg.exon_count_distributions[t] = {k: v / total for k, v in counts.items()}
        else:
            logger.info("dbfit: only %d exon-count examples for type=%s, keeping base config",
                        len(exon_buckets[t]), t)

        total_splice = sum(splice_counts[t].values())
        if total_splice >= MIN_EXAMPLES_TO_FIT:
            cfg.splice_site_distributions[t] = {k: v / total_splice for k, v in splice_counts[t].items()}
        else:
            logger.info("dbfit: only %d splice-site examples for type=%s, keeping base config",
                        total_splice, t)

    if n_terminal_matched >= MIN_EXAMPLES_TO_FIT:
        cfg.terminal_exon_prob = round((n_uses_first + n_uses_last) / (2 * n_terminal_matched), 4)
        logger.info("dbfit: terminal_exon_prob fit from %d matched records: first=%.1f%%, last=%.1f%%, using=%.4f",
                    n_terminal_matched, 100 * n_uses_first / n_terminal_matched,
                    100 * n_uses_last / n_terminal_matched, cfg.terminal_exon_prob)
    else:
        logger.info("dbfit: only %d records matched to an annotated transcript, keeping base "
                    "config's terminal_exon_prob", n_terminal_matched)

    cfg.validate()
    return cfg
