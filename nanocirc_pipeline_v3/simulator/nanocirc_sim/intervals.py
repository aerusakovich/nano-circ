"""In-memory interval index for region classification (replaces bedtools)."""
from __future__ import annotations

from collections import defaultdict

import numpy as np


class IntervalUnion:
    """A set of disjoint, sorted (start, end) intervals on one (chrom, strand)."""

    __slots__ = ("starts", "ends")

    def __init__(self, starts: np.ndarray, ends: np.ndarray):
        self.starts = starts
        self.ends = ends

    @classmethod
    def from_intervals(cls, intervals) -> "IntervalUnion":
        iv = sorted(intervals)
        if not iv:
            return cls(np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))
        merged = [[iv[0][0], iv[0][1]]]
        for s, e in iv[1:]:
            if s <= merged[-1][1]:
                if e > merged[-1][1]:
                    merged[-1][1] = e
            else:
                merged.append([s, e])
        starts = np.array([m[0] for m in merged], dtype=np.int64)
        ends = np.array([m[1] for m in merged], dtype=np.int64)
        return cls(starts, ends)

    def _candidate_range(self, start: int, end: int):
        if len(self.starts) == 0 or end <= start:
            return 0, 0
        lo = int(np.searchsorted(self.ends, start, side="right"))
        hi = int(np.searchsorted(self.starts, end, side="left"))
        return lo, hi

    def overlaps(self, start: int, end: int) -> bool:
        lo, hi = self._candidate_range(start, end)
        return lo < hi

    def coverage(self, start: int, end: int) -> int:
        lo, hi = self._candidate_range(start, end)
        if lo >= hi:
            return 0
        s = np.maximum(self.starts[lo:hi], start)
        e = np.minimum(self.ends[lo:hi], end)
        return int(np.sum(e - s))

    def fully_covers(self, start: int, end: int) -> bool:
        return self.coverage(start, end) == (end - start)


_EMPTY = IntervalUnion.from_intervals([])


class AnnotationIndex:
    """Per (chrom, strand) gene-union and exon-union interval indexes."""

    def __init__(self, annotation):
        gene_ivs = defaultdict(list)
        exon_ivs = defaultdict(list)

        for gene_id, chrom, strand, start, end in annotation.all_genes():
            gene_ivs[(chrom, strand)].append((start, end))

        for exon in annotation.all_exons():
            exon_ivs[(exon.chrom, exon.strand)].append((exon.start, exon.end))

        self._genes = {key: IntervalUnion.from_intervals(v) for key, v in gene_ivs.items()}
        self._exons = {key: IntervalUnion.from_intervals(v) for key, v in exon_ivs.items()}

        gene_records = defaultdict(list)
        for gene_id, chrom, strand, start, end in annotation.all_genes():
            gene_records[(chrom, strand)].append((start, end, gene_id))
        self._gene_arrays = {}
        for key, recs in gene_records.items():
            recs.sort()
            self._gene_arrays[key] = (
                np.array([r[0] for r in recs], dtype=np.int64),
                np.array([r[1] for r in recs], dtype=np.int64),
                [r[2] for r in recs],
            )

    def genes(self, chrom: str, strand: str) -> IntervalUnion:
        return self._genes.get((chrom, strand), _EMPTY)

    def exons(self, chrom: str, strand: str) -> IntervalUnion:
        return self._exons.get((chrom, strand), _EMPTY)

    def overlapping_genes(self, chrom: str, strand: str, start: int, end: int) -> list:
        """gene_ids on (chrom, strand) overlapping [start, end). Only used by dbfit.py."""
        arrs = self._gene_arrays.get((chrom, strand))
        if arrs is None:
            return []
        starts, ends, ids = arrs
        hi = int(np.searchsorted(starts, end, side="left"))
        return [ids[i] for i in range(hi) if ends[i] > start]

    def is_intergenic(self, chrom: str, start: int, end: int) -> bool:
        """True if the region overlaps no annotated gene on either strand."""
        return not self.genes(chrom, "+").overlaps(start, end) and not self.genes(chrom, "-").overlaps(start, end)

    def classify_region(self, chrom: str, strand: str, start: int, end: int) -> str:
        """Classify a region as ecirna/eicirna/cirna/intergenic."""
        if not self.genes(chrom, strand).overlaps(start, end):
            return "intergenic"
        exon_cov = self.exons(chrom, strand).coverage(start, end)
        if exon_cov == 0:
            return "cirna"
        if exon_cov == (end - start):
            return "ecirna"
        return "eicirna"

    def classify_region_by_blocks(self, chrom: str, strand: str, blocks: list) -> str:
        """Same as classify_region, but for real BED12 records (coverage summed per block)."""
        span_start, span_end = blocks[0][0], blocks[-1][1]
        if not self.genes(chrom, strand).overlaps(span_start, span_end):
            return "intergenic"
        exon_union = self.exons(chrom, strand)
        total_len = sum(e - s for s, e in blocks)
        covered = sum(exon_union.coverage(s, e) for s, e in blocks)
        if covered == 0:
            return "cirna"
        if covered == total_len:
            return "ecirna"
        return "eicirna"

    def classify_final_type(self, chrom: str, strand: str, blocks: list, base_type: str = None) -> str:
        """The type a circRNA's *final* coordinates would get from the real
        pipeline's classifier (bin/circrna_clean.py::classify_types()).
        `base_type` skips recomputing classify_region_by_blocks's result if the
        caller already has it."""
        if base_type is None:
            base_type = self.classify_region_by_blocks(chrom, strand, blocks)
        if base_type != "intergenic":
            return base_type
        # antisense: overlaps one or more exons of a gene on the OPPOSITE strand.
        opposite_strand = "-" if strand == "+" else "+"
        exon_union_opposite = self.exons(chrom, opposite_strand)
        covered_opposite = sum(exon_union_opposite.coverage(s, e) for s, e in blocks)
        return "antisense" if covered_opposite > 0 else "intergenic"
