"""Indexed genome access via pysam.FastaFile."""
from __future__ import annotations

import logging

import pysam

logger = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


class GenomeLoader:
    def __init__(self, fasta_path: str):
        try:
            self._fasta = pysam.FastaFile(fasta_path)
        except OSError as exc:
            raise OSError(
                f"Could not open/index genome FASTA {fasta_path!r} "
                f"(pysam needs write access to create a .fai index if one doesn't exist): {exc}"
            ) from exc
        self.chrom_lengths = dict(zip(self._fasta.references, self._fasta.lengths))
        logger.info("Loaded genome index %s: %d sequences", fasta_path, len(self.chrom_lengths))

    def fetch(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        """0-based half-open [start, end) fetch, clamped to chromosome bounds."""
        chrom_len = self.chrom_lengths.get(chrom)
        if chrom_len is None:
            logger.warning("Chromosome %s not found in genome", chrom)
            return ""
        start = max(0, start)
        end = min(chrom_len, end)
        if end <= start:
            return ""
        seq = self._fasta.fetch(chrom, start, end).upper()
        if strand == "-":
            seq = reverse_complement(seq)
        return seq

    def close(self):
        self._fasta.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
