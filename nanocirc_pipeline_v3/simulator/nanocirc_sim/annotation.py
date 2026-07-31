"""GTF parsing into Gene/Transcript/Exon structures. Coordinates are
converted to 0-based half-open on load."""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


@dataclass(frozen=True)
class Exon:
    chrom: str
    start: int  # 0-based
    end: int  # half-open
    strand: str
    gene_id: str
    transcript_id: str


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    gene_biotype: str
    gene_name: str
    exons: list = field(default_factory=list)  # sorted by start, ascending genomic order

    @property
    def start(self) -> int:
        return self.exons[0].start

    @property
    def end(self) -> int:
        return self.exons[-1].end

    @property
    def introns(self):
        """List of (start, end) 0-based half-open intron intervals between consecutive exons."""
        out = []
        for a, b in zip(self.exons, self.exons[1:]):
            if a.end < b.start:
                out.append((a.end, b.start))
        return out


def _parse_attributes(field_str: str) -> dict:
    return dict(_ATTR_RE.findall(field_str))


class Annotation:
    """Parsed GTF: gene_id -> transcript_id -> Transcript, plus gene-level metadata."""

    def __init__(self):
        self.tx_structure: dict[str, dict[str, Transcript]] = {}
        self.gene_biotype: dict[str, str] = {}
        self.gene_name: dict[str, str] = {}
        self.gene_span: dict[str, tuple] = {}  # gene_id -> (chrom, strand, start, end)

    @classmethod
    def from_gtf(cls, gtf_path: str) -> "Annotation":
        ann = cls()
        n_exon_lines = 0
        with open(gtf_path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "exon":
                    continue
                n_exon_lines += 1
                chrom = fields[0]
                start = int(fields[3]) - 1  # GTF 1-based closed -> 0-based half-open
                end = int(fields[4])
                strand = fields[6]
                attrs = _parse_attributes(fields[8])
                gene_id = attrs.get("gene_id")
                transcript_id = attrs.get("transcript_id")
                if not gene_id or not transcript_id:
                    continue
                biotype = attrs.get("gene_biotype") or attrs.get("gene_type") or ""
                gene_name = attrs.get("gene_name", gene_id)

                ann.gene_biotype.setdefault(gene_id, biotype)
                ann.gene_name.setdefault(gene_id, gene_name)

                tx_map = ann.tx_structure.setdefault(gene_id, {})
                tx = tx_map.get(transcript_id)
                if tx is None:
                    tx = Transcript(
                        transcript_id=transcript_id,
                        gene_id=gene_id,
                        chrom=chrom,
                        strand=strand,
                        gene_biotype=biotype,
                        gene_name=gene_name,
                    )
                    tx_map[transcript_id] = tx
                tx.exons.append(Exon(chrom, start, end, strand, gene_id, transcript_id))

        n_genes = 0
        for gene_id, tx_map in ann.tx_structure.items():
            n_genes += 1
            gspan_start = None
            gspan_end = None
            chrom = None
            strand = None
            for tx in tx_map.values():
                tx.exons.sort(key=lambda e: e.start)
                chrom = tx.chrom
                strand = tx.strand
                gspan_start = tx.start if gspan_start is None else min(gspan_start, tx.start)
                gspan_end = tx.end if gspan_end is None else max(gspan_end, tx.end)
            ann.gene_span[gene_id] = (chrom, strand, gspan_start, gspan_end)

        logger.info(
            "Parsed GTF %s: %d exon records, %d genes, %d transcripts",
            gtf_path,
            n_exon_lines,
            n_genes,
            sum(len(v) for v in ann.tx_structure.values()),
        )
        return ann

    def is_protein_coding(self, gene_id: str) -> bool:
        return self.gene_biotype.get(gene_id, "") == "protein_coding"

    def is_antisense(self, gene_id: str, biotypes, name_pattern: str | None) -> bool:
        biotype = self.gene_biotype.get(gene_id, "")
        if biotype in biotypes:
            return True
        if name_pattern:
            name = self.gene_name.get(gene_id, "")
            if re.search(name_pattern, name, flags=re.IGNORECASE):
                return True
        return False

    def get_transcripts(self, min_exons: int = 2, biotype_filter=None, name_pattern: str | None = None,
                         protein_coding_only: bool = False):
        """Return list of (gene_id, transcript_id) matching the given filters."""
        out = []
        for gene_id, tx_map in self.tx_structure.items():
            if protein_coding_only and not self.is_protein_coding(gene_id):
                continue
            if biotype_filter is not None:
                if not self.is_antisense(gene_id, biotype_filter, name_pattern):
                    continue
            for tx_id, tx in tx_map.items():
                if len(tx.exons) >= min_exons:
                    out.append((gene_id, tx_id))
        return out

    def all_genes(self):
        """Yield (gene_id, chrom, strand, start, end) for every gene."""
        for gene_id, (chrom, strand, start, end) in self.gene_span.items():
            yield gene_id, chrom, strand, start, end

    def all_exons(self):
        """Yield Exon for every exon in every transcript (may contain duplicates across isoforms)."""
        for tx_map in self.tx_structure.values():
            for tx in tx_map.values():
                yield from tx.exons
