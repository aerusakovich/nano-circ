"""The one record type every generator produces and every writer takes."""
from __future__ import annotations

from dataclasses import dataclass, field

from .blocks import Block, total_length


@dataclass
class Candidate:
    circ_type: str
    chrom: str
    start: int  # 0-based half-open genomic span
    end: int
    strand: str
    blocks: list = field(default_factory=list)  # list[Block], ascending genomic order
    mature_seq: str = ""
    gene_id: str = "NA"
    tx_id: str = "NA"
    locus_id: str = "NA"
    isoform_index: int = 0
    splice_type: str = "other"
    skipped_exons: int = 0
    retained_introns: int = 0
    exon_count: int = 1
    copy_number: int = 1
    pseudo_seq: str = ""
    tpm: float = 0.0
    n_reads_target: int = 0
    n_reads_confirmed: "int | None" = None
    id: str = ""
    # Exon-window key for forcing the same backsplice boundary on a later isoform.
    boundary_key: "tuple | None" = None

    @property
    def mature_length(self) -> int:
        return len(self.mature_seq) if self.mature_seq else total_length(self.blocks)

    @property
    def block_count(self) -> int:
        return len(self.blocks)
