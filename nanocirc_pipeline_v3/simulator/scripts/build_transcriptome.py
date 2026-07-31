#!/usr/bin/env python3
"""Extract a spliced transcriptome FASTA from a GTF + genome."""
import argparse
import sys

sys.path.insert(0, "..")  # allow running as `python scripts/build_transcriptome.py` from simulator/new

from nanocirc_sim.annotation import Annotation
from nanocirc_sim.blocks import Block, assemble_sequence
from nanocirc_sim.genome import GenomeLoader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    annotation = Annotation.from_gtf(args.gtf)
    genome = GenomeLoader(args.genome)

    n_written = 0
    with open(args.out, "w") as fh:
        for tx_map in annotation.tx_structure.values():
            for tx_id, tx in tx_map.items():
                blocks = [Block(e.start, e.end, "exon") for e in tx.exons]
                seq = assemble_sequence(genome, tx.chrom, blocks, tx.strand)
                if not seq:
                    continue
                fh.write(f">{tx_id}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
                n_written += 1
    genome.close()
    print(f"Wrote {n_written} transcript sequences to {args.out}")


if __name__ == "__main__":
    main()
