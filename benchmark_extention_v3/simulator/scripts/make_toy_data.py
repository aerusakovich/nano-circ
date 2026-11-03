#!/usr/bin/env python3
"""Builds a tiny synthetic genome + GTF for run_example.sh. No downloads,
no external data: this makes `simulate` runnable in seconds anywhere."""
import argparse
import random
from pathlib import Path


def make_gene(gtf_lines, seq, gene_id, gene_name, biotype, strand, n_exons, start_1based, tx_id,
              exon_len=100, intron_len=300):
    pos = start_1based
    coords = []
    for _ in range(n_exons):
        coords.append((pos, pos + exon_len - 1))
        pos += exon_len + intron_len
    g_start, g_end = coords[0][0], coords[-1][1]
    gtf_lines.append(f'chr1\ttoy\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t'
                      f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_biotype "{biotype}";')
    gtf_lines.append(f'chr1\ttoy\ttranscript\t{g_start}\t{g_end}\t.\t{strand}\t.\t'
                      f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_biotype "{biotype}";')
    for s, e in coords:
        gtf_lines.append(f'chr1\ttoy\texon\t{s}\t{e}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_biotype "{biotype}";')

    # plant a canonical GT..AG motif at every internal exon-intron junction
    for (_, e_end_1based), (s_next_1based, _) in zip(coords, coords[1:]):
        end0 = e_end_1based
        start0 = s_next_1based - 1
        if strand == "+":
            seq[end0:end0 + 2] = list("GT")
            seq[start0 - 2:start0] = list("AG")
        else:
            seq[end0:end0 + 2] = list("CT")
            seq[start0 - 2:start0] = list("AC")
    return coords


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rnd = random.Random(1)
    chrom_len = 20000
    seq = [rnd.choice("ACGT") for _ in range(chrom_len)]
    gtf_lines = []

    make_gene(gtf_lines, seq, "GENE1", "GENE1", "protein_coding", "+", 8, 1000, "GENE1T1")
    make_gene(gtf_lines, seq, "GENE2", "GENE2", "protein_coding", "-", 6, 12000, "GENE2T1")

    (outdir / "toy_genome.fa").write_text(
        ">chr1\n" + "\n".join("".join(seq[i:i + 70]) for i in range(0, chrom_len, 70)) + "\n"
    )
    (outdir / "toy_annotation.gtf").write_text("\n".join(gtf_lines) + "\n")
    print(f"Wrote {outdir / 'toy_genome.fa'} and {outdir / 'toy_annotation.gtf'}")


if __name__ == "__main__":
    main()
