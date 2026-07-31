import random

import pytest


def _write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


@pytest.fixture
def mid_fixture(tmp_path):
    """One protein-coding gene (chrT, +, 8 exons) and one antisense-biotype
    gene (chrT, -, 4 exons), on a 20kb chromosome."""
    rnd = random.Random(3)
    chrom_len = 20000
    seq = [rnd.choice("ACGT") for _ in range(chrom_len)]

    gtf_lines = []
    exon_len, intron_len = 100, 300

    def make_gene(gene_id, gene_name, biotype, strand, n_exons, start_1based, tx_id):
        pos = start_1based
        coords = []
        for _ in range(n_exons):
            coords.append((pos, pos + exon_len - 1))
            pos += exon_len + intron_len
        g_start, g_end = coords[0][0], coords[-1][1]
        gtf_lines.append(f'chrT\ttest\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_biotype "{biotype}";')
        gtf_lines.append(f'chrT\ttest\ttranscript\t{g_start}\t{g_end}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_biotype "{biotype}";')
        for s, e in coords:
            gtf_lines.append(f'chrT\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t'
                              f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_biotype "{biotype}";')

        # plant a canonical GT..AG motif at every internal exon-intron junction
        for (_, e_end_1based), (s_next_1based, _) in zip(coords, coords[1:]):
            end0 = e_end_1based  # 0-based end of the upstream exon (== 1-based end)
            start0 = s_next_1based - 1  # 0-based start of the downstream exon
            if strand == "+":
                seq[end0:end0 + 2] = list("GT")  # donor: first 2bp of the intron
                seq[start0 - 2:start0] = list("AG")  # acceptor: last 2bp of the intron
            else:
                seq[end0:end0 + 2] = list("CT")  # rc -> AG acceptor (transcript-3' end)
                seq[start0 - 2:start0] = list("AC")  # rc -> GT donor (transcript-5' end)
        return coords

    make_gene("G1", "G1", "protein_coding", "+", 8, 1000, "G1T1")
    make_gene("G2-AS1", "G2-AS1", "lncRNA", "-", 4, 12000, "G2T1")

    genome_path = tmp_path / "mid.fa"
    gtf_path = tmp_path / "mid.gtf"
    _write_fasta(genome_path, {"chrT": "".join(seq)})
    gtf_path.write_text("\n".join(gtf_lines) + "\n")

    return {"genome": str(genome_path), "gtf": str(gtf_path)}


@pytest.fixture
def synthetic_db_bed(tmp_path):
    """A tiny 'real database' BED12 with known composition, matching mid_fixture."""
    lines = []
    for i in range(6):
        # genuine 2-exon ecirna: blocks skip the intron between exon0 and exon1
        lines.append(f"chrT\t999\t1499\tdb_eci_{i}\t0\t+\t999\t1499\t0\t2\t100,100,\t0,400,\n")
    for i in range(6):
        # genuine eicirna: single block spanning exon0+intron+exon1 (intron retained)
        lines.append(f"chrT\t999\t1499\tdb_eic_{i}\t0\t+\t999\t1499\t0\t1\t500,\t0,\n")
    for i in range(6):
        # cirna: fully intronic
        lines.append(f"chrT\t1100\t1399\tdb_cir_{i}\t0\t+\t1100\t1399\t0\t1\t299,\t0,\n")
    for i in range(6):
        # intergenic: far outside any gene
        lines.append(f"chrT\t18000\t18300\tdb_ig_{i}\t0\t+\t18000\t18300\t0\t1\t300,\t0,\n")
    for i in range(6):
        # antisense (literature definition): no gene on this circRNA's own
        # strand ('-') at this locus (only G1, on '+', is there), but the
        # block exactly overlaps G1's exon0 [999,1099), an exon of a gene on
        # the OPPOSITE strand.
        lines.append(f"chrT\t999\t1099\tdb_as_{i}\t0\t-\t999\t1099\t0\t1\t100,\t0,\n")
    path = tmp_path / "db.bed12"
    path.write_text("".join(lines))
    return str(path)
