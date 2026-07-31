import shutil

import pandas as pd
import pytest

from nanocirc_benchmark.type_classification import build_feature_beds, classify_types

BEDTOOLS = shutil.which("bedtools")


def _row(chrom, start, end, name, strand, block_sizes):
    return dict(chrom=chrom, start=start, end=end, name=name, score=0, strand=strand,
                thickStart=start, thickEnd=end, itemRgb=0,
                blockCount=len(block_sizes.rstrip(",").split(",")),
                blockSizes=block_sizes, blockStarts="0,")


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_classify_types_all_five_categories(tmp_path):
    gene_bed = tmp_path / "gene.bed"
    exon_bed = tmp_path / "exon.bed"
    # one gene chr1:1000-5000 (+), with 2 exons: 1000-1200, 4800-5000
    gene_bed.write_text("chr1\t1000\t5000\tgeneA\t0\t+\n")
    exon_bed.write_text("chr1\t1000\t1200\tgeneA\t0\t+\nchr1\t4800\t5000\tgeneA\t0\t+\n")

    rows = [
        _row("chr1", 1000, 1200, "ecirna_test", "+", "200,"),       # fully within one exon
        _row("chr1", 1000, 2000, "eicirna_test", "+", "1000,"),     # partial exon coverage
        _row("chr1", 2000, 3000, "cirna_test", "+", "1000,"),       # inside gene, no exon overlap
        _row("chr1", 1000, 2000, "antisense_test", "-", "1000,"),   # opposite-strand gene overlap
        _row("chr1", 100000, 101000, "intergenic_test", "+", "1000,"),  # no gene overlap
    ]
    df = pd.DataFrame(rows)
    df["row_id"] = range(len(df))
    df["spliced_len"] = df["blockSizes"].str.rstrip(",").apply(lambda s: sum(int(x) for x in s.split(",") if x))

    types = classify_types(df, gene_bed, exon_bed, BEDTOOLS, tmp_path / "work", id_prefix="t")

    expected = ["ecirna", "eicirna", "cirna", "antisense", "intergenic"]
    for row_id, exp in zip(df["row_id"], expected):
        assert types[row_id] == exp, f"row {row_id} ({df.loc[row_id, 'name']}): expected {exp}, got {types[row_id]}"


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_classify_types_empty_df_returns_empty(tmp_path):
    gene_bed = tmp_path / "gene.bed"
    exon_bed = tmp_path / "exon.bed"
    gene_bed.write_text("")
    exon_bed.write_text("")
    df = pd.DataFrame(columns=["row_id", "spliced_len", "chrom", "start", "end", "name", "score", "strand",
                                "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"])
    assert classify_types(df, gene_bed, exon_bed, BEDTOOLS, tmp_path / "work") == {}


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_build_feature_beds_merges_overlapping_exons(tmp_path):
    # Regression test: different transcripts of the same gene contribute
    # overlapping exon intervals (T1's exon1 1000-2000 and T2's exon1
    # 1000-2500 both cover chr1:1000-2000). classify_types() sums bedtools
    # intersect -wo overlap bp per query across every matching exon.bed
    # line. Left un-merged, a query overlapping this region would have its
    # covered bases counted twice (once per overlapping exon.bed line),
    # which can push exon_frac above 1.0 and misclassify an eicirna call as
    # ecirna. build_feature_beds() must merge same-strand overlapping exon
    # intervals before writing exon.bed.
    gtf = tmp_path / "test.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "GENEA";\n'
        'chr1\ttest\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "GENEA"; transcript_id "T1";\n'
        'chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "GENEA"; transcript_id "T1";\n'
        'chr1\ttest\texon\t4000\t5000\t.\t+\t.\tgene_id "GENEA"; transcript_id "T1";\n'
        'chr1\ttest\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "GENEA"; transcript_id "T2";\n'
        'chr1\ttest\texon\t1000\t2500\t.\t+\t.\tgene_id "GENEA"; transcript_id "T2";\n'
        'chr1\ttest\texon\t4000\t5000\t.\t+\t.\tgene_id "GENEA"; transcript_id "T2";\n'
    )
    gene_bed = tmp_path / "gene.bed"
    exon_bed = tmp_path / "exon.bed"
    build_feature_beds(gtf, gene_bed, exon_bed, BEDTOOLS)

    exon_lines = [l.split("\t") for l in exon_bed.read_text().strip().split("\n") if l]
    # T1's 1000-2000 and T2's 1000-2500 must have been merged into ONE
    # interval, not left as two separate overlapping lines. (GTF is 1-based;
    # the awk conversion's $4-1 makes BED start 999, not 1000.)
    starts_ends = sorted((int(l[1]), int(l[2])) for l in exon_lines)
    assert starts_ends == [(999, 2500), (3999, 5000)]

    # A call spanning 1000-3000 with a single unspliced block (2000bp) that
    # only partially overlaps annotated exon (1500bp of the 2000bp block is
    # exonic: 1000-2500) must classify as eicirna (partial), not ecirna.
    # If exon.bed were left un-merged, summed overlap would double-count the
    # 1000-2000 region (500+1500=2000bp instead of the true 1500bp unique
    # coverage), inflating exon_frac to 1.0 and misclassifying it as ecirna.
    row = _row("chr1", 1000, 3000, "eicirna_test", "+", "2000,")
    df = pd.DataFrame([row])
    df["row_id"] = range(len(df))
    df["spliced_len"] = df["blockSizes"].str.rstrip(",").apply(lambda s: sum(int(x) for x in s.split(",") if x))
    types = classify_types(df, gene_bed, exon_bed, BEDTOOLS, tmp_path / "work", id_prefix="t")
    assert types[0] == "eicirna"


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_build_feature_beds_skips_if_already_built(tmp_path):
    gene_bed = tmp_path / "gene.bed"
    exon_bed = tmp_path / "exon.bed"
    gene_bed.write_text("chr1\t0\t100\tSENTINEL\t0\t+\n")
    exon_bed.write_text("chr1\t0\t100\tSENTINEL\t0\t+\n")
    gtf = tmp_path / "test.gtf"
    gtf.write_text('chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "GENEA";\n')

    build_feature_beds(gtf, gene_bed, exon_bed, BEDTOOLS)

    # Pre-existing files must be left untouched, not rebuilt from the GTF.
    assert "SENTINEL" in gene_bed.read_text()
    assert "SENTINEL" in exon_bed.read_text()
