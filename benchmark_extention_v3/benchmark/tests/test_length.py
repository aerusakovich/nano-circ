from nanocirc_benchmark.bed_io import read_bed


def test_spliced_len_sums_blocks(tmp_path):
    bed = tmp_path / "calls.bed12"
    bed.write_text(
        "chr1\t100\t500\tcirc1\t0\t+\t100\t500\t0\t2\t50,80,\t0,320,\n"
        "chr2\t1000\t1200\tcirc2\t0\t-\t1000\t1200\t0\t1\t200,\t0,\n"
    )
    df = read_bed(bed)
    assert list(df["spliced_len"]) == [130, 200]


def test_spliced_len_skips_header_and_blank_lines(tmp_path):
    bed = tmp_path / "calls.bed12"
    bed.write_text(
        "# comment\n"
        "track name=foo\n"
        "\n"
        "chr1\t100\t500\tcirc1\t0\t+\t100\t500\t0\t1\t400,\t0,\n"
    )
    df = read_bed(bed)
    assert list(df["spliced_len"]) == [400]
