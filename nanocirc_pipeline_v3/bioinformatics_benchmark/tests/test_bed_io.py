import pandas as pd

from nanocirc_benchmark.bed_io import BED12_COLS, dedup_by_bsj, read_bed


def _make_df(rows):
    df = pd.DataFrame(rows, columns=BED12_COLS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["row_id"] = range(len(df))
    return df


def test_read_bed_collapses_byte_identical_duplicate_lines(tmp_path):
    # Two literally identical lines (same everything, including block
    # structure). A real duplicate line, not isoform diversity.
    bed = tmp_path / "calls.bed12"
    bed.write_text(
        "chr1\t100\t500\tcirc1\t0\t+\t100\t500\t0\t1\t400,\t0,\n"
        "chr1\t100\t500\tcirc1\t0\t+\t100\t500\t0\t1\t400,\t0,\n"
    )
    df = read_bed(bed)
    assert len(df) == 1


def test_read_bed_preserves_distinct_isoforms_sharing_one_bsj(tmp_path):
    # Same BSJ (chr1:100-500:+), two genuinely different exon structures.
    # Must not be collapsed. The old dedup key was (chrom, start, end,
    # strand) only, which silently dropped one of these.
    bed = tmp_path / "calls.bed12"
    bed.write_text(
        "chr1\t100\t500\tcirc1_main\t0\t+\t100\t500\t0\t2\t50,50,\t0,350,\n"
        "chr1\t100\t500\tcirc1_iso1\t0\t+\t100\t500\t0\t1\t400,\t0,\n"
    )
    df = read_bed(bed)
    assert len(df) == 2
    assert set(df["blockSizes"]) == {"50,50,", "400,"}


def test_read_bed_dedup_false_keeps_all_rows(tmp_path):
    bed = tmp_path / "calls.bed12"
    bed.write_text(
        "chr1\t100\t500\tcirc1\t0\t+\t100\t500\t0\t1\t400,\t0,\n"
        "chr1\t100\t500\tcirc1\t0\t+\t100\t500\t0\t1\t400,\t0,\n"
    )
    df = read_bed(bed, dedup=False)
    assert len(df) == 2


def test_dedup_by_bsj_collapses_isoform_siblings_to_one_row():
    df = _make_df([
        ["chr1", "100", "500", "main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
        ["chr2", "1000", "1200", "other", "0", "-", "1000", "1200", "0", "1", "200,", "0,"],
    ])
    out = dedup_by_bsj(df)
    assert len(out) == 2
    assert set(zip(out["chrom"], out["start"], out["end"], out["strand"])) == {
        ("chr1", 100, 500, "+"), ("chr2", 1000, 1200, "-"),
    }


def test_dedup_by_bsj_keeps_row_id_values_not_renumbered():
    df = _make_df([
        ["chr1", "100", "500", "main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
        ["chr2", "1000", "1200", "other", "0", "-", "1000", "1200", "0", "1", "200,", "0,"],
    ])
    out = dedup_by_bsj(df)
    # row_id 0 (the "main" isoform) should survive as the kept representative;
    # row_id values must reflect the ORIGINAL frame, not be renumbered 0..n-1.
    assert list(out["row_id"]) == [0, 2]


def test_dedup_by_bsj_empty_df_returns_empty():
    df = _make_df([])
    out = dedup_by_bsj(df)
    assert len(out) == 0
