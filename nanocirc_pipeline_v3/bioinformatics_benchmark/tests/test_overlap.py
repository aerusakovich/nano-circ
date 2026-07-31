import pandas as pd

from nanocirc_benchmark.bed_io import BED12_COLS
from nanocirc_benchmark.overlap import _gt_key


def _gt_df_with_isoform_siblings():
    # Same BSJ (chr1:100-500:+), two genuinely distinct exon structures.
    rows = [
        ["chr1", "100", "500", "gt_main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "gt_iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
    ]
    df = pd.DataFrame(rows, columns=BED12_COLS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def test_gt_key_relaxed_bsj_collapses_isoform_siblings_to_one_identity():
    # relaxed_bsj is a BSJ-only measure: isoform siblings sharing one BSJ
    # really are "the same relaxed-BSJ entity" for pooling purposes.
    df = _gt_df_with_isoform_siblings()
    keys = set(_gt_key(df, "relaxed_bsj"))
    assert keys == {"chr1:100-500:+"}


def test_gt_key_exon_based_keeps_isoform_siblings_distinct():
    # exon_based_counts() treats each distinct structure as its own entity
    # (via row_id), so the pooling key must match. Otherwise two callers
    # each matching a different isoform of the same BSJ would wrongly look
    # like they "detected the same thing" in the UpSet matrix.
    df = _gt_df_with_isoform_siblings()
    keys = set(_gt_key(df, "exon_based"))
    assert len(keys) == 2


def test_gt_key_exon_based_same_structure_same_key():
    rows = [
        ["chr1", "100", "500", "a", "0", "+", "100", "500", "0", "1", "400,", "0,"],
        ["chr2", "1000", "1200", "b", "0", "-", "1000", "1200", "0", "1", "200,", "0,"],
    ]
    df = pd.DataFrame(rows, columns=BED12_COLS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    keys = list(_gt_key(df, "exon_based"))
    assert len(set(keys)) == 2
    assert len(keys) == 2
