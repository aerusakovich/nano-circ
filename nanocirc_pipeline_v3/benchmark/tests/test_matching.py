import shutil

import pandas as pd
import pytest

from nanocirc_benchmark.bed_io import BED12_COLS
from nanocirc_benchmark.config import BenchmarkConfig
from nanocirc_benchmark.matching import (
    StructuralMatcher, collapse_near_duplicates, isoform_match_level_params, match_pred_to_pred, match_to_gt,
)


def test_isoform_match_level_params_full_is_exact_and_100pct():
    cfg = BenchmarkConfig.default()
    tol, min_ov = isoform_match_level_params(cfg, "full")
    assert (tol, min_ov) == (0, 1.0)


def test_isoform_match_level_params_relaxed_is_90pct():
    cfg = BenchmarkConfig.default()
    tol, min_ov = isoform_match_level_params(cfg, "relaxed")
    assert tol == cfg.gt_match_tolerance
    assert min_ov == 0.90


def test_isoform_match_level_params_unknown_tag_raises():
    cfg = BenchmarkConfig.default()
    try:
        isoform_match_level_params(cfg, "nonexistent")
        assert False, "expected ValueError"
    except ValueError:
        pass

BEDTOOLS = shutil.which("bedtools")


def test_collapse_near_duplicates_empty():
    assert collapse_near_duplicates([]) == 0


def test_collapse_near_duplicates_takes_max_not_sum():
    # Real observed case: two near-duplicate boundary calls both reporting
    # 187 reads for the same true locus. Summing would double-count it.
    assert collapse_near_duplicates([(187, "a"), (187, "b")]) == 187


def test_collapse_near_duplicates_single():
    assert collapse_near_duplicates([(42, "x")]) == 42


def _gt_df():
    return pd.DataFrame([
        {"bsj_id": "chr1:100-200:+", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"},
        {"bsj_id": "chr1:1000-2000:+", "chrom": "chr1", "start": 1000, "end": 2000, "strand": "+"},
    ])


def test_match_to_gt_within_tolerance():
    call_df = pd.DataFrame([{"chrom": "chr1", "start": 103, "end": 197, "strand": "+", "n_reads": 10}])
    matched = match_to_gt(call_df, _gt_df(), tolerance=5)
    assert matched == {"chr1:100-200:+": [(10, None)]}


def test_match_to_gt_outside_tolerance_no_match():
    call_df = pd.DataFrame([{"chrom": "chr1", "start": 120, "end": 200, "strand": "+", "n_reads": 10}])
    matched = match_to_gt(call_df, _gt_df(), tolerance=5)
    assert matched == {}


def test_match_to_gt_wrong_strand_no_match():
    call_df = pd.DataFrame([{"chrom": "chr1", "start": 100, "end": 200, "strand": "-", "n_reads": 10}])
    matched = match_to_gt(call_df, _gt_df(), tolerance=5)
    assert matched == {}


def test_match_to_gt_unstranded_gt_matches_either_caller_strand():
    # Intergenic GT rows are strand '.' (a coin flip at simulation time, not
    # a real property), so they must match a caller's call on either strand.
    # A real +/- GT row still requires an exact match
    # (see test_match_to_gt_wrong_strand_no_match).
    gt_df = pd.DataFrame([
        {"bsj_id": "chr1:100-200:.", "chrom": "chr1", "start": 100, "end": 200, "strand": "."},
    ])
    for caller_strand in ("+", "-"):
        call_df = pd.DataFrame([{"chrom": "chr1", "start": 103, "end": 197, "strand": caller_strand, "n_reads": 10}])
        matched = match_to_gt(call_df, gt_df, tolerance=5)
        assert matched == {"chr1:100-200:.": [(10, None)]}


def test_match_to_gt_mixed_stranded_and_unstranded_gt():
    # A stranded GT row still only matches its own strand; an unstranded
    # ('.') GT row at a different position matches regardless of caller strand.
    gt_df = pd.DataFrame([
        {"bsj_id": "chr1:100-200:+", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"},
        {"bsj_id": "chr1:1000-2000:.", "chrom": "chr1", "start": 1000, "end": 2000, "strand": "."},
    ])
    call_df = pd.DataFrame([
        {"chrom": "chr1", "start": 100, "end": 200, "strand": "-", "n_reads": 5},   # wrong strand for the stranded row
        {"chrom": "chr1", "start": 1000, "end": 2000, "strand": "-", "n_reads": 8},  # matches the unstranded row regardless
    ])
    matched = match_to_gt(call_df, gt_df, tolerance=5)
    assert matched == {"chr1:1000-2000:.": [(8, None)]}


def test_match_pred_to_pred_unstranded_ref_matches_either_caller_strand():
    ref_df = pd.DataFrame([{"chrom": "chr1", "start": 100, "end": 200, "strand": "."}])
    call_df = pd.DataFrame([{"chrom": "chr1", "start": 101, "end": 199, "strand": "-", "n_reads": 7}])
    pairs = match_pred_to_pred(call_df, ref_df, tolerance=5)
    assert pairs == [(0, 7)]


def test_match_to_gt_nearest_tiebreak():
    gt_df = pd.DataFrame([
        {"bsj_id": "a", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"},
        {"bsj_id": "b", "chrom": "chr1", "start": 102, "end": 200, "strand": "+"},
    ])
    call_df = pd.DataFrame([{"chrom": "chr1", "start": 101, "end": 200, "strand": "+", "n_reads": 5}])
    matched = match_to_gt(call_df, gt_df, tolerance=5)
    assert "a" in matched and "b" not in matched


def test_match_pred_to_pred_returns_ref_row_index():
    ref_df = pd.DataFrame([{"chrom": "chr1", "start": 100, "end": 200, "strand": "+"}])
    call_df = pd.DataFrame([{"chrom": "chr1", "start": 101, "end": 199, "strand": "+", "n_reads": 7}])
    pairs = match_pred_to_pred(call_df, ref_df, tolerance=5)
    assert pairs == [(0, 7)]


# ── StructuralMatcher (needs a real bedtools binary) ────────────────────────


def _bed_df(rows):
    df = pd.DataFrame(rows, columns=BED12_COLS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["row_id"] = range(len(df))
    df["spliced_len"] = df["blockSizes"].str.rstrip(",").apply(lambda s: sum(int(x) for x in s.split(",") if x))
    return df


@pytest.fixture
def matcher(tmp_path):
    genome = tmp_path / "genome.sizes"
    genome.write_text("chr1\t1000000\nchr2\t1000000\n")
    return StructuralMatcher(BEDTOOLS or "bedtools", genome_sizes_file=genome)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_strict_bsj_counts_exact_match(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    tp, fp, fn = matcher.strict_bsj_counts(caller_df, gt_df, tmp_path)
    assert (tp, fp, fn) == (1, 0, 0)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_strict_bsj_counts_boundary_off_by_one_is_fp_and_fn(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "101", "500", "c1", "0", "+", "101", "500", "0", "1", "399,", "0,"]])
    tp, fp, fn = matcher.strict_bsj_counts(caller_df, gt_df, tmp_path)
    assert (tp, fp, fn) == (0, 1, 1)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_relaxed_bsj_counts_within_tolerance_matches(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "103", "497", "c1", "0", "+", "103", "497", "0", "1", "394,", "0,"]])
    tp, fp, fn = matcher.relaxed_bsj_counts(caller_df, gt_df, tmp_path, tolerance=5)
    assert (tp, fp, fn) == (1, 0, 0)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_exon_based_counts_distinguishes_isoform_siblings(matcher, tmp_path):
    # Same BSJ, two distinct structures (as GT). The caller matches only
    # the first one exactly. exon_based must count this as 1 TP + 1 FN (the
    # unmatched isoform sibling), never silently merge them into one.
    gt_df = _bed_df([
        ["chr1", "100", "500", "gt_main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "gt_iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
    ])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"]])
    tp, fp, fn = matcher.exon_based_counts(caller_df, gt_df, tmp_path, min_overlap=0.75)
    assert (tp, fp, fn) == (1, 0, 1)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_isoform_level_counts_matches_when_bsj_and_exon_both_close(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "102", "498", "c1", "0", "+", "102", "498", "0", "1", "396,", "0,"]])
    tp, fp, fn = matcher.isoform_level_counts(caller_df, gt_df, tmp_path, tolerance=5, min_overlap=0.99)
    assert (tp, fp, fn) == (1, 0, 0)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_isoform_level_counts_rejects_close_bsj_but_low_exon_overlap(matcher, tmp_path):
    # BSJ boundaries are within tolerance, but the caller's own exon
    # structure only covers half the GT isoform's spliced length. exon_based
    # alone (a looser threshold) would call this a match; isoform_level's
    # 99% floor must not.
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "+", "100", "500", "0", "1", "200,", "0,"]])
    tp, fp, fn = matcher.isoform_level_counts(caller_df, gt_df, tmp_path, tolerance=5, min_overlap=0.99)
    assert (tp, fp, fn) == (0, 1, 1)
    # but exon_based alone (75% threshold) is unaffected and would still miss
    # this too (200/400=50% < 75%), confirming the fixture itself is sane
    tp2, fp2, fn2 = matcher.exon_based_counts(caller_df, gt_df, tmp_path, min_overlap=0.75)
    assert (tp2, fp2, fn2) == (0, 1, 1)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_isoform_level_counts_rejects_high_exon_overlap_but_far_bsj(matcher, tmp_path):
    # A 50bp shift leaves 350/400=87.5% reciprocal overlap -- well outside
    # BSJ tolerance, but still enough to satisfy a permissive (0.80) exon
    # threshold. exon_based alone doesn't require the outer span to be
    # close and would match at that threshold; isoform_level must not.
    gt_df = _bed_df([["chr1", "1000", "1400", "gt1", "0", "+", "1000", "1400", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "1050", "1450", "c1", "0", "+", "1050", "1450", "0", "1", "400,", "0,"]])
    tp, fp, fn = matcher.isoform_level_counts(caller_df, gt_df, tmp_path, tolerance=5, min_overlap=0.80)
    assert (tp, fp, fn) == (0, 1, 1)
    tp2, fp2, fn2 = matcher.exon_based_counts(caller_df, gt_df, tmp_path, min_overlap=0.80)
    assert (tp2, fp2, fn2) == (1, 0, 0)  # confirms exon_based alone would have matched


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_isoform_level_counts_distinguishes_isoform_siblings(matcher, tmp_path):
    # Same fixture as test_exon_based_counts_distinguishes_isoform_siblings:
    # isoform_level must show the same TP/FN split, not silently merge the
    # two same-BSJ structures the way relaxed_bsj would.
    gt_df = _bed_df([
        ["chr1", "100", "500", "gt_main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "gt_iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
    ])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"]])
    tp, fp, fn = matcher.isoform_level_counts(caller_df, gt_df, tmp_path, tolerance=5, min_overlap=0.99)
    assert (tp, fp, fn) == (1, 0, 1)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_isoform_level_pairs_excludes_far_bsj_pair(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "1000", "1400", "gt1", "0", "+", "1000", "1400", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "1050", "1450", "c1", "0", "+", "1050", "1450", "0", "1", "400,", "0,"]])
    pairs = matcher.isoform_level_pairs(caller_df, gt_df, tmp_path, tolerance=5, min_overlap=0.80)
    assert pairs == []
    # the 50bp shift leaves 350/400=87.5% reciprocal overlap: exon_based
    # alone (a lower threshold) would still match; isoform_level's BSJ
    # tolerance check must reject it regardless.
    exon_pairs = matcher.exon_based_pairs(caller_df, gt_df, tmp_path, min_overlap=0.80)
    assert exon_pairs == [(0, 0)]  # confirms exon_based alone would have matched


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_isoform_level_pairs_includes_close_bsj_and_exon_pair(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", "+", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "102", "498", "c1", "0", "+", "102", "498", "0", "1", "396,", "0,"]])
    pairs = matcher.isoform_level_pairs(caller_df, gt_df, tmp_path, tolerance=5, min_overlap=0.99)
    assert pairs == [(0, 0)]


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_relaxed_bsj_counts_return_ids_expands_to_all_isoform_siblings(matcher, tmp_path):
    # Regression test: relaxed_bsj_counts collapses GT isoform siblings
    # sharing one BSJ to a single entity internally (a BSJ-only measure),
    # but return_ids=True must still report every original GT row sharing
    # that matched BSJ as matched, not just whichever isoform the internal
    # collapse happened to keep. A per-row consumer (sensitivity.py) needs
    # every isoform sibling of a detected BSJ to show as detected.
    gt_df = _bed_df([
        ["chr1", "100", "500", "gt_main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "gt_iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
    ])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"]])
    tp, fp, fn, matched_gt_ids = matcher.relaxed_bsj_counts(caller_df, gt_df, tmp_path, tolerance=5, return_ids=True)
    assert fn == 0
    assert matched_gt_ids == {0, 1}


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_strict_and_relaxed_do_not_multiply_count_isoform_siblings(matcher, tmp_path):
    # A GT locus with 2 isoform-sibling rows sharing one BSJ must still only
    # count as one required detection for a BSJ-only measure (strict or
    # relaxed), not two, since neither measure looks at internal structure.
    gt_df = _bed_df([
        ["chr1", "100", "500", "gt_main", "0", "+", "100", "500", "0", "2", "50,50,", "0,350,"],
        ["chr1", "100", "500", "gt_iso1", "0", "+", "100", "500", "0", "1", "400,", "0,"],
    ])
    empty_caller_df = _bed_df([])
    tp, fp, fn = matcher.strict_bsj_counts(empty_caller_df, gt_df, tmp_path)
    assert fn == 1  # one missed BSJ, not two
    tp, fp, fn = matcher.relaxed_bsj_counts(empty_caller_df, gt_df, tmp_path, tolerance=5)
    assert fn == 1


# Unstranded ('.') GT rows are intergenic circRNAs. Strand is random at
# simulation time, not a real property, so they must match a caller call on
# either strand. Without this, every intergenic locus scored 0 true
# positives, for every caller and tier.


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_strict_bsj_counts_unstranded_gt_matches_either_caller_strand(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", ".", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "-", "100", "500", "0", "1", "400,", "0,"]])
    tp, fp, fn = matcher.strict_bsj_counts(caller_df, gt_df, tmp_path)
    assert (tp, fp, fn) == (1, 0, 0)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_relaxed_bsj_counts_unstranded_gt_matches_either_caller_strand(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", ".", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "103", "497", "c1", "0", "+", "103", "497", "0", "1", "394,", "0,"]])
    tp, fp, fn = matcher.relaxed_bsj_counts(caller_df, gt_df, tmp_path, tolerance=5)
    assert (tp, fp, fn) == (1, 0, 0)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_exon_based_counts_unstranded_gt_matches_either_caller_strand(matcher, tmp_path):
    gt_df = _bed_df([["chr1", "100", "500", "gt1", "0", ".", "100", "500", "0", "1", "400,", "0,"]])
    caller_df = _bed_df([["chr1", "100", "500", "c1", "0", "-", "100", "500", "0", "1", "400,", "0,"]])
    tp, fp, fn = matcher.exon_based_counts(caller_df, gt_df, tmp_path, min_overlap=0.9)
    assert (tp, fp, fn) == (1, 0, 0)


@pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not installed")
def test_strict_bsj_counts_mixed_stranded_and_unstranded_gt(matcher, tmp_path):
    # A stranded GT row still requires an exact strand match; a separate
    # unstranded GT row matches regardless of caller strand. Both must be
    # counted correctly in the same run (no cross-contamination between the
    # two matching passes).
    gt_df = _bed_df([
        ["chr1", "100", "500", "gt_stranded", "0", "+", "100", "500", "0", "1", "400,", "0,"],
        ["chr2", "100", "500", "gt_unstranded", "0", ".", "100", "500", "0", "1", "400,", "0,"],
    ])
    caller_df = _bed_df([
        ["chr1", "100", "500", "c_wrong_strand", "0", "-", "100", "500", "0", "1", "400,", "0,"],
        ["chr2", "100", "500", "c_matches_unstranded", "0", "-", "100", "500", "0", "1", "400,", "0,"],
    ])
    tp, fp, fn = matcher.strict_bsj_counts(caller_df, gt_df, tmp_path)
    # c_wrong_strand: FP (doesn't match gt_stranded, wrong strand)
    # c_matches_unstranded: TP (matches gt_unstranded regardless of strand)
    # gt_stranded: FN (never matched)
    assert (tp, fp, fn) == (1, 1, 1)
