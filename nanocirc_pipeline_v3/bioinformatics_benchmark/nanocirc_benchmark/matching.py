"""Structural circRNA matching (3 strategies), plus coordinate-tolerance
matching for read-count comparisons. Uses bedtools intersect approach.

Strategies:
  strict_bsj   exact whole-span match (`bedtools intersect -s -f 1.0 -r`).
  relaxed_bsj  both boundaries within `tolerance` bp of each other.
  exon_based   reciprocal exon overlap >= min_overlap (`-wo -split`).

GT rows with strand '.' are intergenic circRNAs (strand is random at
simulation time, so we can't judge based upon it). These match a caller's call on
either strand. Every function below splits GT into a stranded subset
(matched with -s) and an unstranded subset (matched with no strand flag),
then combines both results.

Usage as a library:
    from nanocirc_benchmark.matching import StructuralMatcher
    m = StructuralMatcher(bedtools_bin="/path/to/bedtools")
    tp, fp, fn = m.relaxed_bsj_counts(caller_df, gt_df, work_dir, tolerance=10)
"""
import subprocess
from pathlib import Path

import pandas as pd

from .bed_io import BED12_COLS, write_bed_with_id, dedup_by_bsj, Bedtools

_STRANDED = ("+", "-")

STRATEGY_TITLES = {"strict_bsj": "Strict BSJ", "relaxed_bsj": "Relaxed BSJ", "exon_based": "Exon"}


def strategy_footnote(strategy: str, config) -> str:
    """Short, plot-footnote-length description of what a strategy's match
    threshold actually is, for whichever strategy is plotted."""
    if strategy == "strict_bsj":
        return "exact coordinate match"
    if strategy == "exon_based":
        return f">={config.tp_min_overlap:.0%} exon overlap matched"
    return f"coordinate-tolerance matched, +/-{config.gt_match_tolerance}bp"


def _expand_bsj_match(matched_ids_prefixed, id_prefix, deduped_df, orig_df):
    """Expands a BSJ match (one row per BSJ) back to every original row
    sharing that BSJ, so every isoform of a matched BSJ counts as matched."""
    matched_row_ids = {int(x[len(id_prefix):]) for x in matched_ids_prefixed}
    matched_keys = set(
        deduped_df.loc[deduped_df["row_id"].isin(matched_row_ids), ["chrom", "start", "end", "strand"]]
        .itertuples(index=False, name=None)
    )
    orig_keys = list(orig_df[["chrom", "start", "end", "strand"]].itertuples(index=False, name=None))
    return {rid for rid, key in zip(orig_df["row_id"], orig_keys) if key in matched_keys}


def _split_by_strand(df):
    """Returns (stranded_subset, unstranded_subset). Unstranded rows (strand
    '.') need strand-agnostic matching, see module docstring."""
    is_stranded = df["strand"].isin(_STRANDED)
    return df[is_stranded], df[~is_stranded]


def _strand_passes(gt_df):
    """Yields (gt_subset, strand_flag_or_None, tag) for the two matching
    passes every counting function needs: stranded GT rows matched with -s,
    unstranded GT rows ('.') matched with no strand flag at all."""
    stranded, unstranded = _split_by_strand(gt_df)
    yield stranded, "-s", "stranded"
    yield unstranded, None, "unstranded"


class StructuralMatcher:
    """genome_sizes_file: a bedtools-style chrom-sizes file (chrom<TAB>length
    per line), required by relaxed_bsj_counts()'s use of `bedtools slop`."""

    def __init__(self, bedtools_bin="bedtools", genome_sizes_file=None):
        self.bt = Bedtools(bedtools_bin)
        self.genome_sizes_file = genome_sizes_file

    def strict_bsj_counts(self, caller_df, gt_df, work: Path, return_ids=False):
        # BSJ-only: collapse isoforms sharing one BSJ to a single one first,
        # so isoform diversity can't multiply-count one real locus.
        # exon_based_counts must not do this, it needs every structure counted.
        orig_caller_df, orig_gt_df = caller_df, gt_df
        caller_df = dedup_by_bsj(caller_df)
        gt_df = dedup_by_bsj(gt_df)
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")

        tp_ids = set()
        fn_matched_gt = set()
        for gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []

            tp_out = self.bt.intersect("-a", str(caller_bed), "-b", str(gt_bed), *flags, "-f", "1.0", "-r", "-wa", "-u")
            tp_ids |= {line.split("\t")[3] for line in tp_out.strip().split("\n") if line}

            fn_out = self.bt.intersect("-a", str(gt_bed), "-b", str(caller_bed), *flags, "-f", "1.0", "-r", "-wa", "-u")
            fn_matched_gt |= {line.split("\t")[3] for line in fn_out.strip().split("\n") if line}

        tp = len(tp_ids)
        fp = len(caller_df) - tp
        fn = len(gt_df) - len(fn_matched_gt)
        if return_ids == "both":
            matched_caller_row_ids = _expand_bsj_match(tp_ids, "c", caller_df, orig_caller_df)
            matched_gt_row_ids = _expand_bsj_match(fn_matched_gt, "g", gt_df, orig_gt_df)
            return tp, fp, fn, matched_gt_row_ids, matched_caller_row_ids
        if return_ids:
            matched_gt_row_ids = _expand_bsj_match(fn_matched_gt, "g", gt_df, orig_gt_df)
            return tp, fp, fn, matched_gt_row_ids
        return tp, fp, fn

    def relaxed_bsj_counts(self, caller_df, gt_df, work: Path, tolerance: int, return_ids=False):
        # Same BSJ-only collapse as strict_bsj_counts, for the same reason.
        orig_caller_df, orig_gt_df = caller_df, gt_df
        caller_df = dedup_by_bsj(caller_df)
        gt_df = dedup_by_bsj(gt_df)
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")

        def boundary_points(bed_path, out_path, which):
            df = pd.read_csv(bed_path, sep="\t", header=None, names=BED12_COLS)
            if which == "start":
                df["p0"], df["p1"] = df["start"], df["start"] + 1
            else:
                df["p0"], df["p1"] = df["end"] - 1, df["end"]
            df[["chrom", "p0", "p1", "name", "score", "strand"]].to_csv(out_path, sep="\t", header=False, index=False)

        caller_start = work / "caller_start.bed"
        caller_end = work / "caller_end.bed"
        boundary_points(caller_bed, caller_start, "start")
        boundary_points(caller_bed, caller_end, "end")

        def pair_set(a_bed, b_slop_bed, strand_flag):
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(a_bed), "-b", str(b_slop_bed), *flags, "-wa", "-wb")
            pairs = set()
            for line in out.strip().split("\n"):
                if not line:
                    continue
                cols = line.split("\t")
                pairs.add((cols[3], cols[9]))
            return pairs

        genome = str(self.genome_sizes_file)
        matched_pairs = set()
        for gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            gt_start = work / f"gt_{tag}_start.bed"
            gt_end = work / f"gt_{tag}_end.bed"
            boundary_points(gt_bed, gt_start, "start")
            boundary_points(gt_bed, gt_end, "end")

            gt_start_slop = work / f"gt_{tag}_start_slop.bed"
            gt_end_slop = work / f"gt_{tag}_end_slop.bed"
            subprocess.run([self.bt.bedtools_bin, "slop", "-i", str(gt_start), "-g", genome, "-b", str(tolerance)],
                           stdout=open(gt_start_slop, "w"), check=True)
            subprocess.run([self.bt.bedtools_bin, "slop", "-i", str(gt_end), "-g", genome, "-b", str(tolerance)],
                           stdout=open(gt_end_slop, "w"), check=True)

            start_pairs = pair_set(caller_start, gt_start_slop, strand_flag)
            end_pairs = pair_set(caller_end, gt_end_slop, strand_flag)
            matched_pairs |= (start_pairs & end_pairs)

        matched_caller = {p[0] for p in matched_pairs}
        matched_gt = {p[1] for p in matched_pairs}

        tp = len(matched_caller)
        fp = len(caller_df) - tp
        fn = len(gt_df) - len(matched_gt)
        # Expand back to every original row sharing a matched BSJ. relaxed_bsj
        # is a BSJ-only measure, so every isoform of a matched BSJ counts too.
        if return_ids == "both":
            matched_caller_row_ids = _expand_bsj_match(matched_caller, "c", caller_df, orig_caller_df)
            matched_gt_row_ids = _expand_bsj_match(matched_gt, "g", gt_df, orig_gt_df)
            return tp, fp, fn, matched_gt_row_ids, matched_caller_row_ids
        if return_ids:
            matched_gt_row_ids = _expand_bsj_match(matched_gt, "g", gt_df, orig_gt_df)
            return tp, fp, fn, matched_gt_row_ids
        return tp, fp, fn

    def exon_based_counts(self, caller_df, gt_df, work: Path, min_overlap: float, return_ids=False):
        """return_ids: False returns (tp, fp, fn). True also returns matched
        GT row_ids. "both" also returns matched caller row_ids, to find
        which caller rows are FP (not in the matched set)."""
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")

        lines = []
        for gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(caller_bed), "-b", str(gt_bed), *flags, "-wo", "-split")
            lines += [l for l in out.strip().split("\n") if l]

        if not lines:
            if return_ids == "both":
                return 0, len(caller_df), len(gt_df), set(), set()
            if return_ids:
                return 0, len(caller_df), len(gt_df), set()
            return 0, len(caller_df), len(gt_df)

        caller_len = dict(zip("c" + caller_df["row_id"].astype(str), caller_df["spliced_len"]))
        gt_len = dict(zip("g" + gt_df["row_id"].astype(str), gt_df["spliced_len"]))

        overlap_bp = {}
        for line in lines:
            cols = line.split("\t")
            c_id, g_id, ov = cols[3], cols[15], int(cols[-1])
            key = (c_id, g_id)
            overlap_bp[key] = overlap_bp.get(key, 0) + ov

        matched_caller, matched_gt = set(), set()
        for (c_id, g_id), ov in overlap_bp.items():
            cl = caller_len.get(c_id, 0)
            gl = gt_len.get(g_id, 0)
            if cl == 0 or gl == 0:
                continue
            exon_score = min(ov / cl, ov / gl)
            if exon_score >= min_overlap:
                matched_caller.add(c_id)
                matched_gt.add(g_id)

        tp = len(matched_caller)
        fp = len(caller_df) - tp
        fn = len(gt_df) - len(matched_gt)
        if return_ids == "both":
            matched_gt_row_ids = {int(g[1:]) for g in matched_gt}
            matched_caller_row_ids = {int(c[1:]) for c in matched_caller}
            return tp, fp, fn, matched_gt_row_ids, matched_caller_row_ids
        if return_ids:
            matched_gt_row_ids = {int(g[1:]) for g in matched_gt}
            return tp, fp, fn, matched_gt_row_ids
        return tp, fp, fn

    def exon_based_pairs(self, caller_df, gt_df, work: Path, min_overlap: float):
        """Same exon-overlap computation as exon_based_counts, but returns
        the actual matched (caller_row_id, gt_row_id) pairs instead of two
        separate sets, so a read count can be looked up per pair. Used for
        exon-level expression/reproducibility matching, not for TP/FP/FN.

        Each caller row keeps only its single best (highest exon_score) GT
        match, the same nearest-wins tie-break match_to_gt uses. A GT row
        can still receive matches from several caller rows, real isoforms
        or tools legitimately supporting one true locus."""
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")

        lines = []
        for gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(caller_bed), "-b", str(gt_bed), *flags, "-wo", "-split")
            lines += [l for l in out.strip().split("\n") if l]

        if not lines:
            return []

        caller_len = dict(zip("c" + caller_df["row_id"].astype(str), caller_df["spliced_len"]))
        gt_len = dict(zip("g" + gt_df["row_id"].astype(str), gt_df["spliced_len"]))

        overlap_bp = {}
        for line in lines:
            cols = line.split("\t")
            c_id, g_id, ov = cols[3], cols[15], int(cols[-1])
            key = (c_id, g_id)
            overlap_bp[key] = overlap_bp.get(key, 0) + ov

        best_for_caller = {}
        for (c_id, g_id), ov in overlap_bp.items():
            cl = caller_len.get(c_id, 0)
            gl = gt_len.get(g_id, 0)
            if cl == 0 or gl == 0:
                continue
            exon_score = min(ov / cl, ov / gl)
            if exon_score < min_overlap:
                continue
            prev = best_for_caller.get(c_id)
            if prev is None or exon_score > prev[1]:
                best_for_caller[c_id] = (g_id, exon_score)

        return [(int(c_id[1:]), int(g_id[1:])) for c_id, (g_id, _) in best_for_caller.items()]

    @staticmethod
    def compute_metrics(tp, fp, fn):
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        return precision, recall, f1


def match_to_gt(call_df, gt_df, tolerance, extra_col=None):
    """Coordinate-tolerance match: same chrom, start/end both within
    `tolerance` bp (nearest by total distance if >1 candidate). Strand must
    match UNLESS the GT row is unstranded ('.', intergenic).
    Returns gt_bsj_id -> list of (n_reads, extra_col value or None), so the
    caller can both sum counts and, for tiered callers, filter by which
    tier's own bed12 that source row belongs to before aggregating."""
    pairs_by_gt = {}
    if len(call_df) == 0:
        return pairs_by_gt
    gt_by_chrom = {k: g for k, g in gt_df.groupby("chrom")}
    for row in call_df.itertuples(index=False):
        cand = gt_by_chrom.get(row.chrom)
        if cand is None:
            continue
        strand_ok = (cand["strand"] == row.strand) | ~cand["strand"].isin(_STRANDED)
        d_start = (cand["start"] - row.start).abs()
        d_end = (cand["end"] - row.end).abs()
        ok = strand_ok & (d_start <= tolerance) & (d_end <= tolerance)
        if not ok.any():
            continue
        total_d = (d_start + d_end)[ok]
        best_idx = total_d.idxmin()
        gt_bid = cand.loc[best_idx, "bsj_id"]
        extra = getattr(row, extra_col) if extra_col else None
        pairs_by_gt.setdefault(gt_bid, []).append((row.n_reads, extra))
    return pairs_by_gt


def match_pairs_by_strategy(strategy, call_df, gt_df, config, matcher=None, call_bed_df=None,
                             gt_bed_df=None, work=None, extra_col=None):
    """Same return shape as match_to_gt (gt_bsj_id -> list of (n_reads,
    extra_col value or None)), for whichever of the 3 structural strategies.

    strict_bsj and relaxed_bsj both reuse match_to_gt, just with a
    different tolerance (0 for strict_bsj, config.gt_match_tolerance for
    relaxed_bsj). exon_based needs the BED12 versions of call_df/gt_df
    (call_bed_df/gt_bed_df, with real block structure), a StructuralMatcher,
    and a work dir, since it goes through the same real exon-overlap check
    performance.py uses, not a boundary-distance check."""
    if strategy != "exon_based":
        tolerance = 0 if strategy == "strict_bsj" else config.gt_match_tolerance
        return match_to_gt(call_df, gt_df, tolerance, extra_col=extra_col)

    pairs = matcher.exon_based_pairs(call_bed_df, gt_bed_df, work, min_overlap=config.tp_min_overlap)
    call_bed_by_id = {r.row_id: r for r in call_bed_df.itertuples()}
    gt_bed_by_id = {r.row_id: r for r in gt_bed_df.itertuples()}
    count_by_coord = {(r.chrom, r.start, r.end, r.strand): r for r in call_df.itertuples()}

    pairs_by_gt = {}
    for c_row_id, g_row_id in pairs:
        c_bed = call_bed_by_id[c_row_id]
        count_row = count_by_coord.get((c_bed.chrom, c_bed.start, c_bed.end, c_bed.strand))
        if count_row is None:
            continue
        g_bed = gt_bed_by_id[g_row_id]
        gt_bid = f"{g_bed.chrom}:{g_bed.start}-{g_bed.end}:{g_bed.strand}"
        extra = getattr(count_row, extra_col) if extra_col else None
        pairs_by_gt.setdefault(gt_bid, []).append((count_row.n_reads, extra))
    return pairs_by_gt


def match_pred_to_pred(call_df, ref_df, tolerance):
    """Same tolerance match, but keyed to the reference call's own row index
    instead of a GT bsj_id. Used for caller-vs-caller comparisons with no
    GT involved. Same unstranded-row handling as match_to_gt. Returns a list
    of (ref_row_idx, n_reads) pairs."""
    if len(call_df) == 0 or len(ref_df) == 0:
        return []
    ref_by_chrom = {k: g for k, g in ref_df.groupby("chrom")}
    pairs = []
    for row in call_df.itertuples(index=False):
        cand = ref_by_chrom.get(row.chrom)
        if cand is None:
            continue
        strand_ok = (cand["strand"] == row.strand) | ~cand["strand"].isin(_STRANDED)
        d_start = (cand["start"] - row.start).abs()
        d_end = (cand["end"] - row.end).abs()
        ok = strand_ok & (d_start <= tolerance) & (d_end <= tolerance)
        if not ok.any():
            continue
        total_d = (d_start + d_end)[ok]
        best_idx = total_d.idxmin()
        pairs.append((cand.loc[best_idx].name, row.n_reads))
    return pairs


def collapse_near_duplicates(candidates):
    """candidates: [(n_reads, id), ...], all already matched to the same GT
    locus by match_to_gt. These are near-duplicate reports of one true
    circRNA, not independent counts, so summing them would double-count it.
    Returns the max read count instead."""
    if not candidates:
        return 0
    return max(n for n, _ in candidates)
