"""Structural circRNA matching (4 strategies), plus coordinate-tolerance
matching for read-count comparisons. Uses bedtools intersect approach.

Strategies:
  strict_bsj     exact whole-span match (`bedtools intersect -s -f 1.0 -r`).
  relaxed_bsj    both boundaries within `tolerance` bp of each other.
  exon_based     reciprocal exon overlap >= min_overlap (`-wo -split`).
  isoform_level  relaxed_bsj AND exon_based (>=99% overlap) both for
                 the same GT row - the "same molecule" test.
  full/strict/relaxed  config.isoform_match_levels' 3 named tiers: full =
                 exact BSJ + 100% overlap, strict = +/-Nbp BSJ + 99%
                 overlap, relaxed = +/-Nbp BSJ + 90% overlap.

Strand '.' means "unresolved", not "no strand" -- a caller that emits '.'
isn't claiming a strand at all (it's telling the researcher to resolve it
independently, e.g. in the lab), and a GT row with '.' means the ground
truth itself never established one (e.g. an intergenic circRNA). Neither
side lying about strand should ever cost a match: a '.' on EITHER side
credits a match against a real '+'/'-' call on the other side. Only two
opposing real strands ('+' vs '-') are a genuine mismatch. Every function
below runs three passes to cover this: real-vs-real (exact match via -s),
caller-'.'-vs-GT-real (any GT strand), and any-caller-vs-GT-'.' (any
caller strand).

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

STRATEGY_TITLES = {"strict_bsj": "Strict BSJ", "relaxed_bsj": "Relaxed BSJ", "exon_based": "Exon",
                    "isoform_level": "Isoform", "full": "Full", "strict": "Strict", "relaxed": "Relaxed"}

ISOFORM_MATCH_TAGS = ("full", "strict", "relaxed")


def isoform_match_level_params(config, tag: str):
    """(bsj_tolerance, min_overlap) for one of config.isoform_match_levels'
    named tiers (full/strict/relaxed)."""
    for level in config.isoform_match_levels:
        if level["tag"] == tag:
            return level["bsj_tolerance"], level["min_overlap"]
    raise ValueError(f"unknown isoform match tier: {tag}")


def _greedy_one_to_one(candidates):
    """candidates: [(c_id, g_id, score), ...], higher score = better match.
    Greedy assignment by descending score: once a caller row or GT row is
    claimed, it can't be claimed again. This is what keeps TP/FP/FN counts
    from inflating when several caller rows all qualify against one true GT
    entry (e.g. consensus_hybrid's minority-BSJ and multi-isoform-recovery
    rows), which otherwise credits one true locus as several true positives.
    COCO-style greedy matching, not exact (Hungarian) optimal assignment,
    same tradeoff standard detection benchmarks make for speed.
    Returns (matched_caller_ids, matched_gt_ids), equal in size by
    construction: TP == len(matched_caller) == len(matched_gt)."""
    matched_caller, matched_gt = set(), set()
    for c_id, g_id, _score in sorted(candidates, key=lambda x: x[2], reverse=True):
        if c_id in matched_caller or g_id in matched_gt:
            continue
        matched_caller.add(c_id)
        matched_gt.add(g_id)
    return matched_caller, matched_gt


def strategy_footnote(strategy: str, config) -> str:
    """Short, plot-footnote-length description of what a strategy's match
    threshold actually is, for whichever strategy is plotted."""
    if strategy == "strict_bsj":
        return "exact coordinate match"
    if strategy == "exon_based":
        return f">={config.tp_min_overlap:.0%} exon overlap matched"
    if strategy == "isoform_level":
        return f"+/-{config.gt_match_tolerance}bp BSJ AND >={config.isoform_level_min_overlap:.0%} exon overlap matched"
    if strategy in ISOFORM_MATCH_TAGS:
        tol, min_ov = isoform_match_level_params(config, strategy)
        bsj_desc = "exact BSJ" if tol == 0 else f"+/-{tol}bp BSJ"
        return f"{bsj_desc} AND >={min_ov:.0%} exon overlap matched"
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
    """Yields (caller_selector, gt_subset, strand_flag_or_None, tag) for the
    three matching passes every counting function needs:
      1. 'all' caller rows x stranded GT, exact match required (-s). The -s
         flag alone already excludes caller '.' rows here, since '.' never
         equals '+' or '-'.
      2. 'unstranded' caller rows only x stranded GT, no strand flag: a
         caller '.' isn't wrong, it declined to claim, so it can match any
         real GT strand.
      3. 'all' caller rows x unstranded ('.') GT, no strand flag: GT never
         claimed a strand either, so any caller strand can match it.
    caller_selector picks which caller bed a counting function should use
    for that pass ('all' or 'unstranded'), since only the GT side varies
    by content here."""
    stranded, unstranded = _split_by_strand(gt_df)
    yield "all", stranded, "-s", "stranded"
    yield "unstranded", stranded, None, "caller_unstranded"
    yield "all", unstranded, None, "gt_unstranded"


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
        caller_unstranded_df = caller_df[~caller_df["strand"].isin(_STRANDED)]
        caller_unstranded_bed = work / "caller_unstranded.bed"
        if len(caller_unstranded_df) > 0:
            write_bed_with_id(caller_unstranded_df, caller_unstranded_bed, "c")

        tp_ids = set()
        fn_matched_gt = set()
        for caller_sel, gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            if caller_sel == "unstranded":
                if len(caller_unstranded_df) == 0:
                    continue
                a_bed = caller_unstranded_bed
            else:
                a_bed = caller_bed
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []

            tp_out = self.bt.intersect("-a", str(a_bed), "-b", str(gt_bed), *flags, "-f", "1.0", "-r", "-wa", "-u")
            tp_ids |= {line.split("\t")[3] for line in tp_out.strip().split("\n") if line}

            fn_out = self.bt.intersect("-a", str(gt_bed), "-b", str(a_bed), *flags, "-f", "1.0", "-r", "-wa", "-u")
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
        # Dedup alone doesn't guarantee one-to-one here: two distinct (still
        # post-dedup) caller BSJs can each independently fall within
        # tolerance of the same GT BSJ, so this still needs the same greedy
        # one-to-one assignment exon_based_counts uses, scored by nearest
        # boundary distance instead of exon overlap.
        orig_caller_df, orig_gt_df = caller_df, gt_df
        caller_df = dedup_by_bsj(caller_df)
        gt_df = dedup_by_bsj(gt_df)
        caller_bounds = {f"c{r.row_id}": (r.start, r.end) for r in caller_df.itertuples()}
        gt_bounds = {f"g{r.row_id}": (r.start, r.end) for r in gt_df.itertuples()}
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")
        caller_unstranded_df = caller_df[~caller_df["strand"].isin(_STRANDED)]
        caller_unstranded_bed = work / "caller_unstranded.bed"
        if len(caller_unstranded_df) > 0:
            write_bed_with_id(caller_unstranded_df, caller_unstranded_bed, "c")

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
        caller_unstranded_start = work / "caller_unstranded_start.bed"
        caller_unstranded_end = work / "caller_unstranded_end.bed"
        if len(caller_unstranded_df) > 0:
            boundary_points(caller_unstranded_bed, caller_unstranded_start, "start")
            boundary_points(caller_unstranded_bed, caller_unstranded_end, "end")

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
        for caller_sel, gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            if caller_sel == "unstranded":
                if len(caller_unstranded_df) == 0:
                    continue
                a_start, a_end = caller_unstranded_start, caller_unstranded_end
            else:
                a_start, a_end = caller_start, caller_end
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

            start_pairs = pair_set(a_start, gt_start_slop, strand_flag)
            end_pairs = pair_set(a_end, gt_end_slop, strand_flag)
            matched_pairs |= (start_pairs & end_pairs)

        candidates = []
        for c_id, g_id in matched_pairs:
            c_start, c_end = caller_bounds[c_id]
            g_start, g_end = gt_bounds[g_id]
            distance = abs(c_start - g_start) + abs(c_end - g_end)
            candidates.append((c_id, g_id, -distance))  # nearest wins: less distance = higher score

        matched_caller, matched_gt = _greedy_one_to_one(candidates)

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
        which caller rows are FP (not in the matched set).

        One caller row can claim at most one GT row and vice versa (greedy
        one-to-one, best exon_score first): several caller rows qualifying
        against the same true GT entry no longer inflate TP, they compete
        for it, and the losers are FPs. See _greedy_one_to_one."""
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")
        caller_unstranded_df = caller_df[~caller_df["strand"].isin(_STRANDED)]
        caller_unstranded_bed = work / "caller_unstranded.bed"
        if len(caller_unstranded_df) > 0:
            write_bed_with_id(caller_unstranded_df, caller_unstranded_bed, "c")

        lines = []
        for caller_sel, gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            if caller_sel == "unstranded":
                if len(caller_unstranded_df) == 0:
                    continue
                a_bed = caller_unstranded_bed
            else:
                a_bed = caller_bed
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(a_bed), "-b", str(gt_bed), *flags, "-wo", "-split")
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

        candidates = []
        for (c_id, g_id), ov in overlap_bp.items():
            cl = caller_len.get(c_id, 0)
            gl = gt_len.get(g_id, 0)
            if cl == 0 or gl == 0:
                continue
            exon_score = min(ov / cl, ov / gl)
            if exon_score >= min_overlap:
                candidates.append((c_id, g_id, exon_score))

        matched_caller, matched_gt = _greedy_one_to_one(candidates)

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
        caller_unstranded_df = caller_df[~caller_df["strand"].isin(_STRANDED)]
        caller_unstranded_bed = work / "caller_unstranded.bed"
        if len(caller_unstranded_df) > 0:
            write_bed_with_id(caller_unstranded_df, caller_unstranded_bed, "c")

        lines = []
        for caller_sel, gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            if caller_sel == "unstranded":
                if len(caller_unstranded_df) == 0:
                    continue
                a_bed = caller_unstranded_bed
            else:
                a_bed = caller_bed
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(a_bed), "-b", str(gt_bed), *flags, "-wo", "-split")
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

    def isoform_level_counts(self, caller_df, gt_df, work: Path, tolerance: int, min_overlap: float,
                              return_ids=False):
        """A caller row counts as matched only if BOTH its own BSJ
        boundaries are within `tolerance` bp of a GT locus's own boundaries
        AND its exon structure reciprocally overlaps that same GT isoform
        by >= min_overlap. Same overlap computation as exon_based_counts, with an
        added per-pair boundary check. Same greedy one-to-one assignment as
        exon_based_counts too, for the same reason (see _greedy_one_to_one)."""
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")
        caller_unstranded_df = caller_df[~caller_df["strand"].isin(_STRANDED)]
        caller_unstranded_bed = work / "caller_unstranded.bed"
        if len(caller_unstranded_df) > 0:
            write_bed_with_id(caller_unstranded_df, caller_unstranded_bed, "c")

        lines = []
        for caller_sel, gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            if caller_sel == "unstranded":
                if len(caller_unstranded_df) == 0:
                    continue
                a_bed = caller_unstranded_bed
            else:
                a_bed = caller_bed
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(a_bed), "-b", str(gt_bed), *flags, "-wo", "-split")
            lines += [l for l in out.strip().split("\n") if l]

        if not lines:
            if return_ids == "both":
                return 0, len(caller_df), len(gt_df), set(), set()
            if return_ids:
                return 0, len(caller_df), len(gt_df), set()
            return 0, len(caller_df), len(gt_df)

        caller_len = dict(zip("c" + caller_df["row_id"].astype(str), caller_df["spliced_len"]))
        gt_len = dict(zip("g" + gt_df["row_id"].astype(str), gt_df["spliced_len"]))
        caller_bounds = {f"c{r.row_id}": (r.start, r.end) for r in caller_df.itertuples()}
        gt_bounds = {f"g{r.row_id}": (r.start, r.end) for r in gt_df.itertuples()}

        overlap_bp = {}
        for line in lines:
            cols = line.split("\t")
            c_id, g_id, ov = cols[3], cols[15], int(cols[-1])
            key = (c_id, g_id)
            overlap_bp[key] = overlap_bp.get(key, 0) + ov

        candidates = []
        for (c_id, g_id), ov in overlap_bp.items():
            cl = caller_len.get(c_id, 0)
            gl = gt_len.get(g_id, 0)
            if cl == 0 or gl == 0:
                continue
            exon_score = min(ov / cl, ov / gl)
            if exon_score < min_overlap:
                continue
            c_start, c_end = caller_bounds[c_id]
            g_start, g_end = gt_bounds[g_id]
            if abs(c_start - g_start) > tolerance or abs(c_end - g_end) > tolerance:
                continue
            candidates.append((c_id, g_id, exon_score))

        matched_caller, matched_gt = _greedy_one_to_one(candidates)

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

    def isoform_level_pairs(self, caller_df, gt_df, work: Path, tolerance: int, min_overlap: float,
                             return_scores=False):
        """Same computation as isoform_level_counts, but returns the actual
        matched (caller_row_id, gt_row_id) pairs instead of two sets,
        mirroring exon_based_pairs's relationship to exon_based_counts.
        Used for isoform-level expression/reproducibility matching, not for
        TP/FP/FN.

        return_scores: also return each pair's exon_score, for callers that
        need to pick the single best-scoring pair per GT row (e.g. a
        best-match read count) instead of just knowing a pair matched."""
        caller_bed = work / "caller.bed"
        write_bed_with_id(caller_df, caller_bed, "c")
        caller_unstranded_df = caller_df[~caller_df["strand"].isin(_STRANDED)]
        caller_unstranded_bed = work / "caller_unstranded.bed"
        if len(caller_unstranded_df) > 0:
            write_bed_with_id(caller_unstranded_df, caller_unstranded_bed, "c")

        lines = []
        for caller_sel, gt_subset, strand_flag, tag in _strand_passes(gt_df):
            if len(gt_subset) == 0:
                continue
            if caller_sel == "unstranded":
                if len(caller_unstranded_df) == 0:
                    continue
                a_bed = caller_unstranded_bed
            else:
                a_bed = caller_bed
            gt_bed = work / f"gt_{tag}.bed"
            write_bed_with_id(gt_subset, gt_bed, "g")
            flags = [strand_flag] if strand_flag else []
            out = self.bt.intersect("-a", str(a_bed), "-b", str(gt_bed), *flags, "-wo", "-split")
            lines += [l for l in out.strip().split("\n") if l]

        if not lines:
            return []

        caller_len = dict(zip("c" + caller_df["row_id"].astype(str), caller_df["spliced_len"]))
        gt_len = dict(zip("g" + gt_df["row_id"].astype(str), gt_df["spliced_len"]))
        caller_bounds = {f"c{r.row_id}": (r.start, r.end) for r in caller_df.itertuples()}
        gt_bounds = {f"g{r.row_id}": (r.start, r.end) for r in gt_df.itertuples()}

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
            c_start, c_end = caller_bounds[c_id]
            g_start, g_end = gt_bounds[g_id]
            if abs(c_start - g_start) > tolerance or abs(c_end - g_end) > tolerance:
                continue
            prev = best_for_caller.get(c_id)
            if prev is None or exon_score > prev[1]:
                best_for_caller[c_id] = (g_id, exon_score)

        if return_scores:
            return [(int(c_id[1:]), int(g_id[1:]), score) for c_id, (g_id, score) in best_for_caller.items()]
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
    match UNLESS the GT row is unstranded ('.', intergenic) OR the caller
    row is ('.', declined to claim a strand rather than lying about one).
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
        strand_ok = (cand["strand"] == row.strand) | ~cand["strand"].isin(_STRANDED) \
            | (row.strand not in _STRANDED)
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


def _bed12_pairs(strategy, call_bed_df, gt_bed_df, config, matcher, work, return_scores=False):
    if strategy == "exon_based":
        return matcher.exon_based_pairs(call_bed_df, gt_bed_df, work, min_overlap=config.tp_min_overlap)
    if strategy == "isoform_level":
        return matcher.isoform_level_pairs(call_bed_df, gt_bed_df, work, tolerance=config.gt_match_tolerance,
                                            min_overlap=config.isoform_level_min_overlap,
                                            return_scores=return_scores)
    tol, min_ov = isoform_match_level_params(config, strategy)
    return matcher.isoform_level_pairs(call_bed_df, gt_bed_df, work, tolerance=tol, min_overlap=min_ov,
                                        return_scores=return_scores)


def match_pairs_by_strategy(strategy, call_df, gt_df, config, matcher=None, call_bed_df=None,
                             gt_bed_df=None, work=None, extra_col=None):
    """gt_bsj_id -> list of (n_reads, extra_col value or None). BED12
    strategies read n_reads/extra_col from call_bed_df's matched row, and
    pool matches of sibling isoforms sharing one BSJ into one key -- use
    match_pairs_by_strategy_isoform_resolved instead when isoforms sharing a
    BSJ must stay distinguishable (ISOFORM_MATCH_TAGS)."""
    bed12_strategies = ("exon_based", "isoform_level") + ISOFORM_MATCH_TAGS
    if strategy not in bed12_strategies:
        tolerance = 0 if strategy == "strict_bsj" else config.gt_match_tolerance
        return match_to_gt(call_df, gt_df, tolerance, extra_col=extra_col)

    pairs = _bed12_pairs(strategy, call_bed_df, gt_bed_df, config, matcher, work)
    call_bed_by_id = {r.row_id: r for r in call_bed_df.itertuples()}
    gt_bed_by_id = {r.row_id: r for r in gt_bed_df.itertuples()}

    pairs_by_gt = {}
    for c_row_id, g_row_id in pairs:
        c_bed = call_bed_by_id[c_row_id]
        g_bed = gt_bed_by_id[g_row_id]
        gt_bid = f"{g_bed.chrom}:{g_bed.start}-{g_bed.end}:{g_bed.strand}"
        extra = c_bed.name.split("|")[0] if extra_col == "base_bsj_id" else (
            getattr(c_bed, extra_col) if extra_col else None)
        pairs_by_gt.setdefault(gt_bid, []).append((c_bed.n_reads, extra))
    return pairs_by_gt


def match_pairs_by_strategy_isoform_resolved(strategy, call_bed_df, gt_bed_df, config, matcher, work,
                                              extra_col=None):
    """Like match_pairs_by_strategy, but keyed by gt_bed_df's own row_id
    (one real isoform) instead of pooled BSJ coordinate, so sibling
    isoforms sharing a BSJ don't collapse onto the same key. Only valid for
    ISOFORM_MATCH_TAGS (full/strict/relaxed), which is where distinguishing
    isoforms sharing a BSJ actually matters.
    Returns gt_row_id -> list of (n_reads, extra_col value or None)."""
    if strategy not in ISOFORM_MATCH_TAGS:
        raise ValueError(f"isoform-resolved matching only supports {ISOFORM_MATCH_TAGS}, got {strategy!r}")
    pairs = _bed12_pairs(strategy, call_bed_df, gt_bed_df, config, matcher, work)
    call_bed_by_id = {r.row_id: r for r in call_bed_df.itertuples()}

    pairs_by_gt_row = {}
    for c_row_id, g_row_id in pairs:
        c_bed = call_bed_by_id[c_row_id]
        extra = c_bed.name.split("|")[0] if extra_col == "base_bsj_id" else (
            getattr(c_bed, extra_col) if extra_col else None)
        pairs_by_gt_row.setdefault(g_row_id, []).append((c_bed.n_reads, extra))
    return pairs_by_gt_row


def match_pairs_by_strategy_isoform_resolved_scored(strategy, call_bed_df, gt_bed_df, config, matcher, work,
                                                      extra_col=None):
    """Like match_pairs_by_strategy_isoform_resolved, but keeps each pair's
    exon_score too, for callers that need to filter candidates (e.g. by
    NanoCirc tier membership) before picking the best-scoring one from
    whatever's left, rather than picking the single global best up front.
    Returns gt_row_id -> list of (n_reads, extra_col value or None, exon_score)."""
    if strategy not in ISOFORM_MATCH_TAGS:
        raise ValueError(f"isoform-resolved matching only supports {ISOFORM_MATCH_TAGS}, got {strategy!r}")
    pairs = _bed12_pairs(strategy, call_bed_df, gt_bed_df, config, matcher, work, return_scores=True)
    call_bed_by_id = {r.row_id: r for r in call_bed_df.itertuples()}

    pairs_by_gt_row = {}
    for c_row_id, g_row_id, score in pairs:
        c_bed = call_bed_by_id[c_row_id]
        extra = c_bed.name.split("|")[0] if extra_col == "base_bsj_id" else (
            getattr(c_bed, extra_col) if extra_col else None)
        pairs_by_gt_row.setdefault(g_row_id, []).append((c_bed.n_reads, extra, score))
    return pairs_by_gt_row


def best_match_by_strategy_isoform_resolved(strategy, call_bed_df, gt_bed_df, config, matcher, work,
                                             extra_col=None):
    """Like match_pairs_by_strategy_isoform_resolved, but keeps only the
    single best-scoring caller row per GT row (highest exon_score), instead
    of pooling every qualifying caller row's reads into one number. Several
    near-duplicate caller isoforms matching one true isoform are read as
    "which one call is this", not summed -- summing would collapse
    isoform-level scoring back into locus-level pooling, the exact thing
    isoform-resolved matching exists to avoid.
    Returns gt_row_id -> (n_reads, extra_col value or None)."""
    scored = match_pairs_by_strategy_isoform_resolved_scored(strategy, call_bed_df, gt_bed_df, config, matcher,
                                                               work, extra_col=extra_col)
    return {g: max(candidates, key=lambda c: c[2])[:2] for g, candidates in scored.items()}


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
        strand_ok = (cand["strand"] == row.strand) | ~cand["strand"].isin(_STRANDED) \
            | (row.strand not in _STRANDED)
        d_start = (cand["start"] - row.start).abs()
        d_end = (cand["end"] - row.end).abs()
        ok = strand_ok & (d_start <= tolerance) & (d_end <= tolerance)
        if not ok.any():
            continue
        total_d = (d_start + d_end)[ok]
        best_idx = total_d.idxmin()
        pairs.append((cand.loc[best_idx].name, row.n_reads))
    return pairs


