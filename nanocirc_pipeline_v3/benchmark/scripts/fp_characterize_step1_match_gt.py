#!/usr/bin/env python3
"""Step 1 of FP characterization: matches one caller's bed12 against ground
truth (exon-based) per run. --caller can be a nanocirc tier (TP/FP joined
onto that run's clean_with_counts.tsv) or a raw tool run standalone (every
FP is single-tool support by definition). Output feeds
fp_characterize_step2_trace_reads.py.

Usage:
    python fp_characterize_step1_match_gt.py \\
        --dataset-root /path/to/final_benchmark \\
        --species human --runs human_run1 human_run2 human_run3 \\
        --caller discovery \\
        --bedtools-bin /path/to/bedtools \\
        --out-csv fp_characterize_raw.csv \\
        --work-dir ./fp_characterize_work
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

TIERS = ["discovery", "balanced_recall", "balanced_precision", "high_confidence"]
TOOLS = ["cirilong", "circfl", "isocirc", "circnick"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset-root", required=True, type=Path,
                   help="Benchmark dataset root (fastq/, ground_truth/, results/, benchmark/)")
    p.add_argument("--species", required=True)
    p.add_argument("--runs", nargs="+", required=True, help="Run names, e.g. human_run1 human_run2 human_run3")
    p.add_argument("--caller", default="discovery", choices=TIERS + TOOLS,
                   help="Nanocirc tier or raw tool to characterize")
    p.add_argument("--bedtools-bin", default="bedtools")
    p.add_argument("--min-overlap", type=float, default=0.75, help="Exon-based reciprocal overlap threshold for TP")
    p.add_argument("--out-csv", required=True, type=Path)
    p.add_argument("--work-dir", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    # nanocirc_benchmark package lives one directory up from scripts/.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from nanocirc_benchmark.bed_io import read_bed
    from nanocirc_benchmark.matching import StructuralMatcher
    from nanocirc_benchmark.paths import Paths

    args.work_dir.mkdir(parents=True, exist_ok=True)
    matcher = StructuralMatcher(args.bedtools_bin)
    paths = Paths(args.dataset_root, args.species, args.runs)
    is_tier = args.caller in TIERS

    rows = []
    for run in args.runs:
        bed_path = paths.caller_bed_path(run, args.caller, TIERS)
        gt_path = args.dataset_root / "benchmark" / "gt_prefiltered" / f"{run}_gt_prefiltered.bed12"

        caller_df = read_bed(bed_path)
        gt_df = read_bed(gt_path)

        wdir = args.work_dir / args.species / run
        wdir.mkdir(parents=True, exist_ok=True)
        tp, fp, fn, matched_gt_ids, matched_caller_ids = matcher.exon_based_counts(
            caller_df, gt_df, wdir, min_overlap=args.min_overlap, return_ids="both"
        )
        caller_df = caller_df.assign(is_tp=caller_df["row_id"].isin(matched_caller_ids))

        if is_tier:
            # Join on (chrom, start, end, strand), not bsj_id string equality:
            # some tiers keep several isoform rows per BSJ in the tsv
            # ("|iso1" suffix) but one collapsed row in the bed12.
            tsv_path = args.dataset_root / "results" / args.species / "circrna" / run / "merged" / "clean" / f"{run}_{args.caller}_clean_with_counts.tsv"
            tsv = pd.read_csv(tsv_path, sep="\t").rename(columns={"#chrom": "chrom"})

            def union_tools(series):
                seen = []
                for tools in series:
                    for t in str(tools).split(","):
                        if t and t not in seen:
                            seen.append(t)
                return ",".join(seen)

            tsv_by_locus = tsv.groupby(["chrom", "start", "end", "strand"], as_index=False).agg(
                supporting_tools=("supporting_tools", union_tools),
                type=("type", "first"),
                nanocirc_quant_reads=("nanocirc_quant_reads", "max"),
                sel_block_count=("sel_block_count", "first"),
            )
            merged = caller_df[["chrom", "start", "end", "strand", "is_tp"]].merge(
                tsv_by_locus, on=["chrom", "start", "end", "strand"], how="left")
            n_unmatched = merged["supporting_tools"].isna().sum()
            if n_unmatched:
                print(f"[{args.caller}] {run}: {n_unmatched}/{len(merged)} bed12 loci have no clean_with_counts.tsv match, "
                      f"left as unknown support (bed12 vs tsv are not fully 1:1 for this tier)")
            assert len(merged) == len(caller_df), f"{run}: bed12 {len(caller_df)} rows, joined {len(merged)} rows"

            for _, row in merged.iterrows():
                tools = row.get("supporting_tools")
                tools = "" if pd.isna(tools) else str(tools)
                rows.append({
                    "species": args.species, "run": run, "caller": args.caller, "is_tp": bool(row["is_tp"]),
                    "reads": row["nanocirc_quant_reads"], "type": row.get("type"),
                    "supporting_tools": tools,
                    "n_tools": len([t for t in tools.split(",") if t]),
                    "sel_block_count": row.get("sel_block_count"),
                    "chrom": row.get("chrom"), "start": row.get("start"), "end": row.get("end"),
                    "strand": row.get("strand"),
                })
        else:
            # A raw tool's own output has no clean_with_counts.tsv (that's
            # a nanocirc merge artefact) and no cross-tool voting: every
            # locus it calls is, by construction, single-tool support.
            for row in caller_df.itertuples():
                rows.append({
                    "species": args.species, "run": run, "caller": args.caller, "is_tp": bool(row.is_tp),
                    "reads": None, "type": None,
                    "supporting_tools": args.caller, "n_tools": 1,
                    "sel_block_count": row.blockCount,
                    "chrom": row.chrom, "start": row.start, "end": row.end, "strand": row.strand,
                })

    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)
    print(f"[{args.caller}] {len(df)} rows ({(~df['is_tp']).sum()} FP, {df['is_tp'].sum()} TP) -> {args.out_csv}")


if __name__ == "__main__":
    main()
