#!/usr/bin/env python3
"""
Step 1 of the discovery-tier false-positive characterization behind
docs/benchmark.md's "What are these false positives, really?" table.

For each run, matches the discovery tier's own merged BED12 against ground
truth (exon-based, same StructuralMatcher performance.py uses) and joins
the TP/FP label back onto that run's own clean_with_counts.tsv (read
counts, type, supporting_tools). Output feeds fp_characterize_step2_trace_reads.py.

Usage:
    python fp_characterize_step1_match_gt.py \\
        --dataset-root /path/to/final_benchmark \\
        --species human --runs human_run1 human_run2 human_run3 \\
        --bedtools-bin /path/to/bedtools \\
        --out-csv fp_characterize_raw.csv \\
        --work-dir ./fp_characterize_work
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset-root", required=True, type=Path,
                   help="Benchmark dataset root (fastq/, ground_truth/, results/, benchmark/)")
    p.add_argument("--species", required=True)
    p.add_argument("--runs", nargs="+", required=True, help="Run names, e.g. human_run1 human_run2 human_run3")
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

    args.work_dir.mkdir(parents=True, exist_ok=True)
    matcher = StructuralMatcher(args.bedtools_bin)

    rows = []
    for run in args.runs:
        bed_path = args.dataset_root / "results" / args.species / "circrna" / run / "merged" / f"{run}_discovery.bed12"
        tsv_path = args.dataset_root / "results" / args.species / "circrna" / run / "merged" / "clean" / f"{run}_discovery_clean_with_counts.tsv"
        gt_path = args.dataset_root / "benchmark" / "gt_prefiltered" / f"{run}_gt_prefiltered.bed12"

        caller_df = read_bed(bed_path)
        gt_df = read_bed(gt_path)
        tsv = pd.read_csv(tsv_path, sep="\t").rename(columns={"#chrom": "chrom"})

        wdir = args.work_dir / args.species / run
        wdir.mkdir(parents=True, exist_ok=True)
        tp, fp, fn, matched_gt_ids, matched_caller_ids = matcher.exon_based_counts(
            caller_df, gt_df, wdir, min_overlap=args.min_overlap, return_ids="both"
        )

        # Join bed12 rows to tsv rows on bsj_id/name, not position, in case
        # the two files aren't in lockstep row order.
        caller_df = caller_df.assign(is_tp=caller_df["row_id"].isin(matched_caller_ids))
        merged = tsv.merge(caller_df[["name", "is_tp"]], left_on="bsj_id", right_on="name", how="inner")
        assert len(merged) == len(tsv), f"{run}: tsv {len(tsv)} rows, matched {len(merged)} rows on bsj_id"

        for _, row in merged.iterrows():
            tools = str(row.get("supporting_tools", ""))
            rows.append({
                "species": args.species, "run": run, "is_tp": bool(row["is_tp"]),
                "reads": row["nanocirc_quant_reads"], "type": row.get("type"),
                "supporting_tools": row.get("supporting_tools"),
                "n_tools": len([t for t in tools.split(",") if t]),
                "sel_block_count": row.get("sel_block_count"),
                "chrom": row.get("chrom"), "start": row.get("start"), "end": row.get("end"),
                "strand": row.get("strand"),
            })

    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)
    print(f"{len(df)} rows ({(~df['is_tp']).sum()} FP, {df['is_tp'].sum()} TP) -> {args.out_csv}")


if __name__ == "__main__":
    main()
