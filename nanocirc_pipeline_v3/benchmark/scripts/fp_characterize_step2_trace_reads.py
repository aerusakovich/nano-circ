#!/usr/bin/env python3
"""
Step 2 of the false-positive characterization behind docs/benchmark.md's
"What are these false positives, really?" tables. Consumes
fp_characterize_step1_match_gt.py's output CSV for the same --caller.

Traces every FP locus's supporting reads back to the nanocirc_sim
simulator's own ground truth (circRNA_metadata.tsv) via each tool's native
per-locus read-evidence output, and classifies every supporting read into
exactly one bucket:
  - not from a circle at all (linear background read, no "_circ" suffix)
  - real circle, correct boundary (<=10bp of true source) and correct strand
    (pure boundary-precision noise)
  - real circle, correct boundary (<=10bp) but WRONG strand (e.g. an
    antisense-type true circle called on the host gene's sense strand)
  - real circle, called at a different/distant locus (>50bp from true source,
    or 11-50bp, both folded into "different locus" for the published table)
Also reports what fraction of each locus set (not reads) overlaps a UCSC
RepeatMasker repeat.

Three modes, set by --mode:
  - tool:   --caller is a raw tool (cirilong/circfl/isocirc/circnick) run
            standalone. Traces ALL of that tool's own FPs (every FP locus
            is single-tool support by construction, see step1).
  - sliced: --caller is a nanocirc tier. Traces only the tier's FP loci
            supported by exactly one tool (n_tools==1), one row per
            supporting tool per run, same as this script's original
            discovery-only behaviour.
  - pooled: --caller is a nanocirc tier. Traces ALL of the tier's FP loci
            regardless of n_tools, one row per run (not sliced by tool).
            A multi-tool-supported locus has its reads traced through the
            UNION of every one of its supporting tools' native read
            evidence at that locus, deduped by read_id first (the same
            physical read can appear in more than one tool's own report of
            the same locus).

Requires each run's native per-tool output (isocirc.out, CIRI-long's
.reads, CircFL-seq's circFL_Normal_pass.txt, CircNick-LRS's own BAM via
`bedtools bamtobed`) and the simulator's circRNA_metadata.tsv, which encodes
each true circle's own chrom/start/end/strand plus circ_type. Read IDs
carry their true source ID directly in the header
(`<circ_type><NNNNNNN>_..._circ` for simulated circles,
`ENST...` for linear background, no further guessing needed).

Usage:
    python fp_characterize_step2_trace_reads.py \\
        --dataset-root /path/to/final_benchmark \\
        --species human --runs human_run1 human_run2 human_run3 \\
        --caller discovery --mode pooled \\
        --sim-metadata /path/to/nanocirc_sim/datasets/human_dataset3/circRNA_metadata.tsv \\
        --repeats-bed /path/to/hg38_repeats.bed \\
        --bedtools-bin /path/to/bedtools \\
        --fp-csv fp_characterize_raw.csv \\
        --out-prefix fp_characterize
"""
import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd

TIERS = ["discovery", "balanced_recall", "balanced_precision", "high_confidence"]
TOOLS = ["cirilong", "circfl", "isocirc", "circnick"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset-root", required=True, type=Path)
    p.add_argument("--species", required=True)
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--caller", default="discovery", choices=TIERS + TOOLS)
    p.add_argument("--mode", default="sliced", choices=["sliced", "pooled", "tool"])
    p.add_argument("--sim-metadata", required=True, type=Path,
                   help="nanocirc_sim's circRNA_metadata.tsv for this species/dataset")
    p.add_argument("--repeats-bed", required=True, type=Path,
                   help="UCSC RepeatMasker BED (chrom, start, end, repName, repClass, repFamily)")
    p.add_argument("--bedtools-bin", default="bedtools")
    p.add_argument("--fp-csv", required=True, type=Path, help="Output of step1 for the same --caller")
    p.add_argument("--out-prefix", required=True, type=Path)
    return p.parse_args()


def run_dir(dataset_root, species, run):
    return dataset_root / "results" / species / "circrna" / run


def load_cirilong_reads(dataset_root, species, run):
    df = pd.read_csv(run_dir(dataset_root, species, run) / f"ciri_long/{run}_cirilong/{run}.reads", sep="\t")
    parsed = df["circ_id"].str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    df = pd.concat([df[["read_id"]], parsed], axis=1)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def load_circfl_reads(dataset_root, species, run):
    df = pd.read_csv(run_dir(dataset_root, species, run) / f"circfl_seq/{run}_circfl/RG/circFL_Normal_pass.txt", sep="\t")
    parsed = df["circID"].str.extract(r"^(?P<chrom>[^|]+)\|(?P<start>\d+)\|(?P<end>\d+)$")
    df = pd.concat([df[["readID"]], parsed], axis=1).dropna(subset=["readID"])
    df["read_id"] = df["readID"].str.split(",")
    df = df.explode("read_id")
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df[["read_id", "chrom", "start", "end"]]


def load_isocirc_reads(dataset_root, species, run):
    # line 1 of isocirc.out is a "#isocirc\t<version>" stamp; header is line 2.
    df = pd.read_csv(run_dir(dataset_root, species, run) / "isocirc/isocirc_output/isocirc.out", sep="\t", skiprows=1)
    df = df.rename(columns={"startCoor0based": "start", "endCoor": "end", "readIDs": "read_id"})
    df = df.dropna(subset=["read_id"])
    df["read_id"] = df["read_id"].str.split(",")
    df = df.explode("read_id")
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df[["read_id", "chrom", "start", "end"]]


def load_circnick_reads(dataset_root, species, run, bedtools_bin):
    # circnick's own text output has no per-locus read-ID list, only an
    # aggregate BSJ_reads count. Pull actual read names from its BAM
    # instead. bedtools bamtobed avoids needing a separate samtools container.
    bam = run_dir(dataset_root, species, run) / f"circnick_lrs/{run}_circnick/{run}/{run}.circRNA.sort.bam"
    out = subprocess.run([bedtools_bin, "bamtobed", "-i", str(bam)], capture_output=True, text=True, check=True).stdout
    rows = []
    for line in out.strip().split("\n"):
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        rows.append({"chrom": cols[0], "start": int(cols[1]), "end": int(cols[2]), "read_id": cols[3].split("~")[0]})
    return pd.DataFrame(rows)


LOADERS = {
    "cirilong": load_cirilong_reads,
    "circfl": load_circfl_reads,
    "isocirc": load_isocirc_reads,
    "circnick": load_circnick_reads,
}


def get_strand_lookup(dataset_root, species, run, caller, tiers):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from nanocirc_benchmark.paths import Paths
    bed_path = Paths(dataset_root, species, [run]).caller_bed_path(run, caller, tiers)
    bed = pd.read_csv(bed_path, sep="\t", header=None,
                       names=["chrom", "start", "end", "name", "score", "strand",
                              "ts", "te", "rgb", "bc", "bsizes", "bstarts"])
    return bed.set_index(["chrom", "start", "end"])["strand"].to_dict()


def repeat_overlap_fraction(fp_loci, repeats_bed, bedtools_bin, tmp_dir, tag):
    if not fp_loci:
        return 0.0
    bed_path = tmp_dir / f"_tmp_fp_loci_{tag}.bed"
    with open(bed_path, "w") as fh:
        for i, (chrom, start, end) in enumerate(sorted(fp_loci)):
            fh.write(f"{chrom}\t{start}\t{end}\t{i}\n")
    overlap_path = tmp_dir / f"_tmp_fp_overlap_{tag}.bed"
    subprocess.run(f"{bedtools_bin} intersect -a {bed_path} -b {repeats_bed} -wo > {overlap_path}", shell=True, check=True)
    hit_ids = set()
    with open(overlap_path) as fh:
        for line in fh:
            hit_ids.add(line.split("\t")[3])
    return len(hit_ids) / len(fp_loci)


def classify_reads(read_rows, chrom, start, end, called_strand, meta_by_id, buckets):
    """read_rows: iterable of read_id strings already matched to this locus
    (deduped by caller). Updates buckets in place."""
    for rid in read_rows:
        if not rid.endswith("_circ"):
            buckets["linear"] += 1
            continue
        m = re.match(r"^([a-zA-Z]+\d+)_", rid)
        true = meta_by_id.get(m.group(1)) if m else None
        if true is None:
            continue
        dist = None if true["chrom"] != chrom else max(abs(true["start"] - start), abs(true["end"] - end))
        if dist is None or dist > 50:
            buckets["far"] += 1
        elif dist <= 10:
            if true["strand"] == called_strand:
                buckets["near_right_strand"] += 1
            else:
                buckets["near_wrong_strand"] += 1
        else:
            buckets["mid"] += 1


def reads_at_locus(reads_df, tool, chrom, start, end):
    """Every matching read_id row, NOT deduped: circnick's BAM can carry
    more than one alignment record for the same physical read (split/
    secondary alignments), and the original single-tool methodology counts
    each row. Callers that need cross-tool deduplication (pooled mode) do
    it themselves after collecting reads from every contributing tool."""
    overlap_mode = (tool == "circnick")
    csub = reads_df[reads_df["chrom"] == chrom]
    if overlap_mode:
        sub = csub[(csub["start"] < end) & (csub["end"] > start)]
    else:
        sub = csub[(csub["start"] - start).abs().le(1) & (csub["end"] - end).abs().le(1)]
    return list(sub["read_id"])


def main():
    args = parse_args()
    tmp_dir = args.out_prefix.parent
    tmp_dir.mkdir(parents=True, exist_ok=True)

    fp_all = pd.read_csv(args.fp_csv)
    fp_all = fp_all[(fp_all["species"] == args.species) & (~fp_all["is_tp"])]

    meta = pd.read_csv(args.sim_metadata, sep="\t")
    meta_by_id = meta.set_index("circRNA_id")[["chrom", "start", "end", "strand"]].to_dict("index")

    reads_cache = {}  # (run, tool) -> reads_df, loaded once per run

    def get_reads(run, tool):
        key = (run, tool)
        if key not in reads_cache:
            reads_cache[key] = load_circnick_reads(args.dataset_root, args.species, run, args.bedtools_bin) \
                if tool == "circnick" else LOADERS[tool](args.dataset_root, args.species, run)
        return reads_cache[key]

    per_run_rows = []

    if args.mode == "tool":
        for run in args.runs:
            fp_run = fp_all[fp_all["run"] == run]
            strand_lookup = get_strand_lookup(args.dataset_root, args.species, run, args.caller, TIERS)
            fp_loci = set(zip(fp_run["chrom"], fp_run["start"], fp_run["end"]))
            if not fp_loci:
                continue
            reads_df = get_reads(run, args.caller)
            buckets = {"linear": 0, "far": 0, "near_wrong_strand": 0, "near_right_strand": 0, "mid": 0}
            for chrom, start, end in fp_loci:
                called_strand = strand_lookup.get((chrom, start, end))
                rids = reads_at_locus(reads_df, args.caller, chrom, start, end)
                classify_reads(rids, chrom, start, end, called_strand, meta_by_id, buckets)
            total = sum(buckets.values())
            repeat_frac = repeat_overlap_fraction(fp_loci, args.repeats_bed, args.bedtools_bin, tmp_dir, f"{run}_{args.caller}")
            per_run_rows.append({**buckets, "total_reads": total, "n_loci": len(fp_loci),
                                  "repeat_frac_loci": repeat_frac, "run": run, "caller": args.caller})

    elif args.mode == "sliced":
        fp_all = fp_all[fp_all["n_tools"] == 1]
        for run in args.runs:
            fp_run = fp_all[fp_all["run"] == run]
            strand_lookup = get_strand_lookup(args.dataset_root, args.species, run, args.caller, TIERS)
            for tool in TOOLS:
                tool_fp = fp_run[fp_run["supporting_tools"] == tool]
                fp_loci = set(zip(tool_fp["chrom"], tool_fp["start"], tool_fp["end"]))
                if not fp_loci:
                    continue
                reads_df = get_reads(run, tool)
                buckets = {"linear": 0, "far": 0, "near_wrong_strand": 0, "near_right_strand": 0, "mid": 0}
                for chrom, start, end in fp_loci:
                    called_strand = strand_lookup.get((chrom, start, end))
                    rids = reads_at_locus(reads_df, tool, chrom, start, end)
                    classify_reads(rids, chrom, start, end, called_strand, meta_by_id, buckets)
                total = sum(buckets.values())
                repeat_frac = repeat_overlap_fraction(fp_loci, args.repeats_bed, args.bedtools_bin, tmp_dir, f"{run}_{args.caller}_{tool}")
                per_run_rows.append({**buckets, "total_reads": total, "n_loci": len(fp_loci),
                                      "repeat_frac_loci": repeat_frac, "run": run, "caller": f"{args.caller}:{tool}"})

    else:  # pooled
        for run in args.runs:
            fp_run = fp_all[fp_all["run"] == run]
            strand_lookup = get_strand_lookup(args.dataset_root, args.species, run, args.caller, TIERS)
            fp_loci_tools = list(zip(fp_run["chrom"], fp_run["start"], fp_run["end"], fp_run["supporting_tools"]))
            if not fp_loci_tools:
                continue
            fp_loci = set((c, s, e) for c, s, e, _ in fp_loci_tools)
            buckets = {"linear": 0, "far": 0, "near_wrong_strand": 0, "near_right_strand": 0, "mid": 0}
            for chrom, start, end, tools_str in fp_loci_tools:
                called_strand = strand_lookup.get((chrom, start, end))
                supporting = [t for t in str(tools_str).split(",") if t and t in LOADERS]
                rids = set()
                for tool in supporting:
                    reads_df = get_reads(run, tool)
                    rids.update(reads_at_locus(reads_df, tool, chrom, start, end))
                classify_reads(rids, chrom, start, end, called_strand, meta_by_id, buckets)
            total = sum(buckets.values())
            repeat_frac = repeat_overlap_fraction(fp_loci, args.repeats_bed, args.bedtools_bin, tmp_dir, f"{run}_{args.caller}_pooled")
            per_run_rows.append({**buckets, "total_reads": total, "n_loci": len(fp_loci),
                                  "repeat_frac_loci": repeat_frac, "run": run, "caller": args.caller})

    per_run_df = pd.DataFrame(per_run_rows)
    if per_run_df.empty:
        print(f"[{args.caller}/{args.mode}] no FP loci found for {args.species}, nothing to trace")
        per_run_csv = args.out_prefix.with_name(args.out_prefix.name + "_per_run.csv")
        per_run_df.to_csv(per_run_csv, index=False)
        return

    for col in ["linear", "far", "near_wrong_strand", "near_right_strand", "mid"]:
        per_run_df[col + "_pct"] = per_run_df[col] / per_run_df["total_reads"] * 100
    per_run_df["repeat_pct"] = per_run_df["repeat_frac_loci"] * 100
    per_run_df["different_locus_pct"] = per_run_df["far_pct"] + per_run_df["mid_pct"]
    per_run_csv = args.out_prefix.with_name(args.out_prefix.name + "_per_run.csv")
    per_run_df.to_csv(per_run_csv, index=False)

    pct_cols = ["linear_pct", "near_wrong_strand_pct", "near_right_strand_pct", "different_locus_pct", "repeat_pct"]
    agg = per_run_df.groupby("caller")[pct_cols].agg(["mean", "std"])
    n_reads = per_run_df.groupby("caller")["total_reads"].sum()
    n_loci = per_run_df.groupby("caller")["n_loci"].sum()

    print(f"\n-> {per_run_csv}")
    print(f"\n=== mode={args.mode}, mean +/- SD across runs, n pooled across runs ===")
    for caller in agg.index:
        print(f"\n{caller}: n_loci_total={n_loci[caller]}, n_reads_total={n_reads[caller]}")
        for col in pct_cols:
            mean, std = agg.loc[caller, (col, "mean")], agg.loc[caller, (col, "std")]
            print(f"  {col}: {mean:.1f} +/- {std:.1f}")

    agg_csv = args.out_prefix.with_name(args.out_prefix.name + "_agg.csv")
    agg.to_csv(agg_csv)
    print(f"\n-> {agg_csv}")


if __name__ == "__main__":
    main()
