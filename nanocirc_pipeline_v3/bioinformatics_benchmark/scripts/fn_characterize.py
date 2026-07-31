#!/usr/bin/env python3
"""
What does discovery (the union of all 4 tools, max recall) actually miss?
For every ground-truth circle NOT recovered by ANY tool (exon-based match,
same StructuralMatcher performance.py uses), characterizes it against the
simulator's own circRNA_metadata.tsv and compares to recovered (TP) GT
circles on the same axes:
  - mature (spliced) length
  - exon/block count
  - copy_number (rolling-circle/RCRT repeat count, fewer copies means
    weaker concatemer evidence for the RCA-based tools)
  - expression level (TPM, n_reads_target)
  - circRNA biotype (eciRNA/EIciRNA/ciRNA/antisense/intergenic)
  - splice-site canonicality (GT/AG) at the BSJ boundary, via direct FASTA
    lookup (no samtools needed, plain .fai random access)
  - genomic repeat overlap (UCSC RepeatMasker)
  - gene overlap / how many transcripts that gene has (multi-isoform proxy)

GT bed12's own "name" field is the candidate's circRNA_id directly (see
nanocirc_sim's bed12.py), so the join to circRNA_metadata.tsv needs no
coordinate matching at all.

Usage:
    python fn_characterize.py \\
        --dataset-root /path/to/final_benchmark \\
        --species human --runs human_run1 human_run2 human_run3 \\
        --sim-metadata /path/to/nanocirc_sim/datasets/human_dataset3/circRNA_metadata.tsv \\
        --genome-fasta /path/to/hg38.fa --gtf /path/to/gencode.gtf \\
        --repeats-bed /path/to/hg38_repeats.bed \\
        --bedtools-bin /path/to/bedtools \\
        --min-overlap 0.75 \\
        --out-prefix ./fn_characterize
"""
import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


class Faidx:
    """Minimal random-access FASTA reader using a .fai index, no samtools."""

    def __init__(self, fasta_path):
        self.index = {}
        with open(f"{fasta_path}.fai") as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self.fh = open(fasta_path, "rb")

    def fetch(self, chrom, start, end):
        """0-based half-open [start, end)."""
        if chrom not in self.index or start < 0:
            return ""
        length, offset, linebases, linewidth = self.index[chrom]
        end = min(end, length)
        if start >= end:
            return ""
        start_line = start // linebases
        start_byte = offset + start_line * linewidth + (start % linebases)
        n_bases = end - start
        self.fh.seek(start_byte)
        raw = b""
        while len(raw) < n_bases + (n_bases // linebases) + 2:
            chunk = self.fh.read(4096)
            if not chunk:
                break
            raw += chunk
        return raw.replace(b"\n", b"").replace(b"\r", b"")[:n_bases].decode("ascii").upper()


def splice_canonical(faidx, chrom, start, end, strand):
    if strand == "+":
        donor = faidx.fetch(chrom, end, end + 2)
        acceptor = faidx.fetch(chrom, start - 2, start)
    else:
        donor = revcomp(faidx.fetch(chrom, start - 2, start))
        acceptor = revcomp(faidx.fetch(chrom, end, end + 2))
    return donor == "GT" and acceptor == "AG"


def build_gene_features(gtf_path, work_dir):
    gene_bed = work_dir / "genes.bed"
    if gene_bed.exists():
        tx_count = {}
        with open(work_dir / "tx_count.tsv") as fh:
            for line in fh:
                g, n = line.rstrip("\n").split("\t")
                tx_count[g] = int(n)
        return gene_bed, tx_count

    tx_count = defaultdict(set)
    genes = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature, attrs = cols[2], cols[8]
            gene_id = None
            for part in attrs.split(";"):
                part = part.strip()
                if part.startswith("gene_id"):
                    gene_id = part.split(" ")[1].strip('"')
                    break
            if gene_id is None:
                continue
            if feature == "gene":
                genes.append((cols[0], int(cols[3]) - 1, int(cols[4]), cols[6], gene_id))
            elif feature == "transcript":
                tx_id = None
                for part in attrs.split(";"):
                    part = part.strip()
                    if part.startswith("transcript_id"):
                        tx_id = part.split(" ")[1].strip('"')
                        break
                if tx_id:
                    tx_count[gene_id].add(tx_id)

    with open(gene_bed, "w") as fh:
        for chrom, start, end, strand, gene_id in genes:
            fh.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")
    subprocess.run(f"sort -k1,1 -k2,2n {gene_bed} -o {gene_bed}", shell=True, check=True)
    tx_count_flat = {g: len(txs) for g, txs in tx_count.items()}
    with open(work_dir / "tx_count.tsv", "w") as fh:
        for g, n in tx_count_flat.items():
            fh.write(f"{g}\t{n}\n")
    return gene_bed, tx_count_flat


def repeat_overlap_ids(loci, repeats_bed, bedtools_bin, work_dir, tag):
    if not loci:
        return set()
    bed_path = work_dir / f"_fn_loci_{tag}.bed"
    with open(bed_path, "w") as fh:
        for i, (chrom, start, end) in enumerate(sorted(loci)):
            fh.write(f"{chrom}\t{start}\t{end}\t{i}\n")
    overlap_path = work_dir / f"_fn_overlap_{tag}.bed"
    subprocess.run(f"{bedtools_bin} intersect -a {bed_path} -b {repeats_bed} -wo > {overlap_path}", shell=True, check=True)
    hit = set()
    with open(overlap_path) as fh:
        for line in fh:
            hit.add(int(line.split("\t")[3]))
    return hit


def gene_overlap_lookup(loci, gene_bed, bedtools_bin, work_dir, tag):
    if not loci:
        return {}
    bed_path = work_dir / f"_fn_gloci_{tag}.bed"
    with open(bed_path, "w") as fh:
        for i, (chrom, start, end) in enumerate(sorted(loci)):
            fh.write(f"{chrom}\t{start}\t{end}\t{i}\n")
    overlap_path = work_dir / f"_fn_govl_{tag}.bed"
    subprocess.run(f"{bedtools_bin} intersect -a {bed_path} -b {gene_bed} -wo > {overlap_path}", shell=True, check=True)
    hits = defaultdict(set)
    with open(overlap_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            hits[int(cols[3])].add(cols[7])
    return hits


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset-root", required=True, type=Path)
    p.add_argument("--species", required=True)
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--sim-metadata", required=True, type=Path)
    p.add_argument("--genome-fasta", required=True, type=Path, help="Plain FASTA with a matching .fai index")
    p.add_argument("--gtf", required=True, type=Path)
    p.add_argument("--repeats-bed", required=True, type=Path)
    p.add_argument("--bedtools-bin", default="bedtools")
    p.add_argument("--min-overlap", type=float, default=0.75)
    p.add_argument("--out-prefix", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from nanocirc_benchmark.bed_io import read_bed
    from nanocirc_benchmark.matching import StructuralMatcher

    work_dir = args.out_prefix.parent
    work_dir.mkdir(parents=True, exist_ok=True)

    matcher = StructuralMatcher(args.bedtools_bin)
    faidx = Faidx(args.genome_fasta)
    gene_bed, tx_counts = build_gene_features(args.gtf, work_dir)
    meta = pd.read_csv(args.sim_metadata, sep="\t")
    meta_by_id = meta.set_index("circRNA_id")

    per_run_rows = []
    for run in args.runs:
        bed_path = args.dataset_root / "results" / args.species / "circrna" / run / "merged" / f"{run}_discovery.bed12"
        gt_path = args.dataset_root / "benchmark" / "gt_prefiltered" / f"{run}_gt_prefiltered.bed12"
        caller_df = read_bed(bed_path)
        gt_df = read_bed(gt_path)

        wdir = work_dir / args.species / run
        wdir.mkdir(parents=True, exist_ok=True)
        tp, fp, fn, matched_gt_ids = matcher.exon_based_counts(caller_df, gt_df, wdir, min_overlap=args.min_overlap, return_ids=True)

        gt_df = gt_df.assign(is_recovered=gt_df["row_id"].isin(matched_gt_ids))
        loci = list(zip(gt_df["chrom"], gt_df["start"], gt_df["end"]))
        repeat_hits = repeat_overlap_ids(loci, args.repeats_bed, args.bedtools_bin, wdir, run)
        gene_hits = gene_overlap_lookup(loci, gene_bed, args.bedtools_bin, wdir, run)

        for i, row in gt_df.reset_index(drop=True).iterrows():
            circ_id = row["name"]
            m = meta_by_id.loc[circ_id] if circ_id in meta_by_id.index else None
            if m is None:
                continue
            canonical = splice_canonical(faidx, row["chrom"], int(row["start"]), int(row["end"]), row["strand"])
            genes = gene_hits.get(i, set())
            max_tx = max((tx_counts.get(g, 0) for g in genes), default=0)
            per_run_rows.append({
                "run": run, "circ_id": circ_id, "is_recovered": bool(row["is_recovered"]),
                "mature_length": m["mature_length"], "exon_count": m["exon_count"],
                "copy_number": m["copy_number"], "TPM": m["TPM"], "n_reads_target": m["n_reads_target"],
                "n_reads_confirmed": m["n_reads_confirmed"], "type": m["type"],
                "canonical_splice": canonical, "overlaps_repeat": i in repeat_hits,
                "n_genes_overlapped": len(genes), "max_transcripts_in_gene": max_tx,
            })

    df = pd.DataFrame(per_run_rows)
    df.to_csv(args.out_prefix.with_name(args.out_prefix.name + "_per_locus.csv"), index=False)

    def summarize(sub):
        return pd.Series({
            "n": len(sub),
            "mean_length": sub["mature_length"].mean(),
            "median_length": sub["mature_length"].median(),
            "mean_exon_count": sub["exon_count"].mean(),
            "mean_copy_number": sub["copy_number"].mean(),
            "median_copy_number": sub["copy_number"].median(),
            "frac_copy_number_le3": (sub["copy_number"] <= 3).mean(),
            "mean_TPM": sub["TPM"].mean(),
            "median_n_reads_target": sub["n_reads_target"].median(),
            "frac_canonical_splice": sub["canonical_splice"].mean(),
            "frac_overlaps_repeat": sub["overlaps_repeat"].mean(),
            "frac_intergenic_0genes": (sub["n_genes_overlapped"] == 0).mean(),
            "frac_multi_isoform_gene_ge5tx": (sub["max_transcripts_in_gene"] >= 5).mean(),
        })

    rows = []
    for run in args.runs:
        rdf = df[df["run"] == run]
        for label, sub in [("FN (missed)", rdf[~rdf["is_recovered"]]), ("recovered", rdf[rdf["is_recovered"]])]:
            s = summarize(sub)
            s["run"] = run
            s["label"] = label
            rows.append(s)
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(args.out_prefix.with_name(args.out_prefix.name + "_per_run_summary.csv"), index=False)

    print("\n=== per-run FN vs recovered ===")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(summary_df.round(2).to_string(index=False))

    print("\n=== mean +/- SD across runs ===")
    num_cols = [c for c in summary_df.columns if c not in ("run", "label")]
    agg = summary_df.groupby("label")[num_cols].agg(["mean", "std"])
    for label in agg.index:
        print(f"\n{label}:")
        for col in num_cols:
            mean, std = agg.loc[label, (col, "mean")], agg.loc[label, (col, "std")]
            print(f"  {col}: {mean:.2f} +/- {std:.2f}")

    print("\n=== type distribution, FN vs recovered (pooled across runs) ===")
    print(pd.crosstab(df["type"], df["is_recovered"], normalize="columns").round(3).rename(columns={True: "recovered", False: "FN"}))

    agg.to_csv(args.out_prefix.with_name(args.out_prefix.name + "_agg.csv"))
    print(f"\n-> {args.out_prefix.with_name(args.out_prefix.name + '_agg.csv')}")


if __name__ == "__main__":
    main()
