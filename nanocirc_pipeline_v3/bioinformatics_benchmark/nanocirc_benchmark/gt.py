"""Ground-truth loading helpers."""
import gzip

import pandas as pd


def load_gt_loci(paths, run: str) -> pd.DataFrame:
    df = pd.read_csv(paths.ground_truth_tsv(run), sep="\t")
    df["bsj_id"] = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str) + ":" + df["strand"]
    return df.groupby(["bsj_id", "chrom", "start", "end", "strand"], as_index=False)["n_reads_confirmed"].sum()


def count_reads(fastq_gz_path) -> int:
    with gzip.open(fastq_gz_path, "rt") as f:
        return sum(1 for _ in f) // 4
