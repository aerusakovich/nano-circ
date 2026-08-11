"""Shared BED12 I/O and bedtools-wrapper helpers."""
import subprocess

import pandas as pd

BED12_COLS = ["chrom", "start", "end", "name", "score", "strand",
              "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]


def read_bed(path, dedup=True):
    """Read a BED12(+) file, keep only real data rows, add a stable 0..n-1
    row_id and each row's own spliced (mature) length from its BED12 blocks.

    dedup=True collapses only byte-identical duplicate lines (the full row,
    including block structure), not rows sharing just (chrom, start, end,
    strand). Callers and GT both legitimately carry several distinct
    isoforms sharing one BSJ, and exon_based_counts() needs each of those
    structures counted separately. strict_bsj_counts/relaxed_bsj_counts
    collapse to one row per BSJ on their own, see dedup_by_bsj()."""
    rows = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            key = tuple(cols[:12])
            if dedup:
                if key in seen:
                    continue
                seen.add(key)
            rows.append(cols[:12])
    df = pd.DataFrame(rows, columns=BED12_COLS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["row_id"] = range(len(df))
    block_sizes = df["blockSizes"].str.rstrip(",").apply(lambda s: sum(int(x) for x in s.split(",") if x))
    df["spliced_len"] = block_sizes
    return df


def empty_bed_df():
    """Empty DataFrame with read_bed's column set, for a tool that produced
    no output file at all (e.g. disabled via --run_<tool> false) rather than
    an empty one."""
    return pd.DataFrame(columns=BED12_COLS + ["row_id", "spliced_len"])


def dedup_by_bsj(df):
    """Collapse to one row per (chrom, start, end, strand) BSJ, keeping the
    first isoform seen at each. strict_bsj_counts and relaxed_bsj_counts are
    BSJ-only measures, so several isoforms sharing one BSJ count as a single
    detection there. Counting them separately would multiply-count one real
    locus. exon_based_counts() must not use this, it needs every distinct
    structure counted on its own."""
    return df.drop_duplicates(subset=["chrom", "start", "end", "strand"], keep="first").reset_index(drop=True)


def write_bed_with_id(df, path, id_prefix):
    """Write BED12 with row_id encoded into the name field (col4) so
    bedtools -wa/-wb round-trips it."""
    out = df.copy()
    out["name"] = id_prefix + out["row_id"].astype(str)
    out[BED12_COLS].to_csv(path, sep="\t", header=False, index=False)


class Bedtools:
    """Thin subprocess wrapper around `bedtools intersect`/`slop`, bound to
    one bedtools binary path (config-driven, not hardcoded)."""

    def __init__(self, bedtools_bin="bedtools"):
        self.bedtools_bin = bedtools_bin

    def intersect(self, *args):
        result = subprocess.run([self.bedtools_bin, "intersect"] + list(args),
                                 capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bedtools failed: {' '.join(args)}\n{result.stderr}")
        return result.stdout

    def slop(self, *args):
        result = subprocess.run([self.bedtools_bin, "slop"] + list(args),
                                 capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bedtools slop failed: {' '.join(args)}\n{result.stderr}")
        return result.stdout


def bsj_id(chrom, start, end, strand) -> str:
    """Canonical bsj_id string: chrom:start-end:strand."""
    return f"{chrom}:{start}-{end}:{strand}"
