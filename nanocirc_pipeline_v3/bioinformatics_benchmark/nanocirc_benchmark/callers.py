"""Per-run, per-caller locus-level read counts, and NanoCirc tier membership.

Coordinate conventions (0-based BED-style start/end):
  isocirc.out              : startCoor0based/endCoor already match (no offset)
  circFL_Normal_pass.txt   : native start is 1-based -> unified_start = start - 1
  ciri-long .reads (circ_id "chr:start-end") : native start is 1-based -> unified_start = start - 1
  circnick *.annotated.txt : own chrom/start/end are pre-liftover whenever
                             the dataset needed --circnick_liftover_chain
                             (e.g. a human hg19->hg38 run). This can be tens
                             to hundreds of kb off from every other file's
                             coordinates for the same locus. Datasets needing
                             no liftover (e.g. mouse) are unaffected. Its own
                             BED12 output (built after CIRCNICK_LIFTOVER) has
                             the correct, already-lifted coordinates, and
                             shares the same internal_circRNA_name identifier
                             with the native file. So read counts are always
                             joined onto the BED12's coordinates via that
                             name, never onto the native file's own
                             chrom/start/end.
"""
import pandas as pd

from .bed_io import read_bed


def parse_isocirc(f, bed12_path=None) -> dict:
    """readIDs per (chrom, start, end, strand) locus. Groups by coordinates
    only, not by isocirc.out's own geneStrand column: that column is the
    overlapping gene's strand, not the circRNA's own strand, so it is often
    "NA" (no gene overlap) or comma-repeated ("+,+", several same-strand
    genes). The clean strand comes from isocirc.bed instead, the same file
    every other matching stage already treats as authoritative."""
    df = pd.read_csv(f, sep="\t", skiprows=1, dtype=str)
    df.columns = [c.lstrip("#") for c in df.columns]
    df["start"] = df["startCoor0based"].astype(int)
    df["end"] = df["endCoor"].astype(int)

    strand_by_coord = {}
    if bed12_path is not None:
        with open(bed12_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 6:
                    strand_by_coord[(cols[0], int(cols[1]), int(cols[2]))] = cols[5]

    out = {}
    for (chrom, start, end), g in df.groupby(["chrom", "start", "end"]):
        key = (chrom, int(start), int(end))
        strand = strand_by_coord.get(key)
        if strand is None:
            # No BED12 given, or this locus is missing from it: fall back to
            # geneStrand, cleaning up "NA" and comma-repeated values.
            vals = {v.split(",")[0] for v in g["geneStrand"] if v and v != "NA"}
            strand = next(iter(vals), ".")
        reads = set()
        for ids in g["readIDs"].dropna():
            reads.update(x for x in ids.split(",") if x)
        out[(chrom, int(start), int(end), strand)] = reads
    return out


def parse_circfl(f) -> dict:
    df = pd.read_csv(f, sep="\t", dtype=str)
    df["start0"] = df["start"].astype(int) - 1
    df["end"] = df["end"].astype(int)
    out = {}
    for (chrom, start, end, strand), g in df.groupby(["chr", "start0", "end", "strand"]):
        reads = set()
        for ids in g["readID"].dropna():
            reads.update(x for x in ids.split(",") if x)
        out[(chrom, int(start), int(end), strand)] = reads
    return out


def parse_cirilong(f) -> dict:
    df = pd.read_csv(f, sep="\t", dtype=str)
    loc = df["circ_id"].str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    df = df.join(loc)
    df["start0"] = df["start"].astype(int) - 1
    df["end"] = df["end"].astype(int)
    out = {}
    for (chrom, start, end, strand), g in df.groupby(["chrom", "start0", "end", "strand"]):
        out[(chrom, int(start), int(end), strand)] = set(g["read_id"])
    return out


def parse_circnick_native_counts(f, bed12_path) -> dict:
    """Read counts come from the native file (keyed by internal_circRNA_name).
    Coordinates come from the BED12 (post-liftover, matching every other
    file). See the module docstring."""
    reads_by_name = {}
    df = pd.read_csv(f, sep="\t", dtype=str)
    df["BSJ_reads"] = df["BSJ_reads"].astype(int)
    for _, row in df.iterrows():
        name = row["internal_circRNA_name"]
        reads_by_name[name] = max(reads_by_name.get(name, 0), row["BSJ_reads"])

    out = {}
    with open(bed12_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            name = cols[3]
            if name not in reads_by_name:
                continue
            key = (cols[0], int(cols[1]), int(cols[2]), cols[5])
            out[key] = max(out.get(key, 0), reads_by_name[name])
    return out


def load_tool_locus_counts(paths, run: str, tool: str, tiers: list = None) -> pd.DataFrame:
    """chrom/start/end/strand/n_reads DataFrame, one row per locus that tool
    called. The column is not named "count": pandas Series and itertuples
    both already have a built-in .count() method, so a column with that
    name would silently return the method, not the value."""
    native = paths.native_tool_paths(run)
    if tool == "isocirc":
        bed12_path = paths.caller_bed_path(run, "isocirc", tiers or [])
        readsets = parse_isocirc(native["isocirc"], bed12_path)
    elif tool == "circfl":
        readsets = parse_circfl(native["circfl"])
    elif tool == "cirilong":
        readsets = parse_cirilong(native["cirilong"])
    elif tool == "circnick":
        bed12_path = paths.caller_bed_path(run, "circnick", tiers or [])
        counts = parse_circnick_native_counts(native["circnick"], bed12_path)
        rows = [{"chrom": c, "start": s, "end": e, "strand": st, "n_reads": v} for (c, s, e, st), v in counts.items()]
        return pd.DataFrame(rows)
    else:
        raise ValueError(f"unknown tool: {tool}")
    rows = [{"chrom": c, "start": s, "end": e, "strand": st, "n_reads": len(v)} for (c, s, e, st), v in readsets.items()]
    return pd.DataFrame(rows)


def load_nanocirc_locus_counts(paths, run: str) -> pd.DataFrame:
    """Includes the nanocirc catalog's own bsj_id (needed downstream to check
    tier membership for whichever GT/tool locus a given row matches)."""
    df = pd.read_csv(paths.nanocirc_clean_counts_tsv(run), sep="\t")
    df = df.rename(columns={"#chrom": "chrom"}) if "#chrom" in df.columns else df
    df["base_bsj_id"] = df["bsj_id"].str.split("|").str[0]
    df = df.sort_values("nanocirc_quant_reads", ascending=False).drop_duplicates("base_bsj_id")
    return df.rename(columns={"nanocirc_quant_reads": "n_reads"})[["base_bsj_id", "chrom", "start", "end", "strand", "n_reads"]]


def nanocirc_catalog_bed_df(paths, run: str, work_dir) -> pd.DataFrame:
    """Full nanocirc catalog (every tier's rows) as a BED12 DataFrame, for
    exon-based matching. Built from the same clean_with_counts.tsv row for
    row as load_nanocirc_locus_counts, using its own block columns
    (sel_block_count/sel_block_sizes/sel_block_starts), so coordinates
    line up exactly with load_nanocirc_locus_counts's output."""
    df = pd.read_csv(paths.nanocirc_clean_counts_tsv(run), sep="\t")
    df = df.rename(columns={"#chrom": "chrom"}) if "#chrom" in df.columns else df
    bed_path = work_dir / "nanocirc_catalog.bed12"
    with open(bed_path, "w") as fh:
        for _, r in df.iterrows():
            fh.write("\t".join([
                str(r["chrom"]), str(r["start"]), str(r["end"]), str(r["bsj_id"]), "0", str(r["strand"]),
                str(r["start"]), str(r["start"]), "0",
                str(r["sel_block_count"]), str(r["sel_block_sizes"]), str(r["sel_block_starts"]),
            ]) + "\n")
    return read_bed(bed_path)


def load_tier_membership(paths, tiers: list) -> dict:
    """Reads the crossrun-merged bed12 (group-level, not per-sample). Must
    match what load_nanocirc_locus_counts reads (the crossrun-merged
    clean_with_counts.tsv). Otherwise, loci that are only in the crossrun
    catalog and not in one sample's own bed12 get wrongly left out of that
    tier's membership set."""
    membership = {}
    for tier in tiers:
        path = paths.crossrun_tier_bed12(tier)
        ids = set()
        if path.exists():
            with open(path) as f:
                for line in f:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) >= 4:
                        ids.add(cols[3].split("|")[0])
        membership[tier] = ids
    return membership
