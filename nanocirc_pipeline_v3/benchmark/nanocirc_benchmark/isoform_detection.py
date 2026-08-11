"""Isoform-level detection within multi-isoform BSJ families.

Long-read protocols exist mainly to resolve internal exon structure, not
just back-splice junction (BSJ) position. This asks a sharper question than
the rest of the suite's BSJ/exon matching: among GT loci where the same BSJ
genuinely has more than one real isoform (different internal exon
structure), how often does each caller recover the right isoform, versus
call something at the right BSJ with the wrong structure, versus miss the
BSJ entirely.

Three matching strategies, same three-way classification either way:
  block_identity          exact block-for-block match (blockSizes/
                           blockStarts). GT isoforms are precise known
                           structures, so this is the strict, no-partial-
                           credit view.
  reciprocal_overlap_99   >=99% reciprocal exon overlap, via the same
                           exon_based approach matching.py uses elsewhere.
                           Allows only boundary noise: in practice
                           this is "the exact isoform, +- several bp".
  reciprocal_overlap_90   >=90% reciprocal exon overlap. Allows mismatch, still
                           requiring substantially the same structure, not
                           a different one.
Both reciprocal_overlap thresholds live in config.isoform_overlap_param_sets
(same {tag, min_overlap} shape as structural_param_sets elsewhere), tighter
than overlap.py's 95% on purpose: this discriminates between real isoforms
at one locus, not just "was a locus found at all".

For the reciprocal_overlap strategies, a GT isoform that isn't a true_iso
match is 'wrong_iso' only if the caller made some call within
config.gt_match_tolerance bp of its own BSJ boundary (not literally the
same boundary): the strategy already forgives small boundary noise on the
inside (the overlap score), so the found/missed split should use the same
tolerance rather than an exact-string test, or a caller with a slightly
shifted BSJ call would wrongly register as having missed the locus
entirely. block_identity keeps the exact-string test throughout, matching
its own no-tolerance philosophy.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .bed_io import read_bed
from .gt import count_reads
from .matching import StructuralMatcher

FAMILY_SIZE_BUCKETS = [("2", lambda n: n == 2), ("3", lambda n: n == 3), ("4+", lambda n: n >= 4)]
_STRANDED = ("+", "-")


def _block_key(row):
    return (row.blockSizes.rstrip(","), row.blockStarts.rstrip(","))


def _bsj_key(row):
    return (row.chrom, row.start, row.end, row.strand)


def _strategies(config) -> list:
    return ["block_identity"] + [p["tag"] for p in config.isoform_overlap_param_sets]


def _strategy_min_overlap(strategy: str, config) -> float:
    for p in config.isoform_overlap_param_sets:
        if p["tag"] == strategy:
            return p["min_overlap"]
    raise ValueError(f"unknown reciprocal_overlap strategy: {strategy!r}")


def _strategy_suffix(strategy: str) -> str:
    return "" if strategy == "block_identity" else f"_{strategy}"


def _strategy_title(strategy: str, config) -> str:
    if strategy == "block_identity":
        return "exact block match"
    return f">={_strategy_min_overlap(strategy, config):.0%} reciprocal exon overlap"


def _locus_called_row_ids(families: pd.DataFrame, caller_df: pd.DataFrame, tolerance: int) -> set:
    """row_ids of `families` where caller_df has a call within `tolerance`
    bp of both boundaries (same strand, or GT unstranded): the caller found
    something at roughly this BSJ, regardless of internal structure."""
    if len(caller_df) == 0 or len(families) == 0:
        return set()
    matched = set()
    for chrom, group in families.groupby("chrom"):
        cand = caller_df[caller_df["chrom"] == chrom]
        if cand.empty:
            continue
        for row in group.itertuples():
            strand_ok = (cand["strand"] == row.strand) | (row.strand not in _STRANDED)
            ok = strand_ok & ((cand["start"] - row.start).abs() <= tolerance) & ((cand["end"] - row.end).abs() <= tolerance)
            if ok.any():
                matched.add(row.row_id)
    return matched


def find_multi_isoform_families(gt_df: pd.DataFrame) -> pd.DataFrame:
    """GT rows belonging to a BSJ shared by >=2 distinct exon structures.
    Adds family_id (the shared BSJ key, stringified) and family_size."""
    gt_df = gt_df.copy()
    gt_df["bsj_key"] = gt_df.apply(_bsj_key, axis=1)
    gt_df["block_key"] = gt_df.apply(_block_key, axis=1)

    keep_rows = []
    for bsj, group in gt_df.groupby("bsj_key"):
        distinct_structures = group["block_key"].nunique()
        if distinct_structures < 2:
            continue
        g = group.copy()
        g["family_id"] = f"{bsj[0]}:{bsj[1]}-{bsj[2]}:{bsj[3]}"
        g["family_size"] = distinct_structures
        keep_rows.append(g)

    if not keep_rows:
        return pd.DataFrame(columns=list(gt_df.columns) + ["family_id", "family_size"])
    return pd.concat(keep_rows, ignore_index=True)


def compute_isoform_detection(paths, config, species: str, strategy: str = "block_identity") -> pd.DataFrame:
    """One row per (run, caller, GT isoform row) in a multi-isoform family.
    classification is 'true_iso' (matched under `strategy`), 'wrong_iso'
    (caller has a call near this locus, structure doesn't match this
    isoform under `strategy`), or 'missed' (no caller call anywhere near
    this locus). See module docstring for how the wrong_iso/missed boundary
    itself differs between block_identity and the reciprocal_overlap
    strategies."""
    valid = _strategies(config)
    if strategy not in valid:
        raise ValueError(f"unknown strategy: {strategy!r}, expected one of {valid}")

    is_reciprocal = strategy != "block_identity"
    matcher = StructuralMatcher(bedtools_bin=config.bedtools_bin) if is_reciprocal else None
    min_overlap = _strategy_min_overlap(strategy, config) if is_reciprocal else None

    rows = []
    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        families = find_multi_isoform_families(gt_df)
        if families.empty:
            continue

        for caller in config.quant_keys:
            caller_path = paths.caller_bed_path(run, caller, config.tiers)
            if not caller_path.exists():
                caller_bsj_keys, true_row_ids, locus_called_row_ids = set(), set(), set()
            else:
                caller_df = read_bed(caller_path)

                if not is_reciprocal:
                    caller_bsj_keys = {_bsj_key(r) for r in caller_df.itertuples()}
                    caller_full_keys = {(_bsj_key(r), _block_key(r)) for r in caller_df.itertuples()}
                    true_row_ids = {r.row_id for r in families.itertuples()
                                    if (r.bsj_key, r.block_key) in caller_full_keys}
                    locus_called_row_ids = None
                else:
                    work = paths.bedtools_work_dir / "isoform_detection" / strategy / run / caller
                    work.mkdir(parents=True, exist_ok=True)
                    _, _, _, true_row_ids = matcher.exon_based_counts(
                        caller_df, families, work, min_overlap=min_overlap, return_ids=True)
                    locus_called_row_ids = _locus_called_row_ids(families, caller_df, config.gt_match_tolerance)
                    caller_bsj_keys = None

            for r in families.itertuples():
                if r.row_id in true_row_ids:
                    classification = "true_iso"
                elif not is_reciprocal:
                    classification = "wrong_iso" if r.bsj_key in caller_bsj_keys else "missed"
                else:
                    classification = "wrong_iso" if r.row_id in locus_called_row_ids else "missed"
                rows.append({
                    "run": run, "caller": caller, "family_id": r.family_id,
                    "family_size": r.family_size, "row_id": r.row_id,
                    "classification": classification,
                })

    return pd.DataFrame(rows)


def plot_isoform_detection_stacked(df: pd.DataFrame, paths, config, species: str,
                                    strategy: str = "block_identity", out_path=None):
    """Same visual convention as plot_gt_coverage_stacked: green (true_iso)
    + grey (missed) stack to 100% of the multi-isoform-family GT population,
    red (wrong_iso) stacks separately on top as its own share of that same
    100%, since a wrong-structure call is neither a correct recovery nor a
    simple miss."""
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    mc = {"true_iso": "#2ECC71", "wrong_iso": "#E74C3C", "missed": "#95A5A6"}

    if df.empty:
        raise ValueError(f"No multi-isoform BSJ families found in GT for species={species}")

    counts = (df.groupby(["run", "caller", "classification"]).size()
              .unstack(fill_value=0).reindex(columns=["true_iso", "wrong_iso", "missed"], fill_value=0))
    pct = counts.div(counts.sum(axis=1), axis=0) * 100
    means = pct.groupby("caller").mean().reindex(quant_keys).fillna(0)
    stds = pct.groupby("caller").std().reindex(quant_keys).fillna(0)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    x = range(len(quant_keys))
    ax.bar(x, means["true_iso"], color=mc["true_iso"], edgecolor="black", linewidth=0.9, label="Correct isoform")
    ax.bar(x, means["missed"], bottom=means["true_iso"], color=mc["missed"], edgecolor="black",
           linewidth=0.9, label="Missed (BSJ not called)")
    ax.bar(x, means["wrong_iso"], bottom=means["true_iso"] + means["missed"], color=mc["wrong_iso"],
           edgecolor="black", linewidth=0.9, hatch="///", alpha=0.85, label="Wrong isoform (right BSJ)")
    ax.errorbar(x, means["true_iso"], yerr=stds["true_iso"], fmt="none", ecolor="black",
                elinewidth=1.1, capsize=3, zorder=10)

    for i, caller in enumerate(quant_keys):
        v = means.loc[caller, "true_iso"]
        if v > 1:
            ax.text(i, v / 2, f"{v:.0f}", ha="center", va="center", fontsize=8, fontweight="bold", color="white")

    ax.set_xticks(list(x))
    ax.set_xticklabels([quant_labels[c].replace("\n", " ") for c in quant_keys], rotation=30, ha="right", fontsize=9)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of GT isoforms in multi-isoform BSJ families", fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.legend(loc="upper right", fontsize=9, frameon=True)

    n_reads = count_reads(paths.run_fastq(paths.runs[0]))
    n_families = df["family_id"].nunique()
    fig.suptitle(f"Isoform detection within shared-BSJ families ({species.capitalize()}, "
                 f"{_strategy_title(strategy, config)}, {len(paths.runs)} runs, {n_reads:,} reads)",
                 fontsize=13, fontweight="bold")
    fig.text(0.5, -0.02,
              f"{n_families} GT BSJ positions with 2+ real isoforms, pooled across runs. "
              "True+missed = 100%; wrong-isoform is a separate share of the same 100%, "
              "since it is neither a correct call nor a simple miss.",
              ha="center", va="top", fontsize=8, color="#888888")
    plt.tight_layout()
    out_path = out_path or paths.plot_dir / f"isoform_detection{_strategy_suffix(strategy)}_{species}.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path


def plot_isoform_detection_by_family_size(df: pd.DataFrame, paths, config, species: str,
                                           strategy: str = "block_identity", out_path=None):
    """Same stacked convention, one panel per family-size bucket (2, 3, 4+
    isoforms sharing one BSJ), to see whether callers get worse as more
    isoforms compete at the same junction."""
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    mc = {"true_iso": "#2ECC71", "wrong_iso": "#E74C3C", "missed": "#95A5A6"}

    if df.empty:
        raise ValueError(f"No multi-isoform BSJ families found in GT for species={species}")

    fig, axes = plt.subplots(1, len(FAMILY_SIZE_BUCKETS), figsize=(7 * len(FAMILY_SIZE_BUCKETS), 6.5), sharey=True)
    for col_idx, (label, pred) in enumerate(FAMILY_SIZE_BUCKETS):
        ax = axes[col_idx]
        sub = df[df["family_size"].apply(pred)]
        if sub.empty:
            ax.set_title(f"{label} isoforms/BSJ (n=0)", fontsize=12, fontweight="bold")
            ax.axis("off")
            continue

        counts = (sub.groupby(["run", "caller", "classification"]).size()
                  .unstack(fill_value=0).reindex(columns=["true_iso", "wrong_iso", "missed"], fill_value=0))
        pct = counts.div(counts.sum(axis=1), axis=0) * 100
        means = pct.groupby("caller").mean().reindex(quant_keys).fillna(0)

        x = range(len(quant_keys))
        ax.bar(x, means["true_iso"], color=mc["true_iso"], edgecolor="black", linewidth=0.9)
        ax.bar(x, means["missed"], bottom=means["true_iso"], color=mc["missed"], edgecolor="black", linewidth=0.9)
        ax.bar(x, means["wrong_iso"], bottom=means["true_iso"] + means["missed"], color=mc["wrong_iso"],
               edgecolor="black", linewidth=0.9, hatch="///", alpha=0.85)

        ax.set_xticks(list(x))
        ax.set_xticklabels([quant_labels[c].replace("\n", " ") for c in quant_keys], rotation=30, ha="right", fontsize=8.5)
        n_fam = sub["family_id"].nunique()
        ax.set_title(f"{label} isoforms/BSJ (n={n_fam} families)", fontsize=12, fontweight="bold")
        if col_idx == 0:
            ax.set_ylabel("% of GT isoforms", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 105)
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=mc["true_iso"], edgecolor="black", label="Correct isoform"),
        plt.Rectangle((0, 0), 1, 1, facecolor=mc["missed"], edgecolor="black", label="Missed (BSJ not called)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=mc["wrong_iso"], edgecolor="black", alpha=0.85, hatch="///", label="Wrong isoform (right BSJ)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.05))
    fig.suptitle(f"Isoform detection by family size ({species.capitalize()}, {_strategy_title(strategy, config)})",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    out_path = out_path or paths.plot_dir / f"isoform_detection_by_family_size{_strategy_suffix(strategy)}_{species}.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path


def plot_isoform_family_completeness(df: pd.DataFrame, paths, config, species: str,
                                      strategy: str = "block_identity", out_path=None):
    """Per family (not per isoform): did this caller get every isoform in
    the family right, some of them, or none? A caller can score well on the
    pooled isoform-level view above by getting one isoform right in many
    families while never fully resolving any single one; this view catches
    that."""
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    mc = {"all": "#2ECC71", "some": "#F1C40F", "none": "#95A5A6"}

    if df.empty:
        raise ValueError(f"No multi-isoform BSJ families found in GT for species={species}")

    fam = (df.groupby(["run", "caller", "family_id"])["classification"]
           .apply(lambda s: "all" if (s == "true_iso").all()
                  else ("none" if (s == "true_iso").sum() == 0 else "some"))
           .reset_index(name="outcome"))

    counts = (fam.groupby(["run", "caller", "outcome"]).size()
              .unstack(fill_value=0).reindex(columns=["all", "some", "none"], fill_value=0))
    pct = counts.div(counts.sum(axis=1), axis=0) * 100
    means = pct.groupby("caller").mean().reindex(quant_keys).fillna(0)
    stds = pct.groupby("caller").std().reindex(quant_keys).fillna(0)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    x = range(len(quant_keys))
    ax.bar(x, means["all"], color=mc["all"], edgecolor="black", linewidth=0.9, label="Family fully resolved")
    ax.bar(x, means["some"], bottom=means["all"], color=mc["some"], edgecolor="black",
           linewidth=0.9, label="Family partially resolved")
    ax.bar(x, means["none"], bottom=means["all"] + means["some"], color=mc["none"], edgecolor="black",
           linewidth=0.9, label="Family unresolved")
    ax.errorbar(x, means["all"], yerr=stds["all"], fmt="none", ecolor="black", elinewidth=1.1, capsize=3, zorder=10)

    for i, caller in enumerate(quant_keys):
        v = means.loc[caller, "all"]
        if v > 1:
            ax.text(i, v / 2, f"{v:.0f}", ha="center", va="center", fontsize=8, fontweight="bold", color="white")

    ax.set_xticks(list(x))
    ax.set_xticklabels([quant_labels[c].replace("\n", " ") for c in quant_keys], rotation=30, ha="right", fontsize=9)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of multi-isoform BSJ families", fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.legend(loc="upper right", fontsize=9, frameon=True)

    n_families = df["family_id"].nunique()
    fig.suptitle(f"Isoform family completeness ({species.capitalize()}, {_strategy_title(strategy, config)}, "
                 f"{len(paths.runs)} runs)", fontsize=13, fontweight="bold")
    fig.text(0.5, -0.02,
              f"{n_families} GT BSJ positions with 2+ real isoforms. A family counts as fully resolved "
              "only if every one of its real isoforms was called with the exact right structure.",
              ha="center", va="top", fontsize=8, color="#888888")
    plt.tight_layout()
    out_path = out_path or paths.plot_dir / f"isoform_family_completeness{_strategy_suffix(strategy)}_{species}.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path


def run_isoform_detection_analysis(paths, config, species: str, strategies=None) -> dict:
    """Runs block_identity and every config.isoform_overlap_param_sets
    strategy (all by default), each producing its own CSV + 3 plots
    (suffixed for non-default strategies). Returns {strategy: df}."""
    strategies = strategies or _strategies(config)
    results = {}
    for strategy in strategies:
        print(f"[{species}] isoform detection, {strategy} ({_strategy_title(strategy, config)}) ...")
        df = compute_isoform_detection(paths, config, species, strategy)
        results[strategy] = df
        if df.empty:
            print(f"[{species}] no multi-isoform BSJ families found in GT, skipping isoform detection plots")
            continue
        suffix = _strategy_suffix(strategy)
        df.to_csv(paths.bench_dir / f"isoform_detection{suffix}_{species}.csv", index=False)
        plot_isoform_detection_stacked(df, paths, config, species, strategy)
        plot_isoform_detection_by_family_size(df, paths, config, species, strategy)
        plot_isoform_family_completeness(df, paths, config, species, strategy)
    return results