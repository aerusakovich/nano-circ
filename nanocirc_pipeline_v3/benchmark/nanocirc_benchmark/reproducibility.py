"""Reproducibility scatter-matrix in the style of the isoCirc paper's
Supplementary Fig. 5 (pairwise BSJ read count, density-colored), but with
callers (ground truth, each tool, each NanoCirc tier) on the matrix axes
instead of replicate runs. Answers "does this caller's read count agree
with that caller's read count for the same locus", pooled across runs.

Every axis is matched back to GT with one of the 3 structural strategies
(matching.match_pairs_by_strategy), the same strategies performance.py
scores TP/FP/FN with, so all callers land on one common per-locus axis (a
caller's own called BSJ boundary is rarely byte-identical to another
caller's).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy.stats import gaussian_kde

from .bed_io import read_bed
from .callers import (load_nanocirc_locus_counts, load_tier_membership,
                       load_tool_locus_counts, nanocirc_catalog_bed_df)
from .gt import load_gt_loci
from .matching import (STRATEGY_TITLES, StructuralMatcher, collapse_near_duplicates,
                        match_pairs_by_strategy, strategy_footnote)


def build_caller_wide_df(paths, config, species: str, strategy: str = "relaxed_bsj") -> pd.DataFrame:
    """One row per (bsj_id, run), one column per caller (+ "ground_truth"),
    read count matched to that GT locus for this strategy (0 if
    unmatched)."""
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file) \
        if strategy == "exon_based" else None
    rows = []
    for run in paths.runs:
        gt_df = load_gt_loci(paths, run)
        gt_map = dict(zip(gt_df["bsj_id"], gt_df["n_reads_confirmed"]))
        gt_bed_df = read_bed(paths.gt_prefiltered_bed12(run)) if strategy == "exon_based" else None
        work = None
        if strategy == "exon_based":
            work = paths.bedtools_work_dir / "reproducibility" / strategy / run
            work.mkdir(parents=True, exist_ok=True)

        per_bsj = {bid: {"bsj_id": bid, "run": run, "ground_truth": expr} for bid, expr in gt_map.items()}

        for tool in config.tools:
            call_df = load_tool_locus_counts(paths, run, tool, config.tiers)
            call_bed_df = read_bed(paths.caller_bed_path(run, tool, config.tiers)) if strategy == "exon_based" else None
            matched = match_pairs_by_strategy(strategy, call_df, gt_df, config, matcher, call_bed_df, gt_bed_df, work)
            for bid in gt_map:
                per_bsj[bid][tool] = sum(n for n, _ in matched.get(bid, []))

        nc_df = load_nanocirc_locus_counts(paths, run)
        nc_bed_df = nanocirc_catalog_bed_df(paths, run, work) if strategy == "exon_based" else None
        nc_matched = match_pairs_by_strategy(strategy, nc_df, gt_df, config, matcher, nc_bed_df, gt_bed_df,
                                              work, extra_col="base_bsj_id")
        tier_members = load_tier_membership(paths, config.tiers)
        for tier in config.tiers:
            members = tier_members[tier]
            for bid in gt_map:
                candidates = [(n, src_bid) for n, src_bid in nc_matched.get(bid, []) if src_bid in members]
                per_bsj[bid][tier] = collapse_near_duplicates(candidates)

        rows.extend(per_bsj.values())
    return pd.DataFrame(rows)


def plot_reproducibility_pairs(paths, config, species: str, out_path=None, strategy: str = "relaxed_bsj"):
    callers = ["ground_truth"] + config.quant_keys
    labels = {"ground_truth": "Ground Truth", **config.tool_labels, **config.tier_labels}
    labels = {c: labels[c].replace("\n", " ") for c in callers}

    df = build_caller_wide_df(paths, config, species, strategy)
    # a locus is included if its count is >=1 for at least one caller
    df = df[(df[callers] >= 1).any(axis=1)]
    print(f"[{species}] loci included (>=1 read for at least one caller): {len(df)}")

    n = len(callers)
    fig, axes = plt.subplots(n, n, figsize=(2.6 * n + 2, 2.4 * n + 1.5))
    cmap = plt.get_cmap("plasma")

    def to_plot_coords(v):
        return np.log1p(v)

    max_val = df[callers].max().max()
    tick_vals = [0, 1, 2, 5, 10, int(max_val)] if max_val > 10 else [0, 1, 2, 5, int(max_val)]
    tick_vals = sorted(set(tick_vals))
    tick_pos = [np.log1p(t) for t in tick_vals]
    axis_limit = np.log1p(max_val * 1.15) if max_val else 1

    rng = np.random.default_rng(0)
    for i, caller_y in enumerate(callers):
        for j, caller_x in enumerate(callers):
            ax = axes[i, j]
            if i == j:
                ax.set_facecolor("#f0f0f0")
                ax.set_xticks([]); ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                continue

            x = to_plot_coords(df[caller_x].values)
            y = to_plot_coords(df[caller_y].values)
            x = x + rng.normal(0, 0.03, size=len(x))
            y = y + rng.normal(0, 0.03, size=len(y))

            if len(x) > 3:
                try:
                    density = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
                except np.linalg.LinAlgError:
                    density = np.ones(len(x))
            else:
                density = np.ones(len(x))
            order = np.argsort(density)
            xo, yo, densityo = x[order], y[order], density[order]

            ax.scatter(xo, yo, c=densityo, cmap=cmap, s=5, alpha=0.7,
                       norm=Normalize(vmin=densityo.min(), vmax=densityo.max()) if len(densityo) else None, linewidths=0)
            ax.plot([0, axis_limit], [0, axis_limit], color="#bbbbbb", lw=0.8, linestyle="--", zorder=0)

            ax.set_xlim(0, axis_limit); ax.set_ylim(0, axis_limit)
            ax.set_xticks(tick_pos); ax.set_yticks(tick_pos)
            if i == n - 1:
                ax.set_xticklabels([f"{t:,}" for t in tick_vals], fontsize=7, rotation=45, ha="right")
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_yticklabels([f"{t:,}" for t in tick_vals], fontsize=7)
            else:
                ax.set_yticklabels([])
            ax.grid(alpha=0.25, linewidth=0.5)

        axes[i, n - 1].text(1.06, 0.5, labels[callers[i]], transform=axes[i, n - 1].transAxes,
                            rotation=-90, va="center", ha="left", fontsize=9, fontweight="bold",
                            bbox=dict(boxstyle="square,pad=0.3", facecolor="#e0e0e0", edgecolor="none"))

    for j in range(n):
        axes[0, j].text(0.5, 1.15, labels[callers[j]], transform=axes[0, j].transAxes,
                        ha="center", va="bottom", fontsize=9, fontweight="bold",
                        bbox=dict(boxstyle="square,pad=0.3", facecolor="#e0e0e0", edgecolor="none"))

    fig.suptitle(f"Read count agreement across callers ({species.capitalize()}, {STRATEGY_TITLES[strategy]} match)",
                 fontsize=15, fontweight="bold", y=0.995)

    cbar_ax = fig.add_axes([0.93, 0.35, 0.015, 0.3])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=1))
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Density", fontsize=10, fontweight="bold")
    cbar.set_ticks([])

    fig.text(0.5, 0.005,
             f"Each panel: read count, caller on x-axis (top label) vs caller on y-axis (right label). "
             f"Matched to GT ({strategy_footnote(strategy, config)}).\n"
             f"Loci with >=1 read for at least one caller included (n={len(df):,}), pooled across {len(paths.runs)} runs. "
             "Axes: log1p-scaled read count.",
             ha="center", va="bottom", fontsize=8, color="#666666")

    plt.tight_layout(rect=[0.02, 0.03, 0.91, 0.95])
    suffix = "" if strategy == "relaxed_bsj" else f"_{strategy}"
    out_path = out_path or paths.plot_dir / f"reproducibility_pairs_{species}{suffix}.png"
    fig.savefig(out_path, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path


def run_reproducibility_analysis(paths, config, species: str, strategies=None) -> list:
    """Runs plot_reproducibility_pairs for each of the 3 structural
    strategies (default: all 3). relaxed_bsj keeps the plain, unsuffixed
    file name (the long-standing default this benchmark has always used).
    strict_bsj and exon_based get a strategy suffix. Returns the list of
    output paths written."""
    strategies = strategies or ["relaxed_bsj", "strict_bsj", "exon_based"]
    return [plot_reproducibility_pairs(paths, config, species, strategy=s) for s in strategies]
