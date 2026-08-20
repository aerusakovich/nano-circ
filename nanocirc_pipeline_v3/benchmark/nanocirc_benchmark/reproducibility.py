"""Pairwise read-count scatter matrix (density-colored), callers (ground
truth, each tool, each NanoCirc tier) on both axes, pooled across runs.
Every axis is matched back to GT with matching.match_pairs_by_strategy.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy.stats import gaussian_kde

from .bed_io import read_bed
from .callers import (load_nanocirc_isoform_bed_counts, load_nanocirc_locus_counts, load_tier_membership,
                       load_tool_isoform_bed_counts, load_tool_locus_counts)
from .gt import load_gt_loci
from .matching import (ISOFORM_MATCH_TAGS, STRATEGY_TITLES, StructuralMatcher, collapse_near_duplicates,
                        match_pairs_by_strategy, match_pairs_by_strategy_isoform_resolved, strategy_footnote)

_BED12_STRATEGIES = ("exon_based", "isoform_level") + ISOFORM_MATCH_TAGS


def build_caller_wide_df(paths, config, species: str, strategy: str = "relaxed_bsj") -> pd.DataFrame:
    is_isoform_resolved = strategy in ISOFORM_MATCH_TAGS
    needs_gt_bed = strategy in _BED12_STRATEGIES
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file) \
        if strategy in _BED12_STRATEGIES else None
    rows = []
    for run in paths.runs:
        gt_df = load_gt_loci(paths, run)
        gt_map = dict(zip(gt_df["bsj_id"], gt_df["n_reads_confirmed"]))
        gt_bed_df = read_bed(paths.gt_prefiltered_bed12(run)) if needs_gt_bed else None
        work = None
        if strategy in _BED12_STRATEGIES:
            work = paths.bedtools_work_dir / "reproducibility" / strategy / run
            work.mkdir(parents=True, exist_ok=True)

        if is_isoform_resolved:
            # bsj_id, not r.name (r.name is the circRNA's own id, e.g. "ecirna0000000",
            # not a coordinate string).
            per_key = {r.row_id: {"bsj_id": f"{r.chrom}:{r.start}-{r.end}:{r.strand}", "run": run,
                                   "ground_truth": int(r.score)}
                       for r in gt_bed_df.itertuples()}
        else:
            per_key = {bid: {"bsj_id": bid, "run": run, "ground_truth": expr} for bid, expr in gt_map.items()}

        for tool in config.tools:
            if strategy in _BED12_STRATEGIES:
                call_df = None
                call_bed_df = load_tool_isoform_bed_counts(paths, run, tool)
            else:
                call_df = load_tool_locus_counts(paths, run, tool, config.tiers)
                call_bed_df = None
            if is_isoform_resolved:
                matched = match_pairs_by_strategy_isoform_resolved(strategy, call_bed_df, gt_bed_df, config, matcher, work)
                for key in per_key:
                    per_key[key][tool] = sum(n for n, _ in matched.get(key, []))
            else:
                matched = match_pairs_by_strategy(strategy, call_df, gt_df, config, matcher, call_bed_df, gt_bed_df, work)
                for key in per_key:
                    per_key[key][tool] = sum(n for n, _ in matched.get(key, []))

        nc_df = load_nanocirc_locus_counts(paths, run)
        nc_bed_df = load_nanocirc_isoform_bed_counts(paths, run, None, work) if strategy in _BED12_STRATEGIES else None
        tier_members = load_tier_membership(paths, config.tiers)
        if is_isoform_resolved:
            nc_matched = match_pairs_by_strategy_isoform_resolved(strategy, nc_bed_df, gt_bed_df, config, matcher,
                                                                    work, extra_col="base_bsj_id")
            for tier in config.tiers:
                members = tier_members[tier]
                for key in per_key:
                    candidates = [(n, src_bid) for n, src_bid in nc_matched.get(key, []) if src_bid in members]
                    per_key[key][tier] = collapse_near_duplicates(candidates)
        else:
            nc_matched = match_pairs_by_strategy(strategy, nc_df, gt_df, config, matcher, nc_bed_df, gt_bed_df,
                                                  work, extra_col="base_bsj_id")
            for tier in config.tiers:
                members = tier_members[tier]
                for key in per_key:
                    candidates = [(n, src_bid) for n, src_bid in nc_matched.get(key, []) if src_bid in members]
                    per_key[key][tier] = collapse_near_duplicates(candidates)

        rows.extend(per_key.values())
    return pd.DataFrame(rows)


def plot_reproducibility_pairs(paths, config, species: str, out_path=None, strategy: str = "relaxed_bsj"):
    callers = ["ground_truth"] + config.quant_keys
    labels = {"ground_truth": "Ground Truth", **config.tool_labels, **config.tier_labels}
    labels = {c: labels[c].replace("\n", " ") for c in callers}

    unit = "isoforms" if strategy in ISOFORM_MATCH_TAGS else "loci"
    df = build_caller_wide_df(paths, config, species, strategy)
    df = df[(df[callers] >= 1).any(axis=1)]
    print(f"[{species}] {unit} included (>=1 read for at least one caller): {len(df)}")

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
             f"{unit.capitalize()} with >=1 read for at least one caller included (n={len(df):,}), pooled across {len(paths.runs)} runs. "
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
    """Runs plot_reproducibility_pairs for relaxed_bsj/strict_bsj/
    exon_based plus the 3 named isoform match tiers (default: all 6).
    relaxed_bsj keeps the plain, unsuffixed file name; every other strategy
    gets a suffix. Returns the list of output paths written."""
    strategies = strategies or ["relaxed_bsj", "strict_bsj", "exon_based", "full", "strict", "relaxed"]
    return [plot_reproducibility_pairs(paths, config, species, strategy=s) for s in strategies]
