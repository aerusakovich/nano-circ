"""Mature (spliced) length distribution by caller, TP calls only. A call's
length only means something once it is known to match a real circRNA.

Also the structural TP/FP/FN contingency breakdown by matching strategy
(reuses performance.py's own benchmark_metrics_<species>_<tag>.csv).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bed_io import read_bed
from .gt import count_reads
from .matching import StructuralMatcher


def plot_length_barplot(paths, config, species: str, out_path=None) -> pd.DataFrame:
    caller_order = ["ground_truth"] + config.quant_keys
    caller_labels = {"ground_truth": "Ground Truth", **config.tool_labels, **config.tier_labels}
    caller_colors = {"ground_truth": "#808080", **config.tool_colors, **config.tier_colors}

    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file)
    work_root = paths.bedtools_work_dir / "length_tp"

    rows = []
    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        for length in gt_df["spliced_len"]:
            rows.append({"caller": "ground_truth", "run": run, "length": length})

        for q in config.quant_keys:
            path = paths.caller_bed_path(run, q, config.tiers)
            if not path.exists():
                continue
            caller_df = read_bed(path)
            if len(caller_df) == 0:
                continue
            work = work_root / run / q
            work.mkdir(parents=True, exist_ok=True)
            _, _, _, _, matched_caller_ids = matcher.exon_based_counts(
                caller_df, gt_df, work, config.tp_min_overlap, return_ids="both")
            tp_lengths = caller_df.loc[caller_df["row_id"].isin(matched_caller_ids), "spliced_len"]
            for length in tp_lengths:
                rows.append({"caller": q, "run": run, "length": length})

    df = pd.DataFrame(rows)
    df.to_csv(paths.bench_dir / f"length_data_{species}.csv", index=False)

    stats = df.groupby("caller")["length"].agg(["mean", "sem", "count"]).reindex(caller_order)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = np.arange(len(caller_order))
    colors = [caller_colors[c] for c in caller_order]
    ax.bar(x, stats["mean"], yerr=stats["sem"], color=colors, edgecolor="black",
           linewidth=1.0, capsize=4, alpha=0.9)
    max_mean = stats["mean"].dropna().max() if stats["mean"].notna().any() else 0
    for xi, (m, se, n) in enumerate(zip(stats["mean"], stats["sem"], stats["count"])):
        if pd.isna(m):
            continue
        ax.text(xi, m + se + max_mean * 0.02, f"{m:.0f}±{se:.1f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")
        ax.text(xi, m / 2, f"n={n:,}", ha="center", va="center", fontsize=9, fontweight="bold", color="white")
    ax.set_xticks(x)
    ax.set_xticklabels([caller_labels[c].replace("\n", " ") for c in caller_order], fontsize=10, rotation=15, ha="right")
    ax.set_ylabel("Mean mature (spliced) length, bp", fontsize=12, fontweight="bold")
    n_reads = count_reads(paths.run_fastq(paths.runs[0]))
    ax.set_title(f"Mean circRNA mature length by caller ({species.capitalize()}, {len(paths.runs)} runs, {n_reads:,} reads)",
                 fontsize=13.5, fontweight="bold")
    fig.text(0.5, -0.02, "TP calls only. Error bars show standard error (SE).",
              ha="center", va="top", fontsize=8, color="#888888")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    out_path = out_path or paths.plot_dir / f"length_barplot_{species}.png"
    fig.savefig(out_path, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return df


def plot_contingency_heatmap(paths, config, species: str, metrics_csv_tag="relaxed_10_75", out_path=None):
    """Faceted by the 3 matching strategies (strict_bsj, relaxed_bsj,
    exon_based), sourced from benchmark_metrics_<species>_<tag>.csv.
    strict_bsj is an exact match with no parameters, so it is identical
    whichever tag it is read from."""
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    strategy_label_map = {"strict_bsj": "Strict BSJ", "relaxed_bsj": "Relaxed BSJ", "exon_based": "Exon-based"}
    strategies = ["strict_bsj", "relaxed_bsj", "exon_based"]

    full_csv = paths.bench_dir / f"benchmark_metrics_{species}_{metrics_csv_tag}.csv"
    fig, axes = plt.subplots(3, 3, figsize=(17, 14), sharey=False)
    mc = {"TP": "#2ECC71", "FP": "#E74C3C", "FN": "#95A5A6"}

    cdf_all = pd.read_csv(full_csv)
    for row_idx, strategy in enumerate(strategies):
        cdf = cdf_all[cdf_all["strategy"] == strategy]
        summed = cdf.groupby("caller")[["TP", "FP", "FN"]].sum().reindex(quant_keys)

        for col_idx, metric in enumerate(["TP", "FP", "FN"]):
            ax = axes[row_idx, col_idx]
            vals = summed[metric].fillna(0)
            ax.bar(range(len(quant_keys)), vals, color=mc[metric], edgecolor="black", linewidth=0.9)
            max_v = max(vals) if len(vals) and vals.max() > 0 else 1
            for i, v in enumerate(vals):
                ax.text(i, v * 1.01 + max_v * 0.01, f"{int(v)}", ha="center", va="bottom", fontsize=9)
            ax.set_xticks(range(len(quant_keys)))
            ax.set_xticklabels([quant_labels[c].replace("\n", " ") for c in quant_keys], rotation=30, ha="right", fontsize=8.5)
            if col_idx == 0:
                ax.set_ylabel(f"{strategy_label_map[strategy]}\ncount", fontsize=10, fontweight="bold")
            if row_idx == 0:
                ax.set_title(metric, fontsize=14, fontweight="bold", color=mc[metric])
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)

    fig.suptitle(f"Structural detection contingency ({species.capitalize()}): TP/FP/FN by matching strategy, summed across runs",
                 fontsize=15, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = out_path or paths.plot_dir / f"tool_contingency_heatmap_{species}.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path


def plot_gt_coverage_stacked(paths, config, species: str, metrics_csv_tag="relaxed_10_75", out_path=None):
    """GT-normalized stacked bar per caller. Green (recovered) plus grey
    (missed) always sum to exactly 100% of ground truth, so the split shows
    that caller's recall directly. Red (FP) stacks on top as a separate
    overshoot above the 100% line, on the same scale, so a caller with many
    false calls visibly shoots past 100.

    Uses FN and each run's own GT total to compute recovered/missed, not
    the CSV's "TP" column. TP counts matched caller rows, not matched GT
    rows, and those two can differ under relaxed matching (one caller row
    can satisfy several GT rows). Using GT-side counts keeps the
    green+grey=100% rule always true.
    """
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    strategy_label_map = {"strict_bsj": "Strict BSJ", "relaxed_bsj": "Relaxed BSJ", "exon_based": "Exon-based"}
    strategies = ["strict_bsj", "relaxed_bsj", "exon_based"]
    mc = {"TP": "#2ECC71", "FP": "#E74C3C", "FN": "#95A5A6"}

    gt_totals = {run: len(read_bed(paths.gt_prefiltered_bed12(run))) for run in paths.runs}

    full_csv = paths.bench_dir / f"benchmark_metrics_{species}_{metrics_csv_tag}.csv"
    cdf_all = pd.read_csv(full_csv)
    cdf_all = cdf_all[cdf_all["caller"].isin(quant_keys)].copy()
    cdf_all["gt_total"] = cdf_all["sample"].map(gt_totals)
    cdf_all["recovered_pct"] = 100 * (cdf_all["gt_total"] - cdf_all["FN"]) / cdf_all["gt_total"]
    cdf_all["missed_pct"] = 100 * cdf_all["FN"] / cdf_all["gt_total"]
    cdf_all["fp_pct"] = 100 * cdf_all["FP"] / cdf_all["gt_total"]

    fig, axes = plt.subplots(1, 3, figsize=(21, 6.5), sharey=True)
    max_top = 100.0
    for row_idx, strategy in enumerate(strategies):
        sdf = cdf_all[cdf_all["strategy"] == strategy]
        means = sdf.groupby("caller")[["recovered_pct", "missed_pct", "fp_pct"]].mean().reindex(quant_keys)
        fp_std = sdf.groupby("caller")["fp_pct"].std().reindex(quant_keys).fillna(0)
        max_top = max(max_top, (100 + means["fp_pct"] + fp_std).max())

    for col_idx, strategy in enumerate(strategies):
        ax = axes[col_idx]
        sdf = cdf_all[cdf_all["strategy"] == strategy]
        means = sdf.groupby("caller")[["recovered_pct", "missed_pct", "fp_pct"]].mean().reindex(quant_keys).fillna(0)
        stds = sdf.groupby("caller")[["recovered_pct", "fp_pct"]].std().reindex(quant_keys).fillna(0)

        x = range(len(quant_keys))
        # TP and FN stack bottom-up to exactly 100 (GT size) by construction.
        # FP stacks separately starting AT the 100 line, so it always reads as
        # pure overshoot past GT's size, never overlapping the FN segment.
        ax.bar(x, means["recovered_pct"], color=mc["TP"], edgecolor="black", linewidth=0.9, label="TP (GT recovered)")
        ax.bar(x, means["missed_pct"], bottom=means["recovered_pct"], color=mc["FN"], edgecolor="black",
               linewidth=0.9, label="FN (GT missed)")
        ax.bar(x, means["fp_pct"], bottom=100.0, color=mc["FP"], edgecolor="black",
               linewidth=0.9, hatch="///", alpha=0.85, label="FP (extra, no GT match)")
        ax.errorbar(x, means["recovered_pct"], yerr=stds["recovered_pct"], fmt="none",
                    ecolor="black", elinewidth=1.1, capsize=3, zorder=10)
        ax.errorbar(x, 100 + means["fp_pct"], yerr=stds["fp_pct"], fmt="none",
                    ecolor="#8B0000", elinewidth=1.1, capsize=3, zorder=10)

        ax.axhline(100, color="black", lw=1.0, linestyle="--", alpha=0.6)
        for i, caller in enumerate(quant_keys):
            rec, fp = means.loc[caller, "recovered_pct"], means.loc[caller, "fp_pct"]
            ax.text(i, rec / 2, f"{rec:.0f}", ha="center", va="center", fontsize=8, fontweight="bold", color="white")
            if fp > 1:
                ax.text(i, 100 + fp / 2, f"+{fp:.0f}", ha="center", va="center", fontsize=7.5,
                        fontweight="bold", color="white")

        ax.set_xticks(list(x))
        ax.set_xticklabels([quant_labels[c].replace("\n", " ") for c in quant_keys], rotation=30, ha="right", fontsize=8.5)
        ax.set_ylim(0, max_top * 1.1)
        ax.set_title(strategy_label_map[strategy], fontsize=13, fontweight="bold")
        if col_idx == 0:
            ax.set_ylabel("% of GT", fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=mc["TP"], edgecolor="black", label="TP (GT recovered)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=mc["FN"], edgecolor="black", label="FN (GT missed)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=mc["FP"], edgecolor="black", alpha=0.85, hatch="///", label="FP (extra, no GT match)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=11, frameon=True, bbox_to_anchor=(0.5, -0.06))
    n_reads = count_reads(paths.run_fastq(paths.runs[0]))
    fig.suptitle(f"GT coverage by caller ({species.capitalize()}, {len(paths.runs)} runs, {n_reads:,} reads)",
                 fontsize=14, fontweight="bold", y=1.03)
    fig.text(0.5, -0.09,
              "Green+grey always = 100% of GT (that caller's recall); red overshoot above 100% = false calls "
              "relative to GT size. Mean +/-1 SD across runs.",
              ha="center", va="top", fontsize=8, color="#888888")
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    out_path = out_path or paths.plot_dir / f"gt_coverage_stacked_{species}.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path
