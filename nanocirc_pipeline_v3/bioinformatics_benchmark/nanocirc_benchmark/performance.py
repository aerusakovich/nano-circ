"""Structural precision/recall/F1 (strict_bsj / relaxed_bsj / exon_based),
and the resulting recall-vs-precision scatter across all 3 strategies.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from .bed_io import read_bed
from .gt import count_reads
from .matching import StructuralMatcher

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "grid.color": "#eeeeee",
    "grid.linewidth": 0.8,
})


def run_structural_matching(paths, config, species: str) -> dict:
    """Runs strict_bsj/relaxed_bsj/exon_based for every caller x run x
    structural param set, writing benchmark_metrics_<species>_<tag>.csv.
    Returns {tag: DataFrame}."""
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file)
    out = {}
    for params in config.structural_param_sets:
        tag, tolerance, min_overlap = params["tag"], params["bsj_tolerance"], params["min_overlap"]
        print(f"\n{'=' * 60}\n{species}, param set: {tag}\n{'=' * 60}")

        rows = []
        for run in paths.runs:
            gt_df = read_bed(paths.gt_prefiltered_bed12(run))
            print(f"\n[{run}] GT: {len(gt_df)} circRNAs")
            work = paths.bedtools_work_dir / tag / run
            work.mkdir(parents=True, exist_ok=True)

            for caller, path in paths.caller_paths(run, config.tools, config.tiers).items():
                if not path.exists():
                    print(f"  [SKIP] {caller}, file not found: {path}")
                    continue
                caller_df = read_bed(path)
                print(f"  {caller:28s} n={len(caller_df):5d}", end="  ")

                strict = matcher.strict_bsj_counts(caller_df, gt_df, work)
                relaxed = matcher.relaxed_bsj_counts(caller_df, gt_df, work, tolerance)
                exon = matcher.exon_based_counts(caller_df, gt_df, work, min_overlap)

                for strategy, (tp, fp, fn) in [("strict_bsj", strict), ("relaxed_bsj", relaxed), ("exon_based", exon)]:
                    prec, rec, f1 = matcher.compute_metrics(tp, fp, fn)
                    rows.append({"sample": run, "caller": caller, "strategy": strategy,
                                 "TP": tp, "FP": fp, "FN": fn,
                                 "precision": round(prec, 6), "recall": round(rec, 6), "f1": round(f1, 6)})
                s_p, s_r, _ = matcher.compute_metrics(*strict)
                r_p, r_r, _ = matcher.compute_metrics(*relaxed)
                e_p, e_r, _ = matcher.compute_metrics(*exon)
                print(f"strict P={s_p:.3f}  relaxed P={r_p:.3f} R={r_r:.3f}  exon P={e_p:.3f} R={e_r:.3f}")

        df = pd.DataFrame(rows)
        out_csv = paths.bench_dir / f"benchmark_metrics_{species}_{tag}.csv"
        df.to_csv(out_csv, index=False)
        print(f"\n-> {out_csv}")
        out[tag] = df
    return out


def plot_pr_scatter(paths, config, species: str, metrics_df: pd.DataFrame, out_path) -> None:
    """PR scatter across the 3 structural strategies, pooled (mean) across
    runs, for every tool + every NanoCirc tier. `metrics_df` should be the
    relaxed-tolerance param set's table (matches the tolerance quoted in the
    plot footnote)."""
    n_reads = count_reads(paths.run_fastq(paths.runs[0]))
    pooled = metrics_df.groupby(["caller", "strategy"])[["precision", "recall"]].mean()

    def pt(caller, strategy):
        row = pooled.loc[(caller, strategy)]
        return (row["recall"], row["precision"])

    strategies = ["exon_based", "relaxed_bsj", "strict_bsj"]
    competitors = {config.tool_labels[t]: {"color": config.tool_colors[t],
                                            "pts": [pt(t, s) for s in strategies]} for t in config.tools}
    nc_modes = {config.tier_labels[t].replace("\n", " "): {"color": config.tier_colors[t],
                                                            "pts": [pt(t, s) for s in strategies]} for t in config.tiers}

    METRICS = ["Exon", "Relaxed BSJ", "Strict BSJ"]
    SIZES_C = [260, 120, 55]
    SIZES_NC = [420, 200, 90]
    ALPHAS = [1.0, 0.72, 0.50]

    fig, ax = plt.subplots(figsize=(12, 7.5))
    ax.grid(True, zorder=0)

    all_recalls = [p[0] for v in list(competitors.values()) + list(nc_modes.values()) for p in v["pts"]]
    r_max = max(all_recalls) * 1.15

    r_iso = np.linspace(0.005, r_max, 500)
    for f1 in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60]:
        denom = 2 * r_iso - f1
        p_iso = np.where(denom > 1e-9, f1 * r_iso / denom, np.nan)
        mask = (p_iso >= 0) & (p_iso <= 1.02)
        ax.plot(r_iso[mask], p_iso[mask], color="#e0e0e0", lw=0.9, zorder=0)
        r_lbl = min(r_max * 0.85, f1 * 0.62)
        dl = 2 * r_lbl - f1
        if dl > 1e-9:
            p_lbl = f1 * r_lbl / dl
            if 0.05 < p_lbl < 1.0:
                ax.text(r_lbl, p_lbl, f"F1={int(f1*100)}%", fontsize=7, color="#c0c0c0",
                        ha="center", va="bottom", rotation=-35, zorder=1)

    FRONTIER_STYLES = [
        dict(lw=2.2, linestyle="--", alpha=0.80, label="NanoCirc frontier (Exon)"),
        dict(lw=1.5, linestyle="-.", alpha=0.55, label="NanoCirc frontier (Relaxed BSJ)"),
        dict(lw=1.0, linestyle=":", alpha=0.45, label="NanoCirc frontier (Strict BSJ)"),
    ]
    for m_idx in range(3):
        pts = sorted([(v["pts"][m_idx][0], v["pts"][m_idx][1]) for v in nc_modes.values()])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#444444", zorder=2, **FRONTIER_STYLES[m_idx])

    for tool in competitors.values():
        rr = [p[0] for p in tool["pts"]]
        pp = [p[1] for p in tool["pts"]]
        ax.plot(rr, pp, "-", color=tool["color"], lw=0.9, alpha=0.28, zorder=3)
        for m_idx, (recall, prec) in enumerate(tool["pts"]):
            ax.scatter(recall, prec, color=tool["color"], s=SIZES_C[m_idx], marker="o",
                       alpha=ALPHAS[m_idx], edgecolors="white", linewidths=1.3, zorder=5)

    for mode in nc_modes.values():
        rr = [p[0] for p in mode["pts"]]
        pp = [p[1] for p in mode["pts"]]
        ax.plot(rr, pp, "-", color=mode["color"], lw=1.0, alpha=0.30, zorder=3)
        for m_idx, (recall, prec) in enumerate(mode["pts"]):
            ax.scatter(recall, prec, color=mode["color"], s=SIZES_NC[m_idx], marker="D",
                       alpha=ALPHAS[m_idx], edgecolors="#222222", linewidths=1.4, zorder=6)

    ax.set_xlabel("Recall", fontsize=14, fontweight="bold", labelpad=10)
    ax.set_ylabel("Precision", fontsize=14, fontweight="bold", labelpad=10)
    ax.set_title(f"Performance Metrics ({species.capitalize()}, {len(paths.runs)} runs, {n_reads:,} reads)",
                 fontsize=16, fontweight="bold", pad=16)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x*100)}%"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{int(y*100)}%"))
    ax.set_xlim(-0.01, r_max)
    ax.set_ylim(-0.02, 1.06)

    handles = []

    def section(label):
        return mpatches.Patch(visible=False, label=label)

    handles.append(section("Competing tools"))
    for name, tool in competitors.items():
        handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=tool["color"],
                              markeredgecolor="white", markersize=11, label=name))

    handles.append(section(""))
    handles.append(section("NanoCirc"))
    for name, mode in nc_modes.items():
        handles.append(Line2D([0], [0], marker="D", color="w", markerfacecolor=mode["color"],
                              markeredgecolor="#333333", markeredgewidth=0.8, markersize=11, label=f"NanoCirc: {name}"))

    handles.append(section(""))
    handles.append(section("Metric (dot size)"))
    for i, metric in enumerate(METRICS):
        sz = [11, 8, 6][i]
        handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#888888", markersize=sz, label=metric))

    handles.append(section(""))
    handles.append(section("NanoCirc frontier"))
    handles.append(Line2D([0], [0], linestyle="--", color="#444444", lw=2.0, label="Exon"))
    handles.append(Line2D([0], [0], linestyle="-.", color="#444444", lw=1.5, label="Relaxed BSJ"))
    handles.append(Line2D([0], [0], linestyle=":", color="#444444", lw=1.2, label="Strict BSJ"))

    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), framealpha=0.97,
              edgecolor="#dddddd", fontsize=9.5, borderpad=0.9, labelspacing=0.45, handlelength=1.8)

    fig.text(0.5, -0.02, "GT filtered to circRNAs with >=1 confirmed-observed read in this run's data.",
              ha="center", va="top", fontsize=8, color="#888888")

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
