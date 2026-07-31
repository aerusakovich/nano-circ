"""Structural discovery sensitivity, split by ground-truth expression level
(tertiles). Answers "does this quantifier even find the circRNA at all",
independent of how well it quantifies it once found.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bed_io import read_bed
from .matching import StructuralMatcher

LEVELS = ["Low", "Medium", "High"]
LEVEL_COLORS = {"Low": "#AED6F1", "Medium": "#2ECC71", "High": "#E74C3C"}


def build_structural_sensitivity(paths, config, species: str) -> pd.DataFrame:
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file)
    tolerance = config.sensitivity_bsj_tolerance
    min_overlap = config.sensitivity_min_overlap

    rows = []
    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        gt_expr = dict(zip(gt_df["row_id"], gt_df["score"].astype(int)))
        work = paths.bedtools_work_dir / "expr_sensitivity" / run
        work.mkdir(parents=True, exist_ok=True)

        for q in config.quant_keys:
            path = paths.caller_bed_path(run, q, config.tiers)
            if not path.exists():
                continue
            caller_df = read_bed(path)
            if len(caller_df) == 0:
                continue
            _, _, _, relaxed_ids = matcher.relaxed_bsj_counts(caller_df, gt_df, work, tolerance=tolerance, return_ids=True)
            _, _, _, exon_ids = matcher.exon_based_counts(caller_df, gt_df, work, min_overlap=min_overlap, return_ids=True)
            for rid, expr in gt_expr.items():
                rows.append({"run": run, "quant": q, "row_id": rid, "gt_expr": expr,
                             "relaxed_detected": int(rid in relaxed_ids), "exon_detected": int(rid in exon_ids)})
    return pd.DataFrame(rows)


def plot_sensitivity_by_expression(struct_sens: pd.DataFrame, config, species: str, out_path) -> None:
    fig, ax = plt.subplots(figsize=(16, 7))
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    quant_colors = {**config.tool_colors, **config.tier_colors}
    n_q = len(quant_keys)
    x = np.arange(len(LEVELS))

    p33, p67 = struct_sens["gt_expr"].quantile([0.333, 0.667])

    def label_tertile(v):
        if v <= p33:
            return "Low"
        if v <= p67:
            return "Medium"
        return "High"

    struct_sens = struct_sens.copy()
    struct_sens["level"] = struct_sens["gt_expr"].apply(label_tertile)
    runs = sorted(struct_sens["run"].unique())

    def per_run_mean_sd(q, strategy_col, lv):
        per_run = [struct_sens[(struct_sens["level"] == lv) & (struct_sens["run"] == run) & (struct_sens["quant"] == q)][strategy_col].mean()
                   for run in runs if len(struct_sens[(struct_sens["level"] == lv) & (struct_sens["run"] == run) & (struct_sens["quant"] == q)])]
        mean = np.mean(per_run) if per_run else 0
        sd = np.std(per_run, ddof=1) if len(per_run) > 1 else 0
        return mean, sd

    w = 0.85 / (n_q * 2)
    for qi, q in enumerate(quant_keys):
        for si, strategy_col in enumerate(["relaxed_detected", "exon_detected"]):
            vals, errs = [], []
            for lv in LEVELS:
                mean, sd = per_run_mean_sd(q, strategy_col, lv)
                vals.append(mean)
                errs.append(sd)
            offset = (qi * 2 + si - n_q + 0.5) * w
            is_exon = strategy_col == "exon_detected"
            ax.bar(x + offset, vals, w * 0.92, yerr=errs, capsize=2.5, error_kw=dict(lw=1.0, capthick=1.0),
                   color=quant_colors[q], edgecolor="black", linewidth=1.0,
                   alpha=0.45 if is_exon else 0.9, hatch="///" if is_exon else None,
                   label=(quant_labels[q].replace("\n", " ") if not is_exon else None))
            for xi, (v, e) in enumerate(zip(vals, errs)):
                if v > 0.01:
                    ax.text(x[xi] + offset, v + e + 0.012, f"{v*100:.0f}", ha="center", va="bottom",
                            fontsize=6, fontweight="bold", rotation=90)

    tick_labels = []
    for xi, lv in enumerate(LEVELS):
        ld = struct_sens[struct_sens["level"] == lv]
        lo, hi = ld["gt_expr"].min(), ld["gt_expr"].max()
        rng = f"{lo:.0f} read" if lo == hi else f"{lo:.0f}-{hi:.0f} reads"
        n_loci = ld["row_id"].nunique()
        tick_labels.append(f"{lv}\n({rng}, n={n_loci})")
        ax.axvspan(xi - 0.48, xi + 0.48, alpha=0.06, color=LEVEL_COLORS[lv], zorder=0)

    ax.set_xticks(x); ax.set_xticklabels(tick_labels, fontsize=12, fontweight="bold")
    ax.set_xlabel("Ground Truth Expression Level (Tertiles)", fontsize=13, fontweight="bold", labelpad=12)
    ax.set_ylabel("Sensitivity (Structural Recall)", fontsize=13, fontweight="bold")
    ax.set_ylim(0, 1.15)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v*100:.0f}%"))
    ax.set_title(f"Quantifier Sensitivity by Expression Level ({species.capitalize()})", fontsize=13.5, fontweight="bold", y=1.09)
    caller_legend = ax.legend(fontsize=9, loc="upper left", ncol=4, framealpha=0.95)
    ax.add_artist(caller_legend)
    strategy_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#555555", edgecolor="black", alpha=0.9,
                       label=f"Relaxed BSJ (+/-{int(struct_sens.attrs.get('tolerance', 10))}bp)"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#555555", edgecolor="black", alpha=0.45, hatch="///",
                       label="Exon-based"),
    ]
    ax.legend(handles=strategy_handles, fontsize=9, loc="upper right", framealpha=0.95)
    ax.grid(axis="y", alpha=0.35, linestyle="--")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run_sensitivity_analysis(paths, config, species: str) -> pd.DataFrame:
    print(f"[{species}] structural sensitivity by GT-expression tertile "
          f"(relaxed_bsj +/-{config.sensitivity_bsj_tolerance}bp, exon >={config.sensitivity_min_overlap:.0%}) ...")
    struct_sens = build_structural_sensitivity(paths, config, species)
    struct_sens.attrs["tolerance"] = config.sensitivity_bsj_tolerance

    csv_path = paths.bench_dir / f"sensitivity_by_expression_{species}.csv"
    struct_sens.to_csv(csv_path, index=False)
    print(f"-> {csv_path}")

    out_path = paths.plot_dir / f"sensitivity_by_expression_{species}.png"
    plot_sensitivity_by_expression(struct_sens, config, species, out_path)
    print(f"-> {out_path}")
    return struct_sens
