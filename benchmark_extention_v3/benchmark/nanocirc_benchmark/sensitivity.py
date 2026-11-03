"""Structural discovery sensitivity, split by ground-truth expression level
(tertiles). Answers "does this quantifier even find the circRNA at all",
independent of how well it quantifies it once found. Two views, matching
the original benchmark paper's own combined-plot convention (one plot,
several strategies grouped per caller) rather than one plot per strategy:

- structural: strict_bsj / relaxed_bsj / exon_based, one plot.
- isoform: full / strict / relaxed (config.isoform_match_levels), one plot.

strict_bsj/relaxed_bsj/exon_based use this module's own dedicated
sensitivity thresholds (config.sensitivity_bsj_tolerance/
sensitivity_min_overlap). full/strict/relaxed use their own fixed tier
definitions (config.isoform_match_levels), the same ones
performance.py/expression.py/reproducibility.py/overlap.py use.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bed_io import read_bed
from .matching import ISOFORM_MATCH_TAGS, StructuralMatcher, isoform_match_level_params

LEVELS = ["Low", "Medium", "High"]
LEVEL_COLORS = {"Low": "#AED6F1", "Medium": "#2ECC71", "High": "#E74C3C"}

STRUCTURAL_STRATEGIES = ["strict_bsj", "relaxed_bsj", "exon_based"]
ISOFORM_STRATEGIES = list(ISOFORM_MATCH_TAGS)  # full, strict, relaxed
VIEWS = {"structural": STRUCTURAL_STRATEGIES, "isoform": ISOFORM_STRATEGIES}

# alpha/hatch per position in a 3-strategy group: tightest first (solid), loosest last
_STYLE_BY_POSITION = [
    {"alpha": 0.95, "hatch": None},
    {"alpha": 0.65, "hatch": "..."},
    {"alpha": 0.40, "hatch": "///"},
]


def _strategy_label(strategy, config):
    if strategy == "strict_bsj":
        return "exact coordinate match"
    if strategy == "relaxed_bsj":
        return f"relaxed BSJ, +/-{config.sensitivity_bsj_tolerance}bp"
    if strategy == "exon_based":
        return f"exon-based, >={config.sensitivity_min_overlap:.0%} overlap"
    tol, min_ov = isoform_match_level_params(config, strategy)
    bsj_desc = "exact BSJ" if tol == 0 else f"+/-{tol}bp BSJ"
    return f"{strategy} isoform match, {bsj_desc}, >={min_ov:.0%} overlap"


def _detected_ids(matcher, caller_df, gt_df, work, strategy, config):
    if strategy == "strict_bsj":
        _, _, _, ids = matcher.strict_bsj_counts(caller_df, gt_df, work, return_ids=True)
    elif strategy == "relaxed_bsj":
        _, _, _, ids = matcher.relaxed_bsj_counts(caller_df, gt_df, work, tolerance=config.sensitivity_bsj_tolerance,
                                                    return_ids=True)
    elif strategy == "exon_based":
        _, _, _, ids = matcher.exon_based_counts(caller_df, gt_df, work, min_overlap=config.sensitivity_min_overlap,
                                                   return_ids=True)
    elif strategy in ISOFORM_MATCH_TAGS:
        tol, min_ov = isoform_match_level_params(config, strategy)
        _, _, _, ids = matcher.isoform_level_counts(caller_df, gt_df, work, tol, min_ov, return_ids=True)
    else:
        raise ValueError(f"unknown strategy: {strategy}")
    return ids


def build_structural_sensitivity(paths, config, species: str, strategies: list) -> pd.DataFrame:
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file)

    rows = []
    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        gt_expr = dict(zip(gt_df["row_id"], gt_df["score"].astype(int)))

        for q in config.quant_keys:
            path = paths.caller_bed_path(run, q, config.tiers)
            if not path.exists():
                continue
            caller_df = read_bed(path)
            if len(caller_df) == 0:
                continue
            ids_by_strategy = {}
            for strategy in strategies:
                work = paths.bedtools_work_dir / "expr_sensitivity" / strategy / run
                work.mkdir(parents=True, exist_ok=True)
                ids_by_strategy[strategy] = _detected_ids(matcher, caller_df, gt_df, work, strategy, config)
            for rid, expr in gt_expr.items():
                for strategy in strategies:
                    rows.append({"run": run, "quant": q, "row_id": rid, "gt_expr": expr,
                                 "strategy": strategy, "detected": int(rid in ids_by_strategy[strategy])})
    return pd.DataFrame(rows)


def plot_sensitivity_by_expression(struct_sens: pd.DataFrame, config, species: str, out_path,
                                    strategies: list) -> None:
    fig, ax = plt.subplots(figsize=(28, 8))
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    quant_colors = {**config.tool_colors, **config.tier_colors}
    n_q = len(quant_keys)
    n_s = len(strategies)
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

    def per_run_mean_sd(q, strategy, lv):
        sub = struct_sens[(struct_sens["level"] == lv) & (struct_sens["quant"] == q) &
                           (struct_sens["strategy"] == strategy)]
        per_run = [sub[sub["run"] == run]["detected"].mean() for run in runs if len(sub[sub["run"] == run])]
        mean = np.mean(per_run) if per_run else 0
        sd = np.std(per_run, ddof=1) if len(per_run) > 1 else 0
        return mean, sd

    w = 0.85 / (n_q * n_s)
    for qi, q in enumerate(quant_keys):
        for si, strategy in enumerate(strategies):
            vals, errs = [], []
            for lv in LEVELS:
                mean, sd = per_run_mean_sd(q, strategy, lv)
                vals.append(mean)
                errs.append(sd)
            g = qi * n_s + si
            offset = (g - (n_q * n_s - 1) / 2) * w
            style = _STYLE_BY_POSITION[min(si, len(_STYLE_BY_POSITION) - 1)]
            ax.bar(x + offset, vals, w * 0.92, yerr=errs, capsize=2.5, error_kw=dict(lw=1.0, capthick=1.0),
                   color=quant_colors[q], edgecolor="black", linewidth=1.0,
                   alpha=style["alpha"], hatch=style["hatch"],
                   label=(quant_labels[q].replace("\n", " ") if si == 0 else None))
            for xi, (v, e) in enumerate(zip(vals, errs)):
                if v > 0.01:
                    ax.text(x[xi] + offset, v + e + 0.015, f"{v*100:.0f}", ha="center", va="bottom",
                            fontsize=10, fontweight="bold", rotation=0)

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
    ax.set_ylabel("Sensitivity (Recall)", fontsize=13, fontweight="bold")
    ax.set_ylim(0, 1.22)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v*100:.0f}%"))
    ax.set_title(f"Quantifier Sensitivity by Expression Level ({species.capitalize()})", fontsize=13.5,
                 fontweight="bold", y=1.09)
    caller_legend = ax.legend(fontsize=9, loc="upper left", ncol=4, framealpha=0.95)
    ax.add_artist(caller_legend)
    strategy_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#555555", edgecolor="black",
                       alpha=_STYLE_BY_POSITION[si]["alpha"], hatch=_STYLE_BY_POSITION[si]["hatch"],
                       label=strategy.replace("_", " "))
        for si, strategy in enumerate(strategies)
    ]
    ax.legend(handles=strategy_handles, fontsize=9, loc="upper right", framealpha=0.95)
    ax.grid(axis="y", alpha=0.35, linestyle="--")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run_sensitivity_analysis(paths, config, species: str, views: dict = None) -> pd.DataFrame:
    """Runs each view (default: structural + isoform, see module docstring),
    one combined plot per view. Returns the structural view's DataFrame."""
    views = views or VIEWS
    default_result = None
    struct_sens = None
    for view_name, strategies in views.items():
        labels = ", ".join(_strategy_label(s, config) for s in strategies)
        print(f"[{species}] {view_name} sensitivity by GT-expression tertile ({labels}) ...")
        struct_sens = build_structural_sensitivity(paths, config, species, strategies)

        suffix = "" if view_name == "structural" else f"_{view_name}"
        csv_path = paths.bench_dir / f"sensitivity_by_expression_{species}{suffix}.csv"
        struct_sens.to_csv(csv_path, index=False)
        print(f"-> {csv_path}")

        out_path = paths.plot_dir / f"sensitivity_by_expression_{species}{suffix}.png"
        plot_sensitivity_by_expression(struct_sens, config, species, out_path, strategies)
        print(f"-> {out_path}")

        if view_name == "structural":
            default_result = struct_sens
    return default_result if default_result is not None else struct_sens
