"""Read-count correlation between each quantifier and ground truth, for the
same strategies performance.py scores TP/FP/FN with (strict_bsj/relaxed_bsj/
exon_based/isoform_level, full/strict/relaxed). See matching.py's module
docstring for what each strategy matches on.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy import stats

from .bed_io import read_bed
from .callers import (load_nanocirc_isoform_bed_counts, load_nanocirc_locus_counts, load_tier_membership,
                       load_tool_isoform_bed_counts, load_tool_locus_counts)
from .gt import count_reads, load_gt_loci
from .matching import (ISOFORM_MATCH_TAGS, STRATEGY_TITLES, StructuralMatcher, collapse_near_duplicates,
                        isoform_match_level_params, match_pairs_by_strategy,
                        match_pairs_by_strategy_isoform_resolved, match_pred_to_pred, strategy_footnote)
from .viz import sig_stars

_BED12_STRATEGIES = ("exon_based", "isoform_level") + ISOFORM_MATCH_TAGS


def build_long_df(paths, config, species: str, strategy: str = "relaxed_bsj") -> pd.DataFrame:
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
            work = paths.bedtools_work_dir / "expression" / strategy / run
            work.mkdir(parents=True, exist_ok=True)

        gt_iso = {r.row_id: (f"{r.chrom}:{r.start}-{r.end}:{r.strand}", int(r.score))
                  for r in gt_bed_df.itertuples()} if needs_gt_bed else None

        for tool in config.tools:
            if strategy in _BED12_STRATEGIES:
                call_df = None
                call_bed_df = load_tool_isoform_bed_counts(paths, run, tool)
            else:
                call_df = load_tool_locus_counts(paths, run, tool, config.tiers)
                call_bed_df = None
            if is_isoform_resolved:
                matched = match_pairs_by_strategy_isoform_resolved(strategy, call_bed_df, gt_bed_df, config, matcher, work)
                for g_row_id, (bid, gt_expr) in gt_iso.items():
                    pred = sum(n for n, _ in matched.get(g_row_id, []))
                    rows.append({"bsj_id": bid, "run": run, "quant": tool, "gt_expr": gt_expr, "pred_expr": pred})
            else:
                matched = match_pairs_by_strategy(strategy, call_df, gt_df, config, matcher, call_bed_df, gt_bed_df, work)
                for bid, gt_expr in gt_map.items():
                    pred = sum(n for n, _ in matched.get(bid, []))
                    rows.append({"bsj_id": bid, "run": run, "quant": tool, "gt_expr": gt_expr, "pred_expr": pred})

        nc_df = load_nanocirc_locus_counts(paths, run)
        nc_bed_df = load_nanocirc_isoform_bed_counts(paths, run, None, work) if strategy in _BED12_STRATEGIES else None
        tier_members = load_tier_membership(paths, config.tiers)
        if is_isoform_resolved:
            nc_matched = match_pairs_by_strategy_isoform_resolved(strategy, nc_bed_df, gt_bed_df, config, matcher,
                                                                    work, extra_col="base_bsj_id")
            for tier in config.tiers:
                members = tier_members[tier]
                for g_row_id, (bid, gt_expr) in gt_iso.items():
                    candidates = [(n, src_bid) for n, src_bid in nc_matched.get(g_row_id, []) if src_bid in members]
                    pred = collapse_near_duplicates(candidates)
                    rows.append({"bsj_id": bid, "run": run, "quant": tier, "gt_expr": gt_expr, "pred_expr": pred})
        else:
            nc_matched = match_pairs_by_strategy(strategy, nc_df, gt_df, config, matcher, nc_bed_df, gt_bed_df,
                                                  work, extra_col="base_bsj_id")
            for tier in config.tiers:
                members = tier_members[tier]
                for bid, gt_expr in gt_map.items():
                    candidates = [(n, src_bid) for n, src_bid in nc_matched.get(bid, []) if src_bid in members]
                    pred = collapse_near_duplicates(candidates)
                    rows.append({"bsj_id": bid, "run": run, "quant": tier, "gt_expr": gt_expr, "pred_expr": pred})
    return pd.DataFrame(rows)


def build_tool_pair_df(paths, config, species: str, strategy: str = "relaxed_bsj") -> pd.DataFrame:
    """One row per (nanocirc-catalog-locus, run) with its own count and each
    tool's matched count (for this strategy), for the NanoCirc-vs-tool
    scatter panels."""
    is_bed12 = strategy in _BED12_STRATEGIES
    tol = 0 if strategy == "strict_bsj" else config.gt_match_tolerance
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file) if is_bed12 else None
    rows = []
    for run in paths.runs:
        if is_bed12:
            work = paths.bedtools_work_dir / "expression" / strategy / run
            work.mkdir(parents=True, exist_ok=True)
            nc_bed_df = load_nanocirc_isoform_bed_counts(paths, run, None, work)
            row_data = {i: {"bsj_id": r.name, "run": run, "nanocirc_count": r.n_reads}
                        for i, r in nc_bed_df.iterrows()}
        else:
            nc_df = load_nanocirc_locus_counts(paths, run).reset_index(drop=True)
            row_data = {i: {"bsj_id": f"{r.chrom}:{r.start}-{r.end}:{r.strand}", "run": run,
                            "nanocirc_count": r.n_reads} for i, r in nc_df.iterrows()}

        for tool in config.tools:
            agg = {}
            if is_bed12:
                call_bed_df = load_tool_isoform_bed_counts(paths, run, tool)
                if strategy == "exon_based":
                    pairs = matcher.exon_based_pairs(call_bed_df, nc_bed_df, work, min_overlap=config.tp_min_overlap)
                elif strategy == "isoform_level":
                    pairs = matcher.isoform_level_pairs(call_bed_df, nc_bed_df, work, tolerance=config.gt_match_tolerance,
                                                         min_overlap=config.isoform_level_min_overlap)
                else:
                    tier_tol, tier_min_ov = isoform_match_level_params(config, strategy)
                    pairs = matcher.isoform_level_pairs(call_bed_df, nc_bed_df, work, tolerance=tier_tol,
                                                         min_overlap=tier_min_ov)
                call_bed_by_id = {r.row_id: r for r in call_bed_df.itertuples()}
                for c_row_id, nc_row_id in pairs:
                    if nc_row_id not in row_data:
                        continue
                    agg[nc_row_id] = agg.get(nc_row_id, 0) + call_bed_by_id[c_row_id].n_reads
            else:
                call_df = load_tool_locus_counts(paths, run, tool, config.tiers)
                pairs = match_pred_to_pred(call_df, nc_df, tol)
                for idx, cnt in pairs:
                    agg[idx] = agg.get(idx, 0) + cnt
            for i in row_data:
                row_data[i][f"{tool}_count"] = agg.get(i, 0)
        rows.extend(row_data.values())
    return pd.DataFrame(rows)


def compute_corr(merged: pd.DataFrame, config) -> pd.DataFrame:
    """Per-run Pearson/Spearman, then averaged across runs (mean, with sd
    for spread), rather than pooling all runs into one sample. Pooling
    inflates n and can mask a quant that is inconsistent run-to-run."""
    rows = []
    for q in config.quant_keys:
        td = merged[merged["quant"] == q]
        n_detected_total = len(td[(td["gt_expr"] > 0) & (td["pred_expr"] > 0)])
        per_run = {"pearson_r": [], "pearson_p": [], "spearman_rho": [], "spearman_p": [], "r_squared": []}
        for run in sorted(td["run"].unique()):
            nz = td[(td["run"] == run) & (td["gt_expr"] > 0) & (td["pred_expr"] > 0)]
            if len(nz) > 10:
                r, p_r = stats.pearsonr(nz["gt_expr"], nz["pred_expr"])
                rho, p_rho = stats.spearmanr(nz["gt_expr"], nz["pred_expr"])
                per_run["pearson_r"].append(r); per_run["pearson_p"].append(p_r)
                per_run["spearman_rho"].append(rho); per_run["spearman_p"].append(p_rho)
                per_run["r_squared"].append(r ** 2)
        if per_run["pearson_r"]:
            rows.append({"quant": q, "n_detected": n_detected_total,
                         "pearson_r": np.mean(per_run["pearson_r"]), "pearson_r_sd": np.std(per_run["pearson_r"], ddof=1) if len(per_run["pearson_r"]) > 1 else 0,
                         "pearson_p": np.mean(per_run["pearson_p"]),
                         "spearman_rho": np.mean(per_run["spearman_rho"]), "spearman_rho_sd": np.std(per_run["spearman_rho"], ddof=1) if len(per_run["spearman_rho"]) > 1 else 0,
                         "spearman_p": np.mean(per_run["spearman_p"]),
                         "r_squared": np.mean(per_run["r_squared"])})
        else:
            rows.append({"quant": q, "n_detected": n_detected_total, "pearson_r": np.nan, "pearson_r_sd": np.nan, "pearson_p": np.nan,
                         "spearman_rho": np.nan, "spearman_rho_sd": np.nan, "spearman_p": np.nan, "r_squared": np.nan})
    return pd.DataFrame(rows)


def _draw_nanocirc_vs_tool_panel(ax, df, tool, label, cmap, rng):
    col = f"{tool}_count"
    sub = df[["bsj_id", "run", "nanocirc_count", col]].copy()
    sub = sub[(sub[col] > 0) & (sub["nanocirc_count"] > 0)]

    x = np.log1p(sub["nanocirc_count"].values)
    y = np.log1p(sub[col].values)
    max_val = max(sub["nanocirc_count"].max(), sub[col].max()) if len(sub) else 1
    candidates = sorted(set([0, 1, 2, 5, 10, int(max_val)]))
    tick_vals = [candidates[0]]
    for t in candidates[1:]:
        if np.log1p(t) - np.log1p(tick_vals[-1]) > 0.25:
            tick_vals.append(t)
        elif t == candidates[-1]:
            tick_vals[-1] = t
    tick_pos = [np.log1p(t) for t in tick_vals]
    axis_limit = np.log1p(max_val * 1.15) if max_val else 1

    xj = x + rng.normal(0, 0.03, size=len(x))
    yj = y + rng.normal(0, 0.03, size=len(y))

    if len(sub) > 3:
        xy = np.vstack([xj, yj])
        try:
            density = stats.gaussian_kde(xy)(xy)
        except np.linalg.LinAlgError:
            density = np.ones(len(xj))
    else:
        density = np.ones(len(xj))
    order = np.argsort(density)
    xj, yj, density = xj[order], yj[order], density[order]

    ax.scatter(xj, yj, c=density, cmap=cmap, s=10, alpha=0.7, linewidths=0,
               norm=Normalize(vmin=density.min(), vmax=density.max()) if len(density) else None)
    ax.plot([0, axis_limit], [0, axis_limit], color="#bbbbbb", lw=1.0, linestyle="--", zorder=0)

    r, p = (stats.pearsonr(sub["nanocirc_count"], sub[col]) if len(sub) > 10 else (np.nan, np.nan))
    if not np.isnan(r):
        ax.text(0.05, 0.95, f"r={r:.3f} ({sig_stars(p)})\nn={len(sub):,}",
                transform=ax.transAxes, va="top", fontsize=9, fontweight="bold",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="black", pad=0.3))

    ax.set_xlim(0, axis_limit); ax.set_ylim(0, axis_limit)
    ax.set_xticks(tick_pos); ax.set_xticklabels([f"{t:,}" for t in tick_vals], fontsize=8)
    ax.set_yticks(tick_pos); ax.set_yticklabels([f"{t:,}" for t in tick_vals], fontsize=8)
    ax.set_xlabel("NanoCirc", fontsize=9, fontweight="bold")
    ax.set_ylabel(label.replace("\n", " "), fontsize=9, fontweight="bold")
    ax.grid(alpha=0.25, linewidth=0.5)


def plot_correlation_summary(tool_pair_df, corr_df, config, species: str, n_reads: int, n_runs: int, out_path,
                              strategy: str = "relaxed_bsj") -> None:
    fig = plt.figure(figsize=(18, 7.5))
    gs = fig.add_gridspec(1, 2, wspace=0.55, top=0.85, width_ratios=[1, 1.15])
    fig.suptitle(f"Expression agreement between tools and modes "
                 f"({species.capitalize()}, {STRATEGY_TITLES[strategy]} match, {n_runs} runs, {n_reads:,} reads)",
                 fontsize=16, fontweight="bold", y=0.98)

    n_tools = len(config.tools)
    ncols = 2
    nrows = -(-n_tools // ncols)
    inner_gs = gs[0, 0].subgridspec(nrows, ncols, wspace=0.55, hspace=0.75)
    cmap = plt.get_cmap("plasma")
    scatter_rng = np.random.default_rng(0)
    for ti, tool in enumerate(config.tools):
        ax_t = fig.add_subplot(inner_gs[ti // ncols, ti % ncols])
        _draw_nanocirc_vs_tool_panel(ax_t, tool_pair_df, tool, config.tool_labels[tool], cmap, scatter_rng)
    top_left_x = gs[0, 0].get_position(fig).x0 + gs[0, 0].get_position(fig).width / 2
    top_left_y = gs[0, 0].get_position(fig).y1
    fig.text(top_left_x, top_left_y + 0.06, "NanoCirc vs each tool", ha="center", va="bottom", fontsize=13, fontweight="bold")
    fig.text(top_left_x, top_left_y + 0.02, f"({strategy_footnote(strategy, config)})",
              ha="center", va="bottom", fontsize=9, color="#888888")

    ax2 = fig.add_subplot(gs[0, 1])
    metrics = ["pearson_r", "spearman_rho", "r_squared"]
    metric_labels = ["Pearson r", "Spearman rho", "R2"]
    p_cols = {"pearson_r": "pearson_p", "spearman_rho": "spearman_p", "r_squared": "pearson_p"}
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    cdf2 = corr_df.set_index("quant").loc[quant_keys]
    hm = cdf2[metrics].values
    im = ax2.imshow(hm, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
    ax2.set_xticks(range(len(metrics))); ax2.set_xticklabels(metric_labels, fontsize=12, fontweight="bold")
    ax2.set_yticks(range(len(quant_keys)))
    ax2.set_yticklabels([f"{quant_labels[q]}\nn={cdf2.loc[q, 'n_detected']:,.0f}" for q in quant_keys],
                         fontsize=12, fontweight="bold")
    sd_cols = {"pearson_r": "pearson_r_sd", "spearman_rho": "spearman_rho_sd", "r_squared": None}
    for i, q in enumerate(quant_keys):
        for j, m in enumerate(metrics):
            v = hm[i, j]
            if np.isnan(v):
                txt = "n/a"
            else:
                sd_col = sd_cols[m]
                sd_txt = f" +/-{cdf2.loc[q, sd_col]:.3f}" if sd_col else ""
                txt = f"{v:.3f}{sd_txt}\n({sig_stars(cdf2.loc[q, p_cols[m]])})"
            ax2.text(j, i, txt, ha="center", va="center", fontsize=10, fontweight="bold")
    ax2.set_title("Correlation vs Ground truth (mean across runs)", fontsize=13, fontweight="bold")
    plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    ax2.text(0.5, -0.10, "*** p<0.001, ** p<0.01, * p<0.05, ns = not significant. Mean +/- sd across runs.",
             transform=ax2.transAxes, ha="center", va="top", fontsize=8, color="#888888")

    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run_expression_analysis(paths, config, species: str, strategies=None) -> pd.DataFrame:
    """Runs the expression correlation analysis for relaxed_bsj, strict_bsj,
    exon_based, full, strict, relaxed (default: all 6). Returns the
    relaxed_bsj correlation table."""
    strategies = strategies or ["relaxed_bsj", "strict_bsj", "exon_based", "full", "strict", "relaxed"]
    n_reads = count_reads(paths.run_fastq(paths.runs[0]))
    default_corr_df = None

    for strategy in strategies:
        suffix = "" if strategy == "relaxed_bsj" else f"_{strategy}"
        print(f"[{species}] {strategy} ({strategy_footnote(strategy, config)}) per-locus counts ...")
        tool_pair_df = build_tool_pair_df(paths, config, species, strategy)
        merged = build_long_df(paths, config, species, strategy)
        corr_df = compute_corr(merged, config)
        corr_df["strategy"] = strategy
        print(corr_df.to_string(index=False))

        corr_csv = paths.bench_dir / f"expression_correlation_{species}{suffix}.csv"
        corr_df.to_csv(corr_csv, index=False)
        print(f"-> {corr_csv}")

        out_path = paths.plot_dir / f"expression_correlation_{species}{suffix}.png"
        plot_correlation_summary(tool_pair_df, corr_df, config, species, n_reads, len(paths.runs), out_path, strategy)
        print(f"-> {out_path}")

        if strategy == "relaxed_bsj":
            default_corr_df = corr_df

    return default_corr_df if default_corr_df is not None else corr_df
