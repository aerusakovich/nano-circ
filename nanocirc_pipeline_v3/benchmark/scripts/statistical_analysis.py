#!/usr/bin/env python3
"""Statistically compares two nanocirc benchmark datasets classified with
the same type classifier (synthetic vs database, protocol vs protocol,
species vs species, etc), reading each side's
classification_results_<strategy>.csv.

- aggregate/per_type/per_tool_regression: linear fit dataset2 = a + b*dataset1
  per metric, plus paired Wilcoxon and Cohen's d.
- tool_ranking: Kendall's tau-b between the two sides' caller rankings
  (macro-averaged over type and run).
- mantel_haenszel: per-caller Mantel-Haenszel odds ratio (dataset2 vs
  dataset1) on raw TP/FP/FN counts stratified by circRNA type, plus
  Breslow-Day for homogeneity across types.

--dataset1-root/--dataset2-root point at the directory directly containing
data/classification/<group>/; <group> is set with --group or the per-side
--dataset1-group/--dataset2-group.

Usage:
    python3 statistical_analysis.py \\
        --dataset1-root ../../publication/species_benchmark \\
        --dataset2-root ../../publication/database_benchmark \\
        --group human --strategy isoform_level \\
        --config ../configs/default.yaml \\
        --out-dir ../../publication/sim_vs_database_human
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nanocirc_benchmark.config import BenchmarkConfig  # noqa: E402
from nanocirc_benchmark.type_classification import ALL_TYPES  # noqa: E402


def load_per_run_metrics(root: str, group: str, strategy: str, callers: list) -> pd.DataFrame:
    """One row per (run_idx, caller, type): tp, fp, fn, gt_n, precision, recall, f1."""
    path = Path(root) / "data" / "classification" / group / f"classification_results_{strategy}.csv"
    df = pd.read_csv(path)
    run_order = sorted(df["run"].unique())
    run_idx = {r: i for i, r in enumerate(run_order)}
    gt = df[df["caller"] == "ground_truth"].set_index("run")

    rows = []
    for run in run_order:
        gt_totals = {t: gt.loc[run, f"count_{t}"] for t in ALL_TYPES}
        for caller in callers:
            sub = df[(df["caller"] == caller) & (df["run"] == run)]
            if sub.empty:
                continue
            for t in ALL_TYPES:
                gtn = gt_totals[t]
                if gtn == 0:
                    continue
                tp = int(sub[f"tp_{t}"].sum())
                fp = int(sub[f"fp_{t}"].sum())
                fn = gtn - tp
                # precision is undefined (not 0) when the caller made zero calls of this
                # type; leaving it as NaN keeps it out of precision averages/regressions
                # instead of falsely counting "made no calls" as "always wrong".
                precision = tp / (tp + fp) if (tp + fp) else float("nan")
                recall = tp / gtn
                f1 = 2 * precision * recall / (precision + recall) if tp else 0.0
                rows.append({"run_idx": run_idx[run], "caller": caller, "type": t, "gt_n": gtn,
                             "tp": tp, "fp": fp, "fn": fn,
                             "precision": precision, "recall": recall, "f1": f1})
    return pd.DataFrame(rows)


def paired_regression_stats(x: np.ndarray, y: np.ndarray) -> dict:
    """dataset2 (y) = intercept + slope*dataset1 (x), plus a paired gap test."""
    diff = y - x
    slope, intercept, r, p_slope, _ = stats.linregress(x, y)
    try:
        _, p_w = stats.wilcoxon(diff)
    except ValueError:
        p_w = float("nan")
    sd = np.std(diff, ddof=1)
    cohens_d = np.mean(diff) / sd if sd > 0 else float("nan")
    return {
        "n": len(x), "mean_dataset1": np.mean(x), "mean_dataset2": np.mean(y),
        "mean_diff": np.mean(diff), "slope": slope, "intercept": intercept,
        "r2": r ** 2, "slope_p": p_slope, "wilcoxon_p": p_w, "cohens_d": cohens_d,
    }


def _xy(sub: pd.DataFrame, metric: str) -> tuple[np.ndarray, np.ndarray]:
    """Paired x,y for a metric, dropping rows where either side is NaN
    (caller made zero calls of that type in that run, so precision is undefined)."""
    cols = sub[[f"{metric}_1", f"{metric}_2"]].dropna()
    return cols[f"{metric}_1"].values, cols[f"{metric}_2"].values


def aggregate_regression(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in ("precision", "recall", "f1"):
        x, y = _xy(merged, metric)
        rows.append({"metric": metric, **paired_regression_stats(x, y)})
    return pd.DataFrame(rows)


def per_type_regression(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for t in ALL_TYPES:
        sub = merged[merged["type"] == t]
        if sub.empty:
            continue
        for metric in ("precision", "recall", "f1"):
            x, y = _xy(sub, metric)
            if len(x) == 0:
                continue
            rows.append({"type": t, "metric": metric, **paired_regression_stats(x, y)})
    return pd.DataFrame(rows)


def per_tool_regression(merged: pd.DataFrame, callers: list) -> pd.DataFrame:
    rows = []
    for c in callers:
        sub = merged[merged["caller"] == c]
        if sub.empty:
            continue
        for metric in ("precision", "recall", "f1"):
            x, y = _xy(sub, metric)
            if len(x) == 0:
                continue
            rows.append({"caller": c, "metric": metric, **paired_regression_stats(x, y)})
    return pd.DataFrame(rows)


def tool_macro_averages(merged: pd.DataFrame, callers: list) -> pd.DataFrame:
    """One row per caller: mean over type and run (macro-average, each type weighted equally)."""
    return (merged.groupby("caller")[[f"{m}_{s}" for m in ("precision", "recall", "f1") for s in ("1", "2")]]
            .mean().reindex(callers).reset_index())


def tool_ranking_table(merged: pd.DataFrame, callers: list) -> pd.DataFrame:
    """Kendall's tau-b between the two rank vectors of caller macro-averages.
    A caller with no data on either side at all (e.g. CIRI-long on isoCirc-protocol
    reads, which makes zero calls) is dropped from that metric's ranking rather than
    propagating NaN into the whole row."""
    macro = tool_macro_averages(merged, callers)
    rows = []
    for metric in ("precision", "recall", "f1"):
        sub = macro[[f"{metric}_1", f"{metric}_2"]].dropna()
        x, y = sub[f"{metric}_1"].values, sub[f"{metric}_2"].values
        tau, p = stats.kendalltau(x, y)
        rho, p_rho = stats.spearmanr(x, y)
        rows.append({"metric": metric, "n_callers": len(x),
                      "kendall_tau_b": tau, "kendall_p": p,
                      "spearman_rho": rho, "spearman_p": p_rho})
    return pd.DataFrame(rows)


def tool_rank_detail(merged: pd.DataFrame, callers: list) -> pd.DataFrame:
    """Per-caller macro-average F1 and its rank (1 = best) on each side, for inspection."""
    macro = tool_macro_averages(merged, callers)
    macro["rank_dataset1"] = macro["f1_1"].rank(ascending=False)
    macro["rank_dataset2"] = macro["f1_2"].rank(ascending=False)
    return macro[["caller", "f1_1", "f1_2", "rank_dataset1", "rank_dataset2"]]


def _mh_pooled_or(tables: list[tuple[float, float, float, float]]) -> tuple[float, float, float, float]:
    """Mantel-Haenszel pooled OR + RBG variance-based CI + MH chi-square test.
    Each table is (a, b, c, d) = (exposed+outcome+, exposed+outcome-, unexposed+outcome+, unexposed+outcome-)."""
    R = S = P_R = P_S_plus_Q_R = Q_S = 0.0
    sum_a = sum_ea = sum_var = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        r_k = a * d / n
        s_k = b * c / n
        p_k = (a + d) / n
        q_k = (b + c) / n
        R += r_k
        S += s_k
        P_R += p_k * r_k
        P_S_plus_Q_R += p_k * s_k + q_k * r_k
        Q_S += q_k * s_k
        n1 = a + b
        n2 = c + d
        m1 = a + c
        sum_a += a
        sum_ea += n1 * m1 / n
        sum_var += n1 * n2 * m1 * (b + d) / (n ** 2 * (n - 1)) if n > 1 else 0.0

    if S == 0 or R == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    or_mh = R / S
    var_log_or = (P_R / (2 * R ** 2)) + (P_S_plus_Q_R / (2 * R * S)) + (Q_S / (2 * S ** 2))
    se_log_or = np.sqrt(var_log_or) if var_log_or > 0 else float("nan")
    ci_lo = np.exp(np.log(or_mh) - 1.96 * se_log_or) if or_mh > 0 else float("nan")
    ci_hi = np.exp(np.log(or_mh) + 1.96 * se_log_or) if or_mh > 0 else float("nan")
    chi2 = (abs(sum_a - sum_ea) - 0.5) ** 2 / sum_var if sum_var > 0 else float("nan")
    p_mh = 1 - stats.chi2.cdf(chi2, df=1) if not np.isnan(chi2) else float("nan")
    return or_mh, ci_lo, ci_hi, p_mh


def _breslow_day(tables: list[tuple[float, float, float, float]], or_common: float) -> tuple[float, float]:
    """Breslow-Day test for homogeneity of the odds ratio across strata."""
    if np.isnan(or_common) or or_common <= 0:
        return float("nan"), float("nan")
    stat = 0.0
    k_used = 0
    for a, b, c, d in tables:
        n1 = a + b
        n2 = c + d
        m1 = a + c
        m2 = b + d
        n = n1 + n2
        if n == 0 or n1 == 0 or n2 == 0 or m1 == 0 or m2 == 0:
            continue
        if abs(or_common - 1) < 1e-12:
            e_a = n1 * m1 / n
        else:
            coef_a = 1 - or_common
            coef_b = n2 - m1 + or_common * (n1 + m1)
            coef_c = -or_common * n1 * m1
            disc = coef_b ** 2 - 4 * coef_a * coef_c
            if disc < 0:
                continue
            root1 = (-coef_b + np.sqrt(disc)) / (2 * coef_a)
            root2 = (-coef_b - np.sqrt(disc)) / (2 * coef_a)
            lo = max(0.0, n1 - m2)
            hi = min(n1, m1)
            candidates = [r for r in (root1, root2) if lo - 1e-6 <= r <= hi + 1e-6]
            e_a = candidates[0] if candidates else n1 * m1 / n
        denom = (1 / e_a + 1 / (n1 - e_a) + 1 / (m1 - e_a) + 1 / (n2 - m1 + e_a)) if e_a not in (0, n1, m1) and (n2 - m1 + e_a) != 0 else None
        if not denom:
            continue
        var_a = 1 / denom
        if var_a <= 0:
            continue
        stat += (a - e_a) ** 2 / var_a
        k_used += 1
    if k_used < 2:
        return float("nan"), float("nan")
    p = 1 - stats.chi2.cdf(stat, df=k_used - 1)
    return stat, p


def mantel_haenszel_table(counts: pd.DataFrame, callers: list) -> pd.DataFrame:
    """Per caller, stratified (by type) MH odds ratio for the precision (TP vs FP)
    and recall (TP vs FN) outcomes: dataset2 vs dataset1."""
    rows = []
    for c in callers:
        sub = counts[counts["caller"] == c]
        if sub.empty:
            continue
        strata = sub.groupby(["type", "source"])[["tp", "fp", "fn"]].sum().reset_index()
        for outcome, neg_col in (("precision", "fp"), ("recall", "fn")):
            tables = []
            for t in ALL_TYPES:
                st = strata[strata["type"] == t]
                d2_row = st[st["source"] == "dataset2"]
                d1_row = st[st["source"] == "dataset1"]
                if d2_row.empty or d1_row.empty:
                    continue
                a = d2_row["tp"].values[0]
                b = d2_row[neg_col].values[0]
                c_ = d1_row["tp"].values[0]
                d = d1_row[neg_col].values[0]
                tables.append((a, b, c_, d))
            if len(tables) < 2:
                continue
            or_mh, ci_lo, ci_hi, p_mh = _mh_pooled_or(tables)
            bd_stat, bd_p = _breslow_day(tables, or_mh)
            rows.append({"caller": c, "outcome": outcome, "k_strata": len(tables),
                         "or_mh": or_mh, "ci_lo": ci_lo, "ci_hi": ci_hi, "mh_p": p_mh,
                         "breslow_day_stat": bd_stat, "breslow_day_p": bd_p})
    return pd.DataFrame(rows)


def plot_slopegraph(merged: pd.DataFrame, callers: list, config: BenchmarkConfig, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {**config.tool_colors, **config.tier_colors}
    macro = tool_macro_averages(merged, callers)
    fig, axes = plt.subplots(1, 3, figsize=(11, 5.5), sharey=False)
    handles, labels = None, None
    for ax, metric in zip(axes, ("precision", "recall", "f1")):
        for _, row in macro.iterrows():
            marker = "D" if row["caller"] in config.tier_colors else "o"
            ax.plot([0, 1], [row[f"{metric}_1"], row[f"{metric}_2"]], marker=marker,
                     color=colors.get(row["caller"], "#888888"), label=row["caller"])
        ax.set_xlim(-0.15, 1.15)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Dataset 1", "Dataset 2"])
        ax.set_title(metric.capitalize())
        ax.set_ylim(0, 1)
        if handles is None:
            handles, labels = ax.get_legend_handles_labels()
    fig.suptitle("Per-caller macro-average: dataset 1 vs. dataset 2")
    fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.05), fontsize=8)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset1-root", required=True)
    p.add_argument("--dataset2-root", required=True)
    p.add_argument("--group", default=None,
                    help="group subfolder name under data/classification/ for both sides; use "
                         "--dataset1-group/--dataset2-group instead if the two sides differ "
                         "(e.g. comparing wet-lab protocols or species against each other, not just "
                         "synthetic vs database)")
    p.add_argument("--dataset1-group", default=None, help="overrides --group for the dataset1-root side")
    p.add_argument("--dataset2-group", default=None, help="overrides --group for the dataset2-root side")
    p.add_argument("--strategy", default="isoform_level",
                    help="matches a classification_results_<strategy>.csv; default is strict isoform-level "
                         "matching (BSJ +/-10bp AND >=99%% exon overlap)")
    p.add_argument("--config", default=None, help="YAML config for caller list (config.quant_keys); default config if omitted")
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    group1 = args.dataset1_group or args.group
    group2 = args.dataset2_group or args.group
    if not group1 or not group2:
        p.error("need --group, or both --dataset1-group and --dataset2-group")

    config = BenchmarkConfig.from_yaml(args.config) if args.config else BenchmarkConfig.default()
    callers = config.quant_keys

    d1 = load_per_run_metrics(args.dataset1_root, group1, args.strategy, callers)
    d2 = load_per_run_metrics(args.dataset2_root, group2, args.strategy, callers)
    merged = d1.merge(d2, on=["run_idx", "caller", "type"], suffixes=("_1", "_2"))

    d1_c = d1.copy(); d1_c["source"] = "dataset1"
    d2_c = d2.copy(); d2_c["source"] = "dataset2"
    counts = pd.concat([d1_c, d2_c], ignore_index=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "aggregate_regression": aggregate_regression(merged),
        "per_type_regression": per_type_regression(merged),
        "per_tool_regression": per_tool_regression(merged, callers),
        "tool_ranking": tool_ranking_table(merged, callers),
        "tool_rank_detail": tool_rank_detail(merged, callers),
        "mantel_haenszel": mantel_haenszel_table(counts, callers),
    }
    for name, table in tables.items():
        out_path = out_dir / f"{name}.csv"
        table.to_csv(out_path, index=False)
        print(f"\n=== {name} (n rows = {len(table)}) ===")
        print(table.to_string(index=False))
        print(f"-> {out_path}")

    plot_path = out_dir / "slopegraph.png"
    plot_slopegraph(merged, callers, config, plot_path)
    print(f"\n-> {plot_path}")


if __name__ == "__main__":
    main()
