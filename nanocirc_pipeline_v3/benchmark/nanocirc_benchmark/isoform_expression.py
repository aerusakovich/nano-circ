"""Read-count correlation at true isoform granularity, within shared-BSJ
multi-isoform families (same GT population as isoform_detection.py).

expression.py collapses both sides to one value per BSJ (gt.py sums reads
per locus, load_nanocirc_locus_counts keeps the highest-count isoform per
BSJ). This module scores each isoform's own count instead, reusing
isoform_detection.py's family definition and matching strategies.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

from .bed_io import read_bed
from .callers import (load_nanocirc_isoform_bed_counts, load_tier_membership,
                       load_tool_isoform_bed_counts)
from .expression import compute_corr
from .gt import count_reads
from .isoform_detection import (_bsj_key, _block_key, _strategies, _strategy_min_overlap,
                                 _strategy_suffix, _strategy_title, find_multi_isoform_families)
from .matching import StructuralMatcher
from .viz import sig_stars


def _block_identity_pred_map(families: pd.DataFrame, caller_df: pd.DataFrame) -> dict:
    """family row_id -> caller n_reads, matched by exact (bsj_key, block_key),
    same exact-string test isoform_detection.py's block_identity strategy
    uses for its true_iso classification."""
    if len(caller_df) == 0:
        return {}
    caller_df = caller_df.copy()
    caller_df["bsj_key"] = caller_df.apply(_bsj_key, axis=1)
    caller_df["block_key"] = caller_df.apply(_block_key, axis=1)
    exact = {(r.bsj_key, r.block_key): r.n_reads for r in caller_df.itertuples()}
    return {r.row_id: exact[(r.bsj_key, r.block_key)]
            for r in families.itertuples() if (r.bsj_key, r.block_key) in exact}


def _reciprocal_pred_map(families: pd.DataFrame, caller_df: pd.DataFrame, matcher, work, min_overlap: float) -> dict:
    """family row_id -> caller n_reads, matched with reciprocal exon overlap
    (matching.py's exon_based_pairs). Several caller isoforms can match one
    family isoform under this looser strategy; collapsed with max, not sum,
    same principle expression.py's best_match_by_strategy_isoform_resolved
    uses: this scores one specific isoform, so a near-duplicate report of
    it should not inflate its count."""
    if len(caller_df) == 0:
        return {}
    pairs = matcher.exon_based_pairs(caller_df, families, work, min_overlap=min_overlap)
    reads_by_caller_row = dict(zip(caller_df["row_id"], caller_df["n_reads"]))
    pred = {}
    for c_row_id, g_row_id in pairs:
        n = reads_by_caller_row.get(c_row_id, 0)
        pred[g_row_id] = max(pred.get(g_row_id, 0), n)
    return pred


def build_isoform_long_df(paths, config, species: str, strategy: str = "block_identity") -> pd.DataFrame:
    """One row per (run, quant, family GT isoform): gt_expr (that isoform's
    own confirmed read count, gt_prefiltered_bed12's score column) +
    pred_expr (this caller's own count for whichever isoform matched under
    `strategy`, 0 if none), restricted to GT loci with 2+ real isoforms.
    Neither side is ever pooled across isoforms."""
    is_reciprocal = strategy != "block_identity"
    matcher = StructuralMatcher(config.bedtools_bin) if is_reciprocal else None
    min_overlap = _strategy_min_overlap(strategy, config) if is_reciprocal else None

    rows = []
    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        families = find_multi_isoform_families(gt_df)
        if families.empty:
            continue
        work = None
        if is_reciprocal:
            work = paths.bedtools_work_dir / "isoform_expression" / strategy / run
            work.mkdir(parents=True, exist_ok=True)

        tier_members = load_tier_membership(paths, config.tiers)
        for quant in config.quant_keys:
            if quant in config.tools:
                caller_df = load_tool_isoform_bed_counts(paths, run, quant)
            else:
                caller_df = load_nanocirc_isoform_bed_counts(paths, run, tier_members[quant], paths.bedtools_work_dir)

            pred_map = (_block_identity_pred_map(families, caller_df) if strategy == "block_identity"
                        else _reciprocal_pred_map(families, caller_df, matcher, work, min_overlap))

            for r in families.itertuples():
                rows.append({"run": run, "quant": quant, "family_id": r.family_id,
                             "family_size": r.family_size, "gt_expr": int(r.score),
                             "pred_expr": pred_map.get(r.row_id, 0)})
    return pd.DataFrame(rows)


def plot_isoform_expression_corr(corr_df: pd.DataFrame, config, species: str, n_reads: int, n_runs: int,
                                  n_families: int, out_path, strategy: str = "block_identity") -> None:
    """Same heatmap convention as expression.py's plot_correlation_summary
    right-hand panel, on its own since a per-tool scatter panel would need
    isoform-to-isoform pairing across tools, not just against GT."""
    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    metrics = ["pearson_r", "spearman_rho", "r_squared"]
    metric_labels = ["Pearson r", "Spearman rho", "R2"]
    p_cols = {"pearson_r": "pearson_p", "spearman_rho": "spearman_p", "r_squared": "pearson_p"}
    sd_cols = {"pearson_r": "pearson_r_sd", "spearman_rho": "spearman_rho_sd", "r_squared": None}

    cdf = corr_df.set_index("quant").reindex(quant_keys)
    hm = cdf[metrics].astype(float).values

    fig, ax = plt.subplots(figsize=(7.5, 0.55 * len(quant_keys) + 2))
    im = ax.imshow(hm, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(metrics))); ax.set_xticklabels(metric_labels, fontsize=12, fontweight="bold")
    ax.set_yticks(range(len(quant_keys)))
    ax.set_yticklabels([quant_labels[q] for q in quant_keys], fontsize=11, fontweight="bold")
    for i, q in enumerate(quant_keys):
        for j, m in enumerate(metrics):
            v = hm[i, j]
            if np.isnan(v):
                txt = "n/a"
            else:
                sd_col = sd_cols[m]
                sd_txt = f" +/-{cdf.loc[q, sd_col]:.3f}" if sd_col else ""
                txt = f"{v:.3f}{sd_txt}\n({sig_stars(cdf.loc[q, p_cols[m]])})"
            ax.text(j, i, txt, ha="center", va="center", fontsize=9, fontweight="bold")
    ax.set_title(f"Per-isoform expression agreement vs ground truth\n"
                 f"({species.capitalize()}, {_strategy_title(strategy, config)}, {n_runs} runs, "
                 f"{n_reads:,} reads, {n_families} multi-isoform BSJ families)",
                 fontsize=11.5, fontweight="bold")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.text(0.5, -0.06, "*** p<0.001, ** p<0.01, * p<0.05, ns = not significant. Mean +/- sd across runs.",
            transform=ax.transAxes, ha="center", va="top", fontsize=8, color="#888888")

    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run_isoform_expression_analysis(paths, config, species: str, strategies=None) -> dict:
    """Runs block_identity and every config.isoform_overlap_param_sets
    strategy (all by default), each producing its own CSV + heatmap plot
    (suffixed for non-default strategies), same convention as
    run_isoform_detection_analysis. Returns {strategy: corr_df}."""
    strategies = strategies or _strategies(config)
    n_reads = count_reads(paths.run_fastq(paths.runs[0]))
    results = {}

    for strategy in strategies:
        print(f"[{species}] isoform expression, {strategy} ({_strategy_title(strategy, config)}) ...")
        long_df = build_isoform_long_df(paths, config, species, strategy)
        if long_df.empty:
            print(f"[{species}] no multi-isoform BSJ families found in GT, skipping isoform expression correlation")
            results[strategy] = None
            continue

        corr_df = compute_corr(long_df, config)
        corr_df["strategy"] = strategy
        print(corr_df.to_string(index=False))
        results[strategy] = corr_df

        suffix = _strategy_suffix(strategy)
        corr_csv = paths.bench_dir / f"isoform_expression_correlation{suffix}_{species}.csv"
        corr_df.to_csv(corr_csv, index=False)
        print(f"-> {corr_csv}")

        out_path = paths.plot_dir / f"isoform_expression_correlation{suffix}_{species}.png"
        n_families = long_df["family_id"].nunique()
        plot_isoform_expression_corr(corr_df, config, species, n_reads, len(paths.runs), n_families, out_path, strategy)
        print(f"-> {out_path}")

    return results
