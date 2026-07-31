"""UpSet-style caller-overlap plot: every configured tool and NanoCirc tier,
all matched against GT, pooled across runs. "Detected by caller X" means
"matched by caller X's calls in at least one run."

Two matching strategies are made, since callers can fail in opposite ways.
Exon-based rewards correct internal exon/intron structure, but penalizes a
caller whose BSJ boundary is right but whose block layout differs from GT's.
Relaxed BSJ rewards accurate genomic positioning, but does not check if the
exon structure inside is correct. A caller can score very differently on
the two, so both are made instead of picking one as "the" answer.

With many sets, only the top N intersections by size are plotted (the
standard UpSet convention). Otherwise the matrix panel becomes unreadable.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

from .bed_io import read_bed
from .matching import StructuralMatcher

TOP_N_COMBOS = 25
STRATEGIES = ["exon_based", "relaxed_bsj"]


def _gt_key(df, strategy):
    """Pooling identity for a GT row. relaxed_bsj is a BSJ-only measure, so a
    BSJ-only key is correct there: isoforms sharing one BSJ really are "the
    same relaxed-BSJ entity". exon_based counts each distinct structure on
    its own, so it needs a key that includes block structure too. A
    BSJ-only key there would wrongly make two callers that matched two
    different isoforms of the same BSJ look like they matched the same
    thing."""
    base = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str) + ":" + df["strand"]
    if strategy == "exon_based":
        return base + ":" + df["blockSizes"] + ":" + df["blockStarts"]
    return base


def matched_gt_keys(matcher, tool_df, gt_df, gt_keys_by_row_id, work, strategy, min_overlap=None, tolerance=None):
    if strategy == "exon_based":
        _, _, _, matched_ids = matcher.exon_based_counts(tool_df, gt_df, work, min_overlap, return_ids=True)
    elif strategy == "relaxed_bsj":
        _, _, _, matched_ids = matcher.relaxed_bsj_counts(tool_df, gt_df, work, tolerance, return_ids=True)
    else:
        raise ValueError(f"unknown strategy: {strategy}")
    return {gt_keys_by_row_id[rid] for rid in matched_ids if rid in gt_keys_by_row_id}


def _strategy_label(strategy, config):
    if strategy == "exon_based":
        return f"exon-based match, ≥{config.overlap_min_overlap:.0%} overlap"
    return f"relaxed BSJ match, ±{config.overlap_bsj_tolerance}bp"


def _compute_detected(paths, config, strategy):
    quant_keys = config.quant_keys
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file)
    work_root = paths.bedtools_work_dir / "upset" / strategy

    detected = {t: set() for t in quant_keys}
    gt_total_keys = set()

    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        gt_keys = _gt_key(gt_df, strategy)
        gt_keys_by_row_id = dict(zip(gt_df["row_id"], gt_keys))
        gt_total_keys |= set(gt_keys)
        work = work_root / run
        work.mkdir(parents=True, exist_ok=True)
        for tool, path in paths.caller_paths(run, config.tools, config.tiers).items():
            if not path.exists():
                continue
            tool_df = read_bed(path)
            detected[tool] |= matched_gt_keys(matcher, tool_df, gt_df, gt_keys_by_row_id, work, strategy,
                                               min_overlap=config.overlap_min_overlap,
                                               tolerance=config.overlap_bsj_tolerance)
    return detected, gt_total_keys


def run_overlap_analysis(paths, config, species: str, strategy: str = "exon_based", out_path=None):
    if strategy not in STRATEGIES:
        raise ValueError(f"strategy must be one of {STRATEGIES}, got {strategy!r}")

    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    quant_colors = {**config.tool_colors, **config.tier_colors}

    print(f"[{species}] caller overlap ({_strategy_label(strategy, config)}) ...")
    detected, gt_total_keys = _compute_detected(paths, config, strategy)

    gt_total = len(gt_total_keys)
    print(f"Pooled unique GT circRNAs across runs: {gt_total}")
    for t in quant_keys:
        print(f"  {quant_labels[t].replace(chr(10), ' ')}: {len(detected[t])} detected")

    all_detected_union = set().union(*detected.values()) if detected else set()
    matrix = {}
    for key in all_detected_union:
        presence = tuple(1 if key in detected[t] else 0 for t in quant_keys)
        matrix[presence] = matrix.get(presence, 0) + 1

    sorted_combos = sorted(matrix.items(), key=lambda x: x[1], reverse=True)
    n_combos_total = len(sorted_combos)
    sorted_combos = sorted_combos[:TOP_N_COMBOS]
    combinations = [c for c, _ in sorted_combos]
    counts = [n for _, n in sorted_combos]
    n_tools, n_combos = len(quant_keys), len(combinations)
    print(f"Showing top {n_combos} of {n_combos_total} observed intersection patterns")

    combo_colors = {}
    for i, t in enumerate(quant_keys):
        key = tuple(1 if j == i else 0 for j in range(n_tools))
        combo_colors[key] = quant_colors[t]
    default_color = "#555555"
    bar_colors = [combo_colors.get(c, default_color) for c in combinations]

    fig = plt.figure(figsize=(max(16, n_combos * 0.9), 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[3, 0.3, 2.6], width_ratios=[1.1, 6], hspace=0.05, wspace=0.05)

    ax_bars = fig.add_subplot(gs[0, 1])
    x_pos = np.arange(n_combos)
    ax_bars.bar(x_pos, counts, color=bar_colors, edgecolor="black", linewidth=0.8)
    for i, cnt in enumerate(counts):
        pct = 100 * cnt / gt_total if gt_total else 0
        ax_bars.text(i, cnt + max(counts) * 0.02, f"{cnt}\n({pct:.1f}%)", ha="center", va="bottom",
                     fontsize=8, fontweight="bold")
    ax_bars.set_ylabel("Intersection size", fontsize=13, fontweight="bold")
    ax_bars.set_ylim(0, max(counts) * 1.18 if counts else 1)
    ax_bars.set_xlim(-0.5, n_combos - 0.5)
    ax_bars.set_xticks([])
    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)
    ax_bars.grid(axis="y", alpha=0.3, linestyle="--")

    ax_sets = fig.add_subplot(gs[2, 0])
    set_sizes = [len(detected[t]) for t in quant_keys]
    y_pos = np.arange(n_tools)
    ax_sets.barh(y_pos, set_sizes, color=[quant_colors[t] for t in quant_keys], edgecolor="black", linewidth=0.8)
    for i, sz in enumerate(set_sizes):
        ax_sets.text(sz + max(set_sizes) * 0.02, i, f"{sz}", ha="left", va="center", fontsize=10, fontweight="bold")
    ax_sets.set_yticks(y_pos)
    ax_sets.set_yticklabels([quant_labels[t].replace("\n", " ") for t in quant_keys], fontsize=11, fontweight="bold")
    ax_sets.set_xlabel("Set size", fontsize=11, fontweight="bold")
    ax_sets.set_xlim(0, max(set_sizes) * 1.2 if set_sizes else 1)
    ax_sets.invert_yaxis()
    ax_sets.spines["top"].set_visible(False)
    ax_sets.spines["right"].set_visible(False)
    ax_sets.grid(axis="x", alpha=0.3, linestyle="--")

    ax_matrix = fig.add_subplot(gs[2, 1])
    for ci, combo in enumerate(combinations):
        active = [i for i, v in enumerate(combo) if v == 1]
        combo_col = bar_colors[ci]
        for ti in range(n_tools):
            if ti in active:
                circ = patches.Circle((ci, ti), 0.14, facecolor=quant_colors[quant_keys[ti]],
                                      edgecolor="black", linewidth=1.1, zorder=3)
            else:
                circ = patches.Circle((ci, ti), 0.09, facecolor="lightgray", edgecolor="gray", linewidth=0.4, zorder=2)
            ax_matrix.add_patch(circ)
        if len(active) > 1:
            ax_matrix.plot([ci] * len(active), active, color=combo_col, lw=3, zorder=1)
    ax_matrix.set_xlim(-0.5, n_combos - 0.5)
    ax_matrix.set_ylim(-0.5, n_tools - 0.5)
    ax_matrix.set_yticks(range(n_tools))
    ax_matrix.set_yticklabels([])
    ax_matrix.set_xticks([])
    ax_matrix.invert_yaxis()
    ax_matrix.axis("off")

    strategy_label = _strategy_label(strategy, config)
    fig.suptitle(f"Caller overlap: detected GT circRNAs ({species.capitalize()})\n"
                 f"({strategy_label}), pooled across replicate runs",
                 fontsize=15, fontweight="bold", y=0.99)
    fig.text(0.5, 0.005, f"Top {n_combos} of {n_combos_total} observed intersection patterns shown. "
             f"Percentages relative to pooled GT (n={gt_total}). "
             "GT filtered to circRNAs with ≥1 BSJ-crossing read in each run's data.",
             ha="center", fontsize=9, color="#888888")
    tag = "exon" if strategy == "exon_based" else "bsj"
    out_path = out_path or paths.plot_dir / f"upset_plot_{species}_{tag}.png"
    fig.savefig(out_path, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_path}")
    return out_path


def run_overlap_analysis_both(paths, config, species: str) -> dict:
    """Produces both the exon-based and relaxed-BSJ UpSet plots, since they
    tell genuinely different stories about a caller's strengths (see module
    docstring)."""
    return {strategy: run_overlap_analysis(paths, config, species, strategy=strategy) for strategy in STRATEGIES}
