"""circRNA type classification, matching nanocirc's own `classify_types()`
(`bin/circrna_clean.py`, `bin/crossrun_annotate.py`): spliced BED12 blocks
vs gene/exon annotation with `bedtools intersect -split -wo`.

  same-strand gene overlap, exon_frac = exon_overlap_bp / spliced_length:
    exon_frac >= 1.0   -> ecirna
    0 < exon_frac < 1  -> eicirna
    exon_frac == 0     -> cirna
  no same-strand overlap, opposite-strand gene overlap -> antisense
  no gene overlap on either strand                     -> intergenic

Only TP calls count toward a type's "% of GT" bar; FP calls are classified
separately and shown as a red bar on top.
"""
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bed_io import read_bed, write_bed_with_id
from .gt import count_reads
from .matching import StructuralMatcher

ALL_TYPES = ["ecirna", "eicirna", "cirna", "antisense", "intergenic"]

# (strategy name, filename tag)
STRATEGIES = [
    ("exon_based", "exon"),
    ("relaxed_bsj", "relaxed_bsj"),
    ("isoform_level", "isoform_level"),
    ("strict_bsj", "strict_bsj"),
]


def _strategy_label(strategy, config):
    if strategy == "exon_based":
        return f"Exon match (reciprocal overlap {config.tp_min_overlap * 100:.0f}%)"
    if strategy == "relaxed_bsj":
        return f"Relaxed BSJ match (+/-{config.tp_bsj_tolerance}nt)"
    if strategy == "isoform_level":
        return (f"Isoform-level match (+/-{config.gt_match_tolerance}nt BSJ AND "
                 f"{config.isoform_level_min_overlap * 100:.0f}% exon overlap)")
    if strategy == "strict_bsj":
        return "Strict BSJ match"
    raise ValueError(f"unknown strategy: {strategy}")

TYPE_COLORS = {
    "ecirna": "#1b9e77", "eicirna": "#d95f02", "cirna": "#7570b3",
    "antisense": "#e6ab02", "intergenic": "#a6761d",
}
FP_COLOR = "#e41a1c"
FP_HATCH = "///"


def build_feature_beds(gtf_path, gene_bed, exon_bed, bedtools_bin="bedtools") -> None:
    """Builds gene.bed and exon.bed: 6-col BED (chrom, start0, end, gene_id,
    0, strand), sorted, the same way as modules/local/gtf_to_feature_bed.nf.

    exon.bed is merged (same-strand) after extraction. Different transcripts
    of the same gene often give overlapping exon intervals. Without merging,
    a locus with several exon-annotated isoforms gets its bases counted more
    than once, which can push exon_frac above 1.0 and misclassify an
    eicirna call as ecirna."""
    if gene_bed.exists() and exon_bed.exists():
        print(f"[skip] feature beds already built: {gene_bed}, {exon_bed}")
        return
    for feature, out_path in [("gene", gene_bed), ("exon", exon_bed)]:
        awk_script = (
            'BEGIN{OFS="\\t"} !/^#/ && $3=="%s" {'
            'n=".";'
            'for(i=9;i<=NF;i++){if($i=="gene_id"){n=$(i+1);gsub(/"/,"",n);break}}'
            'print $1,$4-1,$5,n,0,$7'
            '}' % feature
        )
        raw_path = out_path if feature == "gene" else out_path.with_suffix(".raw.bed")
        cmd = f"awk '{awk_script}' {gtf_path} | sort -k1,1 -k2,2n > {raw_path}"
        subprocess.run(cmd, shell=True, check=True)
    raw_exon = exon_bed.with_suffix(".raw.bed")
    merge_cmd = (
        f"{bedtools_bin} merge -s -i {raw_exon} -c 6 -o distinct "
        f"| awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,\".\",0,$4}}' > {exon_bed}"
    )
    subprocess.run(merge_cmd, shell=True, check=True)
    raw_exon.unlink()
    print(f"-> {gene_bed}\n-> {exon_bed}")


def classify_types(df, gene_bed, exon_bed, bedtools_bin, work_dir, id_prefix="x") -> dict:
    """df: a bed_io.read_bed()-style DataFrame (row_id, spliced_len, BED12
    columns). Returns {row_id: type}."""
    if len(df) == 0:
        return {}
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    bed_path = work_dir / f"{id_prefix}_calls.bed12"
    write_bed_with_id(df, bed_path, id_prefix)
    spliced = dict(zip(id_prefix + df["row_id"].astype(str), df["spliced_len"]))

    def run_intersect(ref_bed, strand_flag):
        out = subprocess.run(
            [bedtools_bin, "intersect", "-a", str(bed_path), "-b", str(ref_bed), strand_flag, "-split", "-wo"],
            capture_output=True, text=True,
        )
        totals = defaultdict(int)
        for line in out.stdout.strip().split("\n"):
            if not line:
                continue
            p = line.split("\t")
            if len(p) >= 19 and p[3]:
                try:
                    totals[p[3]] += int(p[18])
                except ValueError:
                    pass
        return totals

    gene_sense = run_intersect(gene_bed, "-s")
    gene_antisense = run_intersect(gene_bed, "-S")
    exon_sense = run_intersect(exon_bed, "-s")

    types = {}
    for name, slen in spliced.items():
        row_id = int(name[len(id_prefix):])
        in_sense = name in gene_sense
        in_antisense = name in gene_antisense
        if in_sense:
            exon_bp = exon_sense.get(name, 0)
            exon_frac = exon_bp / slen if slen > 0 else 0.0
            if exon_frac >= 1.0:
                t = "ecirna"
            elif exon_frac > 0:
                t = "eicirna"
            else:
                t = "cirna"
        elif in_antisense:
            t = "antisense"
        else:
            t = "intergenic"
        types[row_id] = t
    return types


def _match_ids(matcher, strategy, caller_df, gt_df, work, config):
    """Returns (tp, fp, fn, matched_gt_row_ids, matched_caller_row_ids) for
    any of the 4 structural strategies, using each one's own tolerance from
    config."""
    if strategy == "exon_based":
        return matcher.exon_based_counts(caller_df, gt_df, work, config.tp_min_overlap, return_ids="both")
    if strategy == "relaxed_bsj":
        return matcher.relaxed_bsj_counts(caller_df, gt_df, work, config.tp_bsj_tolerance, return_ids="both")
    if strategy == "isoform_level":
        return matcher.isoform_level_counts(caller_df, gt_df, work, config.gt_match_tolerance,
                                             config.isoform_level_min_overlap, return_ids="both")
    if strategy == "strict_bsj":
        return matcher.strict_bsj_counts(caller_df, gt_df, work, return_ids="both")
    raise ValueError(f"unknown strategy: {strategy}")


def run_type_classification(paths, config, species: str) -> dict:
    """Runs type classification for all 4 structural strategies (exon_based,
    relaxed_bsj, isoform_level, strict_bsj). A call's type label does not depend on match
    strategy, only which calls count as TP/FP does, so it is computed once
    per run and reused across strategies. Returns {strategy: DataFrame},
    and writes one plot and one CSV per strategy."""
    if paths.gtf_path is None:
        raise ValueError("type classification requires paths.gtf_path to be set (--gtf)")

    data_dir = paths.classification_data_dir
    gene_bed = data_dir / "gene.bed"
    exon_bed = data_dir / "exon.bed"
    build_feature_beds(paths.gtf_path, gene_bed, exon_bed, config.bedtools_bin)

    quant_keys = config.quant_keys
    quant_labels = {**config.tool_labels, **config.tier_labels}
    matcher = StructuralMatcher(config.bedtools_bin, genome_sizes_file=paths.genome_sizes_file)
    work_root = paths.bedtools_work_dir / "type_classification"
    n_reads = count_reads(paths.run_fastq(paths.runs[0]))

    rows_by_strategy = {strategy: [] for strategy, _ in STRATEGIES}
    for run in paths.runs:
        gt_df = read_bed(paths.gt_prefiltered_bed12(run))
        gt_types = classify_types(gt_df, gene_bed, exon_bed, config.bedtools_bin, work_root / run / "gt", id_prefix="g")
        gt_hist = Counter(gt_types.values())
        for strategy, _ in STRATEGIES:
            rows_by_strategy[strategy].append(
                {"run": run, "caller": "ground_truth", **{f"count_{t}": gt_hist.get(t, 0) for t in ALL_TYPES}})
        print(f"[{run}] ground_truth (n={len(gt_df)}): {dict(gt_hist)}")

        for q in quant_keys:
            path = paths.caller_bed_path(run, q, config.tiers)
            if not path.exists():
                print(f"  [SKIP] {q}, not found: {path}")
                continue
            caller_df = read_bed(path)
            if len(caller_df) == 0:
                continue
            caller_work = work_root / run / q
            caller_types = classify_types(caller_df, gene_bed, exon_bed, config.bedtools_bin, caller_work, id_prefix="c")

            for strategy, tag in STRATEGIES:
                strat_work = caller_work / tag
                strat_work.mkdir(parents=True, exist_ok=True)
                tp, fp, fn, matched_gt_ids, matched_caller_ids = _match_ids(
                    matcher, strategy, caller_df, gt_df, strat_work, config)
                # Typed by the caller's own structure, not the matched GT row(s)
                # (a caller row can match >1 near-duplicate GT isoform), so
                # sum(tp_<type>) == tp exactly.
                tp_hist = Counter(caller_types[rid] for rid in matched_caller_ids if rid in caller_types)
                fp_ids = set(caller_df["row_id"]) - matched_caller_ids
                fp_hist = Counter(caller_types[rid] for rid in fp_ids if rid in caller_types)

                rows_by_strategy[strategy].append({
                    "run": run, "caller": q,
                    **{f"tp_{t}": tp_hist.get(t, 0) for t in ALL_TYPES},
                    **{f"fp_{t}": fp_hist.get(t, 0) for t in ALL_TYPES},
                })
                print(f"[{run}] {q} [{strategy}] (n={len(caller_df)}, TP={tp}, FP={fp}): "
                      f"TP={dict(tp_hist)} FP={dict(fp_hist)}")

    out = {}
    for strategy, tag in STRATEGIES:
        label = _strategy_label(strategy, config)
        df = pd.DataFrame(rows_by_strategy[strategy])
        csv_path = data_dir / f"classification_results_{tag}.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n-> {csv_path}")

        gt_pooled = {t: df.loc[df["caller"] == "ground_truth", f"count_{t}"].sum() for t in ALL_TYPES}
        caller_df_pooled = df[df["caller"] != "ground_truth"].groupby("caller")[
            [f"tp_{t}" for t in ALL_TYPES] + [f"fp_{t}" for t in ALL_TYPES]].sum()

        out_png = paths.plot_dir / f"type_classification_{species}_{tag}.png"
        _plot_classification(gt_pooled, caller_df_pooled, quant_keys, quant_labels, species, out_png,
                              len(paths.runs), n_reads, label)
        out[strategy] = df
    return out


def _plot_classification(gt_pooled, caller_pooled, quant_keys, quant_labels, species, out_png, n_runs, n_reads,
                          strategy_label):
    n_panels = len(quant_keys) + 1
    ncols = 4
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 6.5 * nrows), subplot_kw={"projection": "polar"})
    axes = np.atleast_1d(axes).flatten()

    # shared radial scale across every non-GT panel, since FP overcounting
    # can push a sector past 100% of the GT count for that type.
    max_frac = 0.0
    for q in quant_keys:
        if q not in caller_pooled.index:
            continue
        for t in ALL_TYPES:
            gt_n = gt_pooled.get(t, 0)
            if gt_n <= 0:
                continue
            total_frac = (caller_pooled.loc[q, f"tp_{t}"] + caller_pooled.loc[q, f"fp_{t}"]) / gt_n
            max_frac = max(max_frac, total_frac)
    max_frac = max(max_frac, 1.0)
    r_axis_max = np.sqrt(max_frac) * 1.12
    tick_fracs = np.arange(0, max_frac + 0.2001, 0.2)
    tick_pos_shared = np.sqrt(tick_fracs)
    tick_labels_shared = [f"{f*100:.0f}%" for f in tick_fracs]

    def draw_gt(ax, vals, title):
        types = ALL_TYPES
        v = [vals.get(t, 0) for t in types]
        cols = [TYPE_COLORS[t] for t in types]
        theta = np.linspace(0, 2 * np.pi, len(types), endpoint=False)
        width = 2 * np.pi / len(types) * 0.75
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.bar(theta, np.ones(len(types)), width=width, bottom=0, color=cols, edgecolor="black", linewidth=1.0, alpha=0.9)
        ax.set_xticks(theta)
        ax.set_xticklabels([t if vv > 0 else "" for t, vv in zip(types, v)], fontsize=8)
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.5, 1.0]); ax.set_yticklabels(["0%", "50%", "100%"], fontsize=7)
        for ang, t, vv in zip(theta, types, v):
            if vv > 0:
                ax.text(ang, 0.5, f"{int(vv)}", ha="center", va="center", fontsize=9, fontweight="bold",
                        color="white", bbox=dict(boxstyle="round,pad=0.25", facecolor="black", alpha=0.8), zorder=11)
        ax.set_title(title, fontsize=11, fontweight="bold", pad=14)

    def draw_caller(ax, q, title):
        types = ALL_TYPES
        theta = np.linspace(0, 2 * np.pi, len(types), endpoint=False)
        width = 2 * np.pi / len(types) * 0.75
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        for ang, t in zip(theta, types):
            gt_n = gt_pooled.get(t, 0)
            tp_n = caller_pooled.loc[q, f"tp_{t}"] if q in caller_pooled.index else 0
            fp_n = caller_pooled.loc[q, f"fp_{t}"] if q in caller_pooled.index else 0
            tp_frac = (tp_n / gt_n) if gt_n > 0 else 0.0
            fp_frac = (fp_n / gt_n) if gt_n > 0 else 0.0
            r_tp = np.sqrt(min(tp_frac, max_frac))
            r_total = np.sqrt(min(tp_frac + fp_frac, max_frac))

            ax.bar([ang], [r_tp], width=width, bottom=0, color=TYPE_COLORS[t],
                   edgecolor="black", linewidth=1.0, alpha=0.9, zorder=5)
            if r_total > r_tp:
                ax.bar([ang], [r_total - r_tp], width=width, bottom=r_tp, color=FP_COLOR,
                       edgecolor="black", linewidth=0.8, alpha=0.75, hatch=FP_HATCH, zorder=6)
            if tp_n > 0:
                ax.text(ang, r_tp / 2, f"{int(tp_n)}", ha="center", va="center", fontsize=8, fontweight="bold",
                        color="white", bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.8), zorder=11)
            if fp_n > 0:
                ax.text(ang, (r_tp + r_total) / 2, f"+{int(fp_n)}", ha="center", va="center", fontsize=7,
                        fontweight="bold", color="white",
                        bbox=dict(boxstyle="round,pad=0.15", facecolor=FP_COLOR, alpha=0.9), zorder=12)

        ax.set_xticks(theta)
        ax.set_xticklabels(types, fontsize=8)
        ax.set_ylim(0, r_axis_max)
        ax.set_yticks(tick_pos_shared); ax.set_yticklabels(tick_labels_shared, fontsize=6.5)
        ax.set_title(title, fontsize=11, fontweight="bold", pad=14)

    gt_n_total = sum(gt_pooled.values())
    draw_gt(axes[0], gt_pooled, f"Ground Truth\n(n={gt_n_total}, pooled)")
    for i, q in enumerate(quant_keys, start=1):
        label = quant_labels[q].replace("\n", " ")
        n = int(sum(caller_pooled.loc[q, f"tp_{t}"] + caller_pooled.loc[q, f"fp_{t}"] for t in ALL_TYPES)) if q in caller_pooled.index else 0
        draw_caller(axes[i], q, f"{label}\n(n={n}, TP % of GT per type)")

    for ax in axes[n_panels:]:
        ax.set_visible(False)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=TYPE_COLORS[t], edgecolor="black", label=t) for t in ALL_TYPES]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=FP_COLOR, edgecolor="black", alpha=0.75, hatch=FP_HATCH, label="FP overcount"))
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=11, frameon=True,
               bbox_to_anchor=(0.5, -0.02), title="circRNA type")
    fig.suptitle(f"circRNA type detection, {strategy_label} ({species.capitalize()}, {n_runs} runs, {n_reads:,} reads)",
                 fontsize=15, fontweight="bold", y=1.02)
    fig.text(0.5, -0.02,
              "TP % of GT per type, pooled across runs. Bedtools-overlap classification, same algorithm as the "
              "production pipeline; FP = unmatched calls, own-coordinate type.",
              ha="center", va="top", fontsize=8, color="#888888")
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(out_png, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> {out_png}")
