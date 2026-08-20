#!/usr/bin/env python3
"""PCA over the benchmark datasets listed in --config: one row per dataset
(mean precision/recall/f1 per caller x strategy x threshold, averaged
across runs), to see how datasets group and how big each factor's effect is.

Also computes, from "paired_distances" in the config, the distance between
dataset pairs that differ in exactly one factor. See pca_datasets.example.json.

Usage:
    python3 pca_analysis.py --config pca_datasets.json --outdir .
    python3 pca_analysis.py --config pca_datasets.json --root-dir /path/to/benchmarks --outdir .
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_PALETTE = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d"]


def load_config(config_path: Path, root_dir: Path | None) -> tuple[dict, set, list, dict]:
    """Reads the JSON config: datasets (name -> root/species_tag/group),
    drop_callers, paired_distances, and an optional group -> color map.
    Dataset "root" paths are resolved relative to root_dir (default: the
    config file's own parent directory), so the config is portable across
    machines as long as root_dir is passed for a different layout."""
    raw = json.loads(config_path.read_text())
    base = root_dir if root_dir is not None else config_path.resolve().parent

    datasets = {}
    for name, spec in raw["datasets"].items():
        datasets[name] = (base / spec["root"], spec["species_tag"], spec["group"])

    drop_callers = set(raw.get("drop_callers", []))
    paired_distances = [tuple(row) for row in raw.get("paired_distances", [])]

    groups = sorted({group for _, _, group in datasets.values()})
    colors = dict(raw.get("group_colors", {}))
    for i, g in enumerate(groups):
        colors.setdefault(g, DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)])

    return datasets, drop_callers, paired_distances, colors


def load_dataset_rows(bench_dir: Path, species_tag: str, dataset_name: str, drop_callers: set) -> pd.DataFrame:
    rows = []
    for tag, strict_tol, relaxed_label in [("strict_5_90", True, "_tol5"), ("relaxed_10_75", False, "_tol10")]:
        path = bench_dir / f"benchmark_metrics_{species_tag}_{tag}.csv"
        df = pd.read_csv(path)
        df = df[~df["caller"].isin(drop_callers)]
        if strict_tol:
            # strict_bsj is tag-independent (identical in both files); take it only here
            keep = df[df["strategy"] == "strict_bsj"].copy()
            keep["feature"] = keep["caller"] + "__strict_bsj"
            rows.append(keep)
            rb = df[df["strategy"] == "relaxed_bsj"].copy()
            rb["feature"] = rb["caller"] + "__relaxed_bsj_tol5"
            rows.append(rb)
            eb = df[df["strategy"] == "exon_based"].copy()
            eb["feature"] = eb["caller"] + "__exon_based_90"
            rows.append(eb)
        else:
            rb = df[df["strategy"] == "relaxed_bsj"].copy()
            rb["feature"] = rb["caller"] + "__relaxed_bsj_tol10"
            rows.append(rb)
            eb = df[df["strategy"] == "exon_based"].copy()
            eb["feature"] = eb["caller"] + "__exon_based_75"
            rows.append(eb)
    combined = pd.concat(rows, ignore_index=True)
    combined["dataset"] = dataset_name
    return combined


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True,
                    help="JSON file listing datasets (root/species_tag/group), drop_callers, "
                         "paired_distances, and optionally group_colors. See pca_datasets.example.json.")
    p.add_argument("--root-dir", default=None,
                    help="Base directory each dataset's 'root' is resolved against. "
                         "Default: the config file's own parent directory.")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()
    outdir = Path(args.outdir)
    root_dir = Path(args.root_dir) if args.root_dir else None

    datasets, drop_callers, paired_distances, colors = load_config(Path(args.config), root_dir)

    all_rows = []
    group_of = {}
    for dataset_name, (bench_dir, species_tag, group) in datasets.items():
        all_rows.append(load_dataset_rows(bench_dir, species_tag, dataset_name, drop_callers))
        group_of[dataset_name] = group
    long_df = pd.concat(all_rows, ignore_index=True)

    # mean precision/recall/f1 across the replicate runs, per dataset x feature
    agg = long_df.groupby(["dataset", "feature"])[["precision", "recall", "f1"]].mean().reset_index()
    melted = agg.melt(id_vars=["dataset", "feature"], value_vars=["precision", "recall", "f1"],
                       var_name="metric", value_name="value")
    melted["colname"] = melted["feature"] + "__" + melted["metric"]
    wide = melted.pivot(index="dataset", columns="colname", values="value")
    wide = wide.dropna(axis=1)  # keep only features present in every dataset
    print(f"Feature matrix: {wide.shape[0]} datasets x {wide.shape[1]} features")

    X = wide.values
    X_std = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    X_std = np.nan_to_num(X_std)  # a zero-variance column (all datasets tied) -> 0, not NaN

    U, S, Vt = np.linalg.svd(X_std, full_matrices=False)
    scores = U * S
    explained = (S ** 2) / np.sum(S ** 2)

    print("\nExplained variance ratio per PC:")
    for i, ev in enumerate(explained[:5]):
        print(f"  PC{i + 1}: {ev:.3f}")

    print("\nPC1/PC2 coordinates:")
    for i, name in enumerate(wide.index):
        print(f"  {name:30s} PC1={scores[i, 0]:8.3f}  PC2={scores[i, 1]:8.3f}  group={group_of[name]}")

    # loadings: which original features drive each PC the most (biggest |weight|)
    loadings = pd.DataFrame(Vt.T, index=wide.columns, columns=[f"PC{i+1}" for i in range(Vt.shape[0])])
    print("\nTop 10 loadings for PC1 (what separates datasets along the x-axis):")
    print(loadings["PC1"].abs().sort_values(ascending=False).head(10).to_string())
    print("\nTop 10 loadings for PC2 (what separates datasets along the y-axis):")
    print(loadings["PC2"].abs().sort_values(ascending=False).head(10).to_string())

    # Effect-size table: full standardized-feature-space Euclidean distance (all PCs,
    # not just PC1/PC2) between pairs of datasets that differ in exactly one factor.
    dist_rows = []
    for a, b, label in paired_distances:
        if a not in wide.index or b not in wide.index:
            continue
        dist = np.linalg.norm(X_std[list(wide.index).index(a)] - X_std[list(wide.index).index(b)])
        dist_rows.append({"dataset_a": a, "dataset_b": b, "effect": label, "distance": dist})
    dist_df = pd.DataFrame(dist_rows).sort_values("distance", ascending=False)
    print("\nEffect size: full standardized-feature-space distance between paired datasets\n"
          "(differ in exactly one factor, sorted biggest effect first):")
    print(dist_df.to_string(index=False))
    dist_df.to_csv(outdir / "pca_pairwise_distances.csv", index=False)
    print(f"-> {outdir / 'pca_pairwise_distances.csv'}")

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for i, name in enumerate(wide.index):
        g = group_of[name]
        ax.scatter(scores[i, 0], scores[i, 1], c=colors[g], s=140, edgecolors="black", zorder=3)
    # declutter labels: for each point, count how many earlier points are within a
    # small radius, and stack its label further along a diagonal offset per collision
    span = max(scores[:, 0].max() - scores[:, 0].min(), scores[:, 1].max() - scores[:, 1].min())
    radius = span * 0.06
    placed = []
    for i, name in enumerate(wide.index):
        x, y = scores[i, 0], scores[i, 1]
        collisions = sum(1 for px, py in placed if abs(px - x) < radius and abs(py - y) < radius)
        placed.append((x, y))
        dx, dy = span * 0.015, span * (0.035 + 0.045 * collisions)
        ax.annotate(name, (x, y), xytext=(x + dx, y + dy), textcoords="data",
                    fontsize=9, ha="left",
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=0.8) if collisions else None)
    ax.set_ylim(scores[:, 1].min() - span * 0.1, scores[:, 1].max() + span * 0.18)
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}% variance)")
    ax.set_title("Benchmark datasets: PCA over pooled precision/recall/f1")
    ax.axhline(0, color="#dddddd", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#dddddd", linewidth=0.8, zorder=1)
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=10, label=g)
               for g, c in colors.items()]
    ax.legend(handles=handles, loc="best")
    fig.tight_layout()
    out_png = outdir / "pca_benchmarks.png"
    fig.savefig(out_png, dpi=150)
    print(f"\n-> {out_png}")

    wide.to_csv(outdir / "pca_feature_matrix.csv")
    print(f"-> {outdir / 'pca_feature_matrix.csv'}")
    loadings.to_csv(outdir / "pca_loadings.csv")
    print(f"-> {outdir / 'pca_loadings.csv'}")


if __name__ == "__main__":
    main()
