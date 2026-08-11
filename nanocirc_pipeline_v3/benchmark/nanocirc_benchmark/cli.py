from __future__ import annotations

import argparse
import logging
import sys

from .config import BenchmarkConfig
from .paths import Paths


def _add_common_args(p):
    p.add_argument("--dataset-root", required=True, help="Benchmark dataset root (fastq/, ground_truth/, results/, benchmark/)")
    p.add_argument("--species", required=True, help="Species/dataset key, e.g. mouse or human")
    p.add_argument("--config", help="YAML config overriding tool/tier/tolerance defaults (configs/default.yaml if omitted)")
    p.add_argument("--runs", nargs="+", default=None, help="Explicit run names (defaults to config.run_name_template x config.n_runs)")


def _add_all_parser(subparsers):
    p = subparsers.add_parser("all", help="Run every metric module, in order")
    _add_common_args(p)
    p.add_argument("--gtf", default=None, help="Gene annotation GTF (enables circFL-style type classification; omit to skip it)")
    return p


def _add_performance_parser(subparsers):
    p = subparsers.add_parser("performance", help="Structural precision/recall/F1 + PR scatter")
    _add_common_args(p)
    return p


def _add_expression_parser(subparsers):
    p = subparsers.add_parser("expression", help="Read-count correlation vs ground truth")
    _add_common_args(p)
    return p


def _add_sensitivity_parser(subparsers):
    p = subparsers.add_parser("sensitivity", help="Structural discovery sensitivity by GT-expression tertile")
    _add_common_args(p)
    return p


def _add_type_classification_parser(subparsers):
    p = subparsers.add_parser("type-classification", help="circRNA type classification (ecirna/eicirna/cirna/antisense/intergenic), TP-only")
    _add_common_args(p)
    p.add_argument("--gtf", required=True, help="Gene annotation GTF")
    return p


def _add_overlap_parser(subparsers):
    p = subparsers.add_parser("overlap", help="UpSet-style caller-overlap plot (exon-based + relaxed-BSJ variants)")
    _add_common_args(p)
    return p


def _add_length_parser(subparsers):
    p = subparsers.add_parser("length", help="Mature-length barplot (TP calls only) + structural TP/FP/FN contingency heatmap")
    _add_common_args(p)
    return p


def _add_reproducibility_parser(subparsers):
    p = subparsers.add_parser("reproducibility", help="Pairwise read-count agreement scatter-matrix across callers (GT + tools + tiers)")
    _add_common_args(p)
    return p


def _add_isoform_detection_parser(subparsers):
    p = subparsers.add_parser("isoform-detection", help="Isoform-level detection accuracy within shared-BSJ multi-isoform families")
    _add_common_args(p)
    return p


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="nanocirc-benchmark", description="circRNA quantifier/caller benchmark suite")
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_all_parser(subparsers)
    _add_performance_parser(subparsers)
    _add_expression_parser(subparsers)
    _add_sensitivity_parser(subparsers)
    _add_type_classification_parser(subparsers)
    _add_overlap_parser(subparsers)
    _add_length_parser(subparsers)
    _add_reproducibility_parser(subparsers)
    _add_isoform_detection_parser(subparsers)
    return parser


def _load_config(args) -> BenchmarkConfig:
    return BenchmarkConfig.from_yaml(args.config) if args.config else BenchmarkConfig.default()


def _build_paths(args, config, gtf_path=None) -> Paths:
    runs = args.runs or config.runs_for(args.species)
    return Paths(args.dataset_root, args.species, runs, gtf_path=gtf_path)


def _cmd_all(args) -> int:
    from .run_benchmark import run_all

    config = _load_config(args)
    config.validate()
    run_all(args.dataset_root, args.species, config=config, runs=args.runs, gtf_path=args.gtf)
    return 0


def _cmd_performance(args) -> int:
    from .performance import plot_pr_scatter, run_structural_matching
    from .run_benchmark import _relaxed_tag

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    metrics_by_tag = run_structural_matching(paths, config, args.species)
    relaxed_tag = _relaxed_tag(config)
    out_path = paths.plot_dir / f"pr_scatter_{args.species}.png"
    plot_pr_scatter(paths, config, args.species, metrics_by_tag[relaxed_tag], out_path)
    print(f"-> {out_path}")
    return 0


def _cmd_expression(args) -> int:
    from .expression import run_expression_analysis

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    run_expression_analysis(paths, config, args.species)
    return 0


def _cmd_sensitivity(args) -> int:
    from .sensitivity import run_sensitivity_analysis

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    run_sensitivity_analysis(paths, config, args.species)
    return 0


def _cmd_type_classification(args) -> int:
    from .type_classification import run_type_classification

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config, gtf_path=args.gtf)
    run_type_classification(paths, config, args.species)
    return 0


def _cmd_overlap(args) -> int:
    from .overlap import run_overlap_analysis_both

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    run_overlap_analysis_both(paths, config, args.species)
    return 0


def _cmd_length(args) -> int:
    from .length import plot_contingency_heatmap, plot_gt_coverage_stacked, plot_length_barplot
    from .run_benchmark import _relaxed_tag

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    plot_length_barplot(paths, config, args.species)
    plot_contingency_heatmap(paths, config, args.species, metrics_csv_tag=_relaxed_tag(config))
    plot_gt_coverage_stacked(paths, config, args.species, metrics_csv_tag=_relaxed_tag(config))
    return 0


def _cmd_reproducibility(args) -> int:
    from .reproducibility import run_reproducibility_analysis

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    run_reproducibility_analysis(paths, config, args.species)
    return 0


def _cmd_isoform_detection(args) -> int:
    from .isoform_detection import run_isoform_detection_analysis

    config = _load_config(args)
    config.validate()
    paths = _build_paths(args, config)
    run_isoform_detection_analysis(paths, config, args.species)
    return 0


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                         format="%(asctime)s - %(levelname)s - %(message)s")

    handlers = {
        "all": _cmd_all, "performance": _cmd_performance,
        "expression": _cmd_expression, "sensitivity": _cmd_sensitivity,
        "type-classification": _cmd_type_classification,
        "overlap": _cmd_overlap, "length": _cmd_length,
        "reproducibility": _cmd_reproducibility,
    }
    return handlers[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
