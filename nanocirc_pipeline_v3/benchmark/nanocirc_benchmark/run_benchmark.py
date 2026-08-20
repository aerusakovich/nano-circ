"""Top-level orchestration: run every metric module for one species and
write all tables + plots under <dataset_root>/benchmark/.
"""
from .config import BenchmarkConfig
from .expression import run_expression_analysis
from .isoform_detection import run_isoform_detection_analysis
from .isoform_expression import run_isoform_expression_analysis
from .length import plot_contingency_heatmap, plot_gt_coverage_stacked, plot_length_barplot
from .overlap import run_overlap_analysis_all
from .paths import Paths
from .performance import (plot_pr_scatter, plot_pr_scatter_isoform, run_isoform_matching,
                           run_structural_matching)
from .reproducibility import run_reproducibility_analysis
from .sensitivity import run_sensitivity_analysis
from .type_classification import run_type_classification


def _relaxed_tag(config: BenchmarkConfig) -> str:
    for p in config.structural_param_sets:
        if p["bsj_tolerance"] == 10 and abs(p["min_overlap"] - 0.75) < 1e-9:
            return p["tag"]
    return config.structural_param_sets[-1]["tag"]


def run_all(dataset_root, species: str, config: BenchmarkConfig = None, runs: list = None,
            gtf_path=None) -> dict:
    """Runs every metric module in sequence. `gtf_path` is optional. When
    omitted, type classification is skipped, since it is the only module
    that needs gene annotation instead of just BED12 calls."""
    config = config or BenchmarkConfig.default()
    config.validate()
    runs = runs or config.runs_for(species)
    paths = Paths(dataset_root, species, runs, gtf_path=gtf_path)

    metrics_by_tag = run_structural_matching(paths, config, species)

    relaxed_tag = _relaxed_tag(config)
    pr_scatter_path = paths.plot_dir / f"pr_scatter_{species}.png"
    plot_pr_scatter(paths, config, species, metrics_by_tag[relaxed_tag], pr_scatter_path)
    print(f"-> {pr_scatter_path}")

    isoform_metrics = run_isoform_matching(paths, config, species)
    pr_scatter_isoform_path = paths.plot_dir / f"pr_scatter_isoform_{species}.png"
    plot_pr_scatter_isoform(paths, config, species, isoform_metrics, pr_scatter_isoform_path)
    print(f"-> {pr_scatter_isoform_path}")

    corr_df = run_expression_analysis(paths, config, species)
    struct_sens = run_sensitivity_analysis(paths, config, species)

    plot_length_barplot(paths, config, species)
    plot_contingency_heatmap(paths, config, species, metrics_csv_tag=relaxed_tag)
    plot_gt_coverage_stacked(paths, config, species, metrics_csv_tag=relaxed_tag)
    run_overlap_analysis_all(paths, config, species)
    run_reproducibility_analysis(paths, config, species)
    isoform_detection_results = run_isoform_detection_analysis(paths, config, species)
    isoform_expression_results = run_isoform_expression_analysis(paths, config, species)

    classification_df = None
    if paths.gtf_path is not None:
        classification_df = run_type_classification(paths, config, species)
    else:
        print(f"[{species}] no --gtf given, skipping type classification")

    return {
        "metrics_by_tag": metrics_by_tag, "isoform_metrics": isoform_metrics, "correlation": corr_df,
        "sensitivity": struct_sens, "classification": classification_df,
        "isoform_detection": isoform_detection_results, "isoform_expression": isoform_expression_results,
    }
