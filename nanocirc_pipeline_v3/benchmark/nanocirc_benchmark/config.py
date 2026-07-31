"""Benchmark parameters (tool/tier definitions, matching tolerances, plot
styling), in one dataclass that round-trips to/from YAML. Dataset location
and species/run selection are CLI concerns instead (see cli.py, paths.py),
the same split nanocirc_sim uses for genome/gtf/outdir."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields

import yaml


def _default_structural_param_sets():
    return [
        {"tag": "strict_5_90", "bsj_tolerance": 5, "min_overlap": 0.90},
        {"tag": "relaxed_10_75", "bsj_tolerance": 10, "min_overlap": 0.75},
    ]


def _default_tools():
    return ["cirilong", "isocirc", "circnick", "circfl"]


def _default_tool_labels():
    return {"cirilong": "CIRI-long", "isocirc": "IsoCirc", "circnick": "CircNick-LRS", "circfl": "circFL-seq"}


def _default_tool_colors():
    return {"cirilong": "#0072B2", "isocirc": "#009E73", "circnick": "#D55E00", "circfl": "#CC79A7"}


def _default_tiers():
    return ["high_confidence", "balanced_precision", "balanced_recall", "discovery"]


def _default_tier_labels():
    return {
        "high_confidence": "NanoCirc\nHigh confidence",
        "balanced_precision": "NanoCirc\nBalanced (Precision)",
        "balanced_recall": "NanoCirc\nBalanced (Recall)",
        "discovery": "NanoCirc\nDiscovery",
    }


def _default_tier_colors():
    return {
        "high_confidence": "#F1C40F",
        "balanced_precision": "#3EB489",
        "balanced_recall": "#87CEEB",
        "discovery": "#E63946",
    }


@dataclass
class BenchmarkConfig:
    bedtools_bin: str = "bedtools"

    # Coordinate-tolerance matching used for read-count correlation (a
    # caller's own called BSJ boundary is rarely byte-identical to GT's).
    gt_match_tolerance: int = 10

    # Structural (bedtools-based) precision/recall param sets: each produces
    # its own benchmark_metrics_<tag>.csv.
    structural_param_sets: list = field(default_factory=_default_structural_param_sets)

    # Expression-tertile sensitivity panel uses these two strategies against
    # each caller's native structural bed12.
    sensitivity_bsj_tolerance: int = 10
    sensitivity_min_overlap: float = 0.75

    # Exon-based overlap defining "TP" (a call structurally matches a GT
    # locus) for the length and type-classification modules. Both are about
    # a call's own structure, so exon-based match fits better than a
    # boundary-only BSJ tolerance.
    tp_min_overlap: float = 0.75

    # relaxed_bsj tolerance (bp) for the type-classification module's BSJ-basis
    # variant, matching the relaxed_10_75 structural_param_sets default.
    tp_bsj_tolerance: int = 10

    # Caller-overlap (UpSet) module makes both an exon-based view (its own,
    # stricter exon-overlap threshold, tighter than the tertile-sensitivity
    # panel's) and a relaxed-BSJ view (boundary tolerance in bp). A caller
    # can score very differently on the two, see overlap.py's module
    # docstring.
    overlap_min_overlap: float = 0.95
    overlap_bsj_tolerance: int = 10

    run_name_template: str = "{species}_run{i}"
    n_runs: int = 3

    tools: list = field(default_factory=_default_tools)
    tool_labels: dict = field(default_factory=_default_tool_labels)
    tool_colors: dict = field(default_factory=_default_tool_colors)

    tiers: list = field(default_factory=_default_tiers)
    tier_labels: dict = field(default_factory=_default_tier_labels)
    tier_colors: dict = field(default_factory=_default_tier_colors)

    @classmethod
    def default(cls) -> "BenchmarkConfig":
        return cls()

    @classmethod
    def from_yaml(cls, path: str) -> "BenchmarkConfig":
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        base = asdict(cls())
        merged = _deep_merge(base, data)
        valid_keys = {f.name for f in fields(cls)}
        merged = {k: v for k, v in merged.items() if k in valid_keys}
        return cls(**merged)

    def to_yaml(self, path: str) -> None:
        with open(path, "w") as fh:
            yaml.safe_dump(asdict(self), fh, sort_keys=False, default_flow_style=False)

    def validate(self) -> None:
        for t in self.tools:
            if t not in self.tool_labels or t not in self.tool_colors:
                raise ValueError(f"config.tools entry '{t}' is missing a label or color")
        for t in self.tiers:
            if t not in self.tier_labels or t not in self.tier_colors:
                raise ValueError(f"config.tiers entry '{t}' is missing a label or color")
        for p in self.structural_param_sets:
            if not {"tag", "bsj_tolerance", "min_overlap"} <= set(p):
                raise ValueError(f"structural_param_sets entry missing required keys: {p}")

    @property
    def quant_keys(self) -> list:
        """Tools followed by tiers, the fixed column order used in every
        correlation/sensitivity table and plot."""
        return list(self.tools) + list(self.tiers)

    def runs_for(self, species: str) -> list:
        return [self.run_name_template.format(species=species, i=i) for i in range(1, self.n_runs + 1)]


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out
