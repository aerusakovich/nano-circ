"""All simulation distributions/parameters, in one dataclass that round-trips
to/from YAML. Can be produced by dbfit.py from a real database, or
overridden via CLI flags."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields

import numpy as np
import yaml

CIRC_TYPES = ["ecirna", "eicirna", "cirna", "antisense", "intergenic"]

# Types generated from an annotated locus (share the isoform machinery),
# as opposed to cirna (single-intron) and intergenic (no locus).
LOCUS_TYPES = ["ecirna", "eicirna", "antisense"]


def _default_type_fractions():
    return {"ecirna": 0.545, "eicirna": 0.355, "cirna": 0.014, "antisense": 0.015, "intergenic": 0.071}


def _default_length_distributions():
    return {
        "ecirna": {"min": 80, "q1": 200, "median": 393, "q3": 700, "max": 3000},
        "eicirna": {"min": 100, "q1": 300, "median": 519.5, "q3": 900, "max": 4000},
        "cirna": {"min": 80, "q1": 150, "median": 314.5, "q3": 600, "max": 2000},
        "antisense": {"min": 80, "q1": 200, "median": 393, "q3": 700, "max": 3000},
        "intergenic": {"min": 80, "q1": 180, "median": 365, "q3": 650, "max": 2500},
    }


def _default_exon_count_distributions():
    return {
        "ecirna": {"1": 0.115, "2": 0.323, "3": 0.445, "4": 0.149, "5+": 0.193},
        "eicirna": {"1": 0.082, "2": 0.179, "3": 0.213, "4": 0.183, "5+": 0.343},
        "cirna": {"1": 0.594, "2": 0.241, "3": 0.094, "4": 0.071, "5+": 0.0},
        "antisense": {"1": 0.115, "2": 0.323, "3": 0.445, "4": 0.149, "5+": 0.193},
        "intergenic": {"1": 0.350, "2": 0.289, "3": 0.197, "4": 0.093, "5+": 0.069},
    }


def _default_splice_site_distributions():
    return {
        "ecirna": {"GT-AG": 1.0, "GC-AG": 0.0, "AT-AC": 0.0, "other": 0.0},
        "eicirna": {"GT-AG": 0.998, "GC-AG": 0.001, "AT-AC": 0.0002, "other": 0.0008},
        "cirna": {"GT-AG": 0.912, "GC-AG": 0.05, "AT-AC": 0.008, "other": 0.03},
        "antisense": {"GT-AG": 0.95, "GC-AG": 0.03, "AT-AC": 0.005, "other": 0.015},
        "intergenic": {"GT-AG": 0.505, "GC-AG": 0.2, "AT-AC": 0.05, "other": 0.245},
    }


def _default_exon_skip_prob():
    return {"ecirna": 0.05, "eicirna": 0.05, "antisense": 0.05}


def _default_intron_retention_prob():
    return {"eicirna": 0.7}


def _default_copy_number_ranges():
    return [[1, 3], [3, 6], [6, 10], [10, 15], [15, 20], [20, 50], [50, 100], [100, 1500]]


def _default_copy_number_probabilities():
    return [0.30, 0.30, 0.20, 0.10, 0.05, 0.03, 0.015, 0.005]


def _default_antisense_biotypes():
    return ["antisense", "antisense_RNA", "sense_antisense"]


def _default_splice_search_max_adjustment():
    # Per-type +/-bp splice search radius. Wider for antisense/intergenic,
    # which don't start from an already-canonical boundary (see README).
    return {"ecirna": 10, "eicirna": 10, "cirna": 10, "antisense": 30, "intergenic": 30}


@dataclass
class SimConfig:
    seed: int = 42
    total_reads: int = 100_000
    isoform_mean: float = 1.0
    isoform_max: int = 3
    min_intron_length: int = 100
    max_attempts_multiplier: int = 50
    shortfall_warn_threshold: float = 0.95

    # circRNA length bounds: too short or too long for any real long read to confirm.
    min_circle_length: int = 50
    max_circle_length: int = 50_000
    # Read-count floor per circRNA, so low-TPM circRNAs aren't left at zero.
    min_reads_per_circrna: int = 2
    # Probability a locus-based circRNA may use the transcript's first/last exon.
    terminal_exon_prob: float = 0.075
    # Probability a new isoform reuses a prior isoform's exact BSJ boundary.
    same_bsj_isoform_prob: float = 0.4

    # Per-type +/-bp window for the splice-boundary search.
    splice_search_max_adjustment: dict = field(default_factory=_default_splice_search_max_adjustment)

    antisense_biotypes: list = field(default_factory=_default_antisense_biotypes)
    antisense_name_pattern: str = r"-AS\d*$"

    type_fractions: dict = field(default_factory=_default_type_fractions)
    length_distributions: dict = field(default_factory=_default_length_distributions)
    exon_count_distributions: dict = field(default_factory=_default_exon_count_distributions)
    splice_site_distributions: dict = field(default_factory=_default_splice_site_distributions)
    exon_skip_prob: dict = field(default_factory=_default_exon_skip_prob)
    intron_retention_prob: dict = field(default_factory=_default_intron_retention_prob)
    copy_number_ranges: list = field(default_factory=_default_copy_number_ranges)
    copy_number_probabilities: list = field(default_factory=_default_copy_number_probabilities)

    @classmethod
    def default(cls) -> "SimConfig":
        return cls()

    @classmethod
    def from_yaml(cls, path: str) -> "SimConfig":
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
        for name, dist in [
            ("length_distributions", self.length_distributions),
            ("exon_count_distributions", self.exon_count_distributions),
            ("splice_site_distributions", self.splice_site_distributions),
            ("splice_search_max_adjustment", self.splice_search_max_adjustment),
        ]:
            for t in CIRC_TYPES:
                if t not in dist:
                    raise ValueError(f"config.{name} is missing an entry for circRNA type '{t}'")
        total_frac = sum(self.type_fractions.get(t, 0.0) for t in CIRC_TYPES)
        if total_frac <= 0:
            raise ValueError("config.type_fractions must sum to a positive value")

    def resolve_counts(self, total: int | None = None, explicit_counts: dict | None = None) -> dict:
        """Resolve exact per-type integer target counts, from explicit
        counts or by apportioning `total` across `type_fractions`."""
        if explicit_counts and any(v is not None for v in explicit_counts.values()):
            return {t: int(explicit_counts.get(t) or 0) for t in CIRC_TYPES}

        if total is None:
            raise ValueError("Must supply either --total or explicit per-type counts")

        fracs = {t: max(0.0, self.type_fractions.get(t, 0.0)) for t in CIRC_TYPES}
        norm = sum(fracs.values())
        if norm <= 0:
            raise ValueError("config.type_fractions must sum to a positive value")
        fracs = {t: f / norm for t, f in fracs.items()}

        raw = {t: total * fracs[t] for t in CIRC_TYPES}
        floors = {t: int(np.floor(v)) for t, v in raw.items()}
        remainder = total - sum(floors.values())
        order = sorted(CIRC_TYPES, key=lambda t: raw[t] - floors[t], reverse=True)
        counts = dict(floors)
        for t in order[:remainder]:
            counts[t] += 1
        return counts


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out
