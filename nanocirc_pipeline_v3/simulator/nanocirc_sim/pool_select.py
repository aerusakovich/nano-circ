"""Builds the final circular read set from a one-shot over-generated "mega-pool"."""
from __future__ import annotations

import csv
import logging
import random

import numpy as np

logger = logging.getLogger(__name__)


def _read_targets(metadata_path: str) -> dict:
    targets = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        id_col = "circRNA_id" if "circRNA_id" in reader.fieldnames else "name"
        for row in reader:
            try:
                targets[row[id_col]] = int(row.get("n_reads_target", 0) or 0)
            except ValueError:
                targets[row[id_col]] = 0
    return targets


def plan_selection(bsj_reads_tsv: str, metadata_path: str) -> dict:
    """Returns {circRNA_id: {"target": int, "available": [read_id, ...]}} for every circRNA."""
    targets = _read_targets(metadata_path)
    plan = {cid: {"target": t, "available": []} for cid, t in targets.items()}

    with open(bsj_reads_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("label") != "circ":
                continue
            if row.get("bsj_confirmed") not in ("1", "True", "true"):
                continue
            cid = row["circRNA_id"]
            if cid in plan:
                plan[cid]["available"].append(row["read_id"])
    return plan


def selected_read_ids(plan: dict) -> set:
    """The read_ids to keep: up to `target` per circRNA."""
    selected = set()
    for info in plan.values():
        cap = info["target"]
        selected.update(info["available"][:cap] if cap > 0 else info["available"][:0])
    return selected


def write_selection_stats(plan: dict, out_path: str) -> dict:
    n_deficient = 0
    with open(out_path, "w") as fh:
        fh.write("circRNA_id\tn_reads_target\tn_confirmed_available\tn_selected\n")
        for cid, info in sorted(plan.items()):
            n_available = len(info["available"])
            n_selected = min(info["target"], n_available)
            if n_available == 0 and info["target"] > 0:
                n_deficient += 1
            fh.write(f"{cid}\t{info['target']}\t{n_available}\t{n_selected}\n")
    return {"n_circrnas": len(plan), "n_deficient": n_deficient}


def _iter_fastq_records(path: str):
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                return
            yield header, seq, plus, qual


def write_selected_fastq(mega_fastq_path: str, selected_ids: set, out_path: str) -> int:
    """Streams the mega-pool once, writing only the selected reads."""
    n_written = 0
    with open(out_path, "w") as out:
        for header, seq, plus, qual in _iter_fastq_records(mega_fastq_path):
            qname = header[1:].strip().split()[0] if header.startswith("@") else header.strip().split()[0]
            if qname in selected_ids:
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                n_written += 1
    return n_written


def deficient_circle_ids(plan: dict) -> list:
    """circRNA_ids with a nonzero target but zero confirmed reads available."""
    return [cid for cid, info in plan.items() if info["target"] > 0 and not info["available"]]


def undersupplied_shortfalls(plan: dict) -> dict:
    """{circRNA_id: target - available} for every circRNA short of its target,
    zero-available circs included. Used to weight a top-up simulation round
    toward the ones that actually need more reads, instead of every circ
    equally."""
    shortfalls = {}
    for cid, info in plan.items():
        short = info["target"] - len(info["available"])
        if info["target"] > 0 and short > 0:
            shortfalls[cid] = short
    return shortfalls


def write_shortfall_profile(shortfalls: dict, out_path: str) -> None:
    """NanoSim expression-profile format (ID, count, TPM), TPM set to each
    circRNA's own shortfall so simulator.py's abundance-proportional sampling
    puts more reads where the gap is bigger."""
    with open(out_path, "w") as fh:
        fh.write("ID\tcount\tTPM\n")
        for cid, short in sorted(shortfalls.items()):
            fh.write(f"{cid}\t0\t{short}\n")


def plan_run_targets(plan: dict, n_target_total: int, seed: int = 0) -> dict:
    """For one replicate run of `n_target_total` reads, draws each eligible
    circRNA's read count from a Poisson distribution centered on its
    proportional share of the run budget. Returns {circRNA_id: n_reads_for_this_run},
    each capped at actual read availability."""
    eligible = {cid: info for cid, info in plan.items() if info["target"] > 0 and info["available"]}
    if not eligible or n_target_total <= 0:
        return {}

    ids = sorted(eligible.keys())  # deterministic iteration order
    total_target = sum(eligible[cid]["target"] for cid in ids)
    rng = np.random.default_rng(seed)

    alloc = {}
    for cid in ids:
        ideal = eligible[cid]["target"] / total_target * n_target_total
        n = int(rng.poisson(ideal))
        if n <= 0:
            continue
        alloc[cid] = min(n, len(eligible[cid]["available"]))
    return alloc


def select_run_reads(plan: dict, run_targets: dict, seed: int) -> dict:
    """For each circRNA in `run_targets`, seeded-shuffles its available reads
    and takes the first `run_targets[cid]` of them. Returns {circRNA_id: [read_ids]}."""
    rng = random.Random(seed)
    picked = {}
    for cid in sorted(run_targets.keys()):
        n = run_targets[cid]
        if n <= 0:
            continue
        available = sorted(plan[cid]["available"])  # sorted first: deterministic shuffle input
        rng.shuffle(available)
        picked[cid] = available[:n]
    return picked


def run_selected_read_ids(run_reads: dict) -> set:
    """Flattens select_run_reads' per-circRNA read lists into one set."""
    selected = set()
    for ids in run_reads.values():
        selected.update(ids)
    return selected


def write_run_bed(ground_truth_path: str, run_targets: dict, out_path: str) -> int:
    """Writes a BED12 file with one row per circRNA present in this run, score
    column set to its picked read count. Returns the number of rows written."""
    n_written = 0
    with open(ground_truth_path, newline="") as fh, open(out_path, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_picked = run_targets.get(row["name"], 0)
            if n_picked <= 0:
                continue
            fields = [
                row["chrom"], row["start"], row["end"], row["name"],
                str(n_picked),
                row["strand"], row["thickStart"], row["thickEnd"], row["itemRgb"],
                row["blockCount"], row["blockSizes"], row["blockStarts"],
            ]
            out.write("\t".join(fields) + "\n")
            n_written += 1
    return n_written


_GROUND_TRUTH_TSV_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
    "type", "locus_id", "isoform_index", "TPM", "n_reads_target", "n_reads_confirmed",
]


def write_run_ground_truth_tsv(ground_truth_path: str, run_targets: dict, out_path: str) -> int:
    """Same per-run filtering as write_run_bed, but keeps every
    ground_truth.tsv column instead of strict BED12."""
    n_written = 0
    with open(ground_truth_path, newline="") as fh, open(out_path, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        out.write("\t".join(_GROUND_TRUTH_TSV_COLUMNS) + "\n")
        for row in reader:
            n_picked = run_targets.get(row["name"], 0)
            if n_picked <= 0:
                continue
            row["score"] = str(n_picked)
            row["n_reads_confirmed"] = str(n_picked)
            out.write("\t".join(row[col] for col in _GROUND_TRUTH_TSV_COLUMNS) + "\n")
            n_written += 1
    return n_written


def filter_plan_to_selected(plan: dict, selected_ids: set) -> dict:
    """Restricts each circRNA's `available` list to just the ids in `selected_ids`."""
    return {
        cid: {"target": info["target"], "available": [r for r in info["available"] if r in selected_ids]}
        for cid, info in plan.items()
    }
