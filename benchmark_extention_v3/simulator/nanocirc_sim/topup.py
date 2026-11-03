"""Finds circRNAs with too little BSJ-confirmed read support and writes a
follow-up expression profile targeting just those."""
from __future__ import annotations

import csv


def select_deficient_circles(ground_truth_reconciled_path: str, min_confirmed: int = 1) -> list:
    """Returns circRNA_ids with fewer than `min_confirmed` BSJ-confirmed reads."""
    deficient = []
    with open(ground_truth_reconciled_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        id_col = "circRNA_id" if "circRNA_id" in reader.fieldnames else "name"
        for row in reader:
            n_confirmed = int(row.get("n_reads_confirmed_observed", 0) or 0)
            if n_confirmed < min_confirmed:
                deficient.append(row[id_col])
    return deficient


def write_topup_expression_profile(deficient_ids: list, out_path: str) -> None:
    """Writes a NanoSim-compatible expression profile listing only the deficient IDs, equal TPM."""
    with open(out_path, "w") as fh:
        fh.write("ID\tcount\tTPM\n")
        for cid in deficient_ids:
            fh.write(f"{cid}\t0\t1.0\n")
