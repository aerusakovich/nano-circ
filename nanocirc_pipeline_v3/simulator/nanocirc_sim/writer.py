"""Writes all output files."""
from __future__ import annotations

import hashlib
import os
import time
from dataclasses import asdict

import yaml

from . import __version__
from .bed12 import bed12_fields

_FASTA_LINE_WIDTH = 60


def assign_ids(candidates: list) -> None:
    # No underscores: NanoSim's own header parsing mangles them.
    for i, c in enumerate(candidates):
        c.id = f"{c.circ_type}{i:07d}"


def write_fasta(candidates: list, path: str) -> None:
    with open(path, "w") as fh:
        for c in candidates:
            fh.write(f">{c.id}\n")
            seq = c.pseudo_seq or c.mature_seq
            for i in range(0, len(seq), _FASTA_LINE_WIDTH):
                fh.write(seq[i:i + _FASTA_LINE_WIDTH] + "\n")


def write_mature_fasta(candidates: list, path: str) -> None:
    """The single-copy mature sequence per circRNA (not the rolling-circle concatemer)."""
    with open(path, "w") as fh:
        for c in candidates:
            fh.write(f">{c.id}\n")
            for i in range(0, len(c.mature_seq), _FASTA_LINE_WIDTH):
                fh.write(c.mature_seq[i:i + _FASTA_LINE_WIDTH] + "\n")


def write_abundances(candidates: list, path: str) -> None:
    with open(path, "w") as fh:
        fh.write("circRNA_id\tmature_sequence\tTPM\n")
        for c in candidates:
            fh.write(f"{c.id}\t{c.mature_seq}\t{c.tpm:.4f}\n")


_METADATA_COLUMNS = [
    "circRNA_id", "chrom", "start", "end", "strand", "type",
    "gene_id", "tx_id", "locus_id", "isoform_index",
    "mature_length", "exon_count", "skipped_exons", "retained_introns",
    "block_count", "copy_number", "splice_type",
    "TPM", "n_reads_target", "n_reads_confirmed",
]


def write_metadata(candidates: list, path: str) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(_METADATA_COLUMNS) + "\n")
        for c in candidates:
            row = [
                c.id, c.chrom, str(c.start), str(c.end), c.strand, c.circ_type,
                c.gene_id, c.tx_id, c.locus_id, str(c.isoform_index),
                str(c.mature_length), str(c.exon_count), str(c.skipped_exons), str(c.retained_introns),
                str(c.block_count), str(c.copy_number), c.splice_type,
                f"{c.tpm:.4f}", str(c.n_reads_target),
                "" if c.n_reads_confirmed is None else str(c.n_reads_confirmed),
            ]
            fh.write("\t".join(row) + "\n")


def write_strict_bed12(candidates: list, path: str) -> None:
    with open(path, "w") as fh:
        for c in candidates:
            fh.write("\t".join(bed12_fields(c)) + "\n")


_GROUND_TRUTH_BED_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
]
_GROUND_TRUTH_EXTRA_COLS = ["type", "locus_id", "isoform_index", "TPM", "n_reads_target", "n_reads_confirmed"]


def write_ground_truth(candidates: list, path: str) -> None:
    """The authoritative table: BED12 + expression + read counts, per isoform."""
    with open(path, "w") as fh:
        fh.write("\t".join(_GROUND_TRUTH_BED_COLS + _GROUND_TRUTH_EXTRA_COLS) + "\n")
        for c in candidates:
            row = bed12_fields(c) + [
                c.circ_type, c.locus_id, str(c.isoform_index),
                f"{c.tpm:.4f}", str(c.n_reads_target),
                "" if c.n_reads_confirmed is None else str(c.n_reads_confirmed),
            ]
            fh.write("\t".join(row) + "\n")


def _file_fingerprint(path: str, hash_max_bytes: int = 50_000_000) -> dict:
    if not path or not os.path.exists(path):
        return {"path": path, "exists": False}
    st = os.stat(path)
    info = {"path": os.path.abspath(path), "exists": True, "size_bytes": st.st_size,
             "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime))}
    if st.st_size <= hash_max_bytes:
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        info["md5"] = h.hexdigest()
    else:
        info["md5"] = None  # too large to hash cheaply; size+mtime is the provenance signal
    return info


def write_run_manifest(path: str, *, config, seed: int, inputs: dict, type_stats: dict,
                        elapsed_seconds: float) -> None:
    manifest = {
        "nanocirc_sim_version": __version__,
        "seed": seed,
        "elapsed_seconds": round(elapsed_seconds, 2),
        "inputs": {name: _file_fingerprint(p) for name, p in inputs.items()},
        "type_stats": type_stats,
        "config": asdict(config),
    }
    with open(path, "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, default_flow_style=False)


def write_all_outputs(candidates: list, outdir: str, *, config, seed: int, inputs: dict,
                       type_stats: dict, elapsed_seconds: float) -> dict:
    os.makedirs(outdir, exist_ok=True)
    assign_ids(candidates)

    paths = {
        "fasta": os.path.join(outdir, "circRNAs.fa"),
        "mature_fasta": os.path.join(outdir, "circRNAs.mature.fa"),
        "abundances": os.path.join(outdir, "abundances.tsv"),
        "metadata": os.path.join(outdir, "circRNA_metadata.tsv"),
        "bed": os.path.join(outdir, "circRNAs.bed"),
        "ground_truth": os.path.join(outdir, "ground_truth.tsv"),
        "manifest": os.path.join(outdir, "run_manifest.yaml"),
    }
    write_fasta(candidates, paths["fasta"])
    write_mature_fasta(candidates, paths["mature_fasta"])
    write_abundances(candidates, paths["abundances"])
    write_metadata(candidates, paths["metadata"])
    write_strict_bed12(candidates, paths["bed"])
    write_ground_truth(candidates, paths["ground_truth"])
    write_run_manifest(
        paths["manifest"], config=config, seed=seed, inputs=inputs,
        type_stats=type_stats, elapsed_seconds=elapsed_seconds,
    )
    return paths
