"""Checks that reads labeled "circular" actually cross the backsplice junction (BSJ)."""
from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from collections import defaultdict

import pysam

logger = logging.getLogger(__name__)


def _read_fasta(path: str) -> dict:
    seqs = {}
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        seqs[header] = "".join(chunks)
    return seqs


def _read_lengths_from_metadata(path: str) -> dict:
    """circRNA_id -> mature_length, from circRNA_metadata.tsv or ground_truth.tsv."""
    lengths = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "circRNA_id" in header:
            id_idx = header.index("circRNA_id")
        elif "name" in header:
            id_idx = header.index("name")
        else:
            return lengths
        length_idx = header.index("mature_length") if "mature_length" in header else None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= id_idx:
                continue
            cid = fields[id_idx]
            if length_idx is not None and len(fields) > length_idx:
                try:
                    lengths[cid] = int(fields[length_idx])
                except ValueError:
                    pass
    return lengths


def build_doubled_reference(mature_seqs: dict, out_path: str) -> None:
    with open(out_path, "w") as fh:
        for cid, seq in mature_seqs.items():
            if not seq:
                continue
            fh.write(f">{cid}\n")
            doubled = seq + seq
            for i in range(0, len(doubled), 60):
                fh.write(doubled[i:i + 60] + "\n")


def _run_minimap2(minimap2_path: str, ref_fasta: str, fastq_path: str, sam_path: str, threads: int) -> None:
    cmd = [minimap2_path, "-a", "-x", "map-ont", "-t", str(threads), "--secondary=no", ref_fasta, fastq_path]
    logger.info("Running: %s", " ".join(cmd))
    try:
        with open(sam_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not run minimap2 (looked for {minimap2_path!r}). "
            f"check-bsj requires minimap2 on PATH or --minimap2 /path/to/minimap2."
        ) from exc
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode("utf-8", errors="replace") if exc.stderr else ""
        raise RuntimeError(f"minimap2 failed (exit {exc.returncode}): {stderr}") from exc


def _read_label(qname: str) -> str:
    """Label from a NanoSim-style `_circ`/`_linear` header suffix."""
    if qname.endswith("_circ"):
        return "circ"
    if qname.endswith("_linear"):
        return "linear"
    return "unknown"


def parse_alignments_for_bsj(sam_path: str, mature_lengths: dict) -> dict:
    """Read a SAM of reads aligned to doubled references, return {qname: {...}} with bsj_confirmed set."""
    per_read = {}
    with pysam.AlignmentFile(sam_path, "r") as sam:
        for read in sam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            cid = sam.get_reference_name(read.reference_id)
            mature_len = mature_lengths.get(cid)
            if not mature_len:
                continue
            ref_start = read.reference_start
            ref_end = read.reference_end
            confirmed = ref_start < mature_len < ref_end
            per_read[read.query_name] = {
                "circRNA_id": cid, "ref_start": ref_start, "ref_end": ref_end,
                "mature_length": mature_len, "bsj_confirmed": confirmed,
                "label": _read_label(read.query_name),
            }
    return per_read


def _write_reads_table(per_read: dict, path: str) -> None:
    with open(path, "w") as fh:
        fh.write("read_id\tcircRNA_id\tlabel\tref_start\tref_end\tmature_length\tbsj_confirmed\n")
        for qname, r in per_read.items():
            fh.write(f"{qname}\t{r['circRNA_id']}\t{r['label']}\t{r['ref_start']}\t{r['ref_end']}\t"
                      f"{r['mature_length']}\t{int(r['bsj_confirmed'])}\n")


def _write_reconciled_metadata(metadata_path: str, per_circ: dict, out_path: str) -> None:
    with open(metadata_path) as fh, open(out_path, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        if "circRNA_id" in header:
            id_idx = header.index("circRNA_id")
        elif "name" in header:
            id_idx = header.index("name")
        else:
            id_idx = 0
        out.write("\t".join(header + ["n_reads_aligned_observed", "n_reads_confirmed_observed"]) + "\n")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            cid = fields[id_idx] if len(fields) > id_idx else ""
            stats = per_circ.get(cid, {"n_reads_aligned": 0, "n_reads_confirmed": 0})
            out.write("\t".join(fields + [str(stats["n_reads_aligned"]), str(stats["n_reads_confirmed"])]) + "\n")


def _read_fastq_records(path: str):
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                return  # truncated trailing record
            yield header, seq, plus, qual


def write_filtered_fastq(fastq_path: str, per_read: dict, out_path: str) -> dict:
    """Writes a FASTQ keeping only real linear reads and BSJ-confirmed circular reads."""
    n_in = n_kept = n_dropped_circ = 0
    with open(out_path, "w") as out:
        for header, seq, plus, qual in _read_fastq_records(fastq_path):
            n_in += 1
            qname = header[1:].strip().split()[0] if header.startswith("@") else header.strip().split()[0]
            info = per_read.get(qname)
            label = info["label"] if info else _read_label(qname)
            confirmed = info["bsj_confirmed"] if info else False
            if label == "linear" or confirmed:
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                n_kept += 1
            elif label == "circ":
                n_dropped_circ += 1
    return {"n_in": n_in, "n_kept": n_kept, "n_dropped_unconfirmed_circ": n_dropped_circ}


def summarize_per_circrna(per_read: dict) -> dict:
    per_circ = defaultdict(lambda: {"n_reads_aligned": 0, "n_reads_confirmed": 0})
    for r in per_read.values():
        per_circ[r["circRNA_id"]]["n_reads_aligned"] += 1
        if r["bsj_confirmed"]:
            per_circ[r["circRNA_id"]]["n_reads_confirmed"] += 1
    return dict(per_circ)


def check_bsj(*, fastq_path: str, mature_fasta_path: str, metadata_path: str, outdir: str,
              minimap2_path: str = "minimap2", threads: int = 4, write_filtered: bool = False) -> dict:
    os.makedirs(outdir, exist_ok=True)
    mature_seqs = _read_fasta(mature_fasta_path)
    mature_lengths = {cid: len(seq) for cid, seq in mature_seqs.items()}
    mature_lengths.update({k: v for k, v in _read_lengths_from_metadata(metadata_path).items() if v})

    with tempfile.TemporaryDirectory() as tmpdir:
        doubled_ref = os.path.join(tmpdir, "doubled_reference.fa")
        build_doubled_reference(mature_seqs, doubled_ref)
        sam_path = os.path.join(tmpdir, "aligned.sam")
        _run_minimap2(minimap2_path, doubled_ref, fastq_path, sam_path, threads)
        per_read = parse_alignments_for_bsj(sam_path, mature_lengths)

    _write_reads_table(per_read, os.path.join(outdir, "bsj_check_reads.tsv"))
    per_circ = summarize_per_circrna(per_read)
    _write_reconciled_metadata(metadata_path, per_circ, os.path.join(outdir, "ground_truth_reconciled.tsv"))

    if write_filtered:
        filtered_path = os.path.join(outdir, "filtered_confirmed_only.fastq")
        filter_stats = write_filtered_fastq(fastq_path, per_read, filtered_path)
        logger.info("Filtered FASTQ: %s", filter_stats)

    n_confirmed = sum(1 for r in per_read.values() if r["bsj_confirmed"])
    n_circ_labeled = sum(1 for r in per_read.values() if r["label"] == "circ")
    n_circ_labeled_confirmed = sum(1 for r in per_read.values() if r["label"] == "circ" and r["bsj_confirmed"])
    has_labels = n_circ_labeled > 0

    summary = {
        "n_reads_aligned": len(per_read),
        "n_confirmed": n_circ_labeled_confirmed if has_labels else n_confirmed,
        "n_circ_labeled": n_circ_labeled if has_labels else len(per_read),
    }
    logger.info("BSJ check summary: %s", summary)
    return summary
