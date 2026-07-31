"""End-to-end check of the `pool-select --subsample-to` CLI wiring."""
import csv

from nanocirc_sim.cli import main


def _write_metadata(path, targets: dict):
    with open(path, "w") as fh:
        fh.write("circRNA_id\tn_reads_target\n")
        for cid, t in targets.items():
            fh.write(f"{cid}\t{t}\n")


def _write_bsj_reads(path, rows):
    with open(path, "w") as fh:
        fh.write("read_id\tcircRNA_id\tlabel\tref_start\tref_end\tmature_length\tbsj_confirmed\n")
        for read_id, cid, label, confirmed in rows:
            fh.write(f"{read_id}\t{cid}\t{label}\t0\t100\t100\t{int(confirmed)}\n")


def _write_mega_fastq(path, read_ids):
    with open(path, "w") as fh:
        for rid in read_ids:
            fh.write(f"@{rid}\nACGT\n+\nIIII\n")


def test_subsample_to_is_seed_reproducible_and_seed_dependent(tmp_path):
    metadata = tmp_path / "meta.tsv"
    _write_metadata(metadata, {"a": 20, "b": 20})
    bsj = tmp_path / "bsj_reads.tsv"
    rows = [(f"a{i}", "a", "circ", True) for i in range(20)] + [(f"b{i}", "b", "circ", True) for i in range(20)]
    _write_bsj_reads(bsj, rows)
    mega = tmp_path / "mega.fastq"
    _write_mega_fastq(mega, [r[0] for r in rows])

    outdir1 = tmp_path / "run1"
    rc = main(["pool-select", "--bsj-reads", str(bsj), "--metadata", str(metadata),
               "--mega-fastq", str(mega), "--outdir", str(outdir1),
               "--subsample-to", "10", "--seed", "1"])
    assert rc == 0
    n1 = sum(1 for _ in open(outdir1 / "selected_reads.fastq")) // 4
    assert 0 < n1 <= 40  # nonzero, and can't exceed total availability

    outdir2 = tmp_path / "run2_same_seed"
    main(["pool-select", "--bsj-reads", str(bsj), "--metadata", str(metadata),
          "--mega-fastq", str(mega), "--outdir", str(outdir2),
          "--subsample-to", "10", "--seed", "1"])
    assert (outdir1 / "selected_reads.fastq").read_text() == (outdir2 / "selected_reads.fastq").read_text()

    outdir3 = tmp_path / "run3_diff_seed"
    main(["pool-select", "--bsj-reads", str(bsj), "--metadata", str(metadata),
          "--mega-fastq", str(mega), "--outdir", str(outdir3),
          "--subsample-to", "10", "--seed", "2"])
    assert (outdir1 / "selected_reads.fastq").read_text() != (outdir3 / "selected_reads.fastq").read_text()

    # deficient_circle_ids.txt reflects the full mega-pool, not the subsample
    assert (outdir1 / "deficient_circle_ids.txt").read_text() == ""


def test_subsample_to_stats_are_self_consistent_with_selected_fastq(tmp_path):
    """n_selected and n_confirmed_available in the stats file must match
    the actual number of reads written."""
    metadata = tmp_path / "meta.tsv"
    _write_metadata(metadata, {"a": 100})
    bsj = tmp_path / "bsj_reads.tsv"
    rows = [(f"a{i}", "a", "circ", True) for i in range(100)]
    _write_bsj_reads(bsj, rows)
    mega = tmp_path / "mega.fastq"
    _write_mega_fastq(mega, [r[0] for r in rows])

    outdir = tmp_path / "run"
    main(["pool-select", "--bsj-reads", str(bsj), "--metadata", str(metadata),
          "--mega-fastq", str(mega), "--outdir", str(outdir),
          "--subsample-to", "30", "--seed", "5"])

    with open(outdir / "selection_stats.tsv") as fh:
        rows_out = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows_out) == 1
    n_fastq = sum(1 for _ in open(outdir / "selected_reads.fastq")) // 4
    assert int(rows_out[0]["n_selected"]) == n_fastq
    assert int(rows_out[0]["n_confirmed_available"]) == n_fastq
    assert 0 < n_fastq <= 100


def test_subsample_to_low_abundance_circle_not_guaranteed_across_seeds(tmp_path):
    """A circRNA with a tiny target next to a dominant one is not
    guaranteed a read in every run."""
    metadata = tmp_path / "meta.tsv"
    _write_metadata(metadata, {"whale": 9999, "minnow": 1})
    bsj = tmp_path / "bsj_reads.tsv"
    rows = [(f"w{i}", "whale", "circ", True) for i in range(500)] + [("m0", "minnow", "circ", True)]
    _write_bsj_reads(bsj, rows)
    mega = tmp_path / "mega.fastq"
    _write_mega_fastq(mega, [r[0] for r in rows])

    presence = []
    for seed in range(10):
        outdir = tmp_path / f"run_seed{seed}"
        main(["pool-select", "--bsj-reads", str(bsj), "--metadata", str(metadata),
              "--mega-fastq", str(mega), "--outdir", str(outdir),
              "--subsample-to", "100", "--seed", str(seed)])
        content = (outdir / "selected_reads.fastq").read_text()
        presence.append("@m0" in content)

    assert any(p is False for p in presence), "minnow was present in every single run, expected some dropout"
    # whale (target=9999 of 10000) should be present in every run regardless
    for seed in range(10):
        content = (tmp_path / f"run_seed{seed}" / "selected_reads.fastq").read_text()
        assert any(f"@w{i}" in content for i in range(500))


def test_subsample_to_with_ground_truth_writes_run_bed(tmp_path):
    metadata = tmp_path / "meta.tsv"
    _write_metadata(metadata, {"a": 50, "b": 50})
    bsj = tmp_path / "bsj_reads.tsv"
    rows = [(f"a{i}", "a", "circ", True) for i in range(50)] + [(f"b{i}", "b", "circ", True) for i in range(50)]
    _write_bsj_reads(bsj, rows)
    mega = tmp_path / "mega.fastq"
    _write_mega_fastq(mega, [r[0] for r in rows])
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\t"
        "blockCount\tblockSizes\tblockStarts\ttype\tlocus_id\tisoform_index\tTPM\t"
        "n_reads_target\tn_reads_confirmed\n"
        "chr1\t0\t100\ta\t0\t+\t0\t100\t0,0,0\t1\t100,\t0,\tecirna\tLa\t0\t5.0000\t50\t\n"
        "chr1\t200\t300\tb\t0\t+\t200\t300\t0,0,0\t1\t100,\t0,\tecirna\tLb\t0\t5.0000\t50\t\n"
    )

    outdir = tmp_path / "run"
    main(["pool-select", "--bsj-reads", str(bsj), "--metadata", str(metadata),
          "--mega-fastq", str(mega), "--outdir", str(outdir),
          "--subsample-to", "40", "--seed", "1", "--ground-truth", str(gt)])

    bed_lines = (outdir / "run_circrnas.bed").read_text().splitlines()
    assert len(bed_lines) == 2
    scores = {line.split("\t")[3]: int(line.split("\t")[4]) for line in bed_lines}
    # equal targets, both with ideal=20 (well above 1) -> both essentially
    # certain to be present, though the exact counts are now random
    assert scores["a"] > 0 and scores["b"] > 0

    tsv_lines = (outdir / "run_ground_truth.tsv").read_text().splitlines()
    assert len(tsv_lines) == 3  # header + a + b
    header = tsv_lines[0].split("\t")
    assert header == [
        "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb",
        "blockCount", "blockSizes", "blockStarts", "type", "locus_id", "isoform_index",
        "TPM", "n_reads_target", "n_reads_confirmed",
    ]
    rows = [dict(zip(header, line.split("\t"))) for line in tsv_lines[1:]]
    tsv_scores = {r["name"]: int(r["score"]) for r in rows}
    assert tsv_scores == scores  # same per-run filtering/scores as the BED12 file
    for r in rows:
        assert r["n_reads_confirmed"] == r["score"]  # not the always-blank source column
