from nanocirc_sim.pool_select import (
    deficient_circle_ids, filter_plan_to_selected, plan_run_targets, plan_selection, run_selected_read_ids,
    select_run_reads, selected_read_ids, write_run_bed, write_run_ground_truth_tsv, write_selected_fastq,
    write_selection_stats,
)


def _write_metadata(path, targets: dict):
    with open(path, "w") as fh:
        fh.write("circRNA_id\tn_reads_target\n")
        for cid, t in targets.items():
            fh.write(f"{cid}\t{t}\n")


def _write_bsj_reads(path, rows):
    """rows: list of (read_id, circRNA_id, label, bsj_confirmed)."""
    with open(path, "w") as fh:
        fh.write("read_id\tcircRNA_id\tlabel\tref_start\tref_end\tmature_length\tbsj_confirmed\n")
        for read_id, cid, label, confirmed in rows:
            fh.write(f"{read_id}\t{cid}\t{label}\t0\t100\t100\t{int(confirmed)}\n")


def test_plan_selection_caps_and_filters_correctly(tmp_path):
    metadata = tmp_path / "meta.tsv"
    _write_metadata(metadata, {"a": 2, "b": 5, "c": 0})

    bsj = tmp_path / "bsj_reads.tsv"
    _write_bsj_reads(bsj, [
        ("r1", "a", "circ", True),
        ("r2", "a", "circ", True),
        ("r3", "a", "circ", True),  # a has 3 confirmed but target is 2 -> only 2 selected
        ("r4", "a", "circ", False),  # not confirmed, excluded
        ("r5", "b", "circ", True),  # b has only 1 confirmed but target is 5 -> deficient-ish but not zero
        ("r6", "b", "linear", True),  # linear-labeled, excluded from circ selection
        ("r7", "d", "circ", True),  # circRNA not even in metadata, ignored
    ])

    plan = plan_selection(str(bsj), str(metadata))
    assert plan["a"]["target"] == 2
    assert plan["a"]["available"] == ["r1", "r2", "r3"]
    assert plan["b"]["available"] == ["r5"]
    assert plan["c"]["available"] == []
    assert "d" not in plan

    selected = selected_read_ids(plan)
    assert selected == {"r1", "r2", "r5"}  # capped at 2 for a, all 1 for b, none for c (target 0)


def test_write_selection_stats_and_deficient_ids(tmp_path):
    metadata = tmp_path / "meta.tsv"
    _write_metadata(metadata, {"a": 2, "b": 5, "c": 3})
    bsj = tmp_path / "bsj_reads.tsv"
    _write_bsj_reads(bsj, [
        ("r1", "a", "circ", True),
        ("r2", "a", "circ", True),
        # b: no confirmed reads at all -> deficient
        # c: no confirmed reads at all -> deficient
    ])
    plan = plan_selection(str(bsj), str(metadata))
    stats_path = tmp_path / "stats.tsv"
    stats = write_selection_stats(plan, str(stats_path))
    assert stats == {"n_circrnas": 3, "n_deficient": 2}
    assert set(deficient_circle_ids(plan)) == {"b", "c"}

    lines = stats_path.read_text().splitlines()
    assert lines[0] == "circRNA_id\tn_reads_target\tn_confirmed_available\tn_selected"
    assert "a\t2\t2\t2" in lines


def test_write_selected_fastq_streams_only_selected_reads(tmp_path):
    mega = tmp_path / "mega.fastq"
    mega.write_text(
        "@r1\nACGT\n+\nIIII\n"
        "@r2\nTTTT\n+\nIIII\n"
        "@r3\nGGGG\n+\nIIII\n"
    )
    out = tmp_path / "selected.fastq"
    n = write_selected_fastq(str(mega), {"r1", "r3"}, str(out))
    assert n == 2
    content = out.read_text()
    assert "@r1" in content and "@r3" in content and "@r2" not in content


def test_plan_run_targets_excludes_untargeted_and_unavailable_circles():
    plan = {
        "big": {"target": 10000, "available": [f"big{i}" for i in range(200)]},
        "zero_target": {"target": 0, "available": ["zt0"]},  # not meant to be expressed -> excluded
        "deficient": {"target": 50, "available": []},  # can't be represented regardless
    }
    targets = plan_run_targets(plan, n_target_total=100, seed=1)
    assert "zero_target" not in targets
    assert "deficient" not in targets
    assert targets.get("big", 0) > 0  # overwhelming share of the budget -> essentially certain to be drawn


def test_plan_run_targets_caps_at_availability():
    plan = {
        "a": {"target": 1000, "available": [f"a{i}" for i in range(5)]},  # only 5 reads exist for it
        "b": {"target": 10, "available": [f"b{i}" for i in range(1000)]},
    }
    targets = plan_run_targets(plan, n_target_total=50, seed=1)
    assert targets["a"] <= 5  # can never exceed what's actually available


def test_plan_run_targets_is_proportional_in_expectation_not_flat():
    """A circRNA's average read count over many draws should track its
    proportional target share."""
    plan = {
        "high": {"target": 900, "available": [f"h{i}" for i in range(100000)]},
        "low": {"target": 100, "available": [f"l{i}" for i in range(100000)]},
    }
    totals = {"high": 0, "low": 0}
    n_trials = 200
    for seed in range(n_trials):
        targets = plan_run_targets(plan, n_target_total=1000, seed=seed)
        totals["high"] += targets.get("high", 0)
        totals["low"] += targets.get("low", 0)
    mean_high = totals["high"] / n_trials
    mean_low = totals["low"] / n_trials
    assert abs(mean_high - 900) < 30  # within a few sqrt(900) std devs
    assert abs(mean_low - 100) < 15


def test_plan_run_targets_high_abundance_circle_present_in_nearly_every_run_low_abundance_varies():
    """A high-abundance circRNA is detected in nearly every run; a
    low-abundance one is sometimes present, sometimes absent."""
    plan = {
        "abundant": {"target": 100000, "available": [f"a{i}" for i in range(100000)]},
        "rare": {"target": 1, "available": ["r0"]},
    }
    n_trials = 200
    n_abundant_present = sum(1 for seed in range(n_trials)
                              if plan_run_targets(plan, n_target_total=1000, seed=seed).get("abundant", 0) > 0)
    n_rare_present = sum(1 for seed in range(n_trials)
                          if plan_run_targets(plan, n_target_total=1000, seed=seed).get("rare", 0) > 0)
    assert n_abundant_present == n_trials  # essentially always detected
    assert 0 < n_rare_present < n_trials  # sometimes present, sometimes absent: not deterministic either way


def test_plan_run_targets_different_seeds_give_different_circle_sets():
    """Independent replicate runs must not all contain the exact same
    circRNA ground truth set."""
    plan = {
        f"c{i}": {"target": 1, "available": [f"c{i}r{j}" for j in range(20)]}
        for i in range(200)
    }
    sets = [frozenset(plan_run_targets(plan, n_target_total=100, seed=seed).keys()) for seed in (1, 2, 3)]
    assert len(set(sets)) == 3, "all three runs had identical circRNA sets"
    # mostly different, only some similar: no pair should be (near-)identical,
    # and no pair should be totally disjoint either (some overlap expected).
    for i in range(3):
        for j in range(i + 1, 3):
            overlap = len(sets[i] & sets[j])
            union = len(sets[i] | sets[j])
            jaccard = overlap / union
            assert jaccard < 0.95, f"runs {i} and {j} are nearly identical (jaccard={jaccard:.3f})"
            assert overlap > 0, f"runs {i} and {j} share no circRNAs at all (expected some overlap)"


def test_plan_run_targets_same_seed_is_reproducible():
    plan = {
        "a": {"target": 500, "available": [f"a{i}" for i in range(1000)]},
        "b": {"target": 5, "available": [f"b{i}" for i in range(1000)]},
    }
    t1 = plan_run_targets(plan, n_target_total=1000, seed=42)
    t2 = plan_run_targets(plan, n_target_total=1000, seed=42)
    assert t1 == t2


def test_select_run_reads_respects_targets_and_seed_changes_which_reads():
    plan = {"a": {"target": 100, "available": [f"a{i}" for i in range(100)]}}
    targets = {"a": 10}
    r1 = select_run_reads(plan, targets, seed=1)
    r2 = select_run_reads(plan, targets, seed=1)
    r3 = select_run_reads(plan, targets, seed=2)
    assert len(r1["a"]) == 10
    assert r1 == r2  # same seed -> identical picks
    assert r1["a"] != r3["a"]  # different seed -> (almost certainly) different picks
    assert set(r1["a"]).issubset(set(plan["a"]["available"]))


def test_run_selected_read_ids_flattens_across_circles():
    run_reads = {"a": ["r1", "r2"], "b": ["r3"]}
    assert run_selected_read_ids(run_reads) == {"r1", "r2", "r3"}


def test_write_run_bed_sets_score_to_picked_read_count(tmp_path):
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\t"
        "blockCount\tblockSizes\tblockStarts\ttype\n"
        "chr1\t100\t200\tcircA\t500\t+\t100\t200\t0,0,0\t1\t100,\t0,\tecirna\n"
        "chr1\t300\t400\tcircB\t500\t+\t300\t400\t0,0,0\t1\t100,\t0,\tecirna\n"
    )
    run_targets = {"circA": 42}  # circB not in this run -> omitted entirely
    out = tmp_path / "run_circrnas.bed"
    n = write_run_bed(str(gt), run_targets, str(out))
    assert n == 1
    lines = out.read_text().splitlines()
    fields_a = lines[0].split("\t")
    assert fields_a[3] == "circA" and fields_a[4] == "42"
    assert len(fields_a) == 12  # strict BED12, no extra columns leaked through


def test_write_run_bed_omits_circles_not_present_in_this_run(tmp_path):
    """A detection tool has zero reads to find an unselected circRNA from
    in this run. Its ground truth must not include (and thus judge
    against) circRNAs that were never actually sequenced in this run."""
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\t"
        "blockCount\tblockSizes\tblockStarts\ttype\n"
        "chr1\t100\t200\tcircA\t500\t+\t100\t200\t0,0,0\t1\t100,\t0,\tecirna\n"
        "chr1\t300\t400\tcircB\t500\t+\t300\t400\t0,0,0\t1\t100,\t0,\tecirna\n"
        "chr1\t500\t600\tcircC\t500\t+\t500\t600\t0,0,0\t1\t100,\t0,\tecirna\n"
    )
    run_targets = {"circA": 3, "circB": 0}  # circB explicitly 0, circC absent from the dict entirely
    out = tmp_path / "run_circrnas.bed"
    n = write_run_bed(str(gt), run_targets, str(out))
    assert n == 1
    names = [line.split("\t")[3] for line in out.read_text().splitlines()]
    assert names == ["circA"]


def test_write_run_ground_truth_tsv_keeps_full_columns_with_same_filtering(tmp_path):
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\t"
        "blockCount\tblockSizes\tblockStarts\ttype\tlocus_id\tisoform_index\tTPM\t"
        "n_reads_target\tn_reads_confirmed\n"
        "chr1\t100\t200\tcircA\t500\t+\t100\t200\t27,158,119\t1\t100,\t0,\tecirna\tL1\t0\t12.3400\t50\t\n"
        "chr1\t300\t400\tcircB\t500\t+\t300\t400\t27,158,119\t1\t100,\t0,\tecirna\tL2\t0\t0.0100\t2\t\n"
    )
    run_targets = {"circA": 7}  # circB not drawn this run -> omitted
    out = tmp_path / "run_ground_truth.tsv"
    n = write_run_ground_truth_tsv(str(gt), run_targets, str(out))
    assert n == 1
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    assert header == [
        "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb",
        "blockCount", "blockSizes", "blockStarts", "type", "locus_id", "isoform_index",
        "TPM", "n_reads_target", "n_reads_confirmed",
    ]
    row = dict(zip(header, lines[1].split("\t")))
    assert row["name"] == "circA"
    assert row["score"] == "7"  # overwritten with the picked count
    assert row["n_reads_confirmed"] == "7"  # same picked count, not the (always-blank) source column
    assert row["type"] == "ecirna" and row["locus_id"] == "L1" and row["TPM"] == "12.3400"


def test_filter_plan_to_selected_restricts_available_and_stats_match(tmp_path):
    plan = {
        "a": {"target": 5, "available": ["r1", "r2", "r3", "r4"]},
        "b": {"target": 5, "available": ["r5", "r6"]},
    }
    selected = {"r1", "r3", "r5"}  # a subsample that dropped r2, r4, r6
    filtered = filter_plan_to_selected(plan, selected)
    assert filtered["a"]["available"] == ["r1", "r3"]
    assert filtered["b"]["available"] == ["r5"]
    assert filtered["a"]["target"] == 5  # target untouched

    stats = write_selection_stats(filtered, str(tmp_path / "stats.tsv"))
    assert stats == {"n_circrnas": 2, "n_deficient": 0}
