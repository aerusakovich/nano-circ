from nanocirc_sim.topup import select_deficient_circles, write_topup_expression_profile


def test_select_deficient_circles(tmp_path):
    path = tmp_path / "ground_truth_reconciled.tsv"
    path.write_text(
        "circRNA_id\tn_reads_confirmed_observed\n"
        "a\t0\n"
        "b\t1\n"
        "c\t5\n"
        "d\t0\n"
    )
    assert select_deficient_circles(str(path), min_confirmed=1) == ["a", "d"]
    assert select_deficient_circles(str(path), min_confirmed=6) == ["a", "b", "c", "d"]


def test_write_topup_expression_profile(tmp_path):
    out = tmp_path / "topup.tsv"
    write_topup_expression_profile(["a", "d"], str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == "ID\tcount\tTPM"
    assert lines[1] == "a\t0\t1.0"
    assert lines[2] == "d\t0\t1.0"
