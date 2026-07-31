from nanocirc_benchmark.config import BenchmarkConfig


def test_default_validates():
    cfg = BenchmarkConfig.default()
    cfg.validate()


def test_quant_keys_order():
    cfg = BenchmarkConfig.default()
    assert cfg.quant_keys == cfg.tools + cfg.tiers


def test_runs_for():
    cfg = BenchmarkConfig.default()
    assert cfg.runs_for("mouse") == ["mouse_run1", "mouse_run2", "mouse_run3"]


def test_yaml_roundtrip(tmp_path):
    cfg = BenchmarkConfig.default()
    cfg.gt_match_tolerance = 42
    path = tmp_path / "config.yaml"
    cfg.to_yaml(path)
    loaded = BenchmarkConfig.from_yaml(path)
    assert loaded.gt_match_tolerance == 42
    assert loaded.tools == cfg.tools


def test_partial_yaml_override_keeps_defaults(tmp_path):
    path = tmp_path / "partial.yaml"
    path.write_text("gt_match_tolerance: 5\n")
    cfg = BenchmarkConfig.from_yaml(path)
    assert cfg.gt_match_tolerance == 5
    assert cfg.tools == BenchmarkConfig.default().tools


def test_validate_rejects_missing_tool_label():
    cfg = BenchmarkConfig.default()
    cfg.tools.append("newtool")
    try:
        cfg.validate()
        assert False, "expected ValueError"
    except ValueError:
        pass
