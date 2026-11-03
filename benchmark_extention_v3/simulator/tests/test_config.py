from nanocirc_sim.config import SimConfig


def test_default_splice_search_max_adjustment_widens_antisense_and_intergenic():
    """antisense/intergenic get a wider search radius than the sense-side types."""
    cfg = SimConfig.default()
    adj = cfg.splice_search_max_adjustment
    assert adj["ecirna"] == 10
    assert adj["eicirna"] == 10
    assert adj["cirna"] == 10
    assert adj["antisense"] > adj["ecirna"]
    assert adj["intergenic"] > adj["ecirna"]


def test_config_validate_requires_splice_search_max_adjustment_for_every_type():
    cfg = SimConfig.default()
    cfg.validate()  # must not raise
    del cfg.splice_search_max_adjustment["antisense"]
    try:
        cfg.validate()
        assert False, "expected ValueError for missing type"
    except ValueError:
        pass


def test_splice_search_max_adjustment_survives_yaml_roundtrip(tmp_path):
    cfg = SimConfig.default()
    cfg.splice_search_max_adjustment["antisense"] = 42
    path = tmp_path / "cfg.yaml"
    cfg.to_yaml(str(path))
    loaded = SimConfig.from_yaml(str(path))
    assert loaded.splice_search_max_adjustment["antisense"] == 42
    assert loaded.splice_search_max_adjustment["ecirna"] == 10


def test_splice_search_max_adjustment_defaults_when_omitted_from_yaml(tmp_path):
    """An older config YAML written before this field existed must still
    load and pick up the field's defaults, not error or silently drop it."""
    path = tmp_path / "old_cfg.yaml"
    path.write_text("seed: 7\ntotal_reads: 500\n")
    loaded = SimConfig.from_yaml(str(path))
    assert loaded.seed == 7
    assert loaded.splice_search_max_adjustment == SimConfig.default().splice_search_max_adjustment
