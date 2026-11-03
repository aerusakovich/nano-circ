import pytest

from nanocirc_benchmark.paths import Paths


@pytest.fixture
def paths(tmp_path):
    return Paths(tmp_path, "mouse", ["mouse_run1", "mouse_run2"])


def test_run_dir(paths, tmp_path):
    assert paths.run_dir("mouse_run1") == tmp_path / "results" / "mouse" / "circrna" / "mouse_run1"


def test_bench_dir_created(paths, tmp_path):
    d = paths.bench_dir
    assert d.exists()
    assert d == tmp_path / "benchmark"


def test_caller_bed_path_tool(paths):
    p = paths.caller_bed_path("mouse_run1", "cirilong", ["high_confidence"])
    assert p.name == "mouse_run1_cirilong.bed12"


def test_caller_bed_path_tier(paths):
    p = paths.caller_bed_path("mouse_run1", "high_confidence", ["high_confidence"])
    assert p.name == "mouse_run1_high_confidence.bed12"
    assert p.parent.name == "merged"


def test_caller_bed_path_unknown_key_raises(paths):
    with pytest.raises(ValueError):
        paths.caller_bed_path("mouse_run1", "not_a_key", ["high_confidence"])


def test_crossrun_tier_bed12(paths, tmp_path):
    p = paths.crossrun_tier_bed12("discovery")
    assert p == tmp_path / "results" / "mouse" / "circrna" / "crossrun" / "mouse" / "discovery" / "mouse_discovery_crossrun.bed12"


def test_genome_sizes_file_missing_raises(paths):
    with pytest.raises(FileNotFoundError):
        paths.genome_sizes_file


def test_classification_data_dir_is_species_specific(tmp_path):
    # Regression test: classification_data_dir used to be shared across
    # species (bench_dir/"data"/"classification", no species component).
    # build_feature_beds() skips rebuilding gene.bed/exon.bed if they already
    # exist, so whichever species ran type-classification FIRST would poison
    # every other species' classification with its own gene/exon annotation
    # (caught via a real case: human classification silently ran against
    # mouse's gene.bed).
    mouse_paths = Paths(tmp_path, "mouse", ["mouse_run1"])
    human_paths = Paths(tmp_path, "human", ["human_run1"])
    assert mouse_paths.classification_data_dir != human_paths.classification_data_dir
    assert "mouse" in str(mouse_paths.classification_data_dir)
    assert "human" in str(human_paths.classification_data_dir)
