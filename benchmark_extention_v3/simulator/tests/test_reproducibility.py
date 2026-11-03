import filecmp

from nanocirc_sim.annotation import Annotation
from nanocirc_sim.config import SimConfig
from nanocirc_sim.intervals import AnnotationIndex
from nanocirc_sim.simulate import run_simulation


def _run(mid_fixture, outdir, seed):
    cfg = SimConfig.default()
    cfg.isoform_mean = 1.5
    return run_simulation(
        genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"], config=cfg, outdir=str(outdir),
        explicit_counts={"ecirna": 10, "eicirna": 8, "cirna": 5, "antisense": 4, "intergenic": 5},
        seed=seed, show_progress=False,
    )


def test_same_seed_gives_identical_outputs(mid_fixture, tmp_path):
    out1, out2 = tmp_path / "run1", tmp_path / "run2"
    r1 = _run(mid_fixture, out1, seed=123)
    r2 = _run(mid_fixture, out2, seed=123)

    assert len(r1["candidates"]) == len(r2["candidates"]) == 32
    for fname in ["circRNAs.fa", "circRNAs.mature.fa", "abundances.tsv", "circRNA_metadata.tsv",
                  "circRNAs.bed", "ground_truth.tsv"]:
        assert filecmp.cmp(out1 / fname, out2 / fname, shallow=False), f"{fname} differs between identical-seed runs"


def test_different_seed_gives_different_output(mid_fixture, tmp_path):
    out1, out2 = tmp_path / "run1", tmp_path / "run2"
    _run(mid_fixture, out1, seed=1)
    _run(mid_fixture, out2, seed=2)
    assert not filecmp.cmp(out1 / "circRNA_metadata.tsv", out2 / "circRNA_metadata.tsv", shallow=False)


def test_exact_per_type_counts_achieved(mid_fixture, tmp_path):
    result = _run(mid_fixture, tmp_path / "run", seed=7)
    stats = result["type_stats"]
    assert stats["ecirna"]["achieved"] == 10
    assert stats["eicirna"]["achieved"] == 8
    assert stats["cirna"]["achieved"] == 5
    assert stats["antisense"]["achieved"] == 4
    assert stats["intergenic"]["achieved"] == 5

    total_target_reads = SimConfig.default().total_reads
    reads_sum = sum(c.n_reads_target for c in result["candidates"])
    assert reads_sum == total_target_reads


def test_output_type_matches_final_coordinate_classification(mid_fixture, tmp_path):
    """Every written candidate's `type` must equal what
    AnnotationIndex.classify_final_type computes from its final coordinates."""
    result = _run(mid_fixture, tmp_path / "run", seed=11)
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    for c in result["candidates"]:
        blocks = [(b.start, b.end) for b in c.blocks]
        assert idx.classify_final_type(c.chrom, c.strand, blocks) == c.circ_type, c.id
