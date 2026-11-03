from nanocirc_sim.annotation import Annotation
from nanocirc_sim.intervals import AnnotationIndex


def test_gtf_coordinates_converted_to_0based_half_open(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    tx = ann.tx_structure["G1"]["G1T1"]
    assert (tx.exons[0].start, tx.exons[0].end) == (999, 1099)  # GTF 1000-1099 (1-based closed)
    assert tx.introns[0] == (1099, 1399)


def test_classify_region(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    assert idx.classify_region("chrT", "+", 999, 1099) == "ecirna"  # fully in exon0
    assert idx.classify_region("chrT", "+", 999, 1399) == "eicirna"  # exon0 + partial intron
    assert idx.classify_region("chrT", "+", 1099, 1399) == "cirna"  # fully intronic
    assert idx.classify_region("chrT", "+", 15000, 15100) == "intergenic"


def test_is_intergenic_checks_both_strands(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    assert idx.is_intergenic("chrT", 15000, 15100) is True
    # antisense gene G2-AS1 occupies this region on '-' -> not intergenic on either strand check
    assert idx.is_intergenic("chrT", 12000, 12099) is False


def test_classify_region_by_blocks_matches_split_semantics(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    # two exons with the intervening intron NOT included in blocks -> ecirna (mirrors bedtools -split)
    blocks = [(999, 1099), (1399, 1499)]
    assert idx.classify_region_by_blocks("chrT", "+", blocks) == "ecirna"
    # a single block spanning the same range (intron included) -> eicirna
    assert idx.classify_region_by_blocks("chrT", "+", [(999, 1499)]) == "eicirna"


def test_antisense_detection(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    biotypes = ["antisense", "antisense_RNA"]
    pattern = r"-AS\d*$"
    assert ann.is_antisense("G2-AS1", biotypes, pattern) is True
    assert ann.is_antisense("G1", biotypes, pattern) is False


def test_overlapping_genes(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    assert idx.overlapping_genes("chrT", "+", 999, 1099) == ["G1"]
    assert idx.overlapping_genes("chrT", "-", 12000, 12099) == ["G2-AS1"]
    assert idx.overlapping_genes("chrT", "+", 15000, 15100) == []


def test_classify_final_type_delegates_to_classify_region_by_blocks(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    assert idx.classify_final_type("chrT", "+", [(999, 1099), (1399, 1499)]) == "ecirna"
    assert idx.classify_final_type("chrT", "+", [(999, 1499)]) == "eicirna"
    assert idx.classify_final_type("chrT", "+", [(1099, 1399)]) == "cirna"


def test_classify_final_type_splits_intergenic_into_antisense(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    # no gene on '-' at this locus, but it overlaps G1's exon0 (999,1099) on '+'
    assert idx.classify_final_type("chrT", "-", [(999, 1099)]) == "antisense"
    # genuinely nothing on either strand
    assert idx.classify_final_type("chrT", "+", [(15000, 15100)]) == "intergenic"


def test_classify_final_type_accepts_precomputed_base_type(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    # a wrong base_type is trusted as given (caller's contract). This confirms
    # the fast path really skips recomputation rather than ignoring it
    assert idx.classify_final_type("chrT", "+", [(999, 1099), (1399, 1499)], base_type="eicirna") == "eicirna"


def test_classify_final_type_boundary_shift_into_intron_flips_ecirna_to_eicirna(mid_fixture):
    """A block shifted a few bp into the neighboring intron must classify
    as eicirna, not ecirna."""
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    # second block starts 3bp before exon1 (1399,1499), into intron0 (1099,1399)
    blocks = [(999, 1099), (1396, 1499)]
    assert idx.classify_final_type("chrT", "+", blocks) == "eicirna"
