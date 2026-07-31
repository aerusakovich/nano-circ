from nanocirc_sim.annotation import Annotation
from nanocirc_sim.blocks import total_length, validate_blocks
from nanocirc_sim.config import SimConfig
from nanocirc_sim.genome import GenomeLoader, reverse_complement
from nanocirc_sim.intervals import AnnotationIndex
from nanocirc_sim.isoforms import generate_with_isoforms
from nanocirc_sim.rng import RNGFactory
from nanocirc_sim.splice_sites import donor_acceptor
from nanocirc_sim.types import antisense, cirna, ecirna, eicirna, intergenic
from nanocirc_sim.types.base import generate_until_satisfied


def _setup(mid_fixture):
    ann = Annotation.from_gtf(mid_fixture["gtf"])
    idx = AnnotationIndex(ann)
    genome = GenomeLoader(mid_fixture["genome"])
    cfg = SimConfig.default()
    cfg.isoform_mean = 1.8
    rng = RNGFactory(42).generator
    return ann, idx, genome, cfg, rng


def _check_all(cands, expect_type):
    for c in cands:
        assert c.circ_type == expect_type
        validate_blocks(c.blocks)
        assert total_length(c.blocks) == len(c.mature_seq)
        assert c.start == c.blocks[0].start and c.end == c.blocks[-1].end


def test_ecirna_exact_count_and_exon_skipping(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = ecirna.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, attempts, rej = generate_with_isoforms(30, rng, loci, attempt, cfg, label="t", show_progress=False)
    assert len(res) == 30
    _check_all(res, "ecirna")
    assert any(c.skipped_exons > 0 for c in res)


def test_ecirna_multi_isoform_loci_are_distinct(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = ecirna.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, _, _ = generate_with_isoforms(25, rng, loci, attempt, cfg, label="t", show_progress=False)
    from collections import defaultdict
    by_locus = defaultdict(list)
    for c in res:
        by_locus[c.locus_id].append(c)
    multi = {k: v for k, v in by_locus.items() if len(v) > 1}
    assert multi, "expected at least one locus to produce multiple isoforms"
    for locus_id, group in multi.items():
        sigs = {tuple((b.start, b.end, b.kind) for b in c.blocks) for c in group}
        assert len(sigs) == len(group), "isoforms from the same locus must be structurally distinct"
        assert sorted(c.isoform_index for c in group) == list(range(len(group)))


def test_eicirna_always_has_retained_intron(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = eicirna.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, _, _ = generate_with_isoforms(20, rng, loci, attempt, cfg, label="t", show_progress=False)
    assert len(res) == 20
    _check_all(res, "eicirna")
    assert all(c.retained_introns >= 1 for c in res)
    assert all(any(b.kind == "intron_retained" for b in c.blocks) for c in res)


def test_cirna_single_intron_only(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = cirna.make_attempt_fn(rng, genome, ann, idx, cfg, loci)
    # the fixture's one gene has only 7 introns, so target=15 can't fully
    # succeed once duplicate-intron draws are (correctly) deduplicated --
    # see test_cirna_dedup_rejects_repeated_intron for that behavior.
    res, _, _ = generate_until_satisfied(5, attempt, cfg.max_attempts_multiplier, label="t", show_progress=False)
    assert len(res) == 5
    _check_all(res, "cirna")


def test_cirna_dedup_rejects_repeated_intron(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = cirna.make_attempt_fn(rng, genome, ann, idx, cfg, loci)
    # only 7 distinct introns exist in the fixture. Requesting far more
    # than that must not silently return duplicate circRNAs to pad the count.
    res, attempts, rejection_rate = generate_until_satisfied(
        50, attempt, cfg.max_attempts_multiplier, label="t", show_progress=False)
    assert len(res) < 50
    signatures = {(c.chrom, c.strand, tuple((b.start, b.end, b.kind) for b in c.blocks)) for c in res}
    assert len(signatures) == len(res), "duplicate circRNAs should never appear in the results"
    assert all(len(c.blocks) == 1 for c in res)


def test_antisense_built_from_sense_loci_on_flipped_strand(mid_fixture):
    """Antisense candidates are built from the same sense_loci pool as
    ecirna/eicirna, but placed on the flipped strand."""
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = antisense.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, _, _ = generate_with_isoforms(10, rng, loci, attempt, cfg, label="t", show_progress=False)
    assert len(res) == 10
    _check_all(res, "antisense")
    assert all(c.gene_id == "G1" for c in res)  # built from the host gene's exon coordinates
    assert all(c.strand == "-" for c in res)  # G1 itself is '+' -> antisense circRNAs sit on '-'


def test_antisense_rejected_if_gene_exists_on_flipped_strand(tmp_path):
    """If a gene already exists on the circRNA's own (flipped) strand,
    the candidate must be rejected."""
    import random as _random

    rnd = _random.Random(5)
    chrom_len = 5000
    seq = [rnd.choice("ACGT") for _ in range(chrom_len)]
    exon_len, intron_len = 100, 300
    gtf_lines = []

    def make_gene(gene_id, strand, n_exons, start_1based, tx_id):
        pos = start_1based
        coords = []
        for _ in range(n_exons):
            coords.append((pos, pos + exon_len - 1))
            pos += exon_len + intron_len
        g_start, g_end = coords[0][0], coords[-1][1]
        gtf_lines.append(f'chrT\ttest\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; gene_name "{gene_id}"; gene_biotype "protein_coding";')
        gtf_lines.append(f'chrT\ttest\ttranscript\t{g_start}\t{g_end}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_biotype "protein_coding";')
        for s, e in coords:
            gtf_lines.append(f'chrT\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\t'
                              f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_biotype "protein_coding";')

    # Two protein-coding genes at the SAME locus, opposite strands: a
    # classic overlapping sense/antisense gene pair.
    make_gene("PLUS", "+", 4, 1000, "PLUST1")
    make_gene("MINUS", "-", 4, 1000, "MINUST1")

    genome_path = tmp_path / "overlap.fa"
    gtf_path = tmp_path / "overlap.gtf"
    with open(genome_path, "w") as fh:
        fh.write(">chrT\n")
        s = "".join(seq)
        for i in range(0, len(s), 70):
            fh.write(s[i:i + 70] + "\n")
    gtf_path.write_text("\n".join(gtf_lines) + "\n")

    ann = Annotation.from_gtf(str(gtf_path))
    idx = AnnotationIndex(ann)
    genome = GenomeLoader(str(genome_path))
    cfg = SimConfig.default()
    rng = RNGFactory(1).generator

    attempt = antisense.make_attempt_fn(rng, genome, ann, idx, cfg)
    # PLUS's own exon coords, flipped to '-', collide with MINUS (a real
    # gene already on '-' at this exact locus). Must be rejected.
    result = attempt("PLUS", "PLUST1")
    assert result is None


def test_antisense_multi_exon_shifts_internal_junction_independently(tmp_path):
    """A multi-exon antisense circRNA's internal junction must be shifted
    to the antisense-strand-canonical motif, not left at the host gene's
    own sense-direction boundary."""
    length = 2000
    seq = list("CCCC" * (length // 4))
    # G1 (+): exon0=[1000,1100), exon1=[1400,1500), real sense-canonical
    # splice site at the unshifted boundary (typical real gene).
    seq[1100:1102] = list("GT")  # host donor (+ direction)
    seq[1398:1400] = list("AG")  # host acceptor (+ direction)
    # Genuinely antisense(-)-canonical motif planted nearby, NOT at G1's
    # own boundary. Only reachable within antisense's widened radius.
    seq[1105:1107] = list(reverse_complement("AG"))
    seq[1394:1396] = list(reverse_complement("GT"))
    # Outer BSJ: also genuinely '-'-canonical, unshifted, so build() only
    # fails (or not) because of the internal junction under test.
    seq[998:1000] = list(reverse_complement("GT"))
    seq[1500:1502] = list(reverse_complement("AG"))

    gtf_lines = [
        'chrT\ttest\tgene\t1001\t1500\t.\t+\t.\tgene_id "G1"; gene_name "G1"; gene_biotype "protein_coding";',
        'chrT\ttest\ttranscript\t1001\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "G1T1";',
        'chrT\ttest\texon\t1001\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "G1T1";',
        'chrT\ttest\texon\t1401\t1500\t.\t+\t.\tgene_id "G1"; transcript_id "G1T1";',
    ]
    genome_path = tmp_path / "g.fa"
    gtf_path = tmp_path / "g.gtf"
    with open(genome_path, "w") as fh:
        fh.write(">chrT\n")
        s = "".join(seq)
        for i in range(0, len(s), 70):
            fh.write(s[i:i + 70] + "\n")
    gtf_path.write_text("\n".join(gtf_lines) + "\n")

    ann = Annotation.from_gtf(str(gtf_path))
    idx = AnnotationIndex(ann)
    genome = GenomeLoader(str(genome_path))
    cfg = SimConfig.default()
    cfg.splice_site_distributions["antisense"] = {"GT-AG": 1.0, "GC-AG": 0.0, "AT-AC": 0.0, "other": 0.0}
    rng = RNGFactory(1).generator

    attempt = antisense.make_attempt_fn(rng, genome, ann, idx, cfg)
    result = attempt.build_with_boundary("G1", "G1T1", (0, 1))
    assert result is not None
    assert result.strand == "-"
    b0, b1 = result.blocks
    # The internal junction must have moved away from G1's own raw exon
    # boundary (1100/1400). Reusing it verbatim is exactly the bug.
    assert (b0.end, b1.start) != (1100, 1400)
    assert (b0.end, b1.start) == (1105, 1396)  # lands on the planted antisense-canonical motif
    donor, acceptor = donor_acceptor(genome, "chrT", b1.start, b0.end, "-")
    assert (donor, acceptor) == ("GT", "AG")  # genuinely canonical on its own (antisense) strand


def test_ecirna_rejects_boundary_shift_that_lands_in_intron(mid_fixture, monkeypatch):
    """A boundary shifted into the neighboring intron must be rejected,
    not emitted mislabeled as ecirna."""
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    from nanocirc_sim.types import ecirna as ecirna_mod

    # Force the "splice-boundary search" to shift the last block's end 3bp
    # past exon1's real end (1499), into intron1, regardless of what
    # desired splice type it's asked for. This reproduces the exact failure
    # mode found on real data (a few bp shortfall at one boundary).
    monkeypatch.setattr(
        ecirna_mod, "find_boundary_for_desired_type",
        lambda genome, chrom, start, end, strand, desired, max_adjustment=10: (start, end + 3, "other"),
    )
    attempt = ecirna_mod.make_attempt_fn(rng, genome, ann, idx, cfg)
    result = attempt.build_with_boundary("G1", "G1T1", (0, 1))  # exon0(999,1099)+exon1(1399,1499)
    assert result is None


def test_eicirna_rejects_retained_intron_thats_exonic_in_another_isoform(tmp_path):
    """A retained intron that's exonic in another isoform of the same gene
    must be rejected, not emitted mislabeled as eicirna."""
    import random as _random

    rnd = _random.Random(9)
    chrom_len = 1000
    seq = "".join(rnd.choice("ACGT") for _ in range(chrom_len))
    gtf_lines = [
        'chrX\ttest\tgene\t81\t520\t.\t+\t.\tgene_id "GX"; gene_name "GX"; gene_biotype "protein_coding";',
        'chrX\ttest\ttranscript\t101\t500\t.\t+\t.\tgene_id "GX"; transcript_id "GXT1";',
        'chrX\ttest\texon\t101\t200\t.\t+\t.\tgene_id "GX"; transcript_id "GXT1";',
        'chrX\ttest\texon\t401\t500\t.\t+\t.\tgene_id "GX"; transcript_id "GXT1";',
        # A second isoform whose single exon spans (with margin for the
        # splice-boundary search) exactly what GXT1's intron would be.
        'chrX\ttest\ttranscript\t81\t520\t.\t+\t.\tgene_id "GX"; transcript_id "GXT2";',
        'chrX\ttest\texon\t81\t520\t.\t+\t.\tgene_id "GX"; transcript_id "GXT2";',
    ]
    genome_path = tmp_path / "cover.fa"
    gtf_path = tmp_path / "cover.gtf"
    with open(genome_path, "w") as fh:
        fh.write(">chrX\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    gtf_path.write_text("\n".join(gtf_lines) + "\n")

    ann = Annotation.from_gtf(str(gtf_path))
    idx = AnnotationIndex(ann)
    genome = GenomeLoader(str(genome_path))
    cfg = SimConfig.default()
    rng = RNGFactory(1).generator

    attempt = eicirna.make_attempt_fn(rng, genome, ann, idx, cfg)
    result = attempt.build_with_boundary("GX", "GXT1", (0, 1))
    assert result is None


def test_intergenic_avoids_all_genes(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    attempt = intergenic.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, _, _ = generate_until_satisfied(10, attempt, cfg.max_attempts_multiplier, label="t", show_progress=False)
    assert len(res) == 10
    _check_all(res, "intergenic")
    for c in res:
        assert idx.is_intergenic(c.chrom, c.start, c.end)


def test_intergenic_strand_is_unknown_not_a_fabricated_pick(mid_fixture):
    """Intergenic circRNA strand must be "." not a randomly picked "+"/"-"."""
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    attempt = intergenic.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, _, _ = generate_until_satisfied(10, attempt, cfg.max_attempts_multiplier, label="t", show_progress=False)
    assert len(res) == 10
    assert all(c.strand == "." for c in res)


def test_circle_length_bounds_enforced(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    cfg.min_circle_length = 150
    cfg.max_circle_length = 250
    loci = ann.get_transcripts(min_exons=2, protein_coding_only=True)
    attempt = ecirna.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, _, _ = generate_with_isoforms(5, rng, loci, attempt, cfg, label="t", show_progress=False)
    assert len(res) > 0, "expected at least some 2-exon (~200bp) candidates within [150,250]"
    for c in res:
        assert 150 <= c.mature_length <= 250


def test_zero_requested_returns_nothing(mid_fixture):
    ann, idx, genome, cfg, rng = _setup(mid_fixture)
    attempt = intergenic.make_attempt_fn(rng, genome, ann, idx, cfg)
    res, attempts, rej = generate_until_satisfied(0, attempt, cfg.max_attempts_multiplier, show_progress=False)
    assert res == [] and attempts == 0
