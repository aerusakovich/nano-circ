from nanocirc_sim.config import SimConfig
from nanocirc_sim.dbfit import fit_config_from_database


def test_dbfit_recovers_planted_composition(mid_fixture, synthetic_db_bed):
    cfg = fit_config_from_database(
        db_bed_path=synthetic_db_bed, genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"],
        base_config=SimConfig.default(),
    )
    # 6 examples of each of 5 types, planted in conftest.synthetic_db_bed
    for t in ["ecirna", "eicirna", "cirna", "antisense", "intergenic"]:
        assert abs(cfg.type_fractions[t] - 0.2) < 1e-9, cfg.type_fractions

    assert cfg.length_distributions["ecirna"]["median"] == 200  # 2x100bp exon blocks, intron excluded
    assert cfg.length_distributions["eicirna"]["median"] == 500  # exon+intron+exon, full span
    assert cfg.length_distributions["cirna"]["median"] == 299
    assert cfg.length_distributions["intergenic"]["median"] == 300


def test_dbfit_fits_terminal_exon_prob(mid_fixture, tmp_path):
    # G1 has 8 exons (0-based): 999-1099, 1399-1499, ..., 3799-3899 (see conftest.mid_fixture).
    lines = []
    for i in range(8):
        # touches first exon (999), not last: exon0+exon1
        lines.append(f"chrT\t999\t1499\tfirst_{i}\t0\t+\t999\t1499\t0\t2\t100,100,\t0,400,\n")
    for i in range(2):
        # touches last exon (3899), not first: exon6+exon7
        lines.append(f"chrT\t3399\t3899\tlast_{i}\t0\t+\t3399\t3899\t0\t2\t100,100,\t0,400,\n")
    for i in range(10):
        # fully internal: exon2+exon3
        lines.append(f"chrT\t1799\t2299\tinternal_{i}\t0\t+\t1799\t2299\t0\t2\t100,100,\t0,400,\n")
    bed = tmp_path / "terminal.bed12"
    bed.write_text("".join(lines))

    cfg = fit_config_from_database(
        db_bed_path=str(bed), genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"],
        base_config=SimConfig.default(),
    )
    # 8 first-touching + 2 last-touching out of 20 matched records, per side:
    # (8 + 2) / (2 * 20) == 0.25
    assert abs(cfg.terminal_exon_prob - 0.25) < 1e-9, cfg.terminal_exon_prob


def test_dbfit_bed6_only_database_does_not_corrupt_exon_counts(mid_fixture, tmp_path):
    # A BED6 database (no block columns at all, like data/mouse/union_bed.bed)
    # has no way to know real exon structure. Every record parses as "1
    # block", but that must NOT be fit into exon_count_distributions as "these
    # circRNAs are single-exon". It should fall back to the base config,
    # the same as having too few examples.
    lines = [f"chrT\t999\t1099\trec_{i}\t0\t+\n" for i in range(20)]  # plain BED6, exon0 span
    bed = tmp_path / "bed6_only.bed"
    bed.write_text("".join(lines))

    base = SimConfig.default()
    fitted = fit_config_from_database(
        db_bed_path=str(bed), genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"], base_config=base,
    )
    assert fitted.exon_count_distributions["ecirna"] == base.exon_count_distributions["ecirna"]
    # length/type fitting should still work fine from a BED6-only database
    assert fitted.type_fractions["ecirna"] == 1.0


def test_dbfit_antisense_requires_opposite_strand_EXON_not_just_intron_overlap(mid_fixture, tmp_path):
    """A region overlapping only an opposite-strand INTRON (no exon) must
    stay intergenic, not antisense."""
    lines = []
    for i in range(6):
        # No gene on '-' at this locus; overlaps G1's exon0 [999,1099) on
        # the opposite '+' strand: genuine antisense.
        lines.append(f"chrT\t999\t1099\tas_{i}\t0\t-\t999\t1099\t0\t1\t100,\t0,\n")
    for i in range(6):
        # No gene on '-' at this locus; overlaps only G1's INTRON
        # [1099,1399) on '+' (no exon touched). Must stay intergenic.
        lines.append(f"chrT\t1150\t1350\tintron_only_{i}\t0\t-\t1150\t1350\t0\t1\t200,\t0,\n")
    bed = tmp_path / "as_vs_intron.bed12"
    bed.write_text("".join(lines))

    cfg = fit_config_from_database(
        db_bed_path=str(bed), genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"],
        base_config=SimConfig.default(),
    )
    assert abs(cfg.type_fractions["antisense"] - 0.5) < 1e-9, cfg.type_fractions
    assert abs(cfg.type_fractions["intergenic"] - 0.5) < 1e-9, cfg.type_fractions


def test_dbfit_antisense_not_triggered_by_samestrand_gene_of_antisense_biotype(mid_fixture, tmp_path):
    """A same-strand exonic record on an 'antisense'-biotype gene must
    classify as ecirna, not antisense."""
    bed = tmp_path / "as_biotype.bed12"
    bed.write_text("".join(
        f"chrT\t12000\t12099\trec_{i}\t0\t-\t12000\t12099\t0\t1\t99,\t0,\n" for i in range(6)
    ))
    cfg = fit_config_from_database(
        db_bed_path=str(bed), genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"],
        base_config=SimConfig.default(),
    )
    assert cfg.type_fractions["ecirna"] == 1.0
    assert cfg.type_fractions["antisense"] == 0.0


def test_dbfit_keeps_base_config_when_too_few_examples(mid_fixture, tmp_path):
    # only 2 examples total -> below MIN_EXAMPLES_TO_FIT, everything should
    # fall back to the base config's values.
    bed = tmp_path / "sparse.bed12"
    bed.write_text(
        "chrT\t999\t1099\tone\t0\t+\t999\t1099\t0\t1\t100,\t0,\n"
        "chrT\t999\t1099\ttwo\t0\t+\t999\t1099\t0\t1\t100,\t0,\n"
    )
    base = SimConfig.default()
    fitted = fit_config_from_database(
        db_bed_path=str(bed), genome_path=mid_fixture["genome"], gtf_path=mid_fixture["gtf"], base_config=base,
    )
    assert fitted.length_distributions["ecirna"] == base.length_distributions["ecirna"]
    assert fitted.exon_count_distributions["ecirna"] == base.exon_count_distributions["ecirna"]
