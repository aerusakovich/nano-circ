"""Regression tests for splice_sites.py's +/- strand donor/acceptor side assignment."""
import numpy as np
import pytest

from nanocirc_sim.genome import GenomeLoader
from nanocirc_sim.rng import RNGFactory
from nanocirc_sim.splice_sites import (
    classify_splice_type, donor_acceptor, find_boundary_for_desired_type, is_canonical,
    shift_all_junctions_independently,
)

_GT_AG_ONLY = {"GT-AG": 1.0, "GC-AG": 0.0, "AT-AC": 0.0, "other": 0.0}


@pytest.fixture
def splice_genome(tmp_path):
    """chrS with two GT-AG-flanked regions: [100,160) on '+' and [300,360) on '-'."""
    length = 500
    seq = list("TTTT" * (length // 4))
    seq[98:100] = list("AG")
    seq[160:162] = list("GT")
    seq[298:300] = list("AC")
    seq[360:362] = list("CT")
    fasta = tmp_path / "chrS.fa"
    fasta.write_text(f">chrS\n{''.join(seq)}\n")
    return GenomeLoader(str(fasta))


def test_donor_acceptor_plus_strand(splice_genome):
    donor, acceptor = donor_acceptor(splice_genome, "chrS", 100, 160, "+")
    assert (donor, acceptor) == ("GT", "AG")
    assert classify_splice_type(donor, acceptor) == "GT-AG"
    assert is_canonical(classify_splice_type(donor, acceptor))


def test_donor_acceptor_minus_strand(splice_genome):
    donor, acceptor = donor_acceptor(splice_genome, "chrS", 300, 360, "-")
    assert (donor, acceptor) == ("GT", "AG")
    assert classify_splice_type(donor, acceptor) == "GT-AG"


def test_donor_acceptor_minus_strand_uses_opposite_sides_from_plus(splice_genome):
    """Minus strand must read donor from the start-side flank, acceptor
    from the end-side flank: the opposite of plus strand."""
    donor_plus_sides, acceptor_plus_sides = donor_acceptor(splice_genome, "chrS", 300, 360, "+")
    donor_minus, acceptor_minus = donor_acceptor(splice_genome, "chrS", 300, 360, "-")
    assert (donor_minus, acceptor_minus) != (donor_plus_sides, acceptor_plus_sides)
    assert donor_minus == "GT" and acceptor_minus == "AG"


def test_find_boundary_already_gt_ag_needs_no_shift(splice_genome):
    new_start, new_end, achieved = find_boundary_for_desired_type(
        splice_genome, "chrS", 100, 160, "+", "GT-AG")
    assert (new_start, new_end) == (100, 160)
    assert achieved == "GT-AG"

    new_start, new_end, achieved = find_boundary_for_desired_type(
        splice_genome, "chrS", 300, 360, "-", "GT-AG")
    assert (new_start, new_end) == (300, 360)
    assert achieved == "GT-AG"


def test_find_boundary_shifts_to_reach_desired_type(splice_genome):
    """A boundary 3bp off from the true GT-AG site should shift back."""
    new_start, new_end, achieved = find_boundary_for_desired_type(
        splice_genome, "chrS", 103, 157, "+", "GT-AG", max_adjustment=10)
    assert (new_start, new_end) == (100, 160)
    assert achieved == "GT-AG"


def test_real_database_boundaries_are_majority_gt_ag_on_both_strands(tmp_path):
    """Sanity check with a hand-built two-gene fixture (one +, one -), each
    with a real internal GT-AG intron."""
    length = 2000
    seq = list("TTTT" * (length // 4))
    # + strand gene: exon[200..300) intron[300..500) GT...AG exon[500..700)
    seq[300:302] = list("GT")
    seq[498:500] = list("AG")
    # - strand gene: exon[1000..1100) intron[1100..1300) exon[1300..1500);
    # transcript direction is high->low coordinate, so donor(GT) sits at
    # the high-coordinate end of the intron and acceptor(AG) at the low end.
    seq[1298:1300] = list("AC")  # rc -> GT (donor, transcript-5' end of intron)
    seq[1100:1102] = list("CT")  # rc -> AG (acceptor, transcript-3' end of intron)
    fasta = tmp_path / "g.fa"
    fasta.write_text(f">chrG\n{''.join(seq)}\n")
    genome = GenomeLoader(str(fasta))

    # + gene: intron[300,500) with donor GT at its 5' (low-coord) end and
    # acceptor AG at its 3' (high-coord) end. A circRNA of exon0=[200,300)
    # alone has a real flanking intron only on its end side (300) -> donor;
    # a circRNA of exon1=[500,700) alone has one only on its start side
    # (500) -> acceptor.
    donor, _ = donor_acceptor(genome, "chrG", 200, 300, "+")
    assert donor == "GT"
    _, acceptor = donor_acceptor(genome, "chrG", 500, 700, "+")
    assert acceptor == "AG"

    # - gene: exon_high=[1300,1500) is transcript-5' (transcription runs
    # high->low coordinate on minus strand), exon_low=[1000,1100) is
    # transcript-3'. The intron [1100,1300) between them has its donor
    # (GT, transcript-5' end) at the high-coordinate edge (1300) and its
    # acceptor (AG, transcript-3' end) at the low-coordinate edge (1100).
    # A circRNA of exon_high alone has a real flanking intron only on its
    # start side (1300) -> should read as donor GT; a circRNA of exon_low
    # alone has a real flanking intron only on its end side (1100) -> should
    # read as acceptor AG.
    donor_high, _ = donor_acceptor(genome, "chrG", 1300, 1500, "-")
    assert donor_high == "GT"
    _, acceptor_low = donor_acceptor(genome, "chrG", 1000, 1100, "-")
    assert acceptor_low == "AG"


def test_shift_all_junctions_single_block_matches_find_boundary_directly(splice_genome):
    """n=1 must behave identically to calling find_boundary_for_desired_type directly."""
    rng = RNGFactory(1).generator
    result = shift_all_junctions_independently(
        splice_genome, "chrS", "+", [100], [160], rng, _GT_AG_ONLY, max_adjustment=10)
    assert result is not None
    starts, ends, achieved = result
    assert (starts, ends) == ([100], [160])
    assert achieved == "GT-AG"


@pytest.fixture
def multi_junction_genome(tmp_path):
    """chrJ: 3 exons e0=[100,200) e1=[250,300) e2=[400,450), outer BSJ and
    both introns already GT-AG-canonical."""
    length = 600
    seq = list("CCCC" * (length // 4))
    seq[98:100] = list("AG")   # outer BSJ acceptor (exon0's start-side)
    seq[200:202] = list("GT")  # i0 donor (exon0's end-side)
    seq[248:250] = list("AG")  # i0 acceptor (exon1's start-side)
    seq[300:302] = list("GT")  # i1 donor (exon1's end-side)
    seq[398:400] = list("AG")  # i1 acceptor (exon2's start-side)
    seq[450:452] = list("GT")  # outer BSJ donor (exon2's end-side)
    fasta = tmp_path / "chrJ.fa"
    fasta.write_text(f">chrJ\n{''.join(seq)}\n")
    return GenomeLoader(str(fasta))


def test_shift_all_junctions_all_already_canonical_needs_no_shift(multi_junction_genome):
    rng = RNGFactory(1).generator
    result = shift_all_junctions_independently(
        multi_junction_genome, "chrJ", "+", [100, 250, 400], [200, 300, 450],
        rng, _GT_AG_ONLY, max_adjustment=10)
    assert result is not None
    starts, ends, achieved = result
    assert starts == [100, 250, 400]
    assert ends == [200, 300, 450]
    assert achieved == "GT-AG"


def test_shift_all_junctions_moves_internal_boundary_independently_to_find_canonical_motif(tmp_path):
    """An internal junction that isn't already canonical gets its own
    independent search and shift, just like the outer BSJ."""
    length = 500
    seq = list("CCCC" * (length // 4))
    # outer BSJ: already canonical at the unshifted position, no shift needed
    seq[98:100] = list("AG")
    seq[300:302] = list("GT")
    # internal junction (between exon0's unshifted end=200 and exon1's
    # unshifted start=250): canonical motifs planted a few bp inside,
    # NOT at the unshifted positions themselves.
    seq[204:206] = list("GT")  # donor, 4bp after exon0's unshifted end
    seq[246:248] = list("AG")  # acceptor, 2bp before exon1's unshifted start
    fasta = tmp_path / "chrK.fa"
    fasta.write_text(f">chrK\n{''.join(seq)}\n")
    genome = GenomeLoader(str(fasta))

    rng = RNGFactory(1).generator
    result = shift_all_junctions_independently(
        genome, "chrK", "+", [100, 250], [200, 300], rng, _GT_AG_ONLY, max_adjustment=10)
    assert result is not None
    starts, ends, achieved = result
    assert starts[0] == 100 and ends[-1] == 300  # outer BSJ: unchanged, already canonical
    assert ends[0] == 204  # internal donor-side: shifted from 200 to the planted GT
    assert starts[1] == 248  # internal acceptor-side: shifted from 250 to the planted AG
    assert achieved == "GT-AG"


def test_shift_all_junctions_returns_none_when_shift_inverts_a_short_internal_intron(tmp_path):
    """An internal intron too short for two independent searches, causing
    an inverted/overlapping shift, must be rejected (None)."""
    length = 400
    seq = list("CCCC" * (length // 4))
    # exon0=[100,200), exon1=[206,300): only a 6bp internal intron gap.
    seq[209:211] = list("GT")  # donor-side match at offset +9 from center 200
    seq[201:203] = list("AG")  # acceptor-side match at offset -3 from center 206
    fasta = tmp_path / "chrL.fa"
    fasta.write_text(f">chrL\n{''.join(seq)}\n")
    genome = GenomeLoader(str(fasta))

    rng = RNGFactory(1).generator
    result = shift_all_junctions_independently(
        genome, "chrL", "+", [100, 206], [200, 300], rng, _GT_AG_ONLY, max_adjustment=10)
    assert result is None
