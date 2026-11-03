import pytest

from nanocirc_sim.bed12 import bed12_fields
from nanocirc_sim.blocks import (
    Block, block_signature, length_in_bounds, total_length, validate_blocks, with_shifted_boundary,
)
from nanocirc_sim.candidate import Candidate


def test_length_in_bounds():
    blocks = [Block(100, 150, "exon")]  # length 50
    assert length_in_bounds(blocks, 50, 100) is True
    assert length_in_bounds(blocks, 51, 100) is False
    assert length_in_bounds(blocks, 10, 49) is False


def test_validate_blocks_rejects_overlap():
    with pytest.raises(ValueError):
        validate_blocks([Block(0, 100, "exon"), Block(50, 150, "exon")])
    validate_blocks([Block(0, 100, "exon"), Block(100, 200, "exon")])  # touching is fine


def test_with_shifted_boundary_single_block():
    blocks = [Block(100, 200, "region")]
    out = with_shifted_boundary(blocks, 95, 205)
    assert out == [Block(95, 205, "region")]
    assert with_shifted_boundary(blocks, 210, 205) is None  # inverted


def test_with_shifted_boundary_multi_block_rejects_invasion():
    blocks = [Block(100, 200, "exon"), Block(300, 400, "exon")]
    out = with_shifted_boundary(blocks, 90, 410)
    assert out[0] == Block(90, 200, "exon")
    assert out[-1] == Block(300, 410, "exon")
    # shifting start past the first block's own end must be rejected
    assert with_shifted_boundary(blocks, 250, 410) is None


def test_block_signature_dedup():
    a = [Block(0, 100, "exon"), Block(200, 300, "exon")]
    b = [Block(0, 100, "exon"), Block(200, 300, "exon")]
    c = [Block(0, 100, "exon"), Block(200, 301, "exon")]
    assert block_signature(a) == block_signature(b)
    assert block_signature(a) != block_signature(c)


def test_bed12_fields_multi_block():
    blocks = [Block(1000, 1100, "exon"), Block(1400, 1500, "intron_retained"), Block(1800, 1900, "exon")]
    c = Candidate(circ_type="eicirna", chrom="chr1", start=1000, end=1900, strand="+", blocks=blocks,
                  mature_seq="A" * total_length(blocks), id="test_0")
    c.tpm = 50.0
    fields = bed12_fields(c)
    assert fields[0:4] == ["chr1", "1000", "1900", "test_0"]
    assert fields[5] == "+"
    assert fields[9] == "3"
    assert fields[10] == "100,100,100,"
    assert fields[11] == "0,400,800,"
