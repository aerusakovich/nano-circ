import pysam

from nanocirc_sim.bsj_check import (
    _read_fasta,
    build_doubled_reference,
    parse_alignments_for_bsj,
    summarize_per_circrna,
    write_filtered_fastq,
)


def _write_sam(path, reads):
    """reads: list of (name, ref_idx, start, length)."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 400, "SN": "short_circ"}, {"LN": 10000, "SN": "long_circ"}]}
    with pysam.AlignmentFile(path, "wh", header=header) as out:
        for name, ref_idx, start, length in reads:
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * length
            a.flag = 0
            a.reference_id = ref_idx
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = [(0, length)]
            out.write(a)


def test_build_doubled_reference(tmp_path):
    out = tmp_path / "doubled.fa"
    build_doubled_reference({"a": "ACGTACGT", "b": "TTTT"}, str(out))
    loaded = _read_fasta(str(out))
    assert loaded == {"a": "ACGTACGTACGTACGT", "b": "TTTTTTTT"}


def test_bsj_confirmation_matches_junction_crossing(tmp_path):
    mature_lengths = {"short_circ": 200, "long_circ": 5000}
    sam_path = tmp_path / "aligned.sam"
    _write_sam(str(sam_path), [
        ("read1_circ", 0, 50, 300),    # short_circ: [50,350) crosses 200 -> confirmed
        ("read2_circ", 0, 0, 150),     # short_circ: [0,150) does not cross 200 -> NOT confirmed
        ("read3_circ", 1, 0, 3000),    # long_circ: [0,3000) does not cross 5000 -> NOT confirmed (the bug scenario)
        ("read4_circ", 1, 4900, 300),  # long_circ: [4900,5200) crosses 5000 -> confirmed
    ])
    per_read = parse_alignments_for_bsj(str(sam_path), mature_lengths)

    assert per_read["read1_circ"]["bsj_confirmed"] is True
    assert per_read["read2_circ"]["bsj_confirmed"] is False
    assert per_read["read3_circ"]["bsj_confirmed"] is False
    assert per_read["read4_circ"]["bsj_confirmed"] is True

    per_circ = summarize_per_circrna(per_read)
    assert per_circ["short_circ"] == {"n_reads_aligned": 2, "n_reads_confirmed": 1}
    assert per_circ["long_circ"] == {"n_reads_aligned": 2, "n_reads_confirmed": 1}
