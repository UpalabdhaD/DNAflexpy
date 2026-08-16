import pytest

from DNAflexpy.io import extract_intervals, read_bed

GENOME = ">chr1\n" + "ACGT" * 25 + "\n>chr2\n" + "AAAACCCCGGGGTTTT" + "\n"


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(GENOME)
    return path


def write_bed(tmp_path, text):
    path = tmp_path / "peaks.bed"
    path.write_text(text)
    return path


def test_read_bed_parses_the_core_columns(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\t20\n")
    assert read_bed(p) == [("chr1", 10, 20, None, "+")]


def test_read_bed_reads_name_and_strand(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\t20\tpeakA\t0\t-\n")
    assert read_bed(p) == [("chr1", 10, 20, "peakA", "-")]


def test_read_bed_skips_track_and_comment_lines(tmp_path):
    p = write_bed(tmp_path, "# a comment\ntrack name=x\nchr1\t10\t20\n\n")
    assert len(read_bed(p)) == 1


def test_read_bed_rejects_a_short_line(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\n")
    with pytest.raises(ValueError, match="at least 3"):
        read_bed(p)


def test_extract_uses_the_interval_as_is_without_width(genome):
    out = extract_intervals([("chr1", 0, 4, None, "+")], genome)
    assert out == [("chr1:0-4", "ACGT")]


def test_extract_recentres_to_a_fixed_width(genome):
    """A 2bp interval widened to 8 stays centred on the same midpoint."""
    out = extract_intervals([("chr1", 10, 12, None, "+")], genome, width=8)
    assert len(out[0][1]) == 8


def test_extract_uses_the_bed_name_as_the_id(genome):
    out = extract_intervals([("chr1", 0, 4, "peakA", "+")], genome)
    assert out[0][0] == "peakA"


def test_minus_strand_is_reverse_complemented(genome):
    plus = extract_intervals([("chr2", 0, 4, None, "+")], genome)[0][1]
    minus = extract_intervals([("chr2", 0, 4, None, "-")], genome)[0][1]
    assert plus == "AAAA"
    assert minus == "TTTT"


def test_all_extracted_sequences_have_equal_length_with_width(genome):
    """Fixed width is what makes downstream position-wise comparison valid."""
    intervals = [("chr1", 10, 12, None, "+"), ("chr1", 40, 60, None, "+")]
    out = extract_intervals(intervals, genome, width=10)
    assert {len(s) for _, s in out} == {10}


def test_on_edge_drop_skips_and_warns(genome):
    """An interval at position 0 cannot be centred in 20bp."""
    intervals = [("chr1", 0, 2, None, "+"), ("chr1", 40, 42, None, "+")]
    with pytest.warns(UserWarning, match="dropped"):
        out = extract_intervals(intervals, genome, width=20, on_edge="drop")
    assert len(out) == 1


def test_on_edge_error_raises_naming_the_interval(genome):
    with pytest.raises(ValueError, match="chr1"):
        extract_intervals([("chr1", 0, 2, None, "+")], genome, width=20,
                          on_edge="error")


def test_on_edge_pad_pads_with_n_and_warns(genome):
    with pytest.warns(UserWarning, match="padded"):
        out = extract_intervals([("chr1", 0, 2, None, "+")], genome, width=20,
                                on_edge="pad")
    assert len(out[0][1]) == 20
    assert "N" in out[0][1]


def test_unknown_on_edge_raises(genome):
    with pytest.raises(ValueError, match="on_edge"):
        extract_intervals([("chr1", 0, 4, None, "+")], genome, on_edge="shrug")


def test_unknown_contig_raises(genome):
    with pytest.raises(ValueError, match="chr99"):
        extract_intervals([("chr99", 0, 4, None, "+")], genome)
