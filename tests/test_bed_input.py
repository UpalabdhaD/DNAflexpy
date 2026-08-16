import pytest

from DNAflexpy import FlexProfiler
from DNAflexpy.core import ProfileSet
from DNAflexpy.io import extract_intervals, read_bed
from DNAflexpy.results import FlexProfile

GENOME = ">chr1\n" + "ACGT" * 25 + "\n>chr2\n" + "AAAACCCCGGGGTTTT" + "\n"


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(GENOME)
    return path


@pytest.fixture
def genome_iupac(tmp_path):
    path = tmp_path / "iupac.fa"
    path.write_text(">chrI\nARGT\n")
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
    """A 2bp interval widened to 8 stays centred on the same midpoint.

    A length-only assertion cannot tell `start = centre - width // 2` (correct)
    apart from a bug like `start = start - width // 2` (shifts the window):
    on this "ACGT"*25 fixture, ref[7:15], ref[6:14] and ref[10:18] are all
    length 8, so pin the actual seqid and bases too.
    """
    out = extract_intervals([("chr1", 10, 12, None, "+")], genome, width=8)
    assert out == [("chr1:7-15", "TACGTACG")]


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


def test_contig_named_like_a_track_line_is_kept(tmp_path):
    """A prefix test would silently drop this real interval."""
    p = write_bed(tmp_path, "trackXYZ\t10\t20\nbrowser1\t30\t40\n")
    assert [r[0] for r in read_bed(p)] == ["trackXYZ", "browser1"]


def test_real_track_and_browser_lines_are_still_skipped(tmp_path):
    p = write_bed(tmp_path, 'track name="peaks"\nbrowser position chr1\nchr1\t10\t20\n')
    assert read_bed(p) == [("chr1", 10, 20, None, "+")]


def test_bad_line_reports_its_real_file_line(tmp_path):
    """The number must match what the user sees in an editor."""
    p = write_bed(tmp_path, "chr1\t10\t20\nchr1\t30\n")
    with pytest.raises(ValueError, match="line 2"):
        read_bed(p)


def test_header_row_raises_naming_the_file_and_line(tmp_path):
    """A column-header line ('chrom\\tstart\\tend') is not skipped -- only
    track/browser/# lines are -- so int('start') must raise something that
    names the file and line, not a bare int() ValueError."""
    p = write_bed(tmp_path, "chrom\tstart\tend\nchr1\t10\t20\n")
    with pytest.raises(ValueError, match="line 1"):
        read_bed(p)


def test_reversed_interval_raises(tmp_path):
    """start > end must not silently yield an empty extracted sequence."""
    p = write_bed(tmp_path, "chr1\t30\t10\n")
    with pytest.raises(ValueError, match="line 1"):
        read_bed(p)


def test_empty_interval_raises(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\t10\n")
    with pytest.raises(ValueError, match="line 1"):
        read_bed(p)


def test_minus_strand_complements_iupac_codes(genome_iupac):
    """R and Y must swap on the minus strand, not pass through unchanged."""
    plus = extract_intervals([("chrI", 0, 4, None, "+")], genome_iupac)[0][1]
    minus = extract_intervals([("chrI", 0, 4, None, "-")], genome_iupac)[0][1]
    assert plus == "ARGT"
    assert minus == "ACYT"
    assert plus != minus          # guards against a self-complementary fixture


def test_from_bed_returns_a_profile(genome, tmp_path):
    bed = write_bed(tmp_path, "chr1\t10\t30\tpeakA\t0\t+\n")
    prof = FlexProfiler("DNaseI", window_size=10).from_bed(bed, genome)
    assert isinstance(prof, FlexProfile)
    assert prof.seqids == ["peakA"]


def test_from_bed_with_width_gives_equal_length_rows(genome, tmp_path):
    """Fixed width is the point: rows line up for position-wise work."""
    bed = write_bed(tmp_path, "chr1\t10\t12\ta\t0\t+\nchr1\t40\t60\tb\t0\t+\n")
    prof = FlexProfiler("DNaseI", window_size=10).from_bed(bed, genome, width=20)
    assert prof.frame.notna().all().all()
    assert prof.frame.shape == (2, 11)


def test_from_bed_has_no_y(genome, tmp_path):
    """BED input carries no label column."""
    bed = write_bed(tmp_path, "chr1\t10\t30\ta\t0\t+\n")
    assert FlexProfiler("DNaseI", window_size=10).from_bed(bed, genome).y is None


def test_from_bed_multi_feature(genome, tmp_path):
    bed = write_bed(tmp_path, "chr1\t10\t30\ta\t0\t+\n")
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).from_bed(bed, genome)
    assert isinstance(ps, ProfileSet)
    assert set(ps) == {"DNaseI", "trx"}


def test_padded_bases_are_masked_not_zeroed(genome, tmp_path):
    """The N padding must show up as masked, never as a real 0 measurement."""
    bed = write_bed(tmp_path, "chr1\t0\t2\ta\t0\t+\n")
    with pytest.warns(UserWarning, match="padded"):
        prof = FlexProfiler("DNaseI", window_size=0).from_bed(
            bed, genome, width=20, on_edge="pad"
        )
    assert prof.n_masked["a"] > 0
