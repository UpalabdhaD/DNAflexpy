"""Tests for the `DNAflexpy` command.

These call `main(argv)` directly rather than shelling out. That keeps them
fast, and it covers the bash examples in the docs, which
`scripts/check_doc_examples.py` cannot run - it only executes ```python
blocks.
"""
import numpy as np
import pytest

from DNAflexpy import FlexProfiler
from DNAflexpy.cli import _features_in, _infer_format, main

SEQ = "ATGCGTACGTAGCTAGCGTAGCTAGT"


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "in.fa"
    path.write_text(f">a\n{SEQ}\n>b\n{SEQ}\n")
    return path


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    return path


# --- format inference -------------------------------------------------------


@pytest.mark.parametrize("name,expected", [
    ("x.fa", "fasta"), ("x.fasta", "fasta"), ("x.fna", "fasta"),
    ("x.FA", "fasta"), ("x.fa.gz", "fasta"),
    ("x.bed", "bed"),
    ("x.tsv", "table"), ("x.csv", "table"), ("x.txt", "table"),
])
def test_format_is_inferred_from_the_extension(name, expected):
    from pathlib import Path

    assert _infer_format(Path(name)) == expected


def test_an_unknown_extension_asks_for_format():
    from pathlib import Path

    with pytest.raises(SystemExit, match="--format"):
        _infer_format(Path("mystery.dat"))


# --- profile ----------------------------------------------------------------


def test_profile_writes_the_archive_named_file(fasta, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["profile", str(fasta), "--feature", "DNaseI",
                 "--window-size", "10"]) == 0
    assert (tmp_path / "in_w10nt_DNaseI.tsv").exists()


def test_the_written_bytes_are_exactly_to_tsv(fasta, tmp_path):
    """The CLI must not have a serialiser of its own: to_tsv is what the
    230-case byte-equality suite covers."""
    out = tmp_path / "cli.tsv"
    assert main(["profile", str(fasta), "--feature", "DNaseI",
                 "--outfile", str(out)]) == 0

    expected = tmp_path / "direct.tsv"
    FlexProfiler("DNaseI", window_size=10).profile_fasta(fasta).to_tsv(expected)
    assert out.read_bytes() == expected.read_bytes()


def test_window_size_reaches_the_profiler(fasta, tmp_path):
    out = tmp_path / "w0.tsv"
    main(["profile", str(fasta), "--window-size", "0", "--outfile", str(out)])
    expected = tmp_path / "direct.tsv"
    FlexProfiler("DNaseI", window_size=0).profile_fasta(fasta).to_tsv(expected)
    assert out.read_bytes() == expected.read_bytes()


def test_several_features_write_one_file_each(fasta, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["profile", str(fasta), "--feature", "DNaseI", "gc",
                 "--window-size", "5"]) == 0
    assert (tmp_path / "in_w5nt_DNaseI.tsv").exists()
    assert (tmp_path / "in_w5nt_gc.tsv").exists()


def test_outfile_with_several_features_is_refused(fasta, tmp_path):
    with pytest.raises(SystemExit, match="one file per feature"):
        main(["profile", str(fasta), "--feature", "DNaseI", "gc",
              "--outfile", str(tmp_path / "both.tsv")])


def test_an_unknown_feature_lists_the_real_ones(fasta):
    with pytest.raises(SystemExit, match="DNaseI"):
        main(["profile", str(fasta), "--feature", "nosuch"])


def test_a_missing_input_file_says_so(tmp_path):
    with pytest.raises(SystemExit, match="not found"):
        main(["profile", str(tmp_path / "absent.fa")])


def test_seq_prints_to_stdout(capsys):
    assert main(["profile", "--seq", SEQ, "--window-size", "10"]) == 0
    fields = capsys.readouterr().out.strip().split("\t")
    expected = FlexProfiler("DNaseI", window_size=10).profile(SEQ)
    assert fields[0] == "sequence"
    assert len(fields) - 1 == len(expected)
    assert float(fields[1]) == expected[0]


def test_stdout_uses_the_same_serialiser_as_the_file(capsys, tmp_path):
    """str(nan) is 'nan' but to_tsv writes an empty field. A hand-rolled
    print here would be a second, untested format for the same values."""
    masked = "ATGNGTACGT"
    assert main(["profile", "--seq", masked, "--window-size", "0"]) == 0
    printed = capsys.readouterr().out

    out = tmp_path / "same.tsv"
    FlexProfiler("DNaseI", window_size=0).profile_seqs(
        {"sequence": masked}).to_tsv(out)
    assert printed == out.read_text()
    assert "nan" not in printed


def test_an_input_file_and_seq_together_are_refused(fasta):
    with pytest.raises(SystemExit, match="not both"):
        main(["profile", str(fasta), "--seq", SEQ])


def test_no_input_at_all_is_refused():
    with pytest.raises(SystemExit, match="needs an input file"):
        main(["profile"])


# --- table input ------------------------------------------------------------


@pytest.fixture
def table(tmp_path):
    path = tmp_path / "affinity.tsv"
    path.write_text(f"{SEQ}\t1.5\n{SEQ}\t2.5\n")
    return path


def test_a_table_without_the_header_flag_is_refused(table):
    """Phase 3 decided header-ness cannot be guessed. The CLI must not
    quietly undo that by defaulting."""
    with pytest.raises(SystemExit, match="--header"):
        main(["profile", str(table)])


def test_no_header_profiles_the_table(table, tmp_path):
    out = tmp_path / "t.tsv"
    assert main(["profile", str(table), "--no-header",
                 "--outfile", str(out)]) == 0
    assert out.read_text().count("\n") == 2


def test_a_header_row_is_honoured(tmp_path):
    path = tmp_path / "labelled.tsv"
    path.write_text(f"sequence\taffinity\n{SEQ}\t1.5\n")
    out = tmp_path / "t.tsv"
    assert main(["profile", str(path), "--header", "--seq-col", "sequence",
                 "--value-col", "affinity", "--outfile", str(out)]) == 0
    assert out.read_text().count("\n") == 1


# --- BED input --------------------------------------------------------------


def test_bed_without_a_genome_is_refused(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t10\t30\n")
    with pytest.raises(SystemExit, match="--genome"):
        main(["profile", str(bed)])


def test_bed_with_a_genome_profiles_the_intervals(tmp_path, genome):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t10\t30\tpeakA\nchr1\t50\t70\tpeakB\n")
    out = tmp_path / "peaks.tsv"
    assert main(["profile", str(bed), "--genome", str(genome),
                 "--width", "20", "--outfile", str(out)]) == 0
    assert out.read_text().startswith("peakA\t")


# --- encode -----------------------------------------------------------------


def test_features_in_derives_the_lookup_features():
    assert _features_in(["1-mer", "1-DNaseI", "2-DNaseI", "3-gc"]) == ["DNaseI", "gc"]
    assert _features_in(["1-mer", "2-mer"]) == []


def test_encode_writes_a_loadable_npz(fasta, tmp_path):
    out = tmp_path / "X.npz"
    assert main(["encode", str(fasta), "--features", "1-mer", "1-DNaseI",
                 "--window-size", "0", "--out", str(out)]) == 0

    loaded = np.load(out, allow_pickle=False)
    expected = FlexProfiler("DNaseI", window_size=0).profile_fasta(fasta).encode(
        ["1-mer", "1-DNaseI"])
    np.testing.assert_allclose(loaded["X"], expected.X)
    assert list(loaded["columns"]) == expected.columns
    assert list(loaded["seqids"]) == ["a", "b"]
    assert "y" not in loaded


def test_encode_of_a_labelled_table_carries_y(table, tmp_path):
    out = tmp_path / "X.npz"
    assert main(["encode", str(table), "--no-header", "--features", "1-DNaseI",
                 "--window-size", "0", "--out", str(out)]) == 0
    loaded = np.load(out)
    np.testing.assert_allclose(loaded["y"], [1.5, 2.5])
    assert loaded["X"].shape[0] == 2


def test_no_normalize_reaches_encode(fasta, tmp_path):
    raw, scaled = tmp_path / "raw.npz", tmp_path / "scaled.npz"
    main(["encode", str(fasta), "--features", "1-wedge", "--window-size", "0",
          "--out", str(raw), "--no-normalize"])
    main(["encode", str(fasta), "--features", "1-wedge", "--window-size", "0",
          "--out", str(scaled)])
    assert np.load(raw)["X"].max() > 1.0
    assert np.load(scaled)["X"].max() == 1.0


def test_a_sequence_only_request_still_works(fasta, tmp_path):
    """There is no feature to profile, but the profile is what holds the
    sequences, so one still has to be built."""
    out = tmp_path / "X.npz"
    assert main(["encode", str(fasta), "--features", "1-mer",
                 "--out", str(out)]) == 0
    assert np.load(out)["X"].shape == (2, 4 * len(SEQ))


def test_a_malformed_encode_request_is_refused(fasta, tmp_path):
    with pytest.raises(SystemExit, match="malformed"):
        main(["encode", str(fasta), "--features", "mer",
              "--out", str(tmp_path / "X.npz")])


# --- top level --------------------------------------------------------------


def test_no_subcommand_prints_help_and_fails(capsys):
    assert main([]) == 1
    assert "profile" in capsys.readouterr().out


def test_citation_prints_a_template_and_flags_it(capsys):
    assert main(["--citation"]) == 0
    captured = capsys.readouterr()
    assert "@software{DNAflexpy" in captured.out
    assert "ADD_AUTHOR_LIST" in captured.out
    assert "placeholder" in captured.err


def test_version_is_the_package_version(capsys):
    from DNAflexpy import __version__

    with pytest.raises(SystemExit):
        main(["--version"])
    assert __version__ in capsys.readouterr().out


def test_there_is_no_plot_subcommand_yet():
    """A subcommand that exists and errors is worse than one that does not."""
    with pytest.raises(SystemExit):
        main(["plot", "out.tsv"])
