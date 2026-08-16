import pytest

from DNAflexpy.io import read_table


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_reads_sequence_and_value(tmp_path):
    p = write(tmp_path, "a.tsv", "ATGC\t1.5\nGGTT\t2.5\n")
    assert read_table(p, header=False) == [("seq_0", "ATGC", 1.5), ("seq_1", "GGTT", 2.5)]


def test_header_can_be_forced_off(tmp_path):
    """A first data row is not lost when header=False is explicit."""
    p = write(tmp_path, "f.tsv", "ATGC\t1.5\nGGTT\t2.5\n")
    assert len(read_table(p, header=False)) == 2


def test_columns_addressable_by_name(tmp_path):
    p = write(tmp_path, "c.tsv", "affinity\tsequence\n1.5\tATGC\n")
    assert read_table(p, seq_col="sequence", value_col="affinity", header=True) == [
        ("seq_0", "ATGC", 1.5)
    ]


def test_id_column_is_used_when_named(tmp_path):
    p = write(tmp_path, "i.tsv", "name\tsequence\taffinity\npeak1\tATGC\t1.5\n")
    assert read_table(
        p, seq_col="sequence", value_col="affinity", id_col="name", header=True
    ) == [("peak1", "ATGC", 1.5)]


def test_csv_via_sep(tmp_path):
    p = write(tmp_path, "a.csv", "ATGC,1.5\n")
    assert read_table(p, sep=",", header=False) == [("seq_0", "ATGC", 1.5)]


def test_non_numeric_value_raises_naming_the_row(tmp_path):
    """Dropping the row instead would silently corrupt a training set."""
    p = write(tmp_path, "b.tsv", "ATGC\t1.5\nGGTT\thigh\n")
    with pytest.raises(ValueError, match="line 2"):
        read_table(p, header=False)


def test_missing_value_raises(tmp_path):
    p = write(tmp_path, "m.tsv", "ATGC\t1.5\nGGTT\t\n")
    with pytest.raises(ValueError, match="line 2"):
        read_table(p, header=False)


def test_out_of_range_column_raises(tmp_path):
    p = write(tmp_path, "o.tsv", "ATGC\t1.5\n")
    with pytest.raises(ValueError, match="value_col"):
        read_table(p, value_col=9, header=False)


def test_unknown_named_column_raises(tmp_path):
    p = write(tmp_path, "u.tsv", "sequence\taffinity\nATGC\t1.5\n")
    with pytest.raises(ValueError, match="nope"):
        read_table(p, seq_col="nope", value_col="affinity", header=True)


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_table(tmp_path / "nope.tsv")


def test_empty_file_raises(tmp_path):
    p = write(tmp_path, "e.tsv", "")
    with pytest.raises(ValueError, match="no data rows"):
        read_table(p, header=False)


def test_error_line_number_accounts_for_the_header(tmp_path):
    """The number must match what the user sees in an editor."""
    p = write(tmp_path, "hl.tsv", "sequence\taffinity\nATGC\t1.5\nGGTT\thigh\n")
    with pytest.raises(ValueError, match="line 3"):
        read_table(p, header=True)


def test_header_only_file_raises(tmp_path):
    """Covers the frame.empty branch, which the 0-byte test does not reach."""
    p = write(tmp_path, "ho.tsv", "sequence\taffinity\n")
    with pytest.raises(ValueError, match="no data rows"):
        read_table(p, header=True)


def test_iupac_first_row_is_data_not_a_header(tmp_path):
    """A real sequence may carry ambiguity codes; dropping it would be silent."""
    p = write(tmp_path, "iu.tsv", "ATGR\t1.5\nGGTT\t2.5\n")
    with pytest.warns(UserWarning, match="non-ACGTN"):
        rows = read_table(p, header=False)
    assert len(rows) == 2
    assert rows[0] == ("seq_0", "ATGR", 1.5)


def test_warns_about_letters_the_tables_cannot_resolve(tmp_path):
    p = write(tmp_path, "w.tsv", "ATGY\t1.5\n")
    with pytest.warns(UserWarning, match="Y"):
        read_table(p, header=False)


def test_plain_n_does_not_warn(tmp_path):
    """N is an ordinary placeholder, not worth a warning."""
    import warnings as _w
    p = write(tmp_path, "n.tsv", "ATGN\t1.5\n")
    with _w.catch_warnings():
        _w.simplefilter("error")
        assert len(read_table(p, header=False)) == 1


def test_clean_sequences_do_not_warn(tmp_path):
    import warnings as _w
    p = write(tmp_path, "c2.tsv", "ATGC\t1.5\n")
    with _w.catch_warnings():
        _w.simplefilter("error")
        assert len(read_table(p, header=False)) == 1


def test_header_must_be_given(tmp_path):
    """Guessing is unreliable: 'dna', 'rna' and 'tag' are all nucleotide letters."""
    p = write(tmp_path, "g.tsv", "ATGC\t1.5\n")
    with pytest.raises(ValueError, match="header=True or header=False"):
        read_table(p)


def test_header_true_skips_the_first_row(tmp_path):
    p = write(tmp_path, "ht.tsv", "sequence\taffinity\nATGC\t1.5\n")
    assert read_table(p, header=True) == [("seq_0", "ATGC", 1.5)]


def test_header_false_keeps_the_first_row(tmp_path):
    """The case that used to be guessed wrong: a header made of DNA letters."""
    p = write(tmp_path, "hf.tsv", "TAG\t1.5\nATGC\t2.5\n")
    assert read_table(p, header=False) == [("seq_0", "TAG", 1.5), ("seq_1", "ATGC", 2.5)]


import numpy as np

from DNAflexpy import FlexProfiler
from DNAflexpy.core import ProfileSet
from DNAflexpy.results import FlexProfile

SEQ_A = "ATGCGTACGTAGCTAGCGTAGCTAGT"
SEQ_B = "CGTAGCTAGTACGATCGTACGTAGCT"


def test_profile_table_returns_a_profile_with_y(tmp_path):
    p = write(tmp_path, "t.tsv", f"{SEQ_A}\t1.5\n{SEQ_B}\t2.5\n")
    prof = FlexProfiler("DNaseI", window_size=10).profile_table(p, header=False)
    assert isinstance(prof, FlexProfile)
    assert prof.y == [1.5, 2.5]
    assert prof.seqids == ["seq_0", "seq_1"]


def test_y_is_aligned_to_seqids(tmp_path):
    p = write(tmp_path, "t.tsv", f"name\tseq\tv\na\t{SEQ_A}\t9.0\nb\t{SEQ_B}\t8.0\n")
    prof = FlexProfiler("DNaseI", window_size=10).profile_table(
        p, seq_col="seq", value_col="v", id_col="name", header=True
    )
    assert dict(zip(prof.seqids, prof.y)) == {"a": 9.0, "b": 8.0}


def test_profile_seqs_has_no_y():
    """y is only meaningful for labelled input."""
    prof = FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ_A])
    assert prof.y is None


def test_profile_fasta_has_no_y(tmp_path):
    """The real FASTA path, serial branch."""
    fa = tmp_path / "small.fa"
    fa.write_text(f">a\n{SEQ_A}\n>b\n{SEQ_B}\n")
    prof = FlexProfiler("DNaseI", window_size=10).profile_fasta(fa, threads=1)
    assert prof.seqids == ["a", "b"]
    assert prof.y is None


def test_profile_fasta_pooled_has_no_y(tmp_path):
    """The pooled branch, which builds its rows separately from _build."""
    fa = tmp_path / "many.fa"
    fa.write_text("".join(f">s{i}\n{SEQ_A}\n" for i in range(70)))
    prof = FlexProfiler("DNaseI", window_size=10).profile_fasta(fa, threads=2)
    assert len(prof.seqids) == 70
    assert prof.y is None


def test_values_match_the_plain_sequence_path(tmp_path):
    """Reading via a table must not change the numbers."""
    p = write(tmp_path, "t.tsv", f"{SEQ_A}\t1.5\n")
    viafile = FlexProfiler("DNaseI", window_size=10).profile_table(p, header=False)
    direct = FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ_A])
    assert np.array_equal(viafile.frame.to_numpy(), direct.frame.to_numpy())


def test_multi_feature_table_gives_every_profile_the_same_y(tmp_path):
    p = write(tmp_path, "t.tsv", f"{SEQ_A}\t1.5\n{SEQ_B}\t2.5\n")
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).profile_table(p, header=False)
    assert isinstance(ps, ProfileSet)
    assert ps["DNaseI"].y == [1.5, 2.5]
    assert ps["trx"].y == [1.5, 2.5]
