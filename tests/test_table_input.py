import pytest

from DNAflexpy.io import read_table


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_reads_sequence_and_value(tmp_path):
    p = write(tmp_path, "a.tsv", "ATGC\t1.5\nGGTT\t2.5\n")
    assert read_table(p) == [("seq_0", "ATGC", 1.5), ("seq_1", "GGTT", 2.5)]


def test_detects_a_header_row(tmp_path):
    p = write(tmp_path, "h.tsv", "sequence\taffinity\nATGC\t1.5\n")
    assert read_table(p) == [("seq_0", "ATGC", 1.5)]


def test_detects_absence_of_a_header(tmp_path):
    p = write(tmp_path, "n.tsv", "ATGC\t1.5\n")
    assert len(read_table(p)) == 1


def test_header_can_be_forced_off(tmp_path):
    """A first data row is not lost when header=False is explicit."""
    p = write(tmp_path, "f.tsv", "ATGC\t1.5\nGGTT\t2.5\n")
    assert len(read_table(p, header=False)) == 2


def test_columns_addressable_by_name(tmp_path):
    p = write(tmp_path, "c.tsv", "affinity\tsequence\n1.5\tATGC\n")
    assert read_table(p, seq_col="sequence", value_col="affinity") == [
        ("seq_0", "ATGC", 1.5)
    ]


def test_id_column_is_used_when_named(tmp_path):
    p = write(tmp_path, "i.tsv", "name\tsequence\taffinity\npeak1\tATGC\t1.5\n")
    assert read_table(p, seq_col="sequence", value_col="affinity", id_col="name") == [
        ("peak1", "ATGC", 1.5)
    ]


def test_csv_via_sep(tmp_path):
    p = write(tmp_path, "a.csv", "ATGC,1.5\n")
    assert read_table(p, sep=",") == [("seq_0", "ATGC", 1.5)]


def test_non_numeric_value_raises_naming_the_row(tmp_path):
    """Dropping the row instead would silently corrupt a training set."""
    p = write(tmp_path, "b.tsv", "ATGC\t1.5\nGGTT\thigh\n")
    with pytest.raises(ValueError, match="line 2"):
        read_table(p)


def test_missing_value_raises(tmp_path):
    p = write(tmp_path, "m.tsv", "ATGC\t1.5\nGGTT\t\n")
    with pytest.raises(ValueError, match="line 2"):
        read_table(p)


def test_out_of_range_column_raises(tmp_path):
    p = write(tmp_path, "o.tsv", "ATGC\t1.5\n")
    with pytest.raises(ValueError, match="value_col"):
        read_table(p, value_col=9)


def test_unknown_named_column_raises(tmp_path):
    p = write(tmp_path, "u.tsv", "sequence\taffinity\nATGC\t1.5\n")
    with pytest.raises(ValueError, match="nope"):
        read_table(p, seq_col="nope", value_col="affinity")


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_table(tmp_path / "nope.tsv")


def test_empty_file_raises(tmp_path):
    p = write(tmp_path, "e.tsv", "")
    with pytest.raises(ValueError, match="no data rows"):
        read_table(p)


def test_numeric_looking_header_is_still_detected(tmp_path):
    """A header whose label column is a number must not become a data row."""
    p = write(tmp_path, "nh.tsv", "sequence\t2024\nATGC\t1.5\n")
    assert read_table(p) == [("seq_0", "ATGC", 1.5)]


def test_error_line_number_accounts_for_the_header(tmp_path):
    """The number must match what the user sees in an editor."""
    p = write(tmp_path, "hl.tsv", "sequence\taffinity\nATGC\t1.5\nGGTT\thigh\n")
    with pytest.raises(ValueError, match="line 3"):
        read_table(p)


def test_header_only_file_raises(tmp_path):
    """Covers the frame.empty branch, which the 0-byte test does not reach."""
    p = write(tmp_path, "ho.tsv", "sequence\taffinity\n")
    with pytest.raises(ValueError, match="no data rows"):
        read_table(p)


def test_iupac_first_row_is_data_not_a_header(tmp_path):
    """A real sequence may carry ambiguity codes; dropping it would be silent."""
    p = write(tmp_path, "iu.tsv", "ATGR\t1.5\nGGTT\t2.5\n")
    with pytest.warns(UserWarning, match="non-ACGTN"):
        rows = read_table(p)
    assert len(rows) == 2
    assert rows[0] == ("seq_0", "ATGR", 1.5)


def test_warns_about_letters_the_tables_cannot_resolve(tmp_path):
    p = write(tmp_path, "w.tsv", "ATGY\t1.5\n")
    with pytest.warns(UserWarning, match="Y"):
        read_table(p)


def test_plain_n_does_not_warn(tmp_path):
    """N is an ordinary placeholder, not worth a warning."""
    import warnings as _w
    p = write(tmp_path, "n.tsv", "ATGN\t1.5\n")
    with _w.catch_warnings():
        _w.simplefilter("error")
        assert len(read_table(p)) == 1


def test_clean_sequences_do_not_warn(tmp_path):
    import warnings as _w
    p = write(tmp_path, "c2.tsv", "ATGC\t1.5\n")
    with _w.catch_warnings():
        _w.simplefilter("error")
        assert len(read_table(p)) == 1
