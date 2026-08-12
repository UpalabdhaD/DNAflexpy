import math
import pathlib

from DNAflexpy.profile import FlexProfile


def make(rows):
    return FlexProfile(rows, feature="DNaseI", window_size=10, kmer_len=3)


def test_ragged_rows_pad_with_trailing_tabs(tmp_path):
    """Matches w10_mixed.tsv: a too-short sequence is ID plus trailing tabs."""
    out = tmp_path / "o.tsv"
    make([["long", 0.1, 0.2, 0.3], ["short"]]).to_tsv(out)
    lines = out.read_text().splitlines()
    assert lines[0] == "long\t0.1\t0.2\t0.3"
    assert lines[1] == "short\t\t\t"


def test_all_id_only_rows_have_no_trailing_tab(tmp_path):
    """Matches w30.tsv: one column means no padding at all."""
    out = tmp_path / "o.tsv"
    make([["a"], ["b"]]).to_tsv(out)
    assert out.read_text().splitlines() == ["a", "b"]


def test_trailing_zeros_are_dropped(tmp_path):
    out = tmp_path / "o.tsv"
    make([["s", round(0.0100, 3)]]).to_tsv(out)
    assert out.read_text().splitlines()[0] == "s\t0.01"


def test_masked_values_serialise_as_empty(tmp_path):
    out = tmp_path / "o.tsv"
    make([["s", 0.1, float("nan"), 0.3]]).to_tsv(out)
    assert out.read_text().splitlines()[0] == "s\t0.1\t\t0.3"


def test_n_masked_counts_masked_positions():
    p = make([["a", 0.1, float("nan")], ["b", float("nan"), float("nan")]])
    assert p.n_masked == {"a": 1, "b": 2}


def test_n_masked_ignores_ragged_padding():
    """Padding is NaN too, but it is not a masked measurement."""
    p = make([["long", 0.1, 0.2, 0.3], ["short"]])
    assert p.n_masked == {"long": 0, "short": 0}


def test_frame_is_indexed_by_seqid():
    p = make([["a", 0.1, 0.2], ["b", 0.3, 0.4]])
    assert list(p.frame.index) == ["a", "b"]
    assert p.frame.shape == (2, 2)


def test_seqids_preserve_input_order():
    p = make([["z", 0.1], ["a", 0.2]])
    assert p.seqids == ["z", "a"]


def test_tidy_frame_is_long():
    p = make([["a", 0.1, 0.2]])
    tidy = p.to_frame(tidy=True)
    assert list(tidy.columns) == ["seqid", "position", "value", "feature"]
    assert len(tidy) == 2
