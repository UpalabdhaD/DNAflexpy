"""Phase 4: turning profiles into a machine-learning design matrix."""
import warnings

import numpy as np
import pytest

from DNAflexpy import FlexProfiler
from DNAflexpy.encode import _one_hot

SEQS = {"s1": "ACGTACGT", "s2": "TTTTTTTT"}


# --- one-hot sequence encoding ---------------------------------------------


def test_one_hot_1mer_places_each_base_in_its_own_column():
    X, cols = _one_hot(["ACGT"], 1)
    assert cols == [
        "seq.1mer.p1.A", "seq.1mer.p1.C", "seq.1mer.p1.G", "seq.1mer.p1.T",
        "seq.1mer.p2.A", "seq.1mer.p2.C", "seq.1mer.p2.G", "seq.1mer.p2.T",
        "seq.1mer.p3.A", "seq.1mer.p3.C", "seq.1mer.p3.G", "seq.1mer.p3.T",
        "seq.1mer.p4.A", "seq.1mer.p4.C", "seq.1mer.p4.G", "seq.1mer.p4.T",
    ]
    assert X.shape == (1, 16)
    # ACGT walks the alphabet, so the hot column moves by 5 each position.
    # A transposed or reversed layout gives a different pattern and fails.
    assert list(X[0]) == [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]


def test_one_hot_rows_follow_the_input_order():
    X, _ = _one_hot(["ACGT", "TGCA"], 1)
    assert X.shape == (2, 16)
    assert list(X[1][:4]) == [0, 0, 0, 1]      # T first
    assert list(X[1][-4:]) == [1, 0, 0, 0]     # A last


def test_one_hot_2mer_uses_overlapping_positions():
    X, cols = _one_hot(["ACGT"], 2)
    assert X.shape == (1, 16 * 3)              # L - k + 1 == 3 positions
    assert cols[0] == "seq.2mer.p1.AA"
    assert cols[-1] == "seq.2mer.p3.TT"
    hot = [cols[i] for i in np.flatnonzero(X[0])]
    assert hot == ["seq.2mer.p1.AC", "seq.2mer.p2.CG", "seq.2mer.p3.GT"]


def test_one_hot_3mer_has_64_columns_per_position():
    X, cols = _one_hot(["ACGTAC"], 3)
    assert X.shape == (1, 64 * 4)
    assert cols[0] == "seq.3mer.p1.AAA"
    hot = [cols[i] for i in np.flatnonzero(X[0])]
    assert hot == ["seq.3mer.p1.ACG", "seq.3mer.p2.CGT",
                   "seq.3mer.p3.GTA", "seq.3mer.p4.TAC"]


def test_one_hot_is_case_insensitive():
    upper, _ = _one_hot(["ACGT"], 1)
    lower, _ = _one_hot(["acgt"], 1)
    assert np.array_equal(upper, lower)


def test_one_hot_is_float_so_concatenation_does_not_promote():
    X, _ = _one_hot(["ACGT"], 1)
    assert X.dtype == np.float64


def test_a_kmer_with_n_becomes_an_all_zero_position():
    with pytest.warns(UserWarning, match="all-zero"):
        X, _ = _one_hot(["ANGT"], 1)
    assert list(X[0][4:8]) == [0, 0, 0, 0]     # position 2 is unknown
    assert list(X[0][0:4]) == [1, 0, 0, 0]     # position 1 is unaffected


def test_an_ambiguous_code_masks_every_kmer_covering_it():
    """At k=2 one bad letter blanks two overlapping positions, not one."""
    with pytest.warns(UserWarning, match="all-zero"):
        X, _ = _one_hot(["ACYT"], 2)
    assert X[0][0:16].sum() == 1.0             # p1 == AC, fine
    assert X[0][16:32].sum() == 0.0            # p2 == CY
    assert X[0][32:48].sum() == 0.0            # p3 == YT


def test_a_clean_sequence_does_not_warn():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _one_hot(["ACGT"], 1)


def test_one_hot_rejects_k_larger_than_the_sequence():
    with pytest.raises(ValueError, match="k=5"):
        _one_hot(["ACGT"], 5)


def test_one_hot_rejects_k_below_one():
    with pytest.raises(ValueError, match="k must be >= 1"):
        _one_hot(["ACGT"], 0)


def test_one_hot_refuses_mixed_lengths_rather_than_truncating():
    with pytest.raises(ValueError, match=r"\[2, 4\]"):
        _one_hot(["ACGT", "AC"], 1)
