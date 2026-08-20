"""Phase 4: turning profiles into a machine-learning design matrix."""
import warnings

import numpy as np
import pytest

from DNAflexpy import FlexProfiler
from DNAflexpy.encode import _feature_block, _minmax, _one_hot

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


# --- nth-order feature blocks ----------------------------------------------
# `wedge` is used throughout because its range (0.9 .. 8.4) is not already
# [0, 1], so a normaliser that quietly did nothing would fail these.


def _wedge(window_size=0):
    return FlexProfiler("wedge", window_size=window_size).profile_seqs(SEQS)


def test_first_order_block_is_the_profile_values():
    X, cols = _feature_block(_wedge(), 1)
    np.testing.assert_allclose(X, [[1.1, 6.7, 1.1, 0.9, 1.1, 6.7, 1.1],
                                   [7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2]])
    assert cols[0] == "wedge.w0.o1.p1"
    assert cols[-1] == "wedge.w0.o1.p7"


def test_second_order_block_multiplies_adjacent_positions():
    X, cols = _feature_block(_wedge(), 2)
    assert X.shape == (2, 6)                   # m - n + 1
    np.testing.assert_allclose(X, [[7.37, 7.37, 0.99, 0.99, 7.37, 7.37],
                                   [51.84] * 6])
    assert cols[0] == "wedge.w0.o2.p1"
    assert cols[-1] == "wedge.w0.o2.p6"


def test_third_order_block_multiplies_three_positions():
    X, _ = _feature_block(_wedge(), 3)
    assert X.shape == (2, 5)
    np.testing.assert_allclose(X[0][0], 1.1 * 6.7 * 1.1)


def test_window_size_appears_in_the_column_name():
    """Two profilers with different windows must not collide on concat."""
    _, wide = _feature_block(_wedge(window_size=4), 1)
    _, narrow = _feature_block(_wedge(window_size=0), 1)
    assert wide[0] == "wedge.w4.o1.p1"
    assert narrow[0] == "wedge.w0.o1.p1"
    assert set(wide).isdisjoint(narrow)


def test_a_masked_value_propagates_through_the_product():
    """One N blanks the two k-mers covering it, and every product using them.

    ACNTACGT -> k-mers AC CN NT TA AC CG GT, so values 2 and 3 are masked
    and the order-2 products p1 (AC*CN), p2 (CN*NT) and p3 (NT*TA) all go
    with them. p4 (TA*AC) is the first that survives.
    """
    prof = FlexProfiler("wedge", window_size=0).profile_seqs(["ACNTACGT"])
    X, _ = _feature_block(prof, 2)
    assert list(np.isnan(X[0])) == [True, True, True, False, False, False]
    np.testing.assert_allclose(X[0][3], 0.9 * 1.1)


def test_order_below_one_raises():
    with pytest.raises(ValueError, match="order must be >= 1"):
        _feature_block(_wedge(), 0)


def test_order_beyond_the_block_width_raises():
    with pytest.raises(ValueError, match="7"):
        _feature_block(_wedge(), 8)


def test_ragged_rows_raise_naming_the_widths():
    from DNAflexpy.results import FlexProfile

    ragged = FlexProfile([["a", 1.0, 2.0], ["b", 3.0]],
                         feature="gc", window_size=0, kmer_len=2)
    with pytest.raises(ValueError, match=r"\[1, 2\]"):
        _feature_block(ragged, 1)


# --- min-max normalisation --------------------------------------------------


def test_minmax_scales_the_whole_block_not_each_column():
    X, _ = _feature_block(_wedge(), 1)
    scaled = _minmax(X)
    assert scaled.min() == 0.0 and scaled.max() == 1.0
    # lo=0.9, hi=7.2 across the WHOLE block. Per-column scaling would put
    # column 1's 1.1 at 0.0, not at 0.031746.
    np.testing.assert_allclose(scaled[0][0], (1.1 - 0.9) / (7.2 - 0.9))
    np.testing.assert_allclose(scaled[1][0], 1.0)


def test_minmax_keeps_nan_as_nan():
    scaled = _minmax(np.array([[1.0, np.nan], [3.0, 5.0]]))
    assert np.isnan(scaled[0][1])
    np.testing.assert_allclose(scaled[0][0], 0.0)
    np.testing.assert_allclose(scaled[1][1], 1.0)


def test_minmax_of_a_constant_block_is_zero_not_a_division_error():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        scaled = _minmax(np.array([[4.0, 4.0], [4.0, 4.0]]))
    np.testing.assert_allclose(scaled, np.zeros((2, 2)))


def test_minmax_of_an_all_nan_block_is_unchanged_and_silent():
    with warnings.catch_warnings():
        warnings.simplefilter("error")         # numpy's nanmin RuntimeWarning fails
        scaled = _minmax(np.full((2, 2), np.nan))
    assert np.isnan(scaled).all()


def test_minmax_does_not_mutate_its_input():
    block = np.array([[1.0, 3.0]])
    _minmax(block)
    np.testing.assert_allclose(block, [[1.0, 3.0]])
