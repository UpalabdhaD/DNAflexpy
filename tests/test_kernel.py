import math

from DNAflexpy.kernel import kmer_values, profile_values, window_means
from DNAflexpy.lookup import default_table

SEQ = "ATGCGTACGTAGCTAGCGTAGCTAGT"  # 26 nt, the canonical test sequence


def dnase():
    return default_table().table("DNaseI")


def test_kmer_values_match_the_archive():
    from rxv.DNAflexpy.utils import load_feature_data, transform_seq_to_feat

    expected = transform_seq_to_feat(SEQ, 3, "DNaseI", load_feature_data())
    assert kmer_values(SEQ, 3, dnase()) == expected


def test_window_means_match_the_archive():
    """window_size == kmer_len must reproduce the per-kmer values exactly."""
    values = kmer_values(SEQ, 3, dnase())
    assert window_means(values, 3, 3, len(SEQ)) == [round(v, 3) for v in values]


def test_window_count_follows_the_formula():
    values = kmer_values(SEQ, 3, dnase())
    assert len(window_means(values, 10, 3, len(SEQ))) == len(SEQ) - 10 + 1


def test_window_smaller_than_kmer_yields_zeros():
    """No k-mer fits, so the archive's `sum(w)/len(w) if w else 0.0` gives 0.0."""
    values = kmer_values(SEQ, 3, dnase())
    out = window_means(values, 2, 3, len(SEQ))
    assert out == [0.0] * (len(SEQ) - 2 + 1)


def test_window_larger_than_sequence_yields_nothing():
    values = kmer_values("ATGC", 3, dnase())
    assert window_means(values, 30, 3, 4) == []


def test_sequence_is_uppercased():
    assert kmer_values("atgcgt", 3, dnase()) == kmer_values("ATGCGT", 3, dnase())


def test_ambiguous_kmers_are_masked_not_zeroed():
    """The archive scores these 0, which is a real value in these tables."""
    values = kmer_values("ATGNGTACGT", 3, dnase())
    assert not math.isnan(values[0])
    assert all(math.isnan(v) for v in values[1:4])
    assert not math.isnan(values[4])


def test_iupac_codes_are_masked_too():
    values = kmer_values("ATGRGTACGT", 3, dnase())
    assert all(math.isnan(v) for v in values[1:4])


def test_partially_masked_window_averages_the_resolved_kmers():
    values = [1.0, 2.0, float("nan"), 4.0]
    # window of 3 kmers starting at 0 -> mean(1.0, 2.0) == 1.5
    assert window_means(values, 5, 3, 7)[0] == 1.5


def test_fully_masked_window_is_nan():
    nan = float("nan")
    assert math.isnan(window_means([nan, nan], 4, 3, 6)[0])


def test_profile_values_window_zero_returns_per_kmer():
    assert profile_values(SEQ, 3, dnase(), 0) == kmer_values(SEQ, 3, dnase())


def test_window_zero_bytes_match_the_archive_for_integer_features():
    """`==` hides 18 == 18.0; only serialisation exposes it."""
    import contextlib, io
    import pandas as pd
    from rxv.DNAflexpy.utils import load_feature_data, seq_to_numeric_profile

    for feature in ("NPP", "stiffness", "trx", "DNaseI"):
        k = default_table().kmer_len(feature)
        with contextlib.redirect_stdout(io.StringIO()):
            legacy = seq_to_numeric_profile("s", SEQ, k, 0, feature, load_feature_data())
        mine = ["s", *profile_values(SEQ, k, default_table().table(feature), 0)]
        to_bytes = lambda row: pd.DataFrame([row]).to_csv(
            index=False, header=False, sep="\t"
        ).encode()
        assert to_bytes(mine) == to_bytes(legacy), feature


def test_builtin_sum_is_compensated_on_this_interpreter():
    """The package's numbers depend on how builtin `sum()` adds floats.

    Python 3.12 gave `sum()` compensated (Neumaier) summation for floats.
    Before that it added naively, left to right, and produced a different
    result for the same input:

        [0.134, 0.076, -0.077, -0.033, 0.025, 0.025, -0.033, -0.033]
        3.11 -> 0.08399999999999999      3.12 -> 0.084

    Divided by 8 and rounded to three places that is 0.01 against 0.011, so a
    published profile would differ in its last digit. The checked-in expected
    outputs are 3.12 values, which is why `requires-python` is `>=3.12`.

    If this ever fails, the interpreter's summation changed again and every
    recorded expected output needs revisiting.
    """
    import math

    values = [0.134, 0.076, -0.077, -0.033, 0.025, 0.025, -0.033, -0.033]

    naive = 0.0
    for v in values:
        naive += v

    assert sum(values) == math.fsum(values), "sum() is no longer compensated"
    assert sum(values) != naive, "sum() now matches naive addition"
    assert round(sum(values) / len(values), 3) == 0.011
