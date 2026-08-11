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
