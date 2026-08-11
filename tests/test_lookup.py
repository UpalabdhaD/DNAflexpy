import warnings

import pytest

from DNAflexpy.lookup import FeatureTable, default_table


def test_infers_kmer_length_from_keys():
    t = default_table()
    assert t.kmer_len("DNaseI") == 3
    assert t.kmer_len("trx") == 2


def test_inference_matches_the_archives_hardcoded_map():
    """The archive hardcodes feature -> k. Inference must agree on all of them."""
    from rxv.DNAflexpy.utils import get_kmer_len

    t = default_table()
    for feature in t.features:
        hardcoded = get_kmer_len(feature)
        if hardcoded is not None:
            assert t.kmer_len(feature) == hardcoded, feature


def test_unlocks_features_the_archive_could_not_reach():
    """freeen, gc and mechen have no entry in the archive's map."""
    t = default_table()
    for feature in ("freeen", "gc", "mechen"):
        assert feature in t.features
        assert t.kmer_len(feature) == 2


def test_packaged_tables_are_complete():
    t = default_table()
    for feature in t.features:
        assert len(t.table(feature)) == 4 ** t.kmer_len(feature), feature


def test_table_is_read_only():
    """The cached singleton must not be corruptible by a caller."""
    t = default_table()
    with pytest.raises(TypeError):
        t.table("DNaseI")["ZZZ"] = 99.0
    assert "ZZZ" not in default_table().table("DNaseI")


def test_unknown_feature_raises_with_available_names():
    t = default_table()
    with pytest.raises(ValueError, match="DNaseI"):
        t.kmer_len("DNasel")  # lowercase L typo


def test_keys_are_uppercased():
    t = FeatureTable.from_dict({"f": {"aa": 1.0, "at": 2.0, "ta": 3.0, "tt": 4.0}})
    assert t.table("f")["AA"] == 1.0


def test_mixed_key_lengths_rejected():
    """This is what makes k inference safe."""
    with pytest.raises(ValueError, match="mixed k-mer lengths"):
        FeatureTable.from_dict({"f": {"AA": 1.0, "AAT": 2.0}})


def test_non_acgt_keys_rejected():
    with pytest.raises(ValueError, match="non-ACGT"):
        FeatureTable.from_dict({"f": {"AN": 1.0}})


def test_non_numeric_values_rejected():
    with pytest.raises(ValueError, match="non-numeric"):
        FeatureTable.from_dict({"f": {"AA": "high"}})


def test_bool_values_rejected():
    """bool is an int subclass; it must not slip through as numeric."""
    with pytest.raises(ValueError, match="non-numeric"):
        FeatureTable.from_dict({"f": {"AA": True, "AT": 1.0, "TA": 1.0, "TT": 1.0}})


def test_empty_table_rejected():
    with pytest.raises(ValueError, match="empty"):
        FeatureTable.from_dict({"f": {}})


def test_incomplete_table_warns_but_loads():
    with pytest.warns(UserWarning, match="incomplete"):
        t = FeatureTable.from_dict({"f": {"AA": 1.0, "AT": 2.0}})
    assert t.kmer_len("f") == 2


def test_default_table_is_memoised():
    """The 8 ms YAML parse must happen once, not per call."""
    assert default_table() is default_table()
