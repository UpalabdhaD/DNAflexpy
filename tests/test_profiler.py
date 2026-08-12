import numpy as np
import pytest

from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.results import FlexProfile

SEQ = "ATGCGTACGTAGCTAGCGTAGCTAGT"


def test_bare_string_returns_an_array():
    out = FlexProfiler("DNaseI", window_size=10).profile(SEQ)
    assert isinstance(out, np.ndarray)
    assert len(out) == len(SEQ) - 10 + 1


def test_no_file_or_seqid_needed():
    """The archive forced a FASTA round-trip to inspect one sequence."""
    assert len(FlexProfiler("trx", window_size=0).profile("ATGC")) == 3


def test_unknown_feature_fails_at_construction():
    """A typo must not survive until after a full FASTA pass."""
    with pytest.raises(ValueError, match="unknown feature"):
        FlexProfiler("DNasel")


def test_features_the_archive_could_not_reach_now_work():
    assert len(FlexProfiler("gc", window_size=0).profile("ATGC")) == 3


def test_profile_rejects_multi_feature():
    p = FlexProfiler(["DNaseI", "trx"])
    with pytest.raises(ValueError, match="single feature"):
        p.profile(SEQ)


def test_profile_seqs_from_list_generates_ids():
    prof = FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ, SEQ])
    assert isinstance(prof, FlexProfile)
    assert prof.seqids == ["seq_0", "seq_1"]


def test_profile_seqs_from_dict_keeps_ids():
    prof = FlexProfiler("DNaseI", window_size=10).profile_seqs({"a": SEQ, "b": SEQ})
    assert prof.seqids == ["a", "b"]


def test_multi_feature_returns_a_profile_set():
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).profile_seqs({"a": SEQ})
    assert isinstance(ps, ProfileSet)
    assert set(ps) == {"DNaseI", "trx"}
    # Different k means different vector lengths for the same sequence.
    assert ps["DNaseI"].kmer_len == 3
    assert ps["trx"].kmer_len == 2


def test_custom_lookup_from_dict():
    """This 4-entry table is deliberately incomplete for k=2 (needs 16);
    the resulting UserWarning is expected, not a leak."""
    with pytest.warns(UserWarning, match="incomplete"):
        p = FlexProfiler("mine", window_size=0, lookup={"mine": {"AA": 1.0, "AT": 2.0,
                                                                "TA": 3.0, "TT": 4.0}})
    assert list(p.profile("ATA")) == [2.0, 3.0]


def test_short_sequence_yields_no_values():
    assert len(FlexProfiler("DNaseI", window_size=30).profile("ATGC")) == 0
