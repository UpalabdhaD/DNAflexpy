import warnings

import numpy as np
import pytest

from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.results import FlexProfile

SEQ = "ATGCGTACGTAGCTAGCGTAGCTAGT"
AMBIGUOUS_SEQ = "ATGRGTACGTAGCTAGCGTAGCTAGT"


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


def test_profile_warns_on_ambiguous_sequence():
    """read_table warns via its own reader; profile() must warn too, since
    read_fasta never calls warn_if_ambiguous itself."""
    with pytest.warns(UserWarning, match="non-ACGTN"):
        FlexProfiler("DNaseI", window_size=10).profile(AMBIGUOUS_SEQ)


def test_profile_does_not_warn_on_a_clean_sequence():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        FlexProfiler("DNaseI", window_size=10).profile(SEQ)


def test_profile_seqs_warns_on_ambiguous_sequence():
    with pytest.warns(UserWarning, match="non-ACGTN"):
        FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ, AMBIGUOUS_SEQ])


def test_profile_seqs_does_not_warn_on_clean_sequences():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ, SEQ])


# --- Phase 4: the input sequences ride along on the profile -----------------
# `encode(["1-mer", ...])` one-hot encodes the bases themselves, which the
# profile values alone cannot reconstruct. `to_tsv` never reads `seqs`, so
# the byte-equality gate cannot notice a path that leaves it None -- hence a
# test per entry point.


def test_profile_seqs_retains_the_input_sequences():
    prof = FlexProfiler("gc", window_size=0).profile_seqs({"a": "ACGT", "b": "TTTT"})
    assert prof.seqs == ["ACGT", "TTTT"]
    assert prof.seqids == ["a", "b"]


def test_seqs_is_a_copy_not_the_stored_list():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(["ACGT"])
    prof.seqs.append("TTTT")
    assert prof.seqs == ["ACGT"]


def test_every_profile_in_a_set_shares_the_same_sequences():
    profiles = FlexProfiler(["gc", "wedge"], window_size=0).profile_seqs(["ACGT"])
    assert profiles["gc"].seqs == ["ACGT"]
    assert profiles["wedge"].seqs == ["ACGT"]


def test_profile_fasta_retains_sequences_serially(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">a\nACGT\n>b\nTTTT\n")
    prof = FlexProfiler("gc", window_size=0).profile_fasta(fa, threads=1)
    assert prof.seqs == ["ACGT", "TTTT"]


def test_profile_fasta_retains_sequences_through_the_pool(tmp_path):
    """The pooled branch builds rows in workers and never sees `records`.

    It is the one path where `seqs` can silently come back None, because it
    calls `_assemble` directly instead of going through `_build`.
    """
    import DNAflexpy.core as core

    n = core._MIN_RECORDS_FOR_POOL + 2
    fa = tmp_path / "many.fa"
    fa.write_text("".join(f">s{i}\nACGTACGT\n" for i in range(n)))
    prof = FlexProfiler("gc", window_size=0).profile_fasta(fa, threads=2)
    assert len(prof.seqs) == n
    assert prof.seqs[0] == "ACGTACGT"
    assert prof.seqids[-1] == f"s{n - 1}"


def test_profile_table_retains_sequences_and_labels(tmp_path):
    tsv = tmp_path / "t.tsv"
    tsv.write_text("ACGT\t1.5\nTTTT\t2.5\n")
    prof = FlexProfiler("gc", window_size=0).profile_table(tsv, header=False)
    assert prof.seqs == ["ACGT", "TTTT"]
    assert prof.y == [1.5, 2.5]


def test_a_hand_built_profile_has_no_sequences():
    """`seqs` is optional: nothing outside encoding needs it."""
    bare = FlexProfile([["a", 1.0]], feature="gc", window_size=0, kmer_len=2)
    assert bare.seqs is None
