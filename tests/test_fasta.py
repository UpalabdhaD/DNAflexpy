import pathlib
import warnings

import pytest

from DNAflexpy.core import FlexProfiler
from DNAflexpy.io import read_fasta

ROOT = pathlib.Path(__file__).resolve().parent.parent
# A copy of the archive's fixture (rxv/DNAflexpy/data/test_fasta.fa), kept in
# sync byte-for-byte so these unit tests of the NEW package don't depend on
# rxv/ still existing. The differential tests legitimately read the archive's
# own copy via conftest.FASTAS -- that's the point of that gate.
FASTA = ROOT / "tests/test_fasta_basic.fa"


def test_read_fasta_yields_id_and_sequence():
    records = list(read_fasta(FASTA))
    assert records[0][0] == "sequence1"
    assert records[0][1] == "ATGCGTACGTAGCTAGCGTAGCTAGT"
    assert len(records) == 2


def test_read_fasta_handles_missing_trailing_newline():
    """test_fasta.fa ends without one."""
    assert len(list(read_fasta(FASTA))[1][1]) == 26


def test_read_fasta_joins_wrapped_lines(tmp_path):
    fa = tmp_path / "w.fa"
    fa.write_text(">s\nATGC\nGTAC\n")
    assert list(read_fasta(fa)) == [("s", "ATGCGTAC")]


def test_missing_file_raises(tmp_path):
    import pytest

    with pytest.raises(FileNotFoundError):
        list(read_fasta(tmp_path / "nope.fa"))


def test_profile_fasta_single_threaded():
    prof = FlexProfiler("DNaseI", window_size=10).profile_fasta(FASTA, threads=1)
    assert prof.seqids == ["sequence1", "sequence2"]
    assert len(prof.frame.columns) == 17


def test_profile_fasta_pooled_matches_single_threaded(monkeypatch):
    """2 records is below the Pool threshold; lower it so threads=2 still
    genuinely spawns a Pool instead of silently falling back to serial."""
    monkeypatch.setattr("DNAflexpy.core._MIN_RECORDS_FOR_POOL", 1)
    p = FlexProfiler("DNaseI", window_size=10)
    assert p.profile_fasta(FASTA, threads=2).frame.equals(
        p.profile_fasta(FASTA, threads=1).frame
    )


def test_pooled_run_preserves_input_order(monkeypatch):
    monkeypatch.setattr("DNAflexpy.core._MIN_RECORDS_FOR_POOL", 1)
    p = FlexProfiler("DNaseI", window_size=10)
    assert p.profile_fasta(FASTA, threads=2).seqids == ["sequence1", "sequence2"]


def test_multi_feature_fasta_returns_profile_set():
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).profile_fasta(FASTA, threads=1)
    assert set(ps) == {"DNaseI", "trx"}


def test_default_threads_stays_serial_for_small_input(monkeypatch):
    """threads=None must mean "decide automatically", not "always spawn a
    Pool": a small FASTA on the default call path must never create one."""
    import multiprocessing

    def no_pool(*args, **kwargs):
        raise AssertionError("Pool should not be created below the threshold")

    monkeypatch.setattr(multiprocessing, "Pool", no_pool)
    prof = FlexProfiler("DNaseI", window_size=10).profile_fasta(FASTA)
    assert prof.seqids == ["sequence1", "sequence2"]


def test_explicit_threads_above_one_still_pools_for_large_input(monkeypatch):
    """A large-enough input with an explicit threads>1 must still use a real
    Pool -- lowering the threshold must not accidentally disable pooling."""
    import multiprocessing as mp

    monkeypatch.setattr("DNAflexpy.core._MIN_RECORDS_FOR_POOL", 1)
    calls = {"n": 0}
    original_pool = mp.Pool

    def spying_pool(*args, **kwargs):
        calls["n"] += 1
        return original_pool(*args, **kwargs)

    monkeypatch.setattr(mp, "Pool", spying_pool)
    prof = FlexProfiler("DNaseI", window_size=10).profile_fasta(FASTA, threads=2)
    assert prof.seqids == ["sequence1", "sequence2"]
    assert calls["n"] == 1


def test_profile_fasta_warns_on_ambiguous_sequence(tmp_path):
    """read_fasta never calls warn_if_ambiguous itself; profile_fasta must,
    so this path warns the same as profile_table already does."""
    fa = tmp_path / "amb.fa"
    fa.write_text(">a\nATGRGTACGTAGCTAGCGTAGCTAGT\n")
    with pytest.warns(UserWarning, match="non-ACGTN"):
        FlexProfiler("DNaseI", window_size=10).profile_fasta(fa, threads=1)


def test_profile_fasta_does_not_warn_on_a_clean_fasta():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        FlexProfiler("DNaseI", window_size=10).profile_fasta(FASTA, threads=1)
