import pathlib

from DNAflexpy.core import FlexProfiler
from DNAflexpy.io import read_fasta

ROOT = pathlib.Path(__file__).resolve().parent.parent
FASTA = ROOT / "rxv/DNAflexpy/data/test_fasta.fa"


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


def test_profile_fasta_pooled_matches_single_threaded():
    p = FlexProfiler("DNaseI", window_size=10)
    assert p.profile_fasta(FASTA, threads=2).frame.equals(
        p.profile_fasta(FASTA, threads=1).frame
    )


def test_pooled_run_preserves_input_order():
    p = FlexProfiler("DNaseI", window_size=10)
    assert p.profile_fasta(FASTA, threads=2).seqids == ["sequence1", "sequence2"]


def test_multi_feature_fasta_returns_profile_set():
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).profile_fasta(FASTA, threads=1)
    assert set(ps) == {"DNaseI", "trx"}
