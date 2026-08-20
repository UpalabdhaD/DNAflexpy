"""The archived package must keep producing exactly what it always produced."""
import pathlib

import pytest

from tests.conftest import ROOT, archive_bytes

# Every checked-in TSV, with the (fasta, feature, window_size) that produces it.
# Verified: all 16 reproduce byte-exactly.
CHECKED_IN = [
    ("tests/test_outputs/output_w0.tsv", "test_fasta", "DNaseI", 0),
    ("tests/test_outputs/output_w3.tsv", "test_fasta", "DNaseI", 3),
    ("tests/test_outputs/output_w10.tsv", "test_fasta", "DNaseI", 10),
    ("tests/test_outputs/output_w15.tsv", "test_fasta", "DNaseI", 15),
    ("tests/edge_case_outputs/exact_kmer_w0.tsv", "exact", "DNaseI", 0),
    ("tests/edge_case_outputs/exact_kmer_w3.tsv", "exact", "DNaseI", 3),
    ("tests/edge_case_outputs/short_w0.tsv", "edge", "DNaseI", 0),
    ("tests/edge_case_outputs/trx_w0.tsv", "edge", "trx", 0),
    ("tests/edge_case_outputs/trx_w1.tsv", "edge", "trx", 1),
    ("tests/edge_case_outputs/trx_w2.tsv", "edge", "trx", 2),
    ("tests/edge_case_outputs/w1.tsv", "edge", "DNaseI", 1),
    ("tests/edge_case_outputs/w2.tsv", "edge", "DNaseI", 2),
    ("tests/edge_case_outputs/w4.tsv", "edge", "DNaseI", 4),
    ("tests/edge_case_outputs/w10_mixed.tsv", "edge", "DNaseI", 10),
    ("tests/edge_case_outputs/w26_exact.tsv", "test_fasta", "DNaseI", 26),
    ("tests/edge_case_outputs/w30.tsv", "edge", "DNaseI", 30),
]


@pytest.mark.parametrize("relpath,fasta_key,feature,window", CHECKED_IN)
def test_archive_reproduces_checked_in_output(relpath, fasta_key, feature, window):
    expected = (ROOT / relpath).read_bytes()
    assert archive_bytes(fasta_key, feature, window) == expected


def test_archive_reads_its_own_lookup_table():
    """Guards the import-path fix: the archive must not read the new package's YAML.

    If rxv/DNAflexpy/utils.py:175 were reverted to `files("DNAflexpy.data")`,
    the archive would silently read the NEW package's lookup table instead
    of its own -- and this failure mode is invisible to every *other* test
    in the suite: DNAflexpy/data/lookupNEW.yaml and
    rxv/DNAflexpy/data/lookupNEW.yaml currently parse to byte-identical
    content, so a reverted archive would still produce byte-identical
    output and the 210-case differential gate in test_differential.py would
    keep reporting every case as passed. It would just be comparing the new
    package against itself, testing nothing. No behavioural or round-trip
    check can catch that; only inspecting the archive's own source for the
    exact call site can. This test is therefore the only thing standing
    between the gate and silent self-comparison -- if it is ever weakened
    or deleted, the gate can pass 210/210 while verifying nothing.
    """
    import importlib.resources

    import rxv.DNAflexpy.utils as archived

    source = pathlib.Path(archived.__file__).read_text()
    # Assert the exact call, not just the package string appearing anywhere
    # in the file: a whole-file substring search for "rxv.DNAflexpy.data"
    # would also be satisfied by that text sitting in a comment or a
    # docstring, without the fixed code path actually being in effect.
    assert 'files("rxv.DNAflexpy.data")' in source

    # Independently confirm the resource genuinely resolves under the
    # repo's rxv/ tree, and to a different location than the new package's
    # own data directory. This is a second, independent signal alongside
    # the source-text check above -- it would catch the two data
    # directories ever being consolidated into one, which the source-text
    # check alone would not.
    # `files()` on a namespace package (no __init__.py in data/) returns a
    # MultiplexedPath, whose str() is a repr, not a filesystem path -- join
    # a real file to get back a concrete PosixPath.
    repo_root = pathlib.Path(__file__).resolve().parent.parent
    archive_data = importlib.resources.files("rxv.DNAflexpy.data").joinpath("lookupNEW.yaml")
    new_data = importlib.resources.files("DNAflexpy.data").joinpath("lookupNEW.yaml")
    archive_path = pathlib.Path(str(archive_data)).resolve()
    new_path = pathlib.Path(str(new_data)).resolve()
    assert archive_path.is_relative_to((repo_root / "rxv").resolve())
    assert archive_path != new_path
