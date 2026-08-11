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
    """Guards the import-path fix: the archive must not read the new package's YAML."""
    import rxv.DNAflexpy.utils as archived

    assert "rxv.DNAflexpy.data" in pathlib.Path(archived.__file__).read_text()
