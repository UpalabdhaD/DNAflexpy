"""Shared test fixtures for verifying the archived DNAflexpy package.

`archive_bytes` reproduces DNAflexpyMP's serialisation without spawning a
Pool (multiprocessing uses `spawn` on macOS, which the Pool needs but tests
don't). It is defined once here, rather than duplicated per test module,
because more than one task's tests need it and two copies could drift --
silently changing what "the archive's output" means.
"""
import contextlib
import io
import pathlib

import pandas as pd

try:  # pragma: no cover - matplotlib is an optional extra
    import matplotlib

    # Must happen before anything imports pyplot: the default backend on
    # macOS opens a window, which hangs a headless test run.
    matplotlib.use("Agg")
except ImportError:
    pass

from rxv.DNAflexpy.utils import (
    get_kmer_len,
    load_feature_data,
    read_fasta,
    seq_to_numeric_profile,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent

FASTAS = {
    "test_fasta": ROOT / "rxv/DNAflexpy/data/test_fasta.fa",
    "edge": ROOT / "tests/test_edge_cases.fasta",
    "exact": ROOT / "tests/test_exact_kmer.fasta",
}


def archive_bytes(fasta_key, feature, window_size):
    """Reproduce DNAflexpyMP's serialisation without spawning a Pool.

    The Pool only parallelises the same per-record call, so bypassing it is
    equivalent and far faster in tests.
    """
    lookup = load_feature_data()
    kmer_len = get_kmer_len(feature)
    rows = []
    # The archive prints warnings to stdout for short sequences; silence them.
    with contextlib.redirect_stdout(io.StringIO()):
        for seqid, seq in read_fasta(str(FASTAS[fasta_key])):
            rows.append(
                seq_to_numeric_profile(seqid, seq, kmer_len, window_size, feature, lookup)
            )
    buf = io.StringIO()
    pd.DataFrame(rows).to_csv(buf, index=False, header=False, sep="\t")
    return buf.getvalue().encode()
