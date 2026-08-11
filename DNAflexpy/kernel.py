"""Pure numeric core of profiling. No I/O, no classes, no pandas.

Arithmetic here is deliberately conservative: builtin `sum` over Python
lists, left to right, matching the archive exactly. Swapping in a numpy
reduction changes the last bits and flips values at `round(v, 3)`
boundaries, which breaks byte equality against the archive.
"""
from __future__ import annotations

import math

MASKED = float("nan")


def kmer_values(sequence: str, kmer_len: int, table: dict[str, float]) -> list[float]:
    """Per-position k-mer values, left to right.

    Any k-mer absent from `table` - because it contains N or an IUPAC code,
    or because a user-supplied table is incomplete - resolves to NaN. The
    archive resolved these to 0, which is indistinguishable from a real
    measurement (gc AA=0.0, trx AA=0).
    """
    seq = sequence.upper()
    return [
        table.get(seq[i : i + kmer_len], MASKED)
        for i in range(len(seq) - kmer_len + 1)
    ]


def window_means(
    values: list[float], window_size: int, kmer_len: int, n_positions: int
) -> list[float]:
    """Mean of the k-mer values inside each sliding window, rounded to 3 dp.

    `n_positions` is the sequence length, which sets the window count
    independently of how many k-mers exist.
    """
    out: list[float] = []
    per_window = window_size - kmer_len + 1
    for start in range(n_positions - window_size + 1):
        window = values[start : start + per_window] if per_window > 0 else []
        out.append(round(_mean(window), 3))
    return out


def _mean(window: list[float]) -> float:
    """Archive-compatible mean, extended to skip masked values.

    With no masked values this is exactly `sum(w) / len(w)`, preserving the
    archive's arithmetic. An empty window gives 0.0, matching the archive's
    `if w else 0.0`. A fully masked window gives NaN, which the archive
    could not express.
    """
    if not window:
        return 0.0
    if not any(math.isnan(v) for v in window):
        return sum(window) / len(window)
    resolved = [v for v in window if not math.isnan(v)]
    if not resolved:
        return MASKED
    return sum(resolved) / len(resolved)


def profile_values(
    sequence: str, kmer_len: int, table: dict[str, float], window_size: int
) -> list[float]:
    """Per-k-mer values when `window_size` is 0, otherwise per-window means."""
    if window_size < 0:
        raise ValueError(f"window_size must be >= 0, got {window_size}")
    values = kmer_values(sequence, kmer_len, table)
    if window_size == 0:
        return values
    return window_means(values, window_size, kmer_len, len(sequence))
