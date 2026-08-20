"""Turn profiles into a machine-learning design matrix.

Modelled on DNAshapeR's `encodeSeqShape`: a request is a list of feature
names, each naming one block of columns, and the blocks are concatenated
left to right in the order asked for.

    "1-mer", "2-mer", "3-mer"   one-hot sequence encoding
    "1-DNaseI", "2-stiffness"   nth-order terms for a profiled feature
"""
from __future__ import annotations

import itertools
import warnings

import numpy as np

_ALPHABET = "ACGT"


def _one_hot(seqs, k: int) -> tuple[np.ndarray, list[str]]:
    """One-hot encode every k-mer position of each sequence.

    Layout, fixed here and relied on by the column names:

    - A sequence of length `L` has `L - k + 1` overlapping k-mer positions.
    - Each position contributes `4**k` binary columns, one per k-mer, in
      lexicographic ACGT order: `AA, AC, AG, AT, CA, ...`.
    - Positions are 1-based, matching `FlexProfile.frame`'s columns.
    - Column names are `seq.{k}mer.p{i}.{KMER}`, e.g. `seq.2mer.p3.GT`.

    The DNAshapeR manual gives the bit strings it uses (A=0001, C=0010,
    G=0100, T=1000) but not the layout of the resulting matrix, so this
    makes no claim to reproduce that package's column order.

    A k-mer containing any letter outside ACGT - `N`, or any other IUPAC
    code - gets an all-zero vector: the standard "unknown" encoding, and
    the direct analogue of the profile masking that same k-mer as NaN.

    Values are floats, not ints or bools, so concatenating this block with
    a float feature block cannot silently promote a dtype.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    lengths = {len(s) for s in seqs}
    if len(lengths) > 1:
        # `encode` checks this first and raises a fuller message. Repeating
        # the guard here means a direct call cannot silently truncate every
        # sequence to the shortest one.
        raise ValueError(f"sequences must all be the same length, got {sorted(lengths)}")
    shortest = min(lengths, default=0)
    if k > shortest:
        raise ValueError(
            f"k={k} is longer than the sequence ({shortest} base(s)); "
            "there is no k-mer to encode"
        )

    kmers = ["".join(t) for t in itertools.product(_ALPHABET, repeat=k)]
    offsets = {kmer: i for i, kmer in enumerate(kmers)}
    width = len(kmers)
    n_positions = shortest - k + 1

    X = np.zeros((len(seqs), width * n_positions), dtype=float)
    unknown_positions, unknown_seqs = 0, 0
    for row, sequence in enumerate(seqs):
        seq = sequence.upper()
        blanks = 0
        for pos in range(n_positions):
            offset = offsets.get(seq[pos : pos + k])
            if offset is None:
                blanks += 1
                continue
            X[row, pos * width + offset] = 1.0
        if blanks:
            unknown_positions += blanks
            unknown_seqs += 1

    if unknown_positions:
        warnings.warn(
            f"{unknown_positions} k-mer position(s) across {unknown_seqs} "
            f"sequence(s) contain a letter outside ACGT and were encoded "
            "as an all-zero column block. The same k-mers are masked as NaN "
            "in the profile values.",
            UserWarning,
            stacklevel=2,
        )

    columns = [
        f"seq.{k}mer.p{pos + 1}.{kmer}"
        for pos in range(n_positions)
        for kmer in kmers
    ]
    return X, columns


def _feature_block(profile, order: int) -> tuple[np.ndarray, list[str]]:
    """The `order`-th order terms of one profiled feature.

    Following DNAshapeR's `encodeNstOrderShape`, the nth-order block holds
    the products of the same feature at n adjacent positions, and nothing
    else. Order 1 is the identity. A profile of `m` values per sequence
    gives `m - n + 1` columns, column `i` being `prod(v[i : i + n])`.

    So `n-<feature>` yields only the nth-order terms; ask for `1-gc` and
    `2-gc` together if you want both.

    NaN propagates through the product, which is right: a product involving
    an unresolvable k-mer is itself unresolvable.

    Column names are `{feature}.w{window_size}.o{order}.p{i}`, 1-based. The
    window size is part of the name because it changes what a position
    means - at `window_size=10`, `p1` is the mean over sequence positions
    1-10; at `window_size=0` it is the k-mer starting at position 1. Without
    it, two matrices built at different window sizes would carry identical
    column names for different quantities.
    """
    if order < 1:
        raise ValueError(f"order must be >= 1, got {order}")

    widths = {len(row) - 1 for row in profile._rows}
    if len(widths) > 1:
        raise ValueError(
            f"feature {profile.feature!r} has ragged rows with widths "
            f"{sorted(widths)}; encoding needs one value per position for "
            "every sequence"
        )
    values = np.array([row[1:] for row in profile._rows], dtype=float)
    m = widths.pop() if widths else 0
    if order > m:
        raise ValueError(
            f"order {order} exceeds the {m} value(s) per sequence in "
            f"feature {profile.feature!r}; there are no {order} adjacent "
            "positions to multiply"
        )

    block = values[:, : m - order + 1].copy()
    for step in range(1, order):
        block *= values[:, step : m - order + 1 + step]

    columns = [
        f"{profile.feature}.w{profile.window_size}.o{order}.p{i}"
        for i in range(1, m - order + 2)
    ]
    return block, columns


def _minmax(block: np.ndarray) -> np.ndarray:
    """Min-max scale a whole block to [0, 1].

    Scalar bounds for the entire block, matching DNAshapeR's
    `normalize(x, max, min)` - not per column. Scaling each column
    separately would erase the between-position differences that make a
    profile a profile.

    NaN stays NaN. A constant block becomes all zeros rather than dividing
    by zero, and an all-NaN block is returned untouched: `np.nanmin` on one
    would emit its own RuntimeWarning and return NaN.
    """
    finite = np.isfinite(block)
    if not finite.any():
        return block.astype(float, copy=True)
    lo = float(block[finite].min())
    hi = float(block[finite].max())
    if hi == lo:
        out = np.zeros_like(block, dtype=float)
        out[~finite] = block[~finite]
        return out
    return (block - lo) / (hi - lo)
