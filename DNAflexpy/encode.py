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
