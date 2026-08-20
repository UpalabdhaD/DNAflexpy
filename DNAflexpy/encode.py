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


class FeatureMatrix:
    """A design matrix and its column names, ready for scikit-learn.

    `.X` is a plain numpy array, not a DataFrame, because these get wide:
    3-mer encoding of 500-position sequences is about 32,000 columns, and a
    DataFrame that wide is slow to build and slow to index. `.to_frame()`
    is there for when a DataFrame is genuinely what you want.
    """

    def __init__(self, X, columns, seqids, y=None, window_size=None,
                 feature_names=None):
        self.X = X
        self.columns = list(columns)
        self.seqids = list(seqids)
        self.y = y
        self.window_size = window_size
        self.feature_names = list(feature_names) if feature_names else []

    @property
    def shape(self) -> tuple[int, int]:
        return self.X.shape

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame(self.X, index=self.seqids, columns=self.columns)

    def __repr__(self) -> str:
        rows, cols = self.X.shape
        # `is not None`, not truthiness: an all-zero label vector is a
        # real label vector.
        return (f"<FeatureMatrix {rows} x {cols}, "
                f"features={self.feature_names}, "
                f"y={'yes' if self.y is not None else 'no'}>")


def encode(profiles, feature_names, normalize: bool = True) -> FeatureMatrix:
    """Build a design matrix from one or more profiles.

    `feature_names` is a list, and each entry names one block of columns:

    - `"1-mer"`, `"2-mer"`, `"3-mer"` - one-hot sequence encoding, 4, 16 or
      64 binary columns per position.
    - `"{n}-{feature}"` - the nth-order terms of a profiled feature, i.e.
      the products of that feature at n adjacent positions. `"1-DNaseI"` is
      the profile values themselves; `"2-DNaseI"` is the adjacent products.

    Blocks are concatenated left to right in the order requested, and the
    column names carry enough to tell them apart afterwards.

    Blocks may differ in width, which is fine. At `window_size > 0` every
    feature yields `L - window_size + 1` values, so all feature blocks
    match. At `window_size = 0` a block is `L - k + 1` wide, so a 2-mer
    feature such as `gc` gives one more column than a 3-mer one such as
    `DNaseI`. A `1-mer` sequence block is `L` wide either way.

    `normalize=True` min-max scales each feature block on its own range -
    features span wildly different units (`stiffness` runs to 5500, `gc` to
    1.0), so one shared scale would flatten the smaller one to nothing.
    One-hot blocks are left alone; they are already 0/1. Note that min-max
    uses **this dataset's** range: if you are going to split train and test,
    pass `normalize=False` and scale inside a scikit-learn pipeline instead,
    so the test set does not leak into the scaling.

    NaN marks a position that could not be resolved - an `N`, an IUPAC code,
    or padding from `from_bed(on_edge="pad")`. It is left as NaN rather than
    filled, so it stays distinguishable from a real measurement of zero;
    reach for `sklearn.impute.SimpleImputer` if your model needs it filled.

    Every sequence must be the same length, since the columns are positions.
    `from_bed(..., width=...)` is the usual way to get that.

    `DNAflexpy.profile()` returns a bare array, so there is nothing to
    encode from the one-liner; start from `profile_seqs`, `profile_fasta`,
    `profile_table` or `from_bed`.
    """
    by_feature = _as_mapping(profiles)
    if isinstance(feature_names, str):
        raise TypeError(
            "feature_names must be a list of names, not a single string; "
            f"did you mean [{feature_names!r}]?"
        )
    requests = list(feature_names)
    if not requests:
        raise ValueError("at least one feature name is required, e.g. ['1-mer']")
    duplicates = {n for n in requests if requests.count(n) > 1}
    if duplicates:
        raise ValueError(
            f"duplicate feature name(s) {sorted(duplicates)}; each name "
            "would produce identically named columns"
        )

    reference = next(iter(by_feature.values()))
    for name, profile in by_feature.items():
        if profile.seqids != reference.seqids:
            raise ValueError(
                f"every profile must cover the same sequences in the same "
                f"order; {name!r} does not match {reference.feature!r}"
            )

    parsed = [_parse(request, by_feature) for request in requests]
    if any(kind == "mer" for kind, _, _ in parsed):
        _require_equal_lengths(reference)

    blocks, columns = [], []
    for kind, order, target in parsed:
        if kind == "mer":
            block, names = _one_hot(reference.seqs, order)
        else:
            block, names = _feature_block(by_feature[target], order)
            if normalize:
                block = _minmax(block)
        blocks.append(block)
        columns.extend(names)

    return FeatureMatrix(
        np.hstack(blocks),
        columns,
        reference.seqids,
        y=reference.y,
        window_size=reference.window_size,
        feature_names=requests,
    )


def _as_mapping(profiles) -> dict:
    """Accept a FlexProfile, a ProfileSet, or a plain {feature: profile}."""
    if isinstance(profiles, dict):
        if not profiles:
            raise ValueError("no profiles to encode")
        return dict(profiles)
    if hasattr(profiles, "feature") and hasattr(profiles, "_rows"):
        return {profiles.feature: profiles}
    raise TypeError(
        "profiles must be a FlexProfile or a ProfileSet, got "
        f"{type(profiles).__name__}"
    )


def _parse(request, by_feature: dict) -> tuple[str, int, str | None]:
    """Split `"{n}-{what}"` into a block kind, an order and a feature."""
    grammar = (
        "expected '<n>-mer' for sequence encoding or '<n>-<feature>' for "
        "an nth-order feature term, with n a positive whole number"
    )
    if not isinstance(request, str) or "-" not in request:
        raise ValueError(f"malformed feature name {request!r}: {grammar}")
    head, _, what = request.partition("-")
    if not head.isdigit() or int(head) < 1:
        raise ValueError(f"malformed feature name {request!r}: {grammar}")
    order = int(head)

    if what == "mer":
        if "mer" in by_feature:
            raise ValueError(
                f"{request!r} is ambiguous: the lookup table has a feature "
                "literally named 'mer', so this could mean sequence one-hot "
                "encoding or that feature. Rename the feature in the table."
            )
        return "mer", order, None
    if what not in by_feature:
        raise ValueError(
            f"unknown feature {what!r} in {request!r}. Profiled features "
            f"are: {', '.join(by_feature)}"
        )
    return "feature", order, what


def _require_equal_lengths(profile) -> None:
    """One-hot encoding needs the bases, and needs them all the same length.

    Checked on the sequences, not the row widths: two sequences shorter
    than the window both produce an empty row, so equal row widths would
    wrongly pass.
    """
    if profile.seqs is None:
        raise ValueError(
            "one-hot sequence encoding needs the input sequences, which this "
            "profile does not carry. Build it with profile_seqs, "
            "profile_fasta, profile_table or from_bed, all of which keep them."
        )
    lengths = {len(s) for s in profile.seqs}
    if len(lengths) > 1:
        raise ValueError(
            f"encoding needs equal-length sequences, found {sorted(lengths)}. "
            "Use from_bed(..., width=N) to cut fixed-width windows, or "
            "trim the sequences yourself."
        )
