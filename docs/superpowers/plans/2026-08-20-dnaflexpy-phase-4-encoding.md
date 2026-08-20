# DNAflexpy Phase 4 — ML Feature Encoding (`encode.py`)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn a `FlexProfile` / `ProfileSet` into a numeric design matrix that scikit-learn can consume directly, modelled on DNAshapeR's `encodeSeqShape`. This is the payoff for `profile_table`'s `y` vector: one file in, `X` and `y` out.

**Architecture:** A new module `DNAflexpy/encode.py` holds a pure function `encode(profiles, feature_names, normalize=True)` and the small result object `FeatureMatrix`. `FlexProfile.encode` and `ProfileSet.encode` are thin wrappers. The only change to existing code outside `encode.py` is threading the input sequences onto `FlexProfile` so one-hot sequence encoding has something to encode.

**Tech Stack:** Python 3.12+, numpy 2.x, pandas 2.x, pytest. **No new dependency and no new optional extra** — numpy is already a hard dependency.

## Global Constraints

- **Spec:** `docs/superpowers/specs/2026-08-12-dnaflex-rewrite-design.md`, the "ML encoding (`encode.py`)" section.
- **Branch:** all work on `dev`. Never commit to `main`.
- **Commit messages:** no `Co-Authored-By` trailers, no AI attribution of any kind.
- **Never modify anything under `rxv/DNAflexpy/`.** It is the frozen archive.
- **Do not break byte equality.** `python -m pytest -q` currently reports **387 passed**, including 230 byte-equality cases. Task 1 edits `_assemble`, which is on the path the byte contract flows through, so Task 1 must re-run the differential gate explicitly:
  `python -m pytest tests/test_differential.py -k byte_identical -q` → 230 passed.
- **Do not touch `DNAflexpy/kernel.py`.** Nothing in this phase needs new arithmetic in the frozen numeric core.
- **`DNAflexpy/results.py` gets exactly one change** — an optional `seqs` argument and its property. `FeatureMatrix` lives in `encode.py`, not here, so the class the byte contract runs through stays as close to unchanged as possible.
- macOS uses the `spawn` start method. Never start a `multiprocessing.Pool` from a `python - <<EOF` stdin heredoc; use a real `.py` file with an `if __name__ == "__main__":` guard, or `threads=1`.

## Decisions already made (do not re-litigate)

These were settled during planning. Each is recorded here with its reason so the implementer does not have to rediscover it.

1. **Min-max normalisation uses the observed block range, not the lookup table's theoretical range.** DNAshapeR's `normalize(x, max, min)` is min-max, and matching it is the default. Theoretical bounds were considered and rejected: for an order-*n* product block the corner-product bound is valid but badly loose (a feature spanning `[-3, 1]` gives an order-3 bound of `[-27, 27]` while real products may only span `[-1, 1]`), which would squeeze every normalised value into a sliver of `[0, 1]` and defeat the point. Every other divergence from DNAshapeR in this project was justified by a demonstrable defect; there is none here.
2. **The train/test leakage question is answered in a docstring, not with a parameter.** Min-max over the dataset's own range does leak when you later split. Say so, and point at `normalize=False` plus scaling inside a scikit-learn pipeline. Do **not** build a `"table" | "data" | False` mode — the spec asked for a boolean.
3. **NaN stays NaN. There is no `fill=` parameter.** Phase 2 deliberately replaced the archive's `0` with `NaN` so an unresolvable k-mer stops looking like a measurement. A `fill=0.0` argument would reinstate exactly that ambiguity one layer up. Point users at `sklearn.impute.SimpleImputer` in a docstring sentence.
4. **The result is a `FeatureMatrix`, not a `DataFrame`.** 3-mer encoding of 500-position sequences is ~32,000 columns; a DataFrame that wide is slow to build and slow to index. `FeatureMatrix` holds a plain `numpy` array plus the column names, and `.to_frame()` is there when a DataFrame is actually wanted.
5. **`window_size` goes into the column names.** At `window_size=10`, `DNaseI…p1` is the mean over sequence positions 1–10; at `window_size=0` it is the k-mer starting at position 1. Two matrices from profilers with different `window_size` would otherwise get identically-named columns meaning different things, and someone will eventually `pd.concat` them.
6. **Do not claim column-order parity with DNAshapeR.** The manual states the bit strings (A=0001, C=0010, G=0100, T=1000) but not the matrix layout. Pick a layout, state it in the docstring, and make no parity claim that cannot be verified from the manual.

## File structure

| File | Responsibility |
|---|---|
| `DNAflexpy/results.py` (modify) | `FlexProfile` gains an optional `seqs` argument and a `.seqs` property. Nothing else changes. |
| `DNAflexpy/core.py` (modify) | `_build` and `_assemble` thread `seqs` through. **Both** `profile_fasta` branches must pass it. `FlexProfile.encode` / `ProfileSet.encode` wrappers. |
| `DNAflexpy/encode.py` (create) | `encode`, `FeatureMatrix`, and the private block builders. |
| `DNAflexpy/__init__.py` (modify) | Re-export `encode` and `FeatureMatrix`. |
| `tests/test_encode.py` (create) | All Phase 4 tests. |
| `tests/test_profiler.py` (modify) | `seqs` is populated by every entry point, including the pooled one. |
| `docs/usage.md`, `docs/api_reference.md` (modify) | Document encoding. Examples must pass `scripts/check_doc_examples.py`. |

---

### Task 1: Thread input sequences onto `FlexProfile`

**Files:**
- Modify: `DNAflexpy/results.py`, `DNAflexpy/core.py`
- Test: `tests/test_profiler.py`

**Interfaces:**
- Produces: `FlexProfile.seqs -> list[str] | None`, aligned position-for-position with `.seqids`. Tasks 2–4 consume it.

**Why this is its own task, and where it can silently half-work.** `to_tsv` serialises `_rows` only, so **the byte-equality gate cannot detect a path that leaves `seqs` as `None`.** The pooled branch of `profile_fasta` calls `_assemble(rows_by_feature)` directly and never sees `records`; the serial branch goes through `_build`. Thread `seqs` into `_build` alone and four of five paths work, all 387 tests still pass, and `profile_fasta(threads=4)` on a file of ≥64 records quietly returns a profile that cannot encode. So `seqs` must be a parameter of **`_assemble`**, and the pooled call site must pass it.

There is no memory objection: `profile_fasta` already does `records = list(read_fasta(path))`, so the sequences are resident either way.

- [ ] **Step 1: Write the failing test**

Add to `tests/test_profiler.py`:

```python
def test_profile_seqs_retains_the_input_sequences():
    prof = FlexProfiler("gc", window_size=0).profile_seqs({"a": "ACGT", "b": "TTTT"})
    assert prof.seqs == ["ACGT", "TTTT"]
    assert prof.seqids == ["a", "b"]


def test_every_profile_in_a_set_shares_the_same_sequences():
    profiles = FlexProfiler(["gc", "wedge"], window_size=0).profile_seqs(["ACGT"])
    assert profiles["gc"].seqs == ["ACGT"]
    assert profiles["wedge"].seqs == ["ACGT"]


def test_profile_fasta_retains_sequences_serially(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">a\nACGT\n>b\nTTTT\n")
    prof = FlexProfiler("gc", window_size=0).profile_fasta(fa, threads=1)
    assert prof.seqs == ["ACGT", "TTTT"]


def test_profile_table_retains_sequences_and_labels(tmp_path):
    tsv = tmp_path / "t.tsv"
    tsv.write_text("ACGT\t1.5\nTTTT\t2.5\n")
    prof = FlexProfiler("gc", window_size=0).profile_table(tsv, header=False)
    assert prof.seqs == ["ACGT", "TTTT"]
    assert prof.y == [1.5, 2.5]


def test_from_bed_retains_sequences(tmp_path, ...):
    # reuse the existing genome fixture from tests/test_bed_input.py
    ...
```

The pooled path cannot be tested by the same trick — a `Pool` under `spawn` re-imports the module. It **can** be tested from pytest normally (pytest runs from a real file, not a heredoc), so write it as an ordinary test with ≥`_MIN_RECORDS_FOR_POOL` records:

```python
def test_profile_fasta_retains_sequences_through_the_pool(tmp_path):
    """The pooled branch builds rows in workers and never sees `records` --
    it is the one path where `seqs` can silently come back None."""
    import DNAflexpy.core as core

    n = core._MIN_RECORDS_FOR_POOL + 2
    fa = tmp_path / "many.fa"
    fa.write_text("".join(f">s{i}\nACGTACGT\n" for i in range(n)))
    prof = FlexProfiler("gc", window_size=0).profile_fasta(fa, threads=2)
    assert len(prof.seqs) == n
    assert prof.seqs[0] == "ACGTACGT"
    assert prof.seqids[-1] == f"s{n - 1}"
```

Run it and confirm it fails with `AttributeError` / `None` before writing any implementation.

- [ ] **Step 2: Implement**

In `results.py`, add `seqs=None` as the last keyword of `FlexProfile.__init__`, store `self._seqs = list(seqs) if seqs is not None else None`, and expose:

```python
    @property
    def seqs(self) -> list[str] | None:
        """The input sequences, aligned to `.seqids`, or None.

        Only sequence one-hot encoding needs these; everything else works
        from the profile values alone.
        """
        return list(self._seqs) if self._seqs is not None else None
```

In `core.py`:
- `_build(self, pairs, y=None)` → pass `seqs=[s for _, s in pairs]` to `_assemble`.
- `_assemble(self, rows_by_feature, y=None, seqs=None)` → pass `seqs=seqs` into every `FlexProfile(...)` it constructs.
- **The pooled branch of `profile_fasta`** → `return self._assemble(rows_by_feature, seqs=[s for _, s in records])`.

- [ ] **Step 3: Verify**

```bash
python -m pytest tests/test_profiler.py -q
python -m pytest tests/test_differential.py -k byte_identical -q   # must be 230 passed
python -m pytest -q
```

The differential run is not optional here: this task changes a signature on the path the byte contract flows through.

---

### Task 2: `encode.py` — the k-mer one-hot block

**Files:**
- Create: `DNAflexpy/encode.py`
- Test: `tests/test_encode.py`

**Interfaces:**
- Produces: `_one_hot(seqs: list[str], k: int) -> tuple[np.ndarray, list[str]]`. Task 4 consumes it.

**Layout (state this in the docstring, verbatim reasoning):**

- A sequence of length `L` has `L - k + 1` k-mer positions.
- Each position contributes `4**k` binary columns, one per k-mer, in **lexicographic ACGT order** (`itertools.product("ACGT", repeat=k)`): `AA, AC, AG, AT, CA, …`.
- Positions are 1-based, matching `FlexProfile.frame`'s columns.
- Column names: `seq.{k}mer.p{i}.{KMER}` — e.g. `seq.1mer.p1.A`, `seq.2mer.p3.GT`.
- Block width is `4**k * (L - k + 1)`.
- The manual gives DNAshapeR's bit strings but not its matrix layout, so **make no claim that this matches DNAshapeR's column order.**

**Behaviour:**

- Sequences are uppercased before matching, as everywhere else in the package.
- A k-mer containing any letter outside `ACGT` (`N` or any other IUPAC code) gets an **all-zero** vector at that position — the standard "unknown" encoding, and the direct analogue of the profile masking that same k-mer as `NaN`. Emit one `UserWarning` naming how many positions across how many sequences came out all-zero, in the same style as the `io.py` warnings.
- `k < 1` raises. `k > L` raises naming `L`, rather than emitting a zero-column block.
- Unequal sequence lengths are **not** checked here — Task 4 owns that check, so its error message can be raised once for the whole request.

- [ ] **Step 1: Write the failing test**

```python
import numpy as np
import pytest

from DNAflexpy.encode import _one_hot


def test_one_hot_1mer_places_each_base_in_its_own_column():
    X, cols = _one_hot(["ACGT"], 1)
    assert cols == [
        "seq.1mer.p1.A", "seq.1mer.p1.C", "seq.1mer.p1.G", "seq.1mer.p1.T",
        "seq.1mer.p2.A", "seq.1mer.p2.C", "seq.1mer.p2.G", "seq.1mer.p2.T",
        "seq.1mer.p3.A", "seq.1mer.p3.C", "seq.1mer.p3.G", "seq.1mer.p3.T",
        "seq.1mer.p4.A", "seq.1mer.p4.C", "seq.1mer.p4.G", "seq.1mer.p4.T",
    ]
    assert X.shape == (1, 16)
    # ACGT is deliberately NOT a palindrome of its own encoding: each
    # position's hot column moves, so a transposed or reversed layout fails.
    assert list(X[0]) == [1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]


def test_one_hot_2mer_uses_overlapping_positions():
    X, cols = _one_hot(["ACGT"], 2)
    assert X.shape == (1, 16 * 3)          # L - k + 1 == 3 positions
    assert cols[0] == "seq.2mer.p1.AA"
    assert cols[-1] == "seq.2mer.p3.TT"
    hot = [cols[i] for i in np.flatnonzero(X[0])]
    assert hot == ["seq.2mer.p1.AC", "seq.2mer.p2.CG", "seq.2mer.p3.GT"]


def test_one_hot_is_case_insensitive():
    upper, _ = _one_hot(["ACGT"], 1)
    lower, _ = _one_hot(["acgt"], 1)
    assert np.array_equal(upper, lower)


def test_a_kmer_with_n_becomes_an_all_zero_position():
    with pytest.warns(UserWarning, match="all-zero"):
        X, cols = _one_hot(["ANGT"], 1)
    assert list(X[0][4:8]) == [0, 0, 0, 0]      # position 2 is unknown
    assert list(X[0][0:4]) == [1, 0, 0, 0]      # position 1 is unaffected


def test_one_hot_rejects_k_larger_than_the_sequence():
    with pytest.raises(ValueError, match="k=5"):
        _one_hot(["ACGT"], 5)


def test_one_hot_rejects_k_below_one():
    with pytest.raises(ValueError, match="k must be >= 1"):
        _one_hot(["ACGT"], 0)
```

- [ ] **Step 2: Implement**

Build the index map once (`{kmer: offset}`) and fill a preallocated `np.zeros((n_seqs, width), dtype=float)`. Use `float`, not `int` or `bool`: the block is concatenated with float feature blocks, and a dtype promotion at concat time would be a second place for a surprise.

- [ ] **Step 3: Verify** — `python -m pytest tests/test_encode.py -q`

---

### Task 3: `encode.py` — nth-order feature blocks and min-max normalisation

**Files:**
- Modify: `DNAflexpy/encode.py`
- Test: `tests/test_encode.py`

**Interfaces:**
- Produces: `_feature_block(profile, order) -> tuple[np.ndarray, list[str]]` and `_minmax(block) -> np.ndarray`. Task 4 consumes both.

**Order-*n* semantics (from DNAshapeR's `encodeNstOrderShape`):** the *n*th-order block holds the **products of the same feature at *n* adjacent positions** — nothing else. Order 1 is the identity. A block of `m` values per sequence yields `m - n + 1` columns, column `i` being `prod(v[i : i + n])`. `n-<feature>` therefore gives **only** the *n*th-order terms; ask for `1-gc` and `2-gc` together if you want both.

- `order < 1` raises.
- `order > m` raises naming `m`, rather than emitting a zero-column block.
- `NaN` propagates through the product, which is correct: a product involving an unresolvable k-mer is itself unresolvable.

**Column names:** `{feature}.w{window_size}.o{order}.p{i}`, positions 1-based — e.g. `gc.w0.o2.p1`, `DNaseI.w10.o1.p1`. The `w` component is decision 5 above; do not drop it.

**Normalisation** (`_minmax`): min-max over **the whole block**, matching DNAshapeR's `normalize(x, max, min)` taking scalar bounds for a matrix — not per column.

```
scaled = (v - lo) / (hi - lo)      where lo = nanmin(block), hi = nanmax(block)
```

Edge cases, all of which must be tested:
- `NaN` entries stay `NaN`.
- `hi == lo` (a constant block): every finite entry becomes `0.0`. Do not divide.
- An all-`NaN` block: return it unchanged; do not let `np.nanmin` raise or emit its own `RuntimeWarning`. Guard with an explicit finite check before calling `nanmin`/`nanmax`.
- One-hot blocks are **never** normalised — they are already 0/1. Task 4 enforces this; `_minmax` itself is agnostic.

- [ ] **Step 1: Write the failing test**

The numbers below are real — computed from the packaged table with `FlexProfiler("wedge", window_size=0)` on these two sequences. `wedge` was chosen because its range (0.9 … 8.4) is not already `[0, 1]`, so a no-op normaliser would fail.

```python
import warnings

import numpy as np
import pytest

from DNAflexpy import FlexProfiler
from DNAflexpy.encode import _feature_block, _minmax

SEQS = {"s1": "ACGTACGT", "s2": "TTTTTTTT"}


def _wedge():
    return FlexProfiler("wedge", window_size=0).profile_seqs(SEQS)


def test_first_order_block_is_the_profile_values():
    X, cols = _feature_block(_wedge(), 1)
    np.testing.assert_allclose(X, [[1.1, 6.7, 1.1, 0.9, 1.1, 6.7, 1.1],
                                   [7.2, 7.2, 7.2, 7.2, 7.2, 7.2, 7.2]])
    assert cols[0] == "wedge.w0.o1.p1"
    assert cols[-1] == "wedge.w0.o1.p7"


def test_second_order_block_multiplies_adjacent_positions():
    X, cols = _feature_block(_wedge(), 2)
    assert X.shape == (2, 6)                       # m - n + 1
    np.testing.assert_allclose(X, [[7.37, 7.37, 0.99, 0.99, 7.37, 7.37],
                                   [51.84] * 6])
    assert cols[0] == "wedge.w0.o2.p1"
    assert cols[-1] == "wedge.w0.o2.p6"


def test_window_size_appears_in_the_column_name():
    prof = FlexProfiler("wedge", window_size=4).profile_seqs(SEQS)
    _, cols = _feature_block(prof, 1)
    assert cols[0] == "wedge.w4.o1.p1"


def test_order_beyond_the_block_width_raises():
    with pytest.raises(ValueError, match="7"):
        _feature_block(_wedge(), 8)


def test_minmax_scales_the_whole_block_not_each_column():
    X, _ = _feature_block(_wedge(), 1)
    scaled = _minmax(X)
    assert scaled.min() == 0.0 and scaled.max() == 1.0
    # lo=0.9, hi=7.2 across the WHOLE block; per-column scaling would put
    # 1.1 at 0.0 in column 1, not at 0.031746.
    np.testing.assert_allclose(scaled[0][0], (1.1 - 0.9) / (7.2 - 0.9))
    np.testing.assert_allclose(scaled[1][0], 1.0)


def test_minmax_keeps_nan_as_nan():
    block = np.array([[1.0, np.nan], [3.0, 5.0]])
    scaled = _minmax(block)
    assert np.isnan(scaled[0][1])
    np.testing.assert_allclose(scaled[0][0], 0.0)
    np.testing.assert_allclose(scaled[1][1], 1.0)


def test_minmax_of_a_constant_block_is_zero_not_a_division_error():
    scaled = _minmax(np.array([[4.0, 4.0], [4.0, 4.0]]))
    np.testing.assert_allclose(scaled, np.zeros((2, 2)))


def test_minmax_of_an_all_nan_block_is_unchanged_and_silent():
    block = np.full((2, 2), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # any RuntimeWarning fails
        scaled = _minmax(block)
    assert np.isnan(scaled).all()
```

- [ ] **Step 2: Implement**

`_feature_block` reads `profile._rows`, builds `V = np.array([r[1:] for r in rows], dtype=float)`, and raises if the rows are not all the same width (message must name the distinct widths). Order *n* is `np.prod(sliding_window_view(V, n, axis=1), axis=-1)` or an explicit loop of `n` multiplications — either is fine; there is no byte contract here.

- [ ] **Step 3: Verify** — `python -m pytest tests/test_encode.py -q`

---

### Task 4: `encode()`, `FeatureMatrix`, the wrappers, and docs

**Files:**
- Modify: `DNAflexpy/encode.py`, `DNAflexpy/core.py`, `DNAflexpy/results.py`, `DNAflexpy/__init__.py`, `docs/usage.md`, `docs/api_reference.md`
- Test: `tests/test_encode.py`

**Interfaces:**
- Produces: `encode(profiles, feature_names, normalize=True) -> FeatureMatrix`, plus `FlexProfile.encode(...)` and `ProfileSet.encode(...)` wrappers.

**Feature-name grammar:** `"{n}-{what}"`.
- `what == "mer"` → sequence one-hot with `k = n` (`1-mer`, `2-mer`, `3-mer`).
- otherwise → `what` is a feature name present in the input, and `n` is the order (`1-DNaseI`, `2-stiffness`).
- Anything else raises, quoting the expected grammar and the offending string.
- If the input profiles contain a feature literally named `mer`, `n-mer` is ambiguous — raise saying so instead of guessing. (No packaged feature is named `mer`; a user-supplied table could be.)
- A requested feature that is not in the input raises **listing what is available**, matching `FeatureTable._require`'s idiom.
- An empty `feature_names` raises.
- A duplicated name raises: duplicate column names would silently break `.to_frame()`.

**Validation, done once for the whole request:**
- Accepts a `FlexProfile`, a `ProfileSet`, or a plain `{feature: FlexProfile}` dict.
- Every profile must carry the **same** `seqids`, in the same order. Raise if not.
- **Equal sequence length is checked on the sequences, not on the row widths.** These are different checks: two sequences both shorter than the window each produce an empty row, so equal row widths do not imply equal sequence lengths. Validate `{len(s) for s in seqs}`, raise naming the distinct lengths found, and point at `from_bed(width=...)` as the spec specifies.
- `seqs` is required **only** when a `k-mer` term is requested. If it is `None` then — the profile came from somewhere that did not retain sequences — raise a message saying which entry points do retain them. When no `k-mer` term is requested, validate equal row widths per feature block instead.

**Block widths differ between blocks, and that is fine.** At `window_size > 0` every feature yields `L - window_size + 1` values, so all feature blocks are the same width. At `window_size = 0` a block is `L - k + 1` wide and therefore differs between a 2-mer feature (`gc`) and a 3-mer one (`DNaseI`). A `1-mer` sequence block is `L` positions wide in either case. Blocks are simply concatenated left to right in the order requested; the column names disambiguate. Document this.

**`FeatureMatrix`** — a plain class, matching `FlexProfile`'s style (the codebase uses no dataclasses):

```python
class FeatureMatrix:
    """A design matrix and its column names, ready for scikit-learn."""

    def __init__(self, X, columns, seqids, y=None, window_size=None, feature_names=None):
        ...

    @property
    def shape(self): return self.X.shape

    def to_frame(self) -> pd.DataFrame:
        """A DataFrame view. Wide for k-mer encodings -- 3-mer over 500
        positions is ~32,000 columns -- which is why `.X` is the default."""
        return pd.DataFrame(self.X, index=self.seqids, columns=self.columns)

    def __repr__(self): ...      # e.g. <FeatureMatrix 12 x 448, features=['1-mer', '1-gc']>
```

`y` is carried straight through from the input profile, so `profile_table(...).encode([...])` gives `X` and `y` together. Every profile in a `ProfileSet` shares one `y`.

**Docstring must state, in one sentence each:**
- Min-max uses this dataset's own range, so pass `normalize=False` and scale inside a scikit-learn pipeline if you are splitting train and test.
- `NaN` marks a position that could not be resolved; use `sklearn.impute.SimpleImputer` if the model needs it filled.
- `DNAflexpy.profile()` returns a bare array, so encoding starts from `profile_seqs`, `profile_fasta`, `profile_table` or `from_bed`.

- [ ] **Step 1: Write the failing test**

```python
def test_encode_concatenates_blocks_in_the_order_requested():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    fm = prof.encode(["1-mer", "1-gc"], normalize=False)
    assert fm.shape == (2, 4 * 8 + 7)
    assert fm.columns[0] == "seq.1mer.p1.A"
    assert fm.columns[32] == "gc.w0.o1.p1"
    assert fm.seqids == ["s1", "s2"]


def test_reversing_the_request_reverses_the_blocks():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    fm = prof.encode(["1-gc", "1-mer"], normalize=False)
    assert fm.columns[0] == "gc.w0.o1.p1"


def test_one_hot_columns_are_never_normalized():
    prof = FlexProfiler("wedge", window_size=0).profile_seqs(SEQS)
    fm = prof.encode(["1-mer", "1-wedge"], normalize=True)
    onehot = fm.X[:, :32]
    assert set(np.unique(onehot)) <= {0.0, 1.0}
    wedge = fm.X[:, 32:]
    assert wedge.min() == 0.0 and wedge.max() == 1.0


def test_encode_across_a_profile_set():
    profiles = FlexProfiler(["gc", "DNaseI"], window_size=4).profile_seqs(SEQS)
    fm = profiles.encode(["1-gc", "2-DNaseI"], normalize=False)
    assert any(c.startswith("gc.w4.") for c in fm.columns)
    assert any(c.startswith("DNaseI.w4.o2.") for c in fm.columns)


def test_encode_carries_the_label_vector(tmp_path):
    tsv = tmp_path / "t.tsv"
    tsv.write_text("ACGTACGT\t1.5\nTTTTTTTT\t2.5\n")
    prof = FlexProfiler("gc", window_size=0).profile_table(tsv, header=False)
    fm = prof.encode(["1-gc"])
    assert fm.y == [1.5, 2.5]
    assert fm.X.shape[0] == len(fm.y)


def test_unequal_sequence_lengths_raise_and_point_at_from_bed():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(["ACGTAC", "ACGT"])
    with pytest.raises(ValueError, match="from_bed"):
        prof.encode(["1-mer"])


def test_two_short_sequences_have_equal_rows_but_unequal_lengths():
    """Both rows are empty, so a row-width check would wrongly pass."""
    prof = FlexProfiler("gc", window_size=10).profile_seqs(["ACGTA", "ACG"])
    assert [len(r) - 1 for r in prof._rows] == [0, 0]
    with pytest.raises(ValueError, match="5|3"):
        prof.encode(["1-mer"])


def test_unknown_feature_lists_what_is_available():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    with pytest.raises(ValueError, match="gc"):
        prof.encode(["1-nosuch"])


def test_malformed_feature_name_raises():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    for bad in ["gc", "mer", "x-gc", "0-gc", "-gc"]:
        with pytest.raises(ValueError):
            prof.encode([bad])


def test_duplicate_feature_names_raise():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    with pytest.raises(ValueError, match="1-gc"):
        prof.encode(["1-gc", "1-gc"])


def test_empty_request_raises():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    with pytest.raises(ValueError):
        prof.encode([])


def test_a_feature_named_mer_is_reported_as_ambiguous():
    table = {"mer": {"AA": 1.0, "AC": 2.0, "AG": 3.0, "AT": 4.0,
                     "CA": 1.0, "CC": 2.0, "CG": 3.0, "CT": 4.0,
                     "GA": 1.0, "GC": 2.0, "GG": 3.0, "GT": 4.0,
                     "TA": 1.0, "TC": 2.0, "TG": 3.0, "TT": 4.0}}
    prof = FlexProfiler("mer", window_size=0, lookup=table).profile_seqs(SEQS)
    with pytest.raises(ValueError, match="ambiguous"):
        prof.encode(["1-mer"])


def test_to_frame_round_trips_the_columns():
    prof = FlexProfiler("gc", window_size=0).profile_seqs(SEQS)
    frame = prof.encode(["1-gc"], normalize=False).to_frame()
    assert list(frame.columns) == ["gc.w0.o1.p" + str(i) for i in range(1, 8)]
    assert list(frame.index) == ["s1", "s2"]


def test_encode_needs_sequences_for_a_kmer_term():
    from DNAflexpy.results import FlexProfile
    bare = FlexProfile([["a", 1.0, 2.0]], feature="gc", window_size=0, kmer_len=2)
    with pytest.raises(ValueError, match="profile_seqs|profile_fasta"):
        bare.encode(["1-mer"])
    # ... but a feature-only request still works
    assert bare.encode(["1-gc"], normalize=False).shape == (1, 2)
```

- [ ] **Step 2: Implement**

Add `encode` and `FeatureMatrix` to `encode.py`; add thin wrappers:

```python
# results.py
    def encode(self, feature_names, normalize: bool = True):
        from DNAflexpy.encode import encode
        return encode(self, feature_names, normalize=normalize)
```

and the same on `ProfileSet` in `core.py`. Import inside the method, matching how `core.py` already defers its `io` imports.

Re-export from `__init__.py`: add `encode` and `FeatureMatrix` to the imports and to `__all__`. Note `DNAflexpy.encode` will then be the **function**, not the submodule — exactly the situation `results.py` was named to avoid for `profile`. Add a one-line comment there saying so, and document `from DNAflexpy.encode import FeatureMatrix` as the way to reach the module.

- [ ] **Step 3: Documentation**

Add an "Encoding for machine learning" section to `docs/usage.md` with a runnable end-to-end example — `profile_table` → `encode` → a scikit-learn-shaped `X`/`y` handoff — written so it **does not import sklearn** (not a dependency); print the shapes instead. Add the `encode` / `FeatureMatrix` entries to `docs/api_reference.md`.

- [ ] **Step 4: Verify**

```bash
python -m pytest -q                                                # all green
python -m pytest tests/test_differential.py -k byte_identical -q   # 230 passed
python scripts/check_doc_examples.py                               # docs examples run
```

---

## Definition of done

- [ ] `FlexProfile.seqs` is populated by `profile_seqs`, `profile_fasta` (serial **and** pooled), `profile_table` and `from_bed`, each with a test.
- [ ] `encode.py` exists with `encode`, `FeatureMatrix`, `_one_hot`, `_feature_block`, `_minmax`.
- [ ] `1-mer` / `2-mer` / `3-mer` one-hot and `n-<feature>` order terms both work, and can be mixed in one request.
- [ ] `normalize=True` min-max scales feature blocks and leaves one-hot blocks alone.
- [ ] `NaN` survives encoding; no `fill=` parameter exists.
- [ ] Unequal sequence lengths raise, pointing at `from_bed(width=...)`.
- [ ] `y` rides through from `profile_table` to `FeatureMatrix.y`.
- [ ] `python -m pytest -q` is green and the 230 byte-equality cases still pass.
- [ ] `scripts/check_doc_examples.py` passes.

## Deferred to later phases

- Plotting (Phase 5) — `heatmap`, `metaprofile`, `trackplot`, and the `plot` extra.
- Dinucleotide-shuffle backgrounds and per-position statistics (Phase 6).
- Feature-redundancy reporting and `ProfileSet.correlate()` (Phase 7). Highly correlated feature blocks are a real concern for a linear model built on this matrix, but the reporting belongs with the rest of the statistics work.
- The CLI (Phase 8) — no `--encode` flag until then.
