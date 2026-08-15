# DNAflexpy Phase 3 — Labelled Table and BED Inputs

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add two more ways to get sequences into `FlexProfiler` — a labelled table (sequence + numeric value, e.g. binding affinity) and a BED file read against a reference genome.

**Architecture:** Both new readers live in `DNAflexpy/io.py` and return plain Python data. `FlexProfiler` gains two thin methods that call a reader and hand the result to the existing `_build`/`_assemble` path, so all four input types converge on the same `FlexProfile`. The labelled table also carries a `y` vector, which `_assemble` threads onto every `FlexProfile` it builds.

**Tech Stack:** Python 3.11+, pandas 2.3, numpy 2.2, pyfaidx 0.9 (optional extra), pytest 9.

## Global Constraints

- **Spec:** `docs/superpowers/specs/2026-08-12-dnaflex-rewrite-design.md`, the "Input formats (`io.py`)" section. This plan covers only `profile_table` and `from_bed`; the rest of Phase 3 already shipped in the previous plan.
- **Branch:** all work on `dev`. Never commit to `main`.
- **Commit messages:** no `Co-Authored-By` trailers, no AI attribution of any kind.
- **Never modify anything under `rxv/DNAflexpy/`.** It is the frozen archive. `tests/test_differential.py` byte-compares against it.
- **Do not break byte equality.** `python -m pytest -q` currently reports 319 passed, including 230 byte-equality cases. Nothing in this plan should change any existing numeric output. If the differential gate goes red, stop and report it.
- **Do not modify `Notebooks/*.ipynb`.**
- **`pyfaidx` is an optional extra.** Never import it at module top level — `DNAflexpy/io.py` is imported by `DNAflexpy/core.py`, so a top-level import would break every install that did not ask for the `bed` extra. Import it inside the function that needs it and raise a clear message if it is missing.
- **One correction to the spec.** The spec says `on_edge="pad"` pads with `N`, "which the lookup zero-fills". That is out of date: Phase 2 changed unresolvable k-mers to resolve to `NaN`, not `0`. Padded positions are therefore **masked**, and show up in `FlexProfile.n_masked`. Say that, not the old wording.
- macOS uses the `spawn` start method. Never start a `multiprocessing.Pool` from a `python - <<EOF` stdin heredoc; use a real `.py` file with an `if __name__ == "__main__":` guard, or `threads=1`.

## File structure

| File | Responsibility |
|---|---|
| `DNAflexpy/io.py` (modify) | All input readers. Gains `read_table`, `read_bed`, `extract_intervals`, and two private helpers. `read_fasta` is unchanged. |
| `DNAflexpy/core.py` (modify) | `FlexProfiler` gains `profile_table` and `from_bed`. `_build` and `_assemble` gain an optional `y` parameter. |
| `DNAflexpy/results.py` | Unchanged. `FlexProfile.y` already exists and is already accepted by `__init__`; it is simply never filled in today. |
| `tests/test_table_input.py` (create) | Tests for `read_table` and `profile_table`. |
| `tests/test_bed_input.py` (create) | Tests for `read_bed`, `extract_intervals`, and `from_bed`. |
| `pyproject.toml` (modify) | Adds the `bed` optional extra. |

---

### Task 1: `read_table` — parse a labelled table

**Files:**
- Modify: `DNAflexpy/io.py`
- Test: `tests/test_table_input.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `read_table(path, seq_col=0, value_col=1, id_col=None, header=None, sep="\t") -> list[tuple[str, str, float]]`, returning `(seqid, sequence, value)` triples in file order. Task 2 consumes it.

Behaviour required by the spec:
- Header row auto-detected: if the first row's `value_col` cannot be parsed as a float, it is a header. Forceable with `header=True|False`.
- Columns addressable by position (`seq_col=0`) or by name (`seq_col="sequence"`). Named columns imply a header.
- Rows get generated IDs `seq_0, seq_1, …` unless `id_col=` names one.
- Non-numeric or missing values **raise**, naming the offending row. Silently dropping labelled rows would corrupt a training set.
- Tab by default; `sep=` allows CSV.

- [ ] **Step 1: Write the failing test**

Create `tests/test_table_input.py`:

```python
import pytest

from DNAflexpy.io import read_table


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_reads_sequence_and_value(tmp_path):
    p = write(tmp_path, "a.tsv", "ATGC\t1.5\nGGTT\t2.5\n")
    assert read_table(p) == [("seq_0", "ATGC", 1.5), ("seq_1", "GGTT", 2.5)]


def test_detects_a_header_row(tmp_path):
    p = write(tmp_path, "h.tsv", "sequence\taffinity\nATGC\t1.5\n")
    assert read_table(p) == [("seq_0", "ATGC", 1.5)]


def test_detects_absence_of_a_header(tmp_path):
    p = write(tmp_path, "n.tsv", "ATGC\t1.5\n")
    assert len(read_table(p)) == 1


def test_header_can_be_forced_off(tmp_path):
    """A first data row is not lost when header=False is explicit."""
    p = write(tmp_path, "f.tsv", "ATGC\t1.5\nGGTT\t2.5\n")
    assert len(read_table(p, header=False)) == 2


def test_columns_addressable_by_name(tmp_path):
    p = write(tmp_path, "c.tsv", "affinity\tsequence\n1.5\tATGC\n")
    assert read_table(p, seq_col="sequence", value_col="affinity") == [
        ("seq_0", "ATGC", 1.5)
    ]


def test_id_column_is_used_when_named(tmp_path):
    p = write(tmp_path, "i.tsv", "name\tsequence\taffinity\npeak1\tATGC\t1.5\n")
    assert read_table(p, seq_col="sequence", value_col="affinity", id_col="name") == [
        ("peak1", "ATGC", 1.5)
    ]


def test_csv_via_sep(tmp_path):
    p = write(tmp_path, "a.csv", "ATGC,1.5\n")
    assert read_table(p, sep=",") == [("seq_0", "ATGC", 1.5)]


def test_non_numeric_value_raises_naming_the_row(tmp_path):
    """Dropping the row instead would silently corrupt a training set."""
    p = write(tmp_path, "b.tsv", "ATGC\t1.5\nGGTT\thigh\n")
    with pytest.raises(ValueError, match="row 1"):
        read_table(p)


def test_missing_value_raises(tmp_path):
    p = write(tmp_path, "m.tsv", "ATGC\t1.5\nGGTT\t\n")
    with pytest.raises(ValueError, match="row 1"):
        read_table(p)


def test_out_of_range_column_raises(tmp_path):
    p = write(tmp_path, "o.tsv", "ATGC\t1.5\n")
    with pytest.raises(ValueError, match="value_col"):
        read_table(p, value_col=9)


def test_unknown_named_column_raises(tmp_path):
    p = write(tmp_path, "u.tsv", "sequence\taffinity\nATGC\t1.5\n")
    with pytest.raises(ValueError, match="nope"):
        read_table(p, seq_col="nope", value_col="affinity")


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_table(tmp_path / "nope.tsv")


def test_empty_file_raises(tmp_path):
    p = write(tmp_path, "e.tsv", "")
    with pytest.raises(ValueError, match="no data rows"):
        read_table(p)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_table_input.py -q`
Expected: FAIL with `ImportError: cannot import name 'read_table' from 'DNAflexpy.io'`.

- [ ] **Step 3: Implement `read_table` in `DNAflexpy/io.py`**

Add these imports at the top of the file, alongside the existing ones:

```python
import math

import pandas as pd
```

Then append:

```python
def read_table(path, seq_col=0, value_col=1, id_col=None, header=None, sep="\t"):
    """Read a labelled table of sequences and numeric values.

    Returns `(seqid, sequence, value)` triples in file order. This is the
    natural input for model fitting, since it carries X and y in one file.

    A row whose value is missing or non-numeric raises rather than being
    dropped: silently discarding labelled rows would corrupt a training set.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"table not found: {path}")

    if header is None:
        header = _looks_like_header(path, seq_col, value_col, id_col, sep)

    frame = pd.read_csv(path, sep=sep, header=0 if header else None, dtype=str)
    if frame.empty:
        raise ValueError(f"table has no data rows: {path}")

    seqs = _pick_column(frame, seq_col, "seq_col")
    values = _pick_column(frame, value_col, "value_col")
    ids = _pick_column(frame, id_col, "id_col") if id_col is not None else None

    out = []
    for i in range(len(frame)):
        raw = values.iloc[i]
        if raw is None or (isinstance(raw, float) and math.isnan(raw)) or str(raw).strip() == "":
            raise ValueError(f"row {i} of {path} has a missing value in value_col")
        try:
            value = float(raw)
        except (TypeError, ValueError):
            raise ValueError(
                f"row {i} of {path} has non-numeric value {raw!r} in value_col"
            ) from None
        seqid = str(ids.iloc[i]) if ids is not None else f"seq_{i}"
        out.append((seqid, str(seqs.iloc[i]), value))
    return out


def _pick_column(frame, col, argname):
    """Select a column by position or by name, with a useful error."""
    if isinstance(col, int):
        if col >= len(frame.columns):
            raise ValueError(
                f"{argname}={col} is out of range; the table has "
                f"{len(frame.columns)} column(s)"
            )
        return frame.iloc[:, col]
    if col not in frame.columns:
        raise ValueError(
            f"{argname}={col!r} not found; columns are {list(frame.columns)}"
        )
    return frame[col]


def _looks_like_header(path, seq_col, value_col, id_col, sep):
    """True when the first non-blank row looks like column names.

    Named columns require a header by definition. Otherwise the test is
    whether the first row's value column parses as a float.
    """
    if any(isinstance(c, str) for c in (seq_col, value_col, id_col) if c is not None):
        return True
    with Path(path).open() as handle:
        for line in handle:
            if line.strip():
                fields = line.rstrip("\n").split(sep)
                break
        else:
            return False
    if value_col >= len(fields):
        return False
    try:
        float(fields[value_col])
    except ValueError:
        return True
    return False
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_table_input.py -q`
Expected: PASS, 13 passed.

- [ ] **Step 5: Run the whole suite**

Run: `python -m pytest -q`
Expected: 332 passed. The 230 byte-equality cases must still pass — this task adds a reader and touches nothing on the numeric path.

- [ ] **Step 6: Commit**

```bash
git add DNAflexpy/io.py tests/test_table_input.py
git commit -m "Add read_table for labelled sequence/value input

Reads a table carrying a sequence and a numeric label per row, which is
the natural input for model fitting since it holds X and y in one file.

Header presence is auto-detected by trying to parse the first row's value
column as a float, and named columns imply a header. Rows with a missing
or non-numeric value raise and name the row rather than being dropped -
silently discarding labelled rows would corrupt a training set."
```

---

### Task 2: `FlexProfiler.profile_table` and the `y` vector

**Files:**
- Modify: `DNAflexpy/core.py`
- Test: `tests/test_table_input.py`

**Interfaces:**
- Consumes: `read_table(...)` from Task 1; the existing `FlexProfiler._build(pairs)` and `FlexProfiler._assemble(rows_by_feature)`; `FlexProfile(rows, feature, window_size, kmer_len, y=None)`.
- Produces: `FlexProfiler.profile_table(path, seq_col=0, value_col=1, id_col=None, header=None, sep="\t") -> FlexProfile | ProfileSet`, whose results carry `.y` as a `list[float]` aligned to `.seqids`. Later phases (encoding, plotting) read `.y`.

`_build` and `_assemble` both gain an optional `y=None`. Existing callers (`profile_seqs`, `profile_fasta`) pass nothing and keep getting `y=None`, which is what FASTA input should produce.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_table_input.py`:

```python
import numpy as np

from DNAflexpy import FlexProfiler
from DNAflexpy.core import ProfileSet
from DNAflexpy.results import FlexProfile

SEQ_A = "ATGCGTACGTAGCTAGCGTAGCTAGT"
SEQ_B = "CGTAGCTAGTACGATCGTACGTAGCT"


def test_profile_table_returns_a_profile_with_y(tmp_path):
    p = write(tmp_path, "t.tsv", f"{SEQ_A}\t1.5\n{SEQ_B}\t2.5\n")
    prof = FlexProfiler("DNaseI", window_size=10).profile_table(p)
    assert isinstance(prof, FlexProfile)
    assert prof.y == [1.5, 2.5]
    assert prof.seqids == ["seq_0", "seq_1"]


def test_y_is_aligned_to_seqids(tmp_path):
    p = write(tmp_path, "t.tsv", f"name\tseq\tv\na\t{SEQ_A}\t9.0\nb\t{SEQ_B}\t8.0\n")
    prof = FlexProfiler("DNaseI", window_size=10).profile_table(
        p, seq_col="seq", value_col="v", id_col="name"
    )
    assert dict(zip(prof.seqids, prof.y)) == {"a": 9.0, "b": 8.0}


def test_fasta_input_still_has_no_y():
    """y is only meaningful for labelled input."""
    prof = FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ_A])
    assert prof.y is None


def test_values_match_the_plain_sequence_path(tmp_path):
    """Reading via a table must not change the numbers."""
    p = write(tmp_path, "t.tsv", f"{SEQ_A}\t1.5\n")
    viafile = FlexProfiler("DNaseI", window_size=10).profile_table(p)
    direct = FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ_A])
    assert np.array_equal(viafile.frame.to_numpy(), direct.frame.to_numpy())


def test_multi_feature_table_gives_every_profile_the_same_y(tmp_path):
    p = write(tmp_path, "t.tsv", f"{SEQ_A}\t1.5\n{SEQ_B}\t2.5\n")
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).profile_table(p)
    assert isinstance(ps, ProfileSet)
    assert ps["DNaseI"].y == [1.5, 2.5]
    assert ps["trx"].y == [1.5, 2.5]
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_table_input.py -q -k profile_table`
Expected: FAIL with `AttributeError: 'FlexProfiler' object has no attribute 'profile_table'`.

- [ ] **Step 3: Thread `y` through `_build` and `_assemble`**

In `DNAflexpy/core.py`, change `_build` and `_assemble` to accept an optional `y`. Replace their signatures and the `FlexProfile(...)` construction:

```python
    def _build(self, pairs, y=None):
        """Turn (seqid, sequence) pairs into a FlexProfile or ProfileSet."""
        rows_by_feature = {
            feature: [[sid, *self._values(feature, seq)] for sid, seq in pairs]
            for feature in self._features
        }
        return self._assemble(rows_by_feature, y=y)

    def _assemble(self, rows_by_feature: dict[str, list[list]], y=None):
        """Turn per-feature row lists into a FlexProfile or ProfileSet.

        `y` is the aligned label vector from labelled input, or None. Every
        profile in a multi-feature set shares the same labels.
        """
        built = {
            feature: FlexProfile(
                rows,
                feature=feature,
                window_size=self.window_size,
                kmer_len=self._table.kmer_len(feature),
                y=y,
            )
            for feature, rows in rows_by_feature.items()
        }
        if len(self._features) == 1:
            return built[self._features[0]]
        return ProfileSet(built)
```

Keep the rest of both methods exactly as they are. `profile_seqs` and `profile_fasta` call these without `y`, so they keep producing `y=None`.

- [ ] **Step 4: Add `profile_table`**

Add this method to `FlexProfiler`, directly after `profile_fasta`:

```python
    def profile_table(self, path, seq_col=0, value_col=1, id_col=None,
                      header=None, sep="\t"):
        """Profile every sequence in a labelled table.

        The table's numeric column becomes `FlexProfile.y`, aligned to
        `.seqids`, so one file supplies both the features and the targets.
        """
        from DNAflexpy.io import read_table

        records = read_table(
            path, seq_col=seq_col, value_col=value_col, id_col=id_col,
            header=header, sep=sep,
        )
        pairs = [(seqid, seq) for seqid, seq, _ in records]
        y = [value for _, _, value in records]
        return self._build(pairs, y=y)
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_table_input.py -q`
Expected: PASS, 18 passed.

- [ ] **Step 6: Run the whole suite**

Run: `python -m pytest -q`
Expected: 337 passed, with the 230 byte-equality cases still green.

- [ ] **Step 7: Commit**

```bash
git add DNAflexpy/core.py tests/test_table_input.py
git commit -m "Add FlexProfiler.profile_table with an aligned label vector

Profiling a labelled table now fills FlexProfile.y with the table's
numeric column, aligned to seqids, so one file supplies both the features
and the targets for model fitting.

_build and _assemble take an optional y and pass it to every FlexProfile
they construct, so a multi-feature run shares one label vector. FASTA and
in-memory input pass nothing and keep y as None, which is correct - a
label only exists for labelled input."
```

---

### Task 3: BED parsing and genome extraction

**Files:**
- Modify: `DNAflexpy/io.py`
- Test: `tests/test_bed_input.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `read_bed(path) -> list[tuple[str, int, int, str | None, str]]` — `(chrom, start, end, name, strand)`, start/end 0-based half-open as BED defines them, strand defaulting to `"+"`.
  - `extract_intervals(intervals, genome, width=None, on_edge="drop") -> list[tuple[str, str]]` — `(seqid, sequence)` pairs ready for `_build`.

`pyfaidx` must be imported **inside** `extract_intervals`, never at module level.

Behaviour required by the spec:
- `width=None` uses each interval as-is; `width=N` re-centres each interval to exactly `N` bases around its midpoint.
- `-` strand intervals are reverse-complemented.
- Contig boundaries are handled by `on_edge`, default `"drop"`: `"drop"` skips and warns with a count, `"error"` raises naming the intervals, `"pad"` pads with `N` and warns. **Padded positions become NaN-masked, not zero** — Phase 2 changed unresolvable k-mers to `NaN`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_bed_input.py`:

```python
import pytest

from DNAflexpy.io import extract_intervals, read_bed

GENOME = ">chr1\n" + "ACGT" * 25 + "\n>chr2\n" + "AAAACCCCGGGGTTTT" + "\n"


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(GENOME)
    return path


def write_bed(tmp_path, text):
    path = tmp_path / "peaks.bed"
    path.write_text(text)
    return path


def test_read_bed_parses_the_core_columns(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\t20\n")
    assert read_bed(p) == [("chr1", 10, 20, None, "+")]


def test_read_bed_reads_name_and_strand(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\t20\tpeakA\t0\t-\n")
    assert read_bed(p) == [("chr1", 10, 20, "peakA", "-")]


def test_read_bed_skips_track_and_comment_lines(tmp_path):
    p = write_bed(tmp_path, "# a comment\ntrack name=x\nchr1\t10\t20\n\n")
    assert len(read_bed(p)) == 1


def test_read_bed_rejects_a_short_line(tmp_path):
    p = write_bed(tmp_path, "chr1\t10\n")
    with pytest.raises(ValueError, match="at least 3"):
        read_bed(p)


def test_extract_uses_the_interval_as_is_without_width(genome):
    out = extract_intervals([("chr1", 0, 4, None, "+")], genome)
    assert out == [("chr1:0-4", "ACGT")]


def test_extract_recentres_to_a_fixed_width(genome):
    """A 2bp interval widened to 8 stays centred on the same midpoint."""
    out = extract_intervals([("chr1", 10, 12, None, "+")], genome, width=8)
    assert len(out[0][1]) == 8


def test_extract_uses_the_bed_name_as_the_id(genome):
    out = extract_intervals([("chr1", 0, 4, "peakA", "+")], genome)
    assert out[0][0] == "peakA"


def test_minus_strand_is_reverse_complemented(genome):
    plus = extract_intervals([("chr2", 0, 4, None, "+")], genome)[0][1]
    minus = extract_intervals([("chr2", 0, 4, None, "-")], genome)[0][1]
    assert plus == "AAAA"
    assert minus == "TTTT"


def test_all_extracted_sequences_have_equal_length_with_width(genome):
    """Fixed width is what makes downstream position-wise comparison valid."""
    intervals = [("chr1", 10, 12, None, "+"), ("chr1", 40, 60, None, "+")]
    out = extract_intervals(intervals, genome, width=10)
    assert {len(s) for _, s in out} == {10}


def test_on_edge_drop_skips_and_warns(genome):
    """An interval at position 0 cannot be centred in 20bp."""
    intervals = [("chr1", 0, 2, None, "+"), ("chr1", 40, 42, None, "+")]
    with pytest.warns(UserWarning, match="dropped"):
        out = extract_intervals(intervals, genome, width=20, on_edge="drop")
    assert len(out) == 1


def test_on_edge_error_raises_naming_the_interval(genome):
    with pytest.raises(ValueError, match="chr1"):
        extract_intervals([("chr1", 0, 2, None, "+")], genome, width=20,
                          on_edge="error")


def test_on_edge_pad_pads_with_n_and_warns(genome):
    with pytest.warns(UserWarning, match="padded"):
        out = extract_intervals([("chr1", 0, 2, None, "+")], genome, width=20,
                                on_edge="pad")
    assert len(out[0][1]) == 20
    assert "N" in out[0][1]


def test_unknown_on_edge_raises(genome):
    with pytest.raises(ValueError, match="on_edge"):
        extract_intervals([("chr1", 0, 4, None, "+")], genome, on_edge="shrug")


def test_unknown_contig_raises(genome):
    with pytest.raises(ValueError, match="chr99"):
        extract_intervals([("chr99", 0, 4, None, "+")], genome)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_bed_input.py -q`
Expected: FAIL with `ImportError: cannot import name 'extract_intervals' from 'DNAflexpy.io'`.

- [ ] **Step 3: Implement `read_bed` and `extract_intervals`**

Append to `DNAflexpy/io.py`. Add `import warnings` to the imports at the top.

```python
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_bed(path):
    """Read a BED file into `(chrom, start, end, name, strand)` tuples.

    Coordinates stay exactly as BED defines them: 0-based, half-open.
    `track`, `browser` and `#` lines are skipped, as are blank lines.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"BED file not found: {path}")

    out = []
    with path.open() as handle:
        for number, line in enumerate(handle):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"line {number} of {path} has {len(fields)} column(s); "
                    "BED needs at least 3 (chrom, start, end)"
                )
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 and fields[3] != "." else None
            strand = fields[5] if len(fields) > 5 and fields[5] in ("+", "-") else "+"
            out.append((chrom, start, end, name, strand))
    return out


def extract_intervals(intervals, genome, width=None, on_edge="drop"):
    """Pull sequences for BED intervals out of a reference genome FASTA.

    With `width`, every interval is re-centred on its midpoint and cut to
    exactly that many bases, so all outputs are the same length. That is
    what makes position-wise comparison downstream meaningful.

    A centred window can run off the end of a contig. `on_edge` decides:
    `"drop"` skips it and warns with a count, `"error"` raises naming the
    intervals, `"pad"` pads with `N`. Padded positions do not resolve to a
    k-mer, so they are masked as NaN and counted in `FlexProfile.n_masked`
    - they are not silently scored as zero.
    """
    if on_edge not in ("drop", "error", "pad"):
        raise ValueError(
            f"on_edge must be 'drop', 'error' or 'pad', got {on_edge!r}"
        )
    try:
        from pyfaidx import Fasta
    except ImportError:
        raise ImportError(
            "reading BED input needs pyfaidx, which is an optional extra. "
            "Install it with: pip install 'DNAflexpy[bed]'"
        ) from None

    fasta = Fasta(str(genome))
    out, dropped, padded = [], [], []

    for chrom, start, end, name, strand in intervals:
        if chrom not in fasta:
            raise ValueError(
                f"contig {chrom!r} is not in {genome}; "
                f"it has {len(fasta.keys())} contig(s)"
            )
        length = len(fasta[chrom])

        if width is not None:
            centre = (start + end) // 2
            start = centre - width // 2
            end = start + width

        seqid = name if name is not None else f"{chrom}:{start}-{end}"

        if start < 0 or end > length:
            if on_edge == "error":
                raise ValueError(
                    f"interval {seqid} runs past the bounds of {chrom} "
                    f"(length {length}); pass on_edge='drop' or 'pad'"
                )
            if on_edge == "drop":
                dropped.append(seqid)
                continue
            left = "N" * max(0, -start)
            right = "N" * max(0, end - length)
            body = str(fasta[chrom][max(0, start):min(end, length)])
            sequence = left + body + right
            padded.append(seqid)
        else:
            sequence = str(fasta[chrom][start:end])

        if strand == "-":
            sequence = sequence.translate(_COMPLEMENT)[::-1]
        out.append((seqid, sequence))

    if dropped:
        warnings.warn(
            f"{len(dropped)} interval(s) dropped for running past a contig "
            f"boundary: {', '.join(dropped[:5])}"
            + (" ..." if len(dropped) > 5 else ""),
            UserWarning,
            stacklevel=2,
        )
    if padded:
        warnings.warn(
            f"{len(padded)} interval(s) padded with N at a contig boundary. "
            "Padded positions do not resolve to a k-mer and are masked as "
            "NaN, so those windows average fewer values.",
            UserWarning,
            stacklevel=2,
        )
    return out
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_bed_input.py -q`
Expected: PASS, 14 passed.

pyfaidx writes a `.fai` index next to the FASTA on first open. The tests use `tmp_path`, so that index lands in a temporary directory.

- [ ] **Step 5: Run the whole suite**

Run: `python -m pytest -q`
Expected: 351 passed, with the 230 byte-equality cases still green.

- [ ] **Step 6: Commit**

```bash
git add DNAflexpy/io.py tests/test_bed_input.py
git commit -m "Add BED parsing and genome interval extraction

read_bed keeps BED's own 0-based half-open coordinates and skips track,
browser and comment lines. extract_intervals pulls the sequences from a
reference FASTA, re-centring to a fixed width when asked, so every output
is the same length - which is what makes position-wise comparison valid.

A centred window can run off the end of a contig, so on_edge makes the
policy explicit rather than silently producing short sequences: drop and
warn, raise, or pad with N. Padded bases do not resolve to a k-mer, so
they are masked as NaN and counted, not scored as zero.

pyfaidx is imported inside the function because it is an optional extra;
a module-level import would break every install that did not ask for it."
```

---

### Task 4: `FlexProfiler.from_bed`, packaging, and docs

**Files:**
- Modify: `DNAflexpy/core.py`
- Modify: `pyproject.toml`
- Modify: `CLAUDE.md`
- Test: `tests/test_bed_input.py`

**Interfaces:**
- Consumes: `read_bed` and `extract_intervals` from Task 3; `FlexProfiler._build(pairs, y=None)` from Task 2.
- Produces: `FlexProfiler.from_bed(bed, genome, width=None, on_edge="drop") -> FlexProfile | ProfileSet`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_bed_input.py`:

```python
from DNAflexpy import FlexProfiler
from DNAflexpy.core import ProfileSet
from DNAflexpy.results import FlexProfile


def test_from_bed_returns_a_profile(genome, tmp_path):
    bed = write_bed(tmp_path, "chr1\t10\t30\tpeakA\t0\t+\n")
    prof = FlexProfiler("DNaseI", window_size=10).from_bed(bed, genome)
    assert isinstance(prof, FlexProfile)
    assert prof.seqids == ["peakA"]


def test_from_bed_with_width_gives_equal_length_rows(genome, tmp_path):
    """Fixed width is the point: rows line up for position-wise work."""
    bed = write_bed(tmp_path, "chr1\t10\t12\ta\t0\t+\nchr1\t40\t60\tb\t0\t+\n")
    prof = FlexProfiler("DNaseI", window_size=10).from_bed(bed, genome, width=20)
    assert prof.frame.notna().all().all()
    assert prof.frame.shape == (2, 11)


def test_from_bed_has_no_y(genome, tmp_path):
    """BED input carries no label column."""
    bed = write_bed(tmp_path, "chr1\t10\t30\ta\t0\t+\n")
    assert FlexProfiler("DNaseI", window_size=10).from_bed(bed, genome).y is None


def test_from_bed_multi_feature(genome, tmp_path):
    bed = write_bed(tmp_path, "chr1\t10\t30\ta\t0\t+\n")
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).from_bed(bed, genome)
    assert isinstance(ps, ProfileSet)
    assert set(ps) == {"DNaseI", "trx"}


def test_padded_bases_are_masked_not_zeroed(genome, tmp_path):
    """The N padding must show up as masked, never as a real 0 measurement."""
    bed = write_bed(tmp_path, "chr1\t0\t2\ta\t0\t+\n")
    with pytest.warns(UserWarning, match="padded"):
        prof = FlexProfiler("DNaseI", window_size=0).from_bed(
            bed, genome, width=20, on_edge="pad"
        )
    assert prof.n_masked["a"] > 0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_bed_input.py -q -k from_bed`
Expected: FAIL with `AttributeError: 'FlexProfiler' object has no attribute 'from_bed'`.

- [ ] **Step 3: Add `from_bed`**

Add this method to `FlexProfiler` in `DNAflexpy/core.py`, directly after `profile_table`:

```python
    def from_bed(self, bed, genome, width=None, on_edge="drop"):
        """Profile intervals from a BED file against a reference genome.

        With `width`, intervals are re-centred and cut to that many bases,
        so every row is the same length. `on_edge` decides what happens to
        an interval whose centred window runs past a contig boundary.
        """
        from DNAflexpy.io import extract_intervals, read_bed

        pairs = extract_intervals(
            read_bed(bed), genome, width=width, on_edge=on_edge
        )
        return self._build(pairs)
```

- [ ] **Step 4: Add the optional extra**

In `pyproject.toml`, replace the `[project.optional-dependencies]` table with:

```toml
[project.optional-dependencies]
dev = ["pytest>=6.0"]
bed = ["pyfaidx>=0.7"]
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_bed_input.py -q`
Expected: PASS, 19 passed.

- [ ] **Step 6: Run the whole suite and reinstall**

Run:
```bash
python -m pip install -e . -q
python -m pytest -q
```
Expected: 356 passed, with the 230 byte-equality cases still green.

- [ ] **Step 7: Update `CLAUDE.md`**

In the architecture section, the `DNAflexpy/core.py` bullet lists the profiler's entry points and the `DNAflexpy/io.py` bullet says only `read_fasta`. Update both so they describe what now exists:

- `core.py`: add `.profile_table(path, seq_col=..., value_col=...)` for a labelled table, which fills `FlexProfile.y`, and `.from_bed(bed, genome, width=..., on_edge=...)` for BED intervals against a reference genome.
- `io.py`: change `read_fasta` to `read_fasta`, `read_table`, `read_bed`, `extract_intervals` — and note that `extract_intervals` imports `pyfaidx` lazily because it is an optional extra (`DNAflexpy[bed]`).

Keep the file's existing tone and length. Do not restate anything already covered elsewhere in it.

- [ ] **Step 8: Commit**

```bash
git add DNAflexpy/core.py pyproject.toml CLAUDE.md tests/test_bed_input.py
git commit -m "Add FlexProfiler.from_bed and the bed optional extra

Profiling a peak set no longer needs a hand-prepared FASTA: from_bed
reads intervals and pulls their sequences straight from a reference
genome, optionally re-centred to a fixed width so the rows line up.

pyfaidx moves into a 'bed' optional extra rather than a hard dependency,
since it is only needed for this one entry point."
```

---

## Definition of done

- [ ] `python -m pytest -q` passes with no failures.
- [ ] The 230 byte-equality cases still pass — this phase must not change any existing numeric output.
- [ ] `FlexProfiler` has all five entry points: `profile`, `profile_seqs`, `profile_fasta`, `profile_table`, `from_bed`.
- [ ] `profile_table` fills `.y`; `profile_fasta`, `profile_seqs` and `from_bed` leave it `None`.
- [ ] `python -c "import DNAflexpy"` works in an environment without `pyfaidx` installed — proving the lazy import holds. Check with: `python -c "import sys; sys.modules['pyfaidx']=None; import DNAflexpy; print('ok')"`.
- [ ] No commit contains a `Co-Authored-By` trailer or AI attribution.
- [ ] `git diff --stat -- rxv/` is empty across the whole phase.

## Deferred to later phases

Phase 4 onward: `encode.py` (ML feature matrices), plotting, dinucleotide shuffle and per-position statistics, streaming input, the provenance sidecar, `ProfileSet.correlate()`, bigWig/BED export, the `DNAflexpy` CLI, and the container. The spec also lists `heatmap(order_rows="y")` and `metaprofile(groupby=("y", threshold))` as things `.y` unlocks — both belong to the plotting phase, not here.
