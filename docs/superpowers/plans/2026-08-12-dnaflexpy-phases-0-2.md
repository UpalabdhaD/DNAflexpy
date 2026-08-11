# DNAflexpy Rewrite — Phases 0–2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Archive the legacy DNAflexpy package to `rxv/DNAflexpy/`, then build a new class-based `DNAflexpy` package whose profiling output is byte-identical to the archive.

**Architecture:** The legacy package is moved with `git mv` and stays importable, so the new implementation can be verified by differential byte comparison rather than against frozen files. The new package separates a pure numeric kernel from a `FeatureTable` (loaded once, fixing an 8 ms-per-sequence YAML reload), a `FlexProfiler` engine, and a `FlexProfile` result object.

**Tech Stack:** Python 3.12, pandas 2.3, numpy 2.2, PyYAML, pytest 9.

## Global Constraints

- **Spec:** `docs/superpowers/specs/2026-08-12-dnaflex-rewrite-design.md`. This plan covers Phases 0–2, **plus four small Phase 3 items** that the class API is incomplete without: the bare-string `profile()` entry point, `profile_seqs()`, `ProfileSet` for multi-feature runs, user-supplied lookup tables, and the memoised module-level one-liner. They are a few lines each and share the same tests as the core; deferring them would mean shipping a `FlexProfiler` that can only read files. The rest of Phase 3 (`profile_table`, `from_bed`) and Phases 4–9 are out of scope.
- **Branch:** all work on `dev`. Never commit to `main`.
- **Commit messages:** no `Co-Authored-By` trailers, no AI attribution of any kind.
- **Never modify the logic of `rxv/DNAflexpy/`.** Its only sanctioned edit is the import-path fix in Task 1.
- **Byte equality is the gate.** Numeric output must match the archive exactly. Use Python's builtin `sum()` on lists, never `numpy.sum`/`numpy.mean`/`math.fsum` — pairwise or compensated summation produces different last bits, which flip values at `round(v, 3)` boundaries.
- **Do not modify** `Notebooks/*.ipynb`.
- Python floats are serialised by `pandas.to_csv`, which uses `repr`. `round(0.0100, 3)` writes `0.01`, never `0.010`.
- macOS uses the `spawn` start method: any test that starts a `multiprocessing.Pool` must live in a real file (pytest satisfies this) and all pool targets must be module-level functions, never closures or lambdas.

---

### Task 1: Archive the legacy package to `rxv/DNAflexpy/`

**Files:**
- Move: `DNAflexpy/` → `rxv/DNAflexpy/` (via `git mv`)
- Create: `rxv/__init__.py`
- Modify: `rxv/DNAflexpy/utils.py:175`
- Delete: `tests/test_utils.py`, `tests/test_core.py`
- Test: `tests/test_archive.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `rxv.DNAflexpy.core.DNAflexpyMP(input_file, window_size, feature, threads, outfile=None)`, `rxv.DNAflexpy.utils.read_fasta(filepath) -> Generator[tuple[str, str]]`, `rxv.DNAflexpy.utils.load_feature_data() -> dict`, `rxv.DNAflexpy.utils.get_kmer_len(feature) -> int | None`, `rxv.DNAflexpy.utils.seq_to_numeric_profile(seqid, sequence, kmer_len, window_size, feature, feature_lookup) -> list`. Task 8 depends on all of these.

`tests/test_utils.py` imports `calculate_window_averages` from a module named `utils`, neither of which exists; it currently breaks pytest collection for the whole repo. `tests/test_core.py` is empty. Both are deleted rather than ported, per the spec.

- [ ] **Step 1: Move the package and create the namespace**

```bash
mkdir -p rxv
git mv DNAflexpy rxv/DNAflexpy
cat > rxv/__init__.py <<'EOF'
"""Namespace holding the archived, frozen DNAflexpy implementation.

`rxv.DNAflexpy` is the original package, preserved verbatim except for the
import path in `utils.load_feature_data`. It is kept importable so the
rewritten `DNAflexpy` package can be verified against it by differential
byte comparison. Do not change its logic.
"""
EOF
git rm -q tests/test_utils.py tests/test_core.py
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_archive.py`:

```python
"""The archived package must keep producing exactly what it always produced."""
import io
import contextlib
import pathlib

import pandas as pd
import pytest

from rxv.DNAflexpy.utils import (
    read_fasta,
    load_feature_data,
    get_kmer_len,
    seq_to_numeric_profile,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent

FASTAS = {
    "test_fasta": ROOT / "rxv/DNAflexpy/data/test_fasta.fa",
    "edge": ROOT / "tests/test_edge_cases.fasta",
    "exact": ROOT / "tests/test_exact_kmer.fasta",
}

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


def archive_tsv_bytes(fasta_key, feature, window_size):
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


@pytest.mark.parametrize("relpath,fasta_key,feature,window", CHECKED_IN)
def test_archive_reproduces_checked_in_output(relpath, fasta_key, feature, window):
    expected = (ROOT / relpath).read_bytes()
    assert archive_tsv_bytes(fasta_key, feature, window) == expected


def test_archive_reads_its_own_lookup_table():
    """Guards the import-path fix: the archive must not read the new package's YAML."""
    import rxv.DNAflexpy.utils as archived

    assert "rxv.DNAflexpy.data" in pathlib.Path(archived.__file__).read_text()
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `python -m pytest tests/test_archive.py -q`
Expected: all tests FAIL with `FileNotFoundError` or a `RuntimeError` wrapping it, because `load_feature_data` still calls `files("DNAflexpy.data")` and the package no longer lives there.

- [ ] **Step 4: Fix the one sanctioned line in the archive**

In `rxv/DNAflexpy/utils.py`, line 175, change only the package string:

```python
        yamlfilepath = files("rxv.DNAflexpy.data").joinpath("lookupNEW.yaml")
```

Leave every other line of the archive untouched. Without this the archive would resolve `DNAflexpy.data` to the *new* package's table once Task 3 lands, silently reading the wrong data.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_archive.py -q`
Expected: PASS, 17 passed.

- [ ] **Step 6: Commit**

```bash
git add -A rxv tests
git commit -m "Archive legacy DNAflexpy to rxv/DNAflexpy

Moved with git mv so history follows. The only logic-adjacent edit is
utils.py:175, where files(\"DNAflexpy.data\") becomes
files(\"rxv.DNAflexpy.data\"); without it the archive would read the new
package's lookup table once that package exists.

Adds tests asserting all 16 checked-in TSVs still reproduce byte-exactly.
This validates the comparison harness itself before anything is verified
against it.

Deletes tests/test_utils.py, which imported a function that no longer
exists and broke collection for the whole repo, and tests/test_core.py,
which was empty."
```

---

### Task 2: Repackage for two side-by-side packages

**Files:**
- Modify: `pyproject.toml`
- Delete: `setup.py`
- Modify: `MANIFEST.in`
- Modify: `scripts/plot_profiles.py:150-160`
- Test: `tests/test_packaging.py`

**Interfaces:**
- Consumes: `rxv.DNAflexpy` from Task 1.
- Produces: an editable install exposing both `rxv.DNAflexpy` and (after Task 3) `DNAflexpy`. No console script yet — `DNAflexpy.cli` does not exist until Phase 8, so declaring one now would break the install.

`setup.py` is deleted because leaving it beside a `[project]` table makes setuptools error on conflicting metadata.

- [ ] **Step 1: Write the failing test**

Create `tests/test_packaging.py`:

```python
"""The archive and the new package must coexist in one install."""
import pathlib
import tomllib

ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_setup_py_is_gone():
    """setup.py beside a [project] table makes setuptools error."""
    assert not (ROOT / "setup.py").exists()


def test_pyproject_declares_packages_explicitly():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    packages = cfg["tool"]["setuptools"]["packages"]
    # Explicit list, never find_packages(): find_packages() would sweep the
    # archive into the distribution while shipping none of its YAML.
    assert "rxv" in packages
    assert "rxv.DNAflexpy" in packages
    assert "DNAflexpy" in packages


def test_pyproject_ships_both_lookup_tables():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    data = cfg["tool"]["setuptools"]["package-data"]
    assert data["rxv.DNAflexpy"] == ["data/*.yaml"]
    assert data["DNAflexpy"] == ["data/*.yaml"]


def test_manifest_points_at_the_archive():
    assert "rxv/DNAflexpy/data/lookupNEW.yaml" in (ROOT / "MANIFEST.in").read_text()


def test_plot_profiles_invokes_the_archived_cli():
    """It shells out to a module path; unfixed it would hit the new CLI."""
    src = (ROOT / "scripts/plot_profiles.py").read_text()
    assert "rxv.DNAflexpy.cli" in src
    assert '"DNAflexpy.cli"' not in src
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_packaging.py -q`
Expected: FAIL — `setup.py` still exists and `pyproject.toml` has no `[tool.setuptools]` table.

- [ ] **Step 3: Rewrite `pyproject.toml`**

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "DNAflexpy"
dynamic = ["version"]
description = "DNA flexibility profiling from sequence using di- and trinucleotide feature tables"
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "pandas>=2.2.3",
    "pyyaml>=6.0.2",
    "numpy>=1.24",
]

[project.optional-dependencies]
dev = ["pytest>=6.0"]

[tool.setuptools]
# Explicit, never find_packages(): that would sweep the archive into the
# distribution while shipping none of its package data.
packages = ["DNAflexpy", "rxv", "rxv.DNAflexpy"]

[tool.setuptools.package-data]
DNAflexpy = ["data/*.yaml"]
"rxv.DNAflexpy" = ["data/*.yaml"]

[tool.setuptools.dynamic]
version = { attr = "DNAflexpy.__version__" }

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 4: Delete `setup.py` and update `MANIFEST.in`**

```bash
git rm -q setup.py
cat > MANIFEST.in <<'EOF'
include DNAflexpy/data/*.yaml
include rxv/DNAflexpy/data/*.yaml
EOF
```

- [ ] **Step 5: Point `scripts/plot_profiles.py` at the archived CLI**

In `scripts/plot_profiles.py`, inside `main()`, the subprocess command list currently reads `"-m", "DNAflexpy.cli",`. Change that one element:

```python
            cmd = [
                sys.executable,
                "-m",
                "rxv.DNAflexpy.cli",
                args.generate_random_fasta,
```

The script passes the archive's flags (`--window-size`, `--feature`, `--outfile`), so it must keep invoking the archive; the new CLI does not exist yet and will not accept the same flags when it does.

- [ ] **Step 6: Create a minimal `DNAflexpy/__init__.py` so the package resolves**

The `[tool.setuptools.dynamic]` version attribute must import. Task 3 fills this file in properly.

```bash
mkdir -p DNAflexpy/data
cat > DNAflexpy/__init__.py <<'EOF'
"""DNAflexpy — DNA flexibility profiling from sequence."""

__version__ = "0.3.0.dev0"
__all__ = ["__version__"]
EOF
cp rxv/DNAflexpy/data/lookupNEW.yaml DNAflexpy/data/lookupNEW.yaml
```

The new package gets its **own copy** of the YAML so the two share no state.

- [ ] **Step 7: Reinstall and run the tests**

Run:
```bash
python -m pip install -e . -q
python -m pytest tests/test_packaging.py tests/test_archive.py -q
```
Expected: PASS, 22 passed.

- [ ] **Step 8: Commit**

```bash
git add -A pyproject.toml MANIFEST.in scripts tests DNAflexpy
git commit -m "Repackage for the archive and the new package side by side

Migrates metadata from setup.py into pyproject.toml with an explicit
package list. find_packages() would have swept rxv into the distribution
while shipping none of its YAML. setup.py is deleted because setuptools
errors when it sits beside a [project] table.

scripts/plot_profiles.py shells out to a module path and now targets
rxv.DNAflexpy.cli; left alone it would invoke the new CLI with the
archive's flags.

Adds a stub DNAflexpy/__init__.py carrying __version__, which the
dynamic version attribute needs to import, plus the new package's own
copy of the lookup table."
```

---

### Task 3: `FeatureTable` — load once, validate, infer k

**Files:**
- Create: `DNAflexpy/lookup.py`
- Modify: `DNAflexpy/__init__.py`
- Test: `tests/test_lookup.py`

**Interfaces:**
- Consumes: `DNAflexpy.__version__` from Task 2.
- Produces:
  - `FeatureTable.from_yaml(path: str | Path | None = None) -> FeatureTable`
  - `FeatureTable.from_dict(mapping: dict[str, dict[str, float]]) -> FeatureTable`
  - `FeatureTable.features -> list[str]`
  - `FeatureTable.kmer_len(feature: str) -> int`
  - `FeatureTable.table(feature: str) -> dict[str, float]`
  - `DNAflexpy.lookup.default_table() -> FeatureTable` (memoised)
  - Raises `ValueError` for unknown features and malformed tables.

Tasks 4, 6, 7, 9 all consume this.

- [ ] **Step 1: Write the failing test**

Create `tests/test_lookup.py`:

```python
import warnings

import pytest

from DNAflexpy.lookup import FeatureTable, default_table


def test_infers_kmer_length_from_keys():
    t = default_table()
    assert t.kmer_len("DNaseI") == 3
    assert t.kmer_len("trx") == 2


def test_inference_matches_the_archives_hardcoded_map():
    """The archive hardcodes feature -> k. Inference must agree on all of them."""
    from rxv.DNAflexpy.utils import get_kmer_len

    t = default_table()
    for feature in t.features:
        hardcoded = get_kmer_len(feature)
        if hardcoded is not None:
            assert t.kmer_len(feature) == hardcoded, feature


def test_unlocks_features_the_archive_could_not_reach():
    """freeen, gc and mechen have no entry in the archive's map."""
    t = default_table()
    for feature in ("freeen", "gc", "mechen"):
        assert feature in t.features
        assert t.kmer_len(feature) == 2


def test_packaged_tables_are_complete():
    t = default_table()
    for feature in t.features:
        assert len(t.table(feature)) == 4 ** t.kmer_len(feature), feature


def test_unknown_feature_raises_with_available_names():
    t = default_table()
    with pytest.raises(ValueError, match="DNaseI"):
        t.kmer_len("DNasel")  # lowercase L typo


def test_keys_are_uppercased():
    t = FeatureTable.from_dict({"f": {"aa": 1.0, "at": 2.0, "ta": 3.0, "tt": 4.0}})
    assert t.table("f")["AA"] == 1.0


def test_mixed_key_lengths_rejected():
    """This is what makes k inference safe."""
    with pytest.raises(ValueError, match="mixed k-mer lengths"):
        FeatureTable.from_dict({"f": {"AA": 1.0, "AAT": 2.0}})


def test_non_acgt_keys_rejected():
    with pytest.raises(ValueError, match="non-ACGT"):
        FeatureTable.from_dict({"f": {"AN": 1.0}})


def test_non_numeric_values_rejected():
    with pytest.raises(ValueError, match="non-numeric"):
        FeatureTable.from_dict({"f": {"AA": "high"}})


def test_empty_table_rejected():
    with pytest.raises(ValueError, match="empty"):
        FeatureTable.from_dict({"f": {}})


def test_incomplete_table_warns_but_loads():
    with pytest.warns(UserWarning, match="incomplete"):
        t = FeatureTable.from_dict({"f": {"AA": 1.0, "AT": 2.0}})
    assert t.kmer_len("f") == 2


def test_default_table_is_memoised():
    """The 8 ms YAML parse must happen once, not per call."""
    assert default_table() is default_table()
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_lookup.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'DNAflexpy.lookup'`.

- [ ] **Step 3: Implement `DNAflexpy/lookup.py`**

```python
"""Feature lookup tables: load once, validate, infer k-mer length.

The archive hardcoded a feature -> k-mer-length map that had to be kept in
sync with the YAML by hand, and silently produced None rows for any feature
missing from it. Here k is inferred from the table's own keys, which is
safe because validation rejects any table with mixed key lengths.
"""
from __future__ import annotations

import functools
import warnings
from importlib.resources import files
from pathlib import Path

import yaml

_ALPHABET = frozenset("ACGT")


class FeatureTable:
    """A validated collection of k-mer -> value tables, keyed by feature."""

    def __init__(self, tables: dict[str, dict[str, float]], kmer_lens: dict[str, int]):
        self._tables = tables
        self._kmer_lens = kmer_lens

    @property
    def features(self) -> list[str]:
        return list(self._tables)

    def kmer_len(self, feature: str) -> int:
        self._require(feature)
        return self._kmer_lens[feature]

    def table(self, feature: str) -> dict[str, float]:
        self._require(feature)
        return self._tables[feature]

    def _require(self, feature: str) -> None:
        if feature not in self._tables:
            raise ValueError(
                f"unknown feature {feature!r}. Available: {', '.join(self.features)}"
            )

    @classmethod
    def from_dict(cls, mapping: dict) -> "FeatureTable":
        tables: dict[str, dict[str, float]] = {}
        kmer_lens: dict[str, int] = {}
        for feature, raw in mapping.items():
            tables[feature], kmer_lens[feature] = _validate(feature, raw)
        return cls(tables, kmer_lens)

    @classmethod
    def from_yaml(cls, path=None) -> "FeatureTable":
        if path is None:
            source = files("DNAflexpy.data").joinpath("lookupNEW.yaml")
            text = source.read_text()
        else:
            text = Path(path).read_text()
        try:
            loaded = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise ValueError(f"could not parse feature table: {exc}") from exc
        if not isinstance(loaded, dict):
            raise ValueError("feature table must map feature names to k-mer tables")
        return cls.from_dict(loaded)


def _validate(feature: str, raw) -> tuple[dict[str, float], int]:
    """Return an uppercased table and its inferred k-mer length."""
    if not isinstance(raw, dict) or not raw:
        raise ValueError(f"feature {feature!r} has an empty or malformed table")

    table: dict[str, float] = {}
    for key, value in raw.items():
        kmer = str(key).upper()
        if set(kmer) - _ALPHABET:
            raise ValueError(
                f"feature {feature!r} has non-ACGT k-mer {kmer!r}; "
                "tables must contain unambiguous k-mers only"
            )
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(
                f"feature {feature!r} has non-numeric value {value!r} for {kmer!r}"
            )
        table[kmer] = float(value)

    lengths = {len(k) for k in table}
    if len(lengths) > 1:
        raise ValueError(
            f"feature {feature!r} has mixed k-mer lengths {sorted(lengths)}; "
            "k-mer length is inferred from the keys and must be uniform"
        )
    kmer_len = lengths.pop()

    expected = 4 ** kmer_len
    if len(table) < expected:
        warnings.warn(
            f"feature {feature!r} table is incomplete: {len(table)} of {expected} "
            f"{kmer_len}-mers. Missing k-mers resolve to NaN and are masked.",
            UserWarning,
            stacklevel=3,
        )
    return table, kmer_len


@functools.lru_cache(maxsize=1)
def default_table() -> FeatureTable:
    """The packaged lookup table, parsed at most once per process."""
    return FeatureTable.from_yaml()
```

- [ ] **Step 4: Export it**

Replace `DNAflexpy/__init__.py`:

```python
"""DNAflexpy — DNA flexibility profiling from sequence."""

from DNAflexpy.lookup import FeatureTable, default_table

__version__ = "0.3.0.dev0"
__all__ = ["FeatureTable", "default_table", "__version__"]
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_lookup.py -q`
Expected: PASS, 12 passed.

- [ ] **Step 6: Commit**

```bash
git add DNAflexpy tests/test_lookup.py
git commit -m "Add FeatureTable with validation and inferred k-mer length

Replaces the archive's hardcoded feature -> k map, which had to be
hand-synced with the YAML and silently yielded None rows for any feature
missing from it. k is now inferred from the table's own keys, which
validation makes safe by rejecting mixed key lengths.

A test asserts inference agrees with the archive's map on every feature
the archive mapped, and that freeen, gc and mechen - which it did not
map - now resolve to k=2.

The packaged table is memoised, so the 8ms parse happens once."
```

---

### Task 4: Numeric kernel with NaN masking

**Files:**
- Create: `DNAflexpy/kernel.py`
- Test: `tests/test_kernel.py`

**Interfaces:**
- Consumes: nothing (pure functions over primitives).
- Produces:
  - `kmer_values(sequence: str, kmer_len: int, table: dict[str, float]) -> list[float]`
  - `window_means(values: list[float], window_size: int, kmer_len: int, n_positions: int) -> list[float]`
  - `profile_values(sequence: str, kmer_len: int, table: dict[str, float], window_size: int) -> list[float]`
  - `MASKED` (module constant, `float("nan")`)

Tasks 6 and 7 consume `profile_values`.

**Why this preserves byte equality.** The archive recomputes k-mers inside every window slice. Slicing a precomputed value array yields the *same floats in the same order*, because the k-mers of `sequence[s : s+w]` are exactly `values[s : s + w - k + 1]`. Summing identical floats left-to-right gives an identical result, so the optimisation is safe — and turns O(n·w) into O(n).

- [ ] **Step 1: Write the failing test**

Create `tests/test_kernel.py`:

```python
import math

from DNAflexpy.kernel import kmer_values, profile_values, window_means
from DNAflexpy.lookup import default_table

SEQ = "ATGCGTACGTAGCTAGCGTAGCTAGT"  # 26 nt, the canonical test sequence


def dnase():
    return default_table().table("DNaseI")


def test_kmer_values_match_the_archive():
    from rxv.DNAflexpy.utils import load_feature_data, transform_seq_to_feat

    expected = transform_seq_to_feat(SEQ, 3, "DNaseI", load_feature_data())
    assert kmer_values(SEQ, 3, dnase()) == expected


def test_window_means_match_the_archive():
    """window_size == kmer_len must reproduce the per-kmer values exactly."""
    values = kmer_values(SEQ, 3, dnase())
    assert window_means(values, 3, 3, len(SEQ)) == [round(v, 3) for v in values]


def test_window_count_follows_the_formula():
    values = kmer_values(SEQ, 3, dnase())
    assert len(window_means(values, 10, 3, len(SEQ))) == len(SEQ) - 10 + 1


def test_window_smaller_than_kmer_yields_zeros():
    """No k-mer fits, so the archive's `sum(w)/len(w) if w else 0.0` gives 0.0."""
    values = kmer_values(SEQ, 3, dnase())
    out = window_means(values, 2, 3, len(SEQ))
    assert out == [0.0] * (len(SEQ) - 2 + 1)


def test_window_larger_than_sequence_yields_nothing():
    values = kmer_values("ATGC", 3, dnase())
    assert window_means(values, 30, 3, 4) == []


def test_sequence_is_uppercased():
    assert kmer_values("atgcgt", 3, dnase()) == kmer_values("ATGCGT", 3, dnase())


def test_ambiguous_kmers_are_masked_not_zeroed():
    """The archive scores these 0, which is a real value in these tables."""
    values = kmer_values("ATGNGTACGT", 3, dnase())
    assert not math.isnan(values[0])
    assert all(math.isnan(v) for v in values[1:4])
    assert not math.isnan(values[4])


def test_iupac_codes_are_masked_too():
    values = kmer_values("ATGRGTACGT", 3, dnase())
    assert all(math.isnan(v) for v in values[1:4])


def test_partially_masked_window_averages_the_resolved_kmers():
    values = [1.0, 2.0, float("nan"), 4.0]
    # window of 3 kmers starting at 0 -> mean(1.0, 2.0) == 1.5
    assert window_means(values, 5, 3, 7)[0] == 1.5


def test_fully_masked_window_is_nan():
    nan = float("nan")
    assert math.isnan(window_means([nan, nan], 4, 3, 6)[0])


def test_profile_values_window_zero_returns_per_kmer():
    assert profile_values(SEQ, 3, dnase(), 0) == kmer_values(SEQ, 3, dnase())
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_kernel.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'DNAflexpy.kernel'`.

- [ ] **Step 3: Implement `DNAflexpy/kernel.py`**

```python
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_kernel.py -q`
Expected: PASS, 11 passed.

- [ ] **Step 5: Commit**

```bash
git add DNAflexpy/kernel.py tests/test_kernel.py
git commit -m "Add numeric kernel with NaN masking for unresolvable k-mers

Pure functions, no I/O. Arithmetic deliberately mirrors the archive:
builtin sum over Python lists, left to right. A numpy reduction would
change the last bits and flip values at round(v, 3) boundaries.

Windows are computed by slicing a precomputed value array rather than
recomputing k-mers per window as the archive does. The floats and their
order are identical, so results are bit-identical, and the cost drops
from O(n*w) to O(n).

Unresolvable k-mers resolve to NaN rather than 0. A window averages
whichever k-mers resolved, and is NaN only if all of them were masked."
```

---

### Task 5: `FlexProfile` and the `to_tsv` byte contract

**Files:**
- Create: `DNAflexpy/profile.py`
- Test: `tests/test_profile.py`

**Interfaces:**
- Consumes: nothing from earlier tasks at runtime.
- Produces:
  - `FlexProfile(rows: list[list], feature: str, window_size: int, kmer_len: int, y=None)` where each row is `[seqid, *values]`
  - `.frame -> pandas.DataFrame` (rows indexed by sequence ID, columns are positions)
  - `.seqids -> list[str]`
  - `.n_masked -> dict[str, int]`
  - `.to_tsv(path) -> None`
  - `.to_frame(tidy: bool = False) -> pandas.DataFrame`

Task 7 constructs these; Task 8 compares their `to_tsv` output.

**The byte contract**, verified against the checked-in files:
- No header, no index, tab-separated, sequence ID first.
- Ragged rows NaN-pad to the widest row; NaN renders as an empty field, so a short sequence emits its ID then trailing tabs.
- When every row is ID-only the frame has one column and there is no padding: a bare ID, no trailing tab.
- `round(v, 3)` then float repr, so `0.01` never `0.010`.
- Masked positions also serialise as empty. Nothing distinguishes them from padding *in the file* — `n_masked` carries that in memory. There is no single `na_rep` that could separate the two, and choosing `"NA"` would relabel padding and break byte equality.

- [ ] **Step 1: Write the failing test**

Create `tests/test_profile.py`:

```python
import math
import pathlib

from DNAflexpy.profile import FlexProfile


def make(rows):
    return FlexProfile(rows, feature="DNaseI", window_size=10, kmer_len=3)


def test_ragged_rows_pad_with_trailing_tabs(tmp_path):
    """Matches w10_mixed.tsv: a too-short sequence is ID plus trailing tabs."""
    out = tmp_path / "o.tsv"
    make([["long", 0.1, 0.2, 0.3], ["short"]]).to_tsv(out)
    lines = out.read_text().splitlines()
    assert lines[0] == "long\t0.1\t0.2\t0.3"
    assert lines[1] == "short\t\t\t"


def test_all_id_only_rows_have_no_trailing_tab(tmp_path):
    """Matches w30.tsv: one column means no padding at all."""
    out = tmp_path / "o.tsv"
    make([["a"], ["b"]]).to_tsv(out)
    assert out.read_text().splitlines() == ["a", "b"]


def test_trailing_zeros_are_dropped(tmp_path):
    out = tmp_path / "o.tsv"
    make([["s", round(0.0100, 3)]]).to_tsv(out)
    assert out.read_text().splitlines()[0] == "s\t0.01"


def test_masked_values_serialise_as_empty(tmp_path):
    out = tmp_path / "o.tsv"
    make([["s", 0.1, float("nan"), 0.3]]).to_tsv(out)
    assert out.read_text().splitlines()[0] == "s\t0.1\t\t0.3"


def test_n_masked_counts_masked_positions():
    p = make([["a", 0.1, float("nan")], ["b", float("nan"), float("nan")]])
    assert p.n_masked == {"a": 1, "b": 2}


def test_n_masked_ignores_ragged_padding():
    """Padding is NaN too, but it is not a masked measurement."""
    p = make([["long", 0.1, 0.2, 0.3], ["short"]])
    assert p.n_masked == {"long": 0, "short": 0}


def test_frame_is_indexed_by_seqid():
    p = make([["a", 0.1, 0.2], ["b", 0.3, 0.4]])
    assert list(p.frame.index) == ["a", "b"]
    assert p.frame.shape == (2, 2)


def test_seqids_preserve_input_order():
    p = make([["z", 0.1], ["a", 0.2]])
    assert p.seqids == ["z", "a"]


def test_tidy_frame_is_long():
    p = make([["a", 0.1, 0.2]])
    tidy = p.to_frame(tidy=True)
    assert list(tidy.columns) == ["seqid", "position", "value", "feature"]
    assert len(tidy) == 2
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_profile.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'DNAflexpy.profile'`.

- [ ] **Step 3: Implement `DNAflexpy/profile.py`**

```python
"""Profiling results and their serialisation."""
from __future__ import annotations

import math

import pandas as pd


class FlexProfile:
    """Profile values for a set of sequences, one row per sequence.

    Rows arrive as `[seqid, *values]` and may be ragged: a sequence shorter
    than the window contributes only its ID.
    """

    def __init__(self, rows, feature: str, window_size: int, kmer_len: int, y=None):
        self._rows = [list(r) for r in rows]
        self.feature = feature
        self.window_size = window_size
        self.kmer_len = kmer_len
        self.y = y

    @property
    def seqids(self) -> list[str]:
        return [r[0] for r in self._rows]

    @property
    def n_masked(self) -> dict[str, int]:
        """Masked measurements per sequence, excluding ragged padding.

        Padding is NaN in the frame too, which is why this counts from the
        original rows rather than from the padded frame.
        """
        return {
            row[0]: sum(1 for v in row[1:] if isinstance(v, float) and math.isnan(v))
            for row in self._rows
        }

    @property
    def frame(self) -> pd.DataFrame:
        """Values as a DataFrame indexed by sequence ID, NaN-padded."""
        widest = max((len(r) - 1 for r in self._rows), default=0)
        data = [r[1:] + [float("nan")] * (widest - (len(r) - 1)) for r in self._rows]
        return pd.DataFrame(data, index=self.seqids, columns=range(1, widest + 1))

    def to_frame(self, tidy: bool = False) -> pd.DataFrame:
        if not tidy:
            return self.frame
        long = (
            self.frame.stack(future_stack=True)
            .rename("value")
            .reset_index()
            .rename(columns={"level_0": "seqid", "level_1": "position"})
        )
        long["feature"] = self.feature
        return long[["seqid", "position", "value", "feature"]]

    def to_tsv(self, path) -> None:
        """Write the archive's exact format.

        Built from ragged lists so pandas pads with NaN and renders each as
        an empty field, reproducing the archive's trailing tabs. Masked
        values are also NaN and therefore also empty: no `na_rep` could
        separate the two without relabelling padding and breaking byte
        equality.
        """
        pd.DataFrame(self._rows).to_csv(path, index=False, header=False, sep="\t")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_profile.py -q`
Expected: PASS, 9 passed.

If `future_stack=True` raises on the installed pandas, use `.stack()` — the argument exists to silence a pandas 2.x FutureWarning and is not load-bearing.

- [ ] **Step 5: Commit**

```bash
git add DNAflexpy/profile.py tests/test_profile.py
git commit -m "Add FlexProfile with the archive's exact TSV byte contract

Output is built from ragged lists so pandas NaN-pads to the widest row
and renders each pad as an empty field, reproducing the archive's
trailing tabs; when every row is ID-only there is one column and no
padding at all.

Masked values serialise as empty too. No na_rep can distinguish them
from padding, since both are NaN in one frame, and choosing NA would
relabel padding and break byte equality. n_masked carries the counts in
memory instead, computed from the original rows so ragged padding is not
miscounted as masking."
```

---

### Task 6: `FlexProfiler` — construction and in-memory profiling

**Files:**
- Create: `DNAflexpy/core.py`
- Modify: `DNAflexpy/__init__.py`
- Test: `tests/test_profiler.py`

**Interfaces:**
- Consumes: `FeatureTable`/`default_table` (Task 3), `profile_values` (Task 4), `FlexProfile` (Task 5).
- Produces:
  - `FlexProfiler(feature: str | list[str], window_size: int = 10, lookup=None)`
  - `.features -> list[str]`, `.kmer_len(feature=None) -> int`
  - `.profile(seq: str) -> numpy.ndarray` (single feature only)
  - `.profile_seqs(seqs: list[str] | dict[str, str]) -> FlexProfile | ProfileSet`
  - `ProfileSet`, a dict-like `{feature: FlexProfile}` with `.to_tidy()`

Task 7 adds `profile_fasta` to this class; Task 9 wraps `profile`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_profiler.py`:

```python
import numpy as np
import pytest

from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.profile import FlexProfile

SEQ = "ATGCGTACGTAGCTAGCGTAGCTAGT"


def test_bare_string_returns_an_array():
    out = FlexProfiler("DNaseI", window_size=10).profile(SEQ)
    assert isinstance(out, np.ndarray)
    assert len(out) == len(SEQ) - 10 + 1


def test_no_file_or_seqid_needed():
    """The archive forced a FASTA round-trip to inspect one sequence."""
    assert len(FlexProfiler("trx", window_size=0).profile("ATGC")) == 3


def test_unknown_feature_fails_at_construction():
    """A typo must not survive until after a full FASTA pass."""
    with pytest.raises(ValueError, match="unknown feature"):
        FlexProfiler("DNasel")


def test_features_the_archive_could_not_reach_now_work():
    assert len(FlexProfiler("gc", window_size=0).profile("ATGC")) == 3


def test_profile_rejects_multi_feature():
    p = FlexProfiler(["DNaseI", "trx"])
    with pytest.raises(ValueError, match="single feature"):
        p.profile(SEQ)


def test_profile_seqs_from_list_generates_ids():
    prof = FlexProfiler("DNaseI", window_size=10).profile_seqs([SEQ, SEQ])
    assert isinstance(prof, FlexProfile)
    assert prof.seqids == ["seq_0", "seq_1"]


def test_profile_seqs_from_dict_keeps_ids():
    prof = FlexProfiler("DNaseI", window_size=10).profile_seqs({"a": SEQ, "b": SEQ})
    assert prof.seqids == ["a", "b"]


def test_multi_feature_returns_a_profile_set():
    ps = FlexProfiler(["DNaseI", "trx"], window_size=10).profile_seqs({"a": SEQ})
    assert isinstance(ps, ProfileSet)
    assert set(ps) == {"DNaseI", "trx"}
    # Different k means different vector lengths for the same sequence.
    assert ps["DNaseI"].kmer_len == 3
    assert ps["trx"].kmer_len == 2


def test_custom_lookup_from_dict():
    p = FlexProfiler("mine", window_size=0, lookup={"mine": {"AA": 1.0, "AT": 2.0,
                                                            "TA": 3.0, "TT": 4.0}})
    assert list(p.profile("ATA")) == [2.0, 3.0]


def test_short_sequence_yields_no_values():
    assert len(FlexProfiler("DNaseI", window_size=30).profile("ATGC")) == 0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_profiler.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'DNAflexpy.core'`.

- [ ] **Step 3: Implement `DNAflexpy/core.py`**

```python
"""The profiling engine."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from DNAflexpy.kernel import profile_values
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.profile import FlexProfile


class ProfileSet(dict):
    """`{feature: FlexProfile}`.

    Multi-feature runs cannot share one matrix: features have different k,
    so their vectors differ in length for the same sequence.
    """

    def to_tidy(self) -> pd.DataFrame:
        return pd.concat([p.to_frame(tidy=True) for p in self.values()],
                         ignore_index=True)


class FlexProfiler:
    """Profiles sequences against one or more flexibility features.

    The lookup table is loaded once, at construction. The archive reloaded
    it on every per-sequence call, which cost 8 ms against 0.1 ms of actual
    computation.
    """

    def __init__(self, feature, window_size: int = 10, lookup=None):
        self._table = _resolve_lookup(lookup)
        self._features = [feature] if isinstance(feature, str) else list(feature)
        if not self._features:
            raise ValueError("at least one feature is required")
        # Validate now so a typo fails immediately, not after a full pass.
        for name in self._features:
            self._table.kmer_len(name)
        if window_size < 0:
            raise ValueError(f"window_size must be >= 0, got {window_size}")
        self.window_size = window_size

    @property
    def features(self) -> list[str]:
        return list(self._features)

    def kmer_len(self, feature: str | None = None) -> int:
        return self._table.kmer_len(feature or self._features[0])

    def profile(self, seq: str) -> np.ndarray:
        """Values for one sequence. Always a 1-D array, never a result object."""
        if len(self._features) != 1:
            raise ValueError(
                "profile() requires a single feature; use profile_seqs() for "
                f"multi-feature profilers (has {len(self._features)})"
            )
        return np.asarray(self._values(self._features[0], seq), dtype=float)

    def profile_seqs(self, seqs):
        """Profile a list of sequences or an {id: sequence} mapping."""
        if isinstance(seqs, dict):
            pairs = list(seqs.items())
        else:
            pairs = [(f"seq_{i}", s) for i, s in enumerate(seqs)]
        return self._build(pairs)

    def _values(self, feature: str, seq: str) -> list[float]:
        return profile_values(
            seq, self._table.kmer_len(feature), self._table.table(feature),
            self.window_size,
        )

    def _build(self, pairs):
        """Turn (seqid, sequence) pairs into a FlexProfile or ProfileSet."""
        built = {}
        for feature in self._features:
            rows = [[sid, *self._values(feature, seq)] for sid, seq in pairs]
            built[feature] = FlexProfile(
                rows, feature=feature, window_size=self.window_size,
                kmer_len=self._table.kmer_len(feature),
            )
        if len(self._features) == 1:
            return built[self._features[0]]
        return ProfileSet(built)


def _resolve_lookup(lookup) -> FeatureTable:
    if lookup is None:
        return default_table()
    if isinstance(lookup, FeatureTable):
        return lookup
    if isinstance(lookup, dict):
        return FeatureTable.from_dict(lookup)
    if isinstance(lookup, (str, Path)):
        return FeatureTable.from_yaml(lookup)
    raise TypeError(
        f"lookup must be None, a path, a dict or a FeatureTable, got {type(lookup).__name__}"
    )
```

- [ ] **Step 4: Export the public API**

Replace `DNAflexpy/__init__.py`:

```python
"""DNAflexpy — DNA flexibility profiling from sequence."""

from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.profile import FlexProfile

__version__ = "0.3.0.dev0"
__all__ = [
    "FlexProfiler",
    "FlexProfile",
    "ProfileSet",
    "FeatureTable",
    "default_table",
    "__version__",
]
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_profiler.py -q`
Expected: PASS, 10 passed.

- [ ] **Step 6: Commit**

```bash
git add DNAflexpy tests/test_profiler.py
git commit -m "Add FlexProfiler with single-sequence and in-memory profiling

The lookup table is loaded once at construction rather than on every
per-sequence call, which is the 8ms defect. Features are validated in
__init__ so a typo fails immediately instead of after a full pass.

profile() takes a bare string and returns a 1-D array: no file, no
sequence ID to invent, no result object to unwrap. Return types stay
predictable - collection entry points always return FlexProfile or
ProfileSet - rather than varying with argument type.

Multi-feature runs return a ProfileSet because differing k means the
vectors cannot share a matrix."
```

---

### Task 7: `profile_fasta` with a pooled worker

**Files:**
- Create: `DNAflexpy/io.py`
- Modify: `DNAflexpy/core.py`
- Test: `tests/test_fasta.py`

**Interfaces:**
- Consumes: everything from Task 6.
- Produces:
  - `DNAflexpy.io.read_fasta(path) -> Iterator[tuple[str, str]]`
  - `FlexProfiler.profile_fasta(path, threads: int | None = None) -> FlexProfile | ProfileSet`

Task 8 compares `profile_fasta` output against the archive.

Workers receive the loaded table through `Pool(initializer=...)`. Passing a *path* for each worker to reload would reintroduce the 8 ms × N cost this whole plan exists to remove. The initializer and worker must be module-level functions: macOS spawns workers by re-importing the module, so closures cannot be pickled.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fasta.py`:

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_fasta.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'DNAflexpy.io'`.

- [ ] **Step 3: Implement `DNAflexpy/io.py`**

```python
"""Input readers."""
from __future__ import annotations

from pathlib import Path
from typing import Iterator


def read_fasta(path) -> Iterator[tuple[str, str]]:
    """Yield `(record_id, sequence)` for each record in a FASTA file.

    Handles wrapped sequence lines and a missing trailing newline.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)
```

- [ ] **Step 4: Add pooling to `DNAflexpy/core.py`**

Append these module-level helpers to `DNAflexpy/core.py`. They must be module-level: macOS spawns workers by re-importing the module, so closures and lambdas cannot be pickled.

```python
_WORKER_STATE: dict = {}


def _init_worker(table: FeatureTable, features: list[str], window_size: int) -> None:
    """Seed each worker once with the already-parsed table.

    Deliberately receives the loaded table, not a path: reloading per worker
    would reintroduce the 8 ms x N cost this design removes.
    """
    _WORKER_STATE["table"] = table
    _WORKER_STATE["features"] = features
    _WORKER_STATE["window_size"] = window_size


def _profile_record(record: tuple[str, str]) -> tuple[str, dict[str, list]]:
    seqid, seq = record
    table = _WORKER_STATE["table"]
    window_size = _WORKER_STATE["window_size"]
    return seqid, {
        feature: profile_values(
            seq, table.kmer_len(feature), table.table(feature), window_size
        )
        for feature in _WORKER_STATE["features"]
    }
```

Then add this method to `FlexProfiler`, directly after `profile_seqs`:

```python
    def profile_fasta(self, path, threads: int | None = None):
        """Profile every record in a FASTA file.

        One read and one pool spawn regardless of how many features are
        requested. With differing k the features cannot share a matrix, so
        the saving is I/O and process startup, not a shared k-mer scan.
        """
        from DNAflexpy.io import read_fasta

        records = list(read_fasta(path))
        if threads is not None and threads <= 1 or len(records) < 2:
            return self._build(records)

        import multiprocessing

        with multiprocessing.Pool(
            processes=threads,
            initializer=_init_worker,
            initargs=(self._table, self._features, self.window_size),
        ) as pool:
            results = pool.map(_profile_record, records)

        # pool.map preserves input order, so output stays deterministic.
        built = {}
        for feature in self._features:
            rows = [[seqid, *values[feature]] for seqid, values in results]
            built[feature] = FlexProfile(
                rows, feature=feature, window_size=self.window_size,
                kmer_len=self._table.kmer_len(feature),
            )
        if len(self._features) == 1:
            return built[self._features[0]]
        return ProfileSet(built)
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_fasta.py -q`
Expected: PASS, 8 passed.

- [ ] **Step 6: Commit**

```bash
git add DNAflexpy tests/test_fasta.py
git commit -m "Add profile_fasta with a pooled worker seeded once

Workers receive the already-parsed table through Pool(initializer=...).
Passing a path for each worker to reload would reintroduce the 8ms x N
cost this design exists to remove.

The initializer and worker are module-level functions because macOS
spawns workers by re-importing the module, so closures cannot be
pickled. pool.map preserves input order, keeping output deterministic
and byte-comparable.

read_fasta handles wrapped lines and a missing trailing newline, which
the packaged test FASTA lacks."
```

---

### Task 8: The differential byte-equality gate

**Files:**
- Create: `tests/test_differential.py`
- Create: `tests/test_ambiguous.fasta`

**Interfaces:**
- Consumes: `rxv.DNAflexpy` (Task 1), `FlexProfiler.profile_fasta` (Task 7), `FlexProfile.to_tsv` (Task 5).
- Produces: nothing — this is the Phase 2 gate.

Matrix: 3 FASTAs × the 10 features the archive can reach × window sizes {0, 2, 3, 10, 15, 26, 30}. Comparison is on **raw bytes**, because parsing would hide exactly the padding and float-repr differences the contract is about.

Allow-listed divergences, asserted explicitly so they cannot mask a regression:
- `gc`, `freeen`, `mechen` — the archive emits `None` rows (verified: `get_kmer_len` returns `None`, so the comparison `window_size < kmer_len` raises `TypeError` and is swallowed), the new package emits numbers.
- Unknown features — the archive emits `None` rows, the new package raises.
- Negative `window_size` — the archive silently returns `None`, the new package raises.
- Ambiguous nucleotides — the archive scores 0, the new package masks to NaN. Not reachable on the matrix above (verified: no test FASTA contains a non-ACGT character), so it gets its own fixture.

- [ ] **Step 1: Create the ambiguity fixture**

```bash
cat > tests/test_ambiguous.fasta <<'EOF'
>clean
ATGCGTACGTAGCTAGCGTAGCTAGT
>has_n
ATGNGTACGTAGCTAGCGTAGCTAGT
>has_iupac
ATGRGTACGTAGCTAGCGTAGCTAGT
>all_n
NNNNNNNNNN
EOF
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_differential.py`:

```python
"""Phase 2 gate: new output must be byte-identical to the archive."""
import contextlib
import io
import math
import pathlib

import pandas as pd
import pytest

from DNAflexpy.core import FlexProfiler
from rxv.DNAflexpy.utils import (
    get_kmer_len,
    load_feature_data,
    read_fasta as archive_read_fasta,
    seq_to_numeric_profile,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent

FASTAS = {
    "test_fasta": ROOT / "rxv/DNAflexpy/data/test_fasta.fa",
    "edge": ROOT / "tests/test_edge_cases.fasta",
    "exact": ROOT / "tests/test_exact_kmer.fasta",
}

# The 10 features the archive's hardcoded map can reach. gc, freeen and
# mechen are excluded on purpose: they are allow-listed divergences below.
ARCHIVE_FEATURES = [
    "DNaseI", "NPP", "bendabilityDNase", "bendabilityConcensus",
    "wedge", "prop", "twistDisp", "stiffness", "bendingStiffness", "trx",
]
WINDOWS = [0, 2, 3, 10, 15, 26, 30]


def archive_bytes(fasta, feature, window):
    lookup = load_feature_data()
    kmer_len = get_kmer_len(feature)
    rows = []
    with contextlib.redirect_stdout(io.StringIO()):
        for seqid, seq in archive_read_fasta(str(fasta)):
            rows.append(
                seq_to_numeric_profile(seqid, seq, kmer_len, window, feature, lookup)
            )
    buf = io.StringIO()
    pd.DataFrame(rows).to_csv(buf, index=False, header=False, sep="\t")
    return buf.getvalue().encode()


def new_bytes(fasta, feature, window, tmp_path):
    out = tmp_path / "new.tsv"
    FlexProfiler(feature, window_size=window).profile_fasta(fasta, threads=1).to_tsv(out)
    return out.read_bytes()


@pytest.mark.parametrize("fasta_key", sorted(FASTAS))
@pytest.mark.parametrize("feature", ARCHIVE_FEATURES)
@pytest.mark.parametrize("window", WINDOWS)
def test_byte_identical_to_archive(fasta_key, feature, window, tmp_path):
    fasta = FASTAS[fasta_key]
    assert new_bytes(fasta, feature, window, tmp_path) == archive_bytes(
        fasta, feature, window
    )


@pytest.mark.parametrize("feature", ["gc", "freeen", "mechen"])
def test_allowlisted_unlocked_features(feature, tmp_path):
    """Archive: None rows, because get_kmer_len has no entry. New: numbers."""
    lookup = load_feature_data()
    assert get_kmer_len(feature) is None
    with contextlib.redirect_stdout(io.StringIO()):
        archived = seq_to_numeric_profile(
            "s", "ATGCGTACGT", get_kmer_len(feature), 10, feature, lookup
        )
    assert archived is None

    values = FlexProfiler(feature, window_size=0).profile("ATGCGTACGT")
    assert len(values) == 9
    assert not any(math.isnan(v) for v in values)

    # And the whole-file path produces a real TSV, not a column of "None".
    text = new_bytes(FASTAS["test_fasta"], feature, 10, tmp_path).decode()
    assert "None" not in text
    assert len(text.splitlines()) == 2


def test_allowlisted_unknown_feature_now_raises():
    with pytest.raises(ValueError, match="unknown feature"):
        FlexProfiler("not_a_feature")


def test_allowlisted_negative_window_now_raises():
    """The archive silently returns None; the new package refuses."""
    lookup = load_feature_data()
    with contextlib.redirect_stdout(io.StringIO()):
        assert seq_to_numeric_profile("s", "ATGCGT", 3, -1, "DNaseI", lookup) is None
    with pytest.raises(ValueError, match="window_size must be >= 0"):
        FlexProfiler("DNaseI", window_size=-1)


def test_allowlisted_ambiguity_masking(tmp_path):
    """Archive scores N as 0; the new package masks it. Off the matrix above."""
    fasta = ROOT / "tests/test_ambiguous.fasta"
    prof = FlexProfiler("DNaseI", window_size=0).profile_fasta(fasta, threads=1)

    assert prof.n_masked["clean"] == 0
    assert prof.n_masked["has_n"] == 3
    assert prof.n_masked["has_iupac"] == 3
    assert prof.n_masked["all_n"] == 8

    archived = archive_bytes(fasta, "DNaseI", 0).decode().splitlines()
    n_row = next(r for r in archived if r.startswith("has_n"))
    assert "\t0\t0\t0\t" in n_row  # the archive's fabricated zeros


def test_no_test_fasta_contains_ambiguity():
    """Guards the claim that masking cannot affect the differential matrix."""
    for fasta in FASTAS.values():
        seqs = "".join(s for _, s in archive_read_fasta(str(fasta))).upper()
        assert not set(seqs) - set("ACGT"), fasta
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `python -m pytest tests/test_differential.py -q`
Expected: FAIL — `tests/test_ambiguous.fasta` may not exist yet, or byte mismatches appear. Investigate any mismatch before proceeding; a failure here means the kernel or the byte contract diverged.

- [ ] **Step 4: Fix any divergence**

Do not weaken the assertions. If a case mismatches, diff the two byte strings and trace it to one of:
- summation order in `kernel._mean` (must be builtin `sum` on a list, left to right)
- `round(v, 3)` applied at the wrong point (the archive rounds per window, not at the end)
- the ragged-row construction in `FlexProfile.to_tsv` (rows must stay ragged; do not pre-pad)
- `window_size == 0` returning window means rather than per-k-mer values

- [ ] **Step 5: Run the full suite**

Run: `python -m pytest -q`
Expected: PASS. 210 differential cases plus every earlier test.

- [ ] **Step 6: Commit**

```bash
git add tests/test_differential.py tests/test_ambiguous.fasta
git commit -m "Add the differential byte-equality gate

Compares the new package against the archive on 3 FASTAs x 10 features x
7 window sizes, asserting on raw bytes: parsing would hide exactly the
NaN-padding and float-repr differences the contract is about.

Every intended divergence is allow-listed and asserted explicitly so it
cannot mask a regression: gc/freeen/mechen become reachable, unknown
features raise instead of yielding None rows, and ambiguous nucleotides
mask instead of scoring 0. A test asserts no test FASTA contains a
non-ACGT character, which is what keeps masking off the matrix."
```

---

### Task 9: Memoised module-level `profile()`

**Files:**
- Modify: `DNAflexpy/__init__.py`
- Test: `tests/test_module_level.py`

**Interfaces:**
- Consumes: `FlexProfiler` (Task 6), `default_table` (Task 3).
- Produces: `DNAflexpy.profile(seq: str, feature: str = "DNaseI", window_size: int = 10) -> numpy.ndarray`

This is the one place the 8 ms defect could return: a naive implementation would build a fresh `FeatureTable` per call. `default_table` is already `lru_cache`d, so the convenience wrapper inherits the memoisation — the test proves it.

- [ ] **Step 1: Write the failing test**

Create `tests/test_module_level.py`:

```python
import numpy as np

import DNAflexpy
from DNAflexpy.lookup import default_table


def test_one_liner_needs_no_object():
    out = DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT", feature="DNaseI",
                            window_size=10)
    assert isinstance(out, np.ndarray)
    assert len(out) == 17


def test_defaults_match_the_documented_cli_defaults():
    assert len(DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT")) == 17


def test_matches_the_class_api():
    seq = "ATGCGTACGTAGCTAGCGTAGCTAGT"
    expected = DNAflexpy.FlexProfiler("trx", window_size=5).profile(seq)
    assert np.array_equal(DNAflexpy.profile(seq, feature="trx", window_size=5), expected)


def test_repeated_calls_parse_the_yaml_once(monkeypatch):
    """The 8ms defect must not sneak back in through the convenience path."""
    default_table.cache_clear()
    calls = {"n": 0}
    original = DNAflexpy.lookup.FeatureTable.from_yaml

    def counting(*args, **kwargs):
        calls["n"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(DNAflexpy.lookup.FeatureTable, "from_yaml", counting)
    for _ in range(1000):
        DNAflexpy.profile("ATGCGTACGT")
    assert calls["n"] == 1
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_module_level.py -q`
Expected: FAIL with `AttributeError: module 'DNAflexpy' has no attribute 'profile'`.

- [ ] **Step 3: Add the wrapper**

Replace `DNAflexpy/__init__.py`:

```python
"""DNAflexpy — DNA flexibility profiling from sequence."""

from DNAflexpy import lookup
from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.profile import FlexProfile

__version__ = "0.3.0.dev0"
__all__ = [
    "FlexProfiler",
    "FlexProfile",
    "ProfileSet",
    "FeatureTable",
    "default_table",
    "profile",
    "lookup",
    "__version__",
]


def profile(seq: str, feature: str = "DNaseI", window_size: int = 10):
    """Profile one sequence without constructing anything.

    The packaged table is memoised by `default_table`, so the first call
    pays the YAML parse and later calls do not.
    """
    return FlexProfiler(feature, window_size=window_size).profile(seq)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_module_level.py -q`
Expected: PASS, 4 passed.

- [ ] **Step 5: Verify the performance claim**

Run:
```bash
python -c "
import time, DNAflexpy
seq='ATGCGTACGTAGCTAGCGTAGCTAGT'*4
DNAflexpy.profile(seq)                      # warm the cache
t=time.perf_counter()
for _ in range(200): DNAflexpy.profile(seq)
print(f'{(time.perf_counter()-t)/200*1000:.3f} ms/seq (archive was 8.0)')"
```
Expected: well under 1 ms per sequence, against the archive's 8.0 ms.

- [ ] **Step 6: Run the full suite and commit**

```bash
python -m pytest -q
git add DNAflexpy/__init__.py tests/test_module_level.py
git commit -m "Add memoised module-level profile()

A one-liner for notebooks that constructs nothing. This is the one place
the 8ms defect could return, since a naive version would build a fresh
FeatureTable per call; a test asserts 1000 calls parse the YAML exactly
once."
```

---

## Definition of done

- [ ] `python -m pytest -q` passes with no failures and no errors.
- [ ] `tests/test_differential.py` passes on all 210 matrix cases.
- [ ] `python -m pip install -e .` succeeds and both `import DNAflexpy` and `import rxv.DNAflexpy` work.
- [ ] `rxv/DNAflexpy/` differs from its pre-move state by exactly one line (`utils.py:175`). Verify: `git diff a9bd674 -- rxv/DNAflexpy | grep -c '^[+-][^+-]'` returns `2`.
- [ ] Per-sequence profiling is under 1 ms, against the archive's 8.0 ms.
- [ ] No commit contains a `Co-Authored-By` trailer.

## Deferred to later phases

Phases 3–9 of the spec are out of scope here: `profile_table`/`FlexProfile.y`, `from_bed`, `encode.py`, plotting, dinucleotide shuffle and per-position statistics, streaming input, the provenance sidecar, `ProfileSet.correlate()`, bigWig/BED export, the `DNAflexpy` CLI, and the container. `FlexProfile.zscale`/`normalize` are declared in the spec's component list but belong with plotting in Phase 5.
