# API reference

Everything importable from `DNAflexpy`.

```python
from DNAflexpy import (
    FlexProfiler, FlexProfile, ProfileSet, FeatureTable, default_table, profile
)
```

---

## `profile(seq, feature="DNaseI", window_size=10)`

Profile one sequence without building anything. Returns a 1-D NumPy array.

The packaged lookup table is cached, so calling this in a loop does not re-read
it from disk.

---

## `FlexProfiler(feature, window_size=10, lookup=None)`

The engine. Reads the lookup table once, when you construct it.

| argument | meaning |
|---|---|
| `feature` | one feature name, or a list of them |
| `window_size` | `0` for per-k-mer values, otherwise the sliding window width |
| `lookup` | `None` for the packaged table, or a path, a dict, or a `FeatureTable` |

An unknown feature name raises straight away, rather than after reading your
whole file. A negative `window_size` raises.

### `.profile(seq)`
One sequence as a string. Returns a 1-D array. Only for single-feature
profilers.

### `.profile_seqs(seqs)`
A list of sequences (named `seq_0`, `seq_1`, …) or a `{name: sequence}` dict.

### `.profile_fasta(path, threads=None)`
Every record in a FASTA file. `threads` spawns processes. Files under 64
records always run serially, whatever you pass.

### `.profile_table(path, seq_col=0, value_col=1, id_col=None, header=None, sep="\t")`
A table of sequences and numeric labels. The labels become `FlexProfile.y`.

`header` must be `True` or `False`. It is not guessed. Columns may be given by
position or by name; names require a header. A missing or non-numeric label
raises, naming the line in the file.

### `.from_bed(bed, genome, width=None, on_edge="drop")`
Intervals from a BED file, with sequences taken from a reference genome FASTA.

`width` re-centres each interval and cuts it to that many bases, so every row
is the same length. `on_edge` handles windows that run past the end of a
chromosome: `"drop"`, `"error"` or `"pad"`.

Needs the `bed` extra (`pip install -e '.[bed]'`).

### `.features`
The feature names this profiler was built with.

### `.kmer_len(feature=None)`
The k-mer length for a feature, defaulting to the first one.

---

## `FlexProfile`

What every profiling call returns.

| member | meaning |
|---|---|
| `.frame` | wide DataFrame: one row per sequence, one column per position |
| `.seqids` | sequence names, in input order |
| `.y` | the labels from a table, or `None` |
| `.seqs` | the input sequences, in the same order, or `None` |
| `.n_masked` | `{name: count}` of positions that could not be scored |
| `.feature`, `.window_size`, `.kmer_len` | what produced it |

### `.to_tsv(path)`
Write the tab-separated format the original package produced, exactly.
Sequences shorter than the window become a name with no values.

### `.to_frame(tidy=False)`
`tidy=False` returns `.frame`. `tidy=True` returns a long table with columns
`seqid`, `position`, `value`, `feature`.

### `.encode(feature_names, normalize=True)`
Build a `FeatureMatrix`. Short for `DNAflexpy.encode.encode(self, ...)`.

---

## `ProfileSet`

A dict of `{feature: FlexProfile}`, returned when you ask for several features
at once. Features are kept apart because different k-mer lengths give different
numbers of values for the same sequence.

### `.to_tidy()`
Every feature stacked into one long DataFrame.

### `.encode(feature_names, normalize=True)`
Build a `FeatureMatrix` spanning several features at once.

---

## `DNAflexpy.encode`

### `encode(profiles, feature_names, normalize=True)`

Turn a `FlexProfile` or a `ProfileSet` into a design matrix. `feature_names` is
a list; each entry names one block of columns and the blocks are joined in the
order given.

| entry | block |
|---|---|
| `"{k}-mer"` | one-hot sequence encoding, `4**k` binary columns per position |
| `"{n}-{feature}"` | the nth-order terms of a profiled feature: products at `n` adjacent positions. `n=1` is the values themselves |

Column names are `seq.{k}mer.p{i}.{KMER}` and `{feature}.w{window}.o{n}.p{i}`,
positions counted from 1. The window size is part of the name because it
changes what a position means.

`normalize=True` min-max scales each feature block to `[0, 1]` on its own
range. One-hot blocks are never scaled. The range comes from the data you
passed, so pass `normalize=False` and scale in a scikit-learn pipeline if you
are splitting train and test.

Raises if the sequences are not all the same length, if a name does not parse,
if a requested feature was not profiled, if a name is repeated, or if a `k-mer`
term is asked for on a profile that did not keep its sequences.

`NaN` is preserved, never filled.

### `FeatureMatrix`

| member | meaning |
|---|---|
| `.X` | the matrix, a NumPy array of floats |
| `.columns` | one name per column |
| `.seqids` | one name per row, in input order |
| `.y` | the labels, if the profile carried any |
| `.window_size` | the window the profile was built with |
| `.feature_names` | the request, as passed |
| `.shape` | `(rows, columns)` |
| `.to_frame()` | a pandas DataFrame with those names and index |

---

## `FeatureTable`

A validated set of lookup tables.

### `FeatureTable.from_yaml(path=None)`
Load from a YAML file. `None` loads the packaged table.

### `FeatureTable.from_dict(mapping)`
Load from `{feature: {kmer: value}}`.

Both check the table as they load it. Mixed k-mer lengths inside one feature
are rejected, because the k-mer length is read from the keys themselves. Keys
outside A, C, G, T are rejected. Non-numeric values are rejected. A table with
fewer k-mers than it should have is allowed but warns.

### `.features`
Every feature name.

### `.kmer_len(feature)`
The k-mer length, worked out from the keys.

### `.table(feature)`
The `{kmer: value}` mapping, read-only.

---

## `default_table()`

The packaged lookup table, parsed at most once per session.

---

## `DNAflexpy.io`

The readers, usable on their own.

### `read_fasta(path)`
Yields `(name, sequence)`. Handles wrapped lines and a missing final newline.

### `read_table(path, seq_col=0, value_col=1, id_col=None, header=None, sep="\t")`
Returns `(name, sequence, value)` triples.

### `read_bed(path)`
Returns `(chrom, start, end, name, strand)` tuples, keeping BED's own
zero-based numbering with the end position excluded. `track`, `browser` and
`#` lines are skipped; a contig whose name merely begins with those words is
not.

### `extract_intervals(intervals, genome, width=None, on_edge="drop")`
Turns those tuples into `(name, sequence)` pairs using a reference genome.
Intervals on the minus strand are reverse-complemented. Imports `pyfaidx` only
when called, so the rest of the package works without it.

### `warn_if_ambiguous(records, source)`
Warns when sequences contain letters outside A, C, G, T, N. Every profiling
entry point calls this once already.
