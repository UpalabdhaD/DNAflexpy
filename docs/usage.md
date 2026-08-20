# Usage

## Installing

```bash
pip install -e .            # the package
pip install -e '.[bed]'     # plus BED / reference-genome input (needs pyfaidx)
pip install -e '.[dev]'     # plus the test suite
```

## The quickest thing

One sequence, no setup:

```python
import DNAflexpy

DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT", feature="DNaseI", window_size=10)
```

That returns a NumPy array of values, one per window.

## Profiling many sequences

Build a `FlexProfiler` once and reuse it. The lookup table is read from disk
only when you construct it, so reusing the profiler is roughly 130x faster than
paying for the table on every call.

```python
from DNAflexpy import FlexProfiler

p = FlexProfiler(feature="DNaseI", window_size=10)
```

There are five ways to give it sequences.

### 1. A single string

```python
values = p.profile("ATGCGTACGTAGCTAGCGTAGCTAGT")
```

Returns a 1-D array. Use this when you just want the numbers.

### 2. A list or a dict

```python
prof = p.profile_seqs(["ATGCGTACGT", "CGTAGCTAGT"])       # ids become seq_0, seq_1
prof = p.profile_seqs({"promoter": "ATGCGTACGT"})          # or name them yourself
```

### 3. A FASTA file

```python
prof = p.profile_fasta("sequences.fa")
prof = p.profile_fasta("sequences.fa", threads=4)          # use several processes
```

`threads` spawns processes, not threads. Below 64 records it always runs
serially anyway, because starting processes costs more than the work saves.

### 4. A table of sequences and labels

For model fitting. Each row holds a sequence and a number, such as a binding
affinity:

```
ATGCGTACGT	1.5
CGTAGCTAGT	2.5
```

```python
prof = p.profile_table("affinity.tsv", header=False)
X, y = prof.frame, prof.y
```

`header` is required. It is not guessed, because words like `dna`, `rna` and
`tag` are spelled with nucleotide letters, so a header row cannot reliably be
told apart from a sequence.

Columns can be named instead of positional. That needs a file with a header
row, like this one:

```
sequence	affinity
ATGCGTACGT	1.5
CGTAGCTAGT	2.5
```

```python
prof = p.profile_table(
    "labelled.tsv", seq_col="sequence", value_col="affinity", header=True
)
```

A row whose label is missing or not a number raises an error naming the line.
It is never skipped: silently dropping labelled rows would corrupt a training
set without telling you.

### 5. A BED file and a reference genome

```python
prof = p.from_bed("peaks.bed", genome="TAIR10.fa", width=200)
```

`width` re-centres every interval on its midpoint and cuts it to exactly that
many bases, so all rows come out the same length. That is what makes
position-by-position comparison meaningful. A zero-length BED entry marking a
TSS works fine with `width`.

Near the start or end of a chromosome a centred window can run off the edge.
`on_edge` decides what happens:

- `"drop"` (the default) skips the interval and warns you how many went
- `"error"` stops and names them
- `"pad"` pads with `N` and warns

Needs `pip install -e '.[bed]'`.

## What you get back

All five return the same object, a `FlexProfile`:

```python
prof.frame          # a pandas DataFrame, one row per sequence
prof.seqids         # the sequence names, in input order
prof.y              # the labels, or None if the input had none
prof.n_masked       # how many positions could not be scored
prof.to_tsv("out.tsv")
```

Ask for several features at once and you get a `ProfileSet`, which behaves
like a dictionary:

```python
ps = FlexProfiler(feature=["DNaseI", "trx"], window_size=10)
result = ps.profile_seqs(["ATGCGTACGT"])
result["DNaseI"].frame
```

Each feature is kept separate because they use different k-mer lengths, so
their value counts differ for the same sequence.

## Features

Thirteen features are available. Four use 3-mers:

`DNaseI`, `NPP`, `bendabilityDNase`, `bendabilityConcensus`

Nine use 2-mers:

`wedge`, `prop`, `freeen`, `gc`, `twistDisp`, `stiffness`, `bendingStiffness`,
`mechen`, `trx`

```python
from DNAflexpy import default_table

default_table().features            # all thirteen names
default_table().kmer_len("DNaseI")  # 3
```

You can supply your own table instead of the packaged one:

```python
p = FlexProfiler(
    feature="mine",
    lookup={"mine": {"AA": 1.0, "AT": 2.0, "TA": 3.0, "TT": 4.0}},
)
```

The k-mer length is worked out from the keys, so there is nothing else to keep
in step.

## Window size

- `window_size=0` gives one value per k-mer, with no averaging.
- `window_size=N` slides a window of `N` bases one step at a time and gives the
  mean of the k-mer values inside it. You get `len(sequence) - N + 1` values.
- If `N` equals the feature's k-mer length, each window holds exactly one
  k-mer, so the output matches `window_size=0`.
- If `N` is larger than the sequence, no windows fit and only the sequence name
  is returned.
- If `N` is smaller than the k-mer length, no k-mer fits and every window is
  `0.0`.

## Letters that are not A, C, G or T

The feature tables only hold ACGT k-mers. Any k-mer covering an `N` or an
ambiguity code such as `R` cannot be looked up, so it is masked rather than
scored as zero — zero is a real value in these tables, so scoring a gap as zero
would be indistinguishable from a measurement.

You get a warning naming the letters found. `prof.n_masked` counts a position
as masked only when *no* k-mer in that window resolved, so at window sizes
above zero it under-reports; a window with some unusable k-mers is simply
averaged over fewer of them.

## Encoding for machine learning

A profile is a table of numbers. To fit a model you usually want a *design
matrix*: one row per sequence, one column per predictor, plus a target vector.
`encode` builds one.

```python
p = FlexProfiler(feature="DNaseI", window_size=0)
prof = p.profile_table("affinity.tsv", header=False)

fm = prof.encode(["1-mer", "1-DNaseI"])
print(fm.shape)          # (2, 48)
print(fm.y)              # [1.5, 2.5] -- the labels, straight from the file
```

`fm.X` is a plain NumPy array and `fm.y` is the label list, so they go into
scikit-learn directly:

```python
# from sklearn.linear_model import Ridge
# Ridge().fit(fm.X, fm.y)
```

### What the names mean

You pass a list. Each entry asks for one block of columns, and the blocks are
joined left to right in the order you wrote them.

| Name | What you get |
|---|---|
| `1-mer` | one-hot sequence, 4 binary columns per position |
| `2-mer` | one-hot dinucleotides, 16 columns per position |
| `3-mer` | one-hot trinucleotides, 64 columns per position |
| `1-DNaseI` | the profile values themselves, one column per position |
| `2-DNaseI` | products of the feature at two adjacent positions |
| `3-DNaseI` | products at three adjacent positions |

Higher orders capture interactions between neighbouring positions that a linear
model would otherwise miss. `2-DNaseI` gives you *only* the second-order terms,
so ask for both if you want both:

```python
fm = prof.encode(["1-DNaseI", "2-DNaseI"])
```

Any feature in the lookup table works, and you can mix features by profiling
several at once:

```python
result = FlexProfiler(["gc", "stiffness"], window_size=0).profile_seqs(
    ["ATGCGTACGT", "CGTAGCTAGT"]
)
fm = result.encode(["1-mer", "1-gc", "2-stiffness"])
```

### Column names

Every column is named, so you can tell afterwards what a coefficient refers to:

```python
fm = prof.encode(["1-DNaseI"])
print(fm.columns[:3])    # ['DNaseI.w0.o1.p1', 'DNaseI.w0.o1.p2', ...]
```

Read it as feature, window size, order, position. The window size is in there
on purpose: at `window_size=10`, position 1 is the average over bases 1-10,
while at `window_size=0` it is the k-mer starting at base 1. Two matrices built
at different window sizes would otherwise carry the same column names for
different things.

`fm.to_frame()` gives a pandas DataFrame with those names, indexed by sequence
ID. It is a convenience, not the default: 3-mer encoding of 500-base sequences
is about 32,000 columns, and a DataFrame that wide is slow.

### Normalising

`normalize=True` is the default. It min-max scales each feature block to
`[0, 1]`, on its own range. Each block separately, because the features are on
wildly different scales -- `stiffness` runs to 5500 and `gc` to 1.0, so a single
shared scale would squash `gc` to nothing. One-hot blocks are left alone; they
are already 0 and 1.

The scaling uses **your dataset's own** minimum and maximum. If you are going to
split into train and test sets, that leaks the test set into the scaling. Pass
`normalize=False` and scale inside a scikit-learn pipeline instead:

```python
fm = prof.encode(["1-DNaseI"], normalize=False)
```

### Two things that will stop you

**Sequences must all be the same length.** Columns are positions, so unequal
lengths have no meaning. `from_bed(..., width=200)` is the usual way to get
fixed-width input.

**Missing values stay missing.** A position covering an `N`, an ambiguity code,
or `on_edge="pad"` padding comes out as `NaN` rather than 0, so it stays
distinguishable from a real measurement of zero. Most scikit-learn estimators
refuse `NaN`; use `sklearn.impute.SimpleImputer`, or a model that handles it
such as `HistGradientBoostingRegressor`.

One last thing: `DNAflexpy.profile()` returns a bare array, so there is nothing
to encode from the one-liner. Start from `profile_seqs`, `profile_fasta`,
`profile_table` or `from_bed`.

## Other things you can do

### A long (tidy) table instead of a wide one

`frame` is wide: one row per sequence, one column per position. For plotting
libraries you usually want it long instead:

```python
prof.to_frame(tidy=True)
```

That gives one row per value, with columns `seqid`, `position`, `value` and
`feature`. For several features at once:

```python
result.to_tidy()
```

which stacks every feature into a single long table.

### CSV instead of tab-separated

```python
prof = p.profile_table("affinity.csv", sep=",", header=True)
```

### Keeping your own row names

If your table has a name column, point at it and those names are used instead
of `seq_0`, `seq_1`:

```
name	sequence	affinity
p1	ATGCGTACGT	1.5
```

```python
prof = p.profile_table(
    "named.tsv",
    id_col="name", seq_col="sequence", value_col="affinity",
    header=True,
)
```

### A lookup table from your own YAML file

```python
from DNAflexpy import FeatureTable

table = FeatureTable.from_yaml("my_features.yaml")
p = FlexProfiler(feature="my_feature", lookup=table)
```

The file has the same shape as the packaged one: each feature name maps to
k-mers, and each k-mer to a number. Mixed k-mer lengths inside one feature are
rejected, because the length is inferred from the keys. An incomplete table is
allowed but warns you, since missing k-mers become masked values.

### Inspecting a table

```python
table = FeatureTable.from_yaml()      # the packaged one
table.features                        # every feature name
table.kmer_len("trx")                 # 2
table.table("trx")["AA"]              # the value for one k-mer
```

`table()` hands back a read-only view, so you cannot corrupt the shared copy
by accident.

### Reading files without profiling them

The readers are useful on their own:

```python
from DNAflexpy.io import read_fasta, read_table, read_bed

list(read_fasta("sequences.fa"))              # [(id, sequence), ...]
read_table("affinity.tsv", header=False)      # [(id, sequence, value), ...]
read_bed("peaks.bed")                         # [(chrom, start, end, name, strand), ...]
```

`read_bed` keeps BED's own numbering: counting starts at zero, and the end
position is not included.

## Command line

The rewritten package has no command-line tool yet. The original one still
works and is kept in the repository:

```bash
python -m rxv.DNAflexpy.cli sequences.fa --window-size 10 --feature DNaseI --outfile out.tsv
```
