# DNAflexpy — Greenfield Rewrite Design

**Date:** 2026-08-12
**Status:** Awaiting user review
**Supersedes:** the old `DNAflexpy/`, which is archived to `rxv/DNAflexpy/`

## Goal

Rewrite DNAflexpy in place. The existing package is archived to `rxv/DNAflexpy/`,
where its source is preserved verbatim (one import-path line excepted, see below)
and stays importable. A new package reusing the `DNAflexpy` name is written from
scratch alongside it.

The rewrite fixes a measured performance defect, replaces the function-based API
with classes, and adds ML feature encoding, plotting, and genomic input modelled
on the R/Bioconductor package DNAshapeR.

## Motivating measurements

Taken from the current package at HEAD (`a9bd674`):

| Path | Cost | Cause |
|---|---|---|
| `load_feature_data()` | 8.1 ms | re-parses a 400-entry YAML |
| `DNAflexpy()` per-sequence API | 8.0 ms/seq | calls `load_feature_data()` on *every* call and discards the `feature_lookup` argument it was passed (`core.py:18`) |
| `seq_to_numeric_profile()` | 0.1 ms/seq | the actual computation |

99% of per-sequence time is YAML re-parsing — an 80× overhead, ~80 s per 10k
sequences. The notebooks hit this directly: they loop
`for id, seq in read_fasta(...)` and call per-sequence functions.

The lookup dict is only 5202 bytes pickled, so multiprocessing argument
serialisation is *not* a bottleneck and needs no optimisation.

## Key finding: k-mer length is inferrable

`get_kmer_len` hardcodes a feature→k map that must be hand-synced with the YAML.
Verified against `lookupNEW.yaml`: every feature's table has **uniform, complete
keys** (64 3-mers or 16 2-mers), no null values, no non-ACGT keys. Inferring k
from key length reproduces the hardcoded map exactly on all 10 mapped features
and correctly yields k=2 for the three unmapped ones.

Therefore the new package infers k from the table and has no hardcoded map. This single
change unlocks the 3 dead features and makes arbitrary user-supplied tables work.

### Consequence for error handling

The approved "fail fast on setup" behaviour was illustrated with
`FlexProfiler(feature="gc")` raising "no k-mer length mapping". Under inference
that error no longer exists: `gc`, `freeen`, and `mechen` become valid k=2
features returning real numbers where the old package returned `None` rows.

This is an intended behaviour change, not a regression, and the differential test
suite partitions it as such. Fail-fast still applies to:

- a feature absent from the lookup table entirely
- a malformed table (mixed key lengths, non-numeric values, non-ACGT keys)

Per-sequence data conditions (window > sequence length, window < k) remain
warnings that let batch processing continue, matching documented behaviour.

## Package layout

The old package is archived to `rxv/DNAflexpy/` via `git mv` (history preserved)
and the new package takes the `DNAflexpy` name and CLI. The name the user cares
about carries forward to the new code; the archive stays importable, which is
what keeps differential testing possible.

```
rxv/
  __init__.py         namespace only
  DNAflexpy/          <- ARCHIVED old package, importable as rxv.DNAflexpy
    core.py  utils.py  cli.py  data/
DNAflexpy/            <- NEW, written from scratch
  __init__.py           FlexProfiler, FlexProfile, ProfileSet, __version__, __all__
  lookup.py             FeatureTable: load + validate + infer k
  io.py                 read_fasta, from_bed
  core.py               FlexProfiler: profiling engine
  profile.py            FlexProfile / ProfileSet: results, TSV, DataFrame
  encode.py             ML feature matrices
  plotting.py           heatmap / metaprofile / trackplot (lazy matplotlib)
  cli.py                `DNAflexpy` command
  data/lookupNEW.yaml   own copy — no shared state with the old package
```

Both packages install side by side; import names never collide. `rxv/DNAflexpy/`
is moved with `git mv` so its history follows it. Its one absolute self-reference
-- `files("DNAflexpy.data")` at `utils.py:175` -- must become
`files("rxv.DNAflexpy.data")`, or the archived code would silently read the NEW
package's lookup table. That single line is the only edit made to old code, and
it is forced by the rename itself.

Version lives only in `DNAflexpy/__init__.py.__version__`, read by packaging
metadata — avoiding the setup.py/requirements.txt drift in the old package.

### Build configuration — needs user sign-off

Root build config must change: `setup.py` calls `find_packages()`, which would
sweep the archive into the distribution while shipping none of its `data/*.yaml`.

Build metadata migrates into the root `pyproject.toml`, declaring one
distribution that ships the new package and the archive:

```toml
[project]
name = "DNAflexpy"
dynamic = ["version"]
scripts = { DNAflexpy = "DNAflexpy.cli:main" }   # new package owns the name

[tool.setuptools]
packages = ["DNAflexpy", "rxv", "rxv.DNAflexpy"]   # explicit, never find_packages()
package-data = { DNAflexpy = ["data/*.yaml"], "rxv.DNAflexpy" = ["data/*.yaml"] }
```

The archive registers no console script; it stays runnable as
`python -m rxv.DNAflexpy.cli` for verification runs. `setup.py` is deleted —
leaving it beside a `[project]` table makes setuptools error on conflicting
metadata. One `pip install -e .` yields both packages, which is what lets the
differential suite import old and new in a single environment. Dropping the
archive later means deleting `rxv/` and its two list entries.

Files the plan modifies that the user already had, in full: `setup.py` (deleted,
superseded), `pyproject.toml`, `MANIFEST.in`, `rxv/DNAflexpy/utils.py:175` (the
forced import-path fix), and `scripts/plot_profiles.py`, which shells out to
`python -m DNAflexpy.cli` and would otherwise silently invoke the new CLI with
old flags.

## Components

### FeatureTable (`lookup.py`)

Loaded exactly once, then reused. This is the fix for the 8 ms defect.

```python
FeatureTable.from_yaml(path=None)   # None -> packaged table
FeatureTable.from_dict(mapping)
.features -> list[str]
.kmer_len(feature) -> int           # inferred, not hardcoded
.table(feature) -> dict[str, float]
```

Validation, on load, for both packaged and user tables:

- reject mixed key lengths within a feature (this is what makes k inference safe)
- uppercase keys; reject non-ACGT characters
- reject non-numeric values
- **warn**, not fail, when a table has fewer than `4**k` entries

An incomplete table is usable but flagged, because any valid ACGT k-mer missing
from it resolves to `NaN` and is masked — the same path as an ambiguous
nucleotide, described below. This is the single rule for every unresolvable
k-mer: the old package's `.get(subseq, 0)` zero-fill is not reproduced anywhere,
since a zero is indistinguishable from a real measurement in these tables.

### FlexProfiler (`core.py`)

```python
FlexProfiler(feature, window_size=10, lookup=None)
    feature: str | list[str]
    lookup:  None | path | dict | FeatureTable
.profile(seq) -> np.ndarray                       # ONE bare string, ~0.1 ms
.profile_seqs(seqs) -> FlexProfile | ProfileSet   # list[str] | dict[id, str]
.profile_fasta(path, threads=None) -> FlexProfile | ProfileSet
.profile_table(path, seq_col=0, value_col=1) -> FlexProfile | ProfileSet
.from_bed(bed, genome, width=None) -> FlexProfile | ProfileSet
```

Features are validated against the table in `__init__`, so a typo fails
immediately rather than after a full FASTA pass.

**No file is required anywhere.** The old API forced a FASTA round-trip even to
inspect a single sequence; here `profile()` takes a bare string and returns the
value array directly — no sequence ID to invent, no result object to unwrap:

```python
p = FlexProfiler("DNaseI", window_size=10)
p.profile("ATGCGTACGTAGCTAGCGTAGCTAGT")     # -> array([0.011, -0.003, ...])
```

Return types stay predictable rather than polymorphic: `profile()` always returns
a 1-D array for one sequence, and the collection entry points always return
`FlexProfile`/`ProfileSet`. A single function returning different types depending
on its argument would cost more in misuse than the extra method name saves.

For notebooks and quick checks, a module-level one-liner constructs nothing:

```python
import DNAflexpy
DNAflexpy.profile("ATGCGTACGT", feature="DNaseI", window_size=10)
```

This is the one place the 8 ms defect could sneak back in — a naive
implementation would build a fresh `FeatureTable` per call. The default table is
therefore memoised at module level (`functools.lru_cache`), so the first call
pays 8 ms and every subsequent call pays ~0.1 ms. A regression test asserts that
1000 calls to `DNAflexpy.profile()` parse the YAML exactly once.

Workers receive the loaded table through `Pool(initializer=...)`. The table is
**never** passed as a path to be reloaded per worker — that would reintroduce the
8 ms × N cost being fixed.

Multi-feature runs return a `ProfileSet`, not one matrix: features have different
k (2 vs 3), so their vectors differ in length for the same sequence and cannot
share a matrix. The saving is one FASTA read and one Pool spawn instead of N —
convenience and I/O, not a shared k-mer scan. Expect modest CPU gains, not N×.

### FlexProfile / ProfileSet (`profile.py`)

`FlexProfile` holds `.frame` (rows = sequence IDs, columns = positions),
`.feature`, `.window_size`, `.kmer_len`, and owns output and analysis:

```python
.to_tsv(path)                   # byte-identical to the old package; contract below
.to_frame(tidy=False)
.zscale(axis="column")          # "column" | "row" | "global"
.normalize(method="minmax")
.encode(...)  .heatmap(...)  .metaprofile(...)  .trackplot(...)
```

`ProfileSet` is a dict-like `{feature: FlexProfile}` with `.to_tidy()`.

### Numerics are frozen

The window loop, `sum(w) / len(w)`, and `round(v, 3)` are reproduced exactly. In
particular the mean is **not** replaced with a numpy reduction: last-bit float
differences would flip values at rounding boundaries and break byte-comparison
against the archive.

The one deliberate divergence is unresolvable k-mers: the archive zero-fills via
`.get(subseq, 0)`, the new package masks to `NaN` (see below). This cannot affect
any compared case, because the packaged tables are complete (`4**k` entries) and
every test FASTA is pure ACGT, so no lookup miss occurs on the differential
matrix. The divergence is exercised by its own N-containing fixture instead.

### `to_tsv` byte contract

The Phase 2 gate is byte equality, so the exact serialisation is part of the
contract — not just the numbers. Verified against the checked-in outputs:

- No header, no index, tab-separated; first field is the sequence ID.
- Ragged rows are **NaN-padded to the widest row**, and NaN renders as an empty
  field. A sequence shorter than the window therefore emits its ID followed by
  trailing tabs — literally `'short_seq\t\t\t…\t'` in `w10_mixed.tsv`.
- When *every* row is ID-only, the frame has one column and there is no padding
  at all: the line is bare `'short_seq'` with no trailing tab (`w30.tsv`).
- Values use Python float repr after `round(v, 3)`, so trailing zeros are
  dropped: `0.01`, never `0.010`.

These are consequences of `pd.DataFrame(ragged_lists)` followed by
`to_csv(header=False, index=False, sep="\t")`. The new package reproduces them by
building output the same way, and the differential tests assert on raw bytes
rather than on parsed values — parsing would hide exactly these differences.

## Verification strategy

Because the old package survives, verification is **differential** rather than
against frozen files — a stronger guarantee. Both implementations run on the
same inputs in the same test process and their TSV output is compared byte for
byte.

Matrix: `DNAflexpy/data/test_fasta.fa`, `tests/test_edge_cases.fasta`,
`tests/test_exact_kmer.fasta` × the 10 working features × window sizes
{0, 2, 3, 10, 15, 26, 30}.

Before relying on the harness, it must first reproduce the TSVs already checked
into `tests/test_outputs/` and `tests/edge_case_outputs/` for overlapping cases.
That validates the harness itself.

Expected, allow-listed differences — asserted explicitly so they cannot hide a
real regression:

- `gc`, `freeen`, `mechen`: old emits `None` rows, new emits numbers
- invalid feature names: old emits `None` rows, new raises

Additional coverage: `FeatureTable` validation units, encoding shape/dtype
assertions, plotting smoke tests on the Agg backend, and a `from_bed`
round-trip. The old package's broken `tests/test_utils.py` is not carried over.

## Ambiguous nucleotides — a correctness fix

The old package resolves any k-mer containing `N` or an IUPAC ambiguity code to
`0` via `.get(subseq, 0)`. Zero is a **legitimate value** in these tables (`gc`
AA=0.0, `trx` AA=0), so a gap is indistinguishable from a measurement:

```
ATGCGTACGT -> [0.134, 0.076, -0.077, -0.033, 0.025, ...]
ATGNGTACGT -> [0.134, 0,      0,      0,     0.025, ...]   <- fabricated
```

One `N` invents `k` zero-valued measurements and drags window means toward zero.
Any peak set overlapping an assembly gap hits this, so `from_bed` on a real
genome would silently produce biased profiles.

The new package resolves unknown k-mers to `NaN`, aggregates with nan-aware
means, and records what was masked:

```python
prof.n_masked        # {'chr1:1200-1400': 3}
prof.frame.mean(skipna=True)
```

A window is `NaN` only if *every* k-mer in it is masked; otherwise it averages
the k-mers that resolved. This keeps partial overlap with a gap usable while
never inventing values.

### Serialisation

Masked positions and ragged padding are both `NaN` in the same frame, so a single
`to_csv(na_rep=...)` cannot distinguish them: `na_rep="NA"` would relabel padding
too and break byte equality against the archive, while `na_rep=""` collapses the
distinction. Rather than build rows as pre-formatted strings to control the two
sources separately, **masked positions serialise as empty fields, exactly like
padding.**

Nothing needs to tell them apart by reading the TSV alone: masked counts are
carried by `prof.n_masked` in memory and by `n_masked_positions` in the
provenance sidecar, and every downstream path (encoding, plotting, statistics)
consumes the in-memory frame where `NaN` is unambiguous. Byte equality against
the archive then falls out for free rather than needing a special case.

Verified that none of the three test FASTAs contains a non-ACGT character, so no
differential case is affected either way. A new N-containing fixture covers the
masking behaviour separately, asserting on `n_masked` and the in-memory frame
rather than on TSV bytes.

## Backgrounds and statistics

### Dinucleotide-preserving shuffle

The metaprofile needs a control, and the obvious one is invalid here. Because the
features *are* dinucleotide and trinucleotide lookups, a mononucleotide shuffle
destroys the exact quantity being measured and manufactures significance.

`prof.shuffle(method="dinuc", n=1000, seed=0)` implements Altschul–Erikson
dinucleotide-preserving shuffling, which holds dinucleotide composition fixed so
the comparison tests for flexibility structure **beyond** base composition — the
actual scientific claim. `method="mononuc"` is available but documented as
inappropriate for these features.

Shuffling is seeded and the seed is recorded in the provenance sidecar.

### Per-position testing

```python
res = prof.test(background=bg, method="mannwhitney")   # or "permutation"
prof.metaprofile(background=bg, mark_significant=True)
```

Compares foreground against background at each position, returning effect size,
p-value, and a Benjamini–Hochberg q-value per position. Positions are correlated
by construction — overlapping windows share k-mers — so the docs state plainly
that BH controls FDR across positions but does not undo that dependence, and the
permutation method is recommended when the correlation matters.

The metaprofile shades regions where `q < 0.05`.

## Scale and provenance

### Streaming input

The old implementation calls `list(read_fasta(input_file))`, materialising every
sequence before dispatch. Input is instead consumed in chunks and fed to the pool
with `imap_unordered`, keeping memory flat and independent of input size. Output
row order is restored to input order before serialisation, so results stay
deterministic and byte-comparable.

### Provenance sidecar

Every `to_tsv` writes `<output>.meta.json` alongside it:

```json
{"DNAflexpy_version": "0.3.0", "feature": "DNaseI", "window_size": 10,
 "kmer_len": 3, "lookup_sha256": "9f2a...", "n_sequences": 4821,
 "n_masked_positions": 12, "shuffle_seed": null}
```

The lookup checksum is the point: it makes a result traceable to the exact
feature table that produced it, which matters when a table is edited or a custom
one is supplied.

### Feature redundancy

`ProfileSet.correlate()` returns a feature-by-feature correlation matrix for
dropping redundant features before model fitting. Measured on 5 kb of random
sequence, only 2 of 36 dinucleotide pairs exceed |r| > 0.9 (`freeen`↔`gc` at
−0.944, `gc`↔`bendingStiffness` at +0.916), so this is a convenience for feature
selection rather than a major result.

### Genome-browser export

`to_bigwig(path, genome=...)` and `to_bed(path)` write profiles as tracks for
IGV/JBrowse. Only meaningful for profiles carrying genomic coordinates, so both
raise on FASTA- or table-derived input, which has no coordinates to write.
Requires the optional `DNAflexpy[bigwig]` extra (pyBigWig).

## ML encoding (`encode.py`)

Modelled on DNAshapeR's `encodeSeqShape`, which is the function that makes that
package useful for statistical learning.

```python
X = prof.encode(["1-mer", "1-flex", "2-DNaseI"], normalize=True)
```

- `1-mer` / `2-mer` / `3-mer` — one-hot sequence encoding, 4/16/64 binary columns
  per position (`encodeKMerSeq`)
- `1-flex` — first-order flexibility values, one column per position
- `n-<feature>` — nth-order terms: products of the same feature at adjacent
  positions, capturing interactions a linear model would otherwise miss
  (`encodeNstOrderShape`)
- `normalize=True` applies min-max scaling (`normalize`)

Returns a 2-D array with named columns, directly consumable by scikit-learn.
Requires equal-length sequences; raises otherwise.

## Plotting (`plotting.py`)

Separate module with lazy matplotlib import, so core profiling stays testable
without a plotting stack. matplotlib is an optional extra, `DNAflexpy[plot]`.

Ragged input **raises** rather than truncating: position-wise comparison is only
meaningful for aligned sequences. The error message names the offending lengths
and points at `from_bed(width=...)` for producing fixed-width input.

### Heatmap

Rows = sequences, columns = positions. Strictly one feature per figure —
`stiffness` (AA=1820) and `DNaseI` (AA=-0.28) cannot share a colour scale.
Diverging colormap centred at 0 for signed features.

- `nbins=` aggregates positions into equal-width bins so large matrices stay
  renderable (`heatShape` nBins)
- `order_rows="cv"` sorts sequences by coefficient of variation, surfacing the
  most variable first instead of input order (`heatShape` ordRow)
- `zscale="column"` — the default

### Metaprofile — and a correction to the approved z-scale choice

Column-wise z-scaling was chosen for both plots. It is correct for the heatmap,
but **degenerate for a single-set metaprofile**: z-scoring each column across
sequences and then averaging down that column yields exactly 0 at every
position, by construction. Verified numerically — the line plot would be flat at
zero.

Pooling foreground and background before column-wise scaling does **not** rescue
it. Verified numerically: after standardizing the pooled matrix per column the
pooled column mean is 0, so `n_fg·mean_fg + n_bg·mean_bg = 0` and therefore
`mean_bg = −(n_fg/n_bg)·mean_fg`. The two curves correlate at exactly −1.0 — the
background is a rescaled mirror of the foreground and carries no independent
information whatsoever.

Resolution, keeping column-wise scaling where it is meaningful:

- **With a background set** (the DNAshapeR `plotShape(background=)` case):
  standardize each column against the **background's** mean and standard
  deviation — a signal-vs-control z-score. The background then sits flat at ~0
  and the foreground's departure from it is a real, independently interpretable
  quantity, in units of background standard deviations. This is the default
  whenever a background is supplied.
- **Without a background set**: default to `zscale="global"`, preserving
  cross-position structure. `zscale="column"` raises with an explanation rather
  than silently drawing a flat line.

Binning is applied **before** z-scoring, so statistics are computed on the values
actually plotted. The two orders give different figures; this one is fixed.

### Trackplot

Single-sequence view with all requested features stacked on a shared x-axis
(`trackShape`).

## Input formats (`io.py`)

Four ways in — bare string, FASTA, labelled table, BED. All converge on the same
`FlexProfile`, so every downstream path is identical regardless of source.

### Labelled TSV — sequence + value

For labelled datasets: first column is the sequence, second a numeric value such
as binding affinity. This is the natural ML input, carrying `X` and `y` in one
file.

```python
prof = FlexProfiler("DNaseI", window_size=10).profile_table(
    "affinity.tsv", seq_col=0, value_col=1,
)
X, y = prof.encode(["1-mer", "1-flex"]), prof.y
```

- Header row auto-detected (first row's `value_col` failing to parse as a float
  means it is a header); forceable with `header=True|False`.
- Columns may be named instead of positional: `seq_col="sequence"`.
- Rows without an ID column get generated IDs `seq_0, seq_1, …`; an existing one
  is named via `id_col=`.
- Non-numeric or missing values **raise**, naming the offending row — silently
  dropping labelled rows would corrupt a training set.
- Tab by default; `sep=` allows CSV.

`FlexProfile.y` holds the aligned value vector (`None` for FASTA input) and
unlocks two things beyond supplying regression targets:

- `heatmap(order_rows="y")` — order sequences by affinity rather than by
  coefficient of variation, making a flexibility gradient across the affinity
  range directly visible.
- `metaprofile(groupby=("y", threshold))` — split into high/low groups and use
  the low group as the background set. That is the signal-vs-control comparison
  described above, with a control already matched to the assay.

### BED + reference genome

```python
FlexProfiler("DNaseI", window_size=10).from_bed("peaks.bed", genome="TAIR10.fa", width=200)
```

Extracts fixed-width, centre-anchored sequences from a reference FASTA via
pyfaidx (optional extra, `DNAflexpy[bed]`), mirroring `getFasta`. Strand-aware:
`-` strand intervals are reverse-complemented.

Fixed width gives alignment by construction **except at contig boundaries**,
where a centred window can run past the start or end of a chromosome. Policy is
explicit via `on_edge=`, defaulting to `"drop"`:

- `"drop"` — discard the interval, reporting how many were dropped
- `"error"` — raise, naming the offending intervals
- `"pad"` — pad with `N`, which the lookup zero-fills; permitted but warned
  about, since padded positions bias profiles toward zero

Without this, "aligned by construction" would quietly be false for any peak set
near a chromosome end.

## CLI

New `DNAflexpy` command with subcommands. The old `DNAflexpy` command is untouched,
so `scripts/plot_profiles.py`, which shells out to `python -m DNAflexpy.cli`,
keeps working unchanged.

```bash
DNAflexpy profile in.fa --feature DNaseI trx --window-size 10 --outfile out.tsv
DNAflexpy profile affinity.tsv --seq-col 0 --value-col 1 --feature DNaseI
DNAflexpy profile peaks.bed --genome TAIR10.fa --width 200 --feature DNaseI
DNAflexpy profile --seq ATGCGTACGT --feature DNaseI       # bare string -> stdout
DNAflexpy encode in.fa --features 1-mer 1-flex --out X.npz
DNAflexpy plot heatmap out.tsv --nbins 20 --order-rows cv --out hm.png
DNAflexpy plot meta out.tsv --background bg.tsv --out meta.png
```

Input format is inferred from extension (`.fa`/`.fasta`/`.fna` → FASTA, `.bed` →
BED, `.tsv`/`.csv` → labelled table) and overridable with `--format`.

## Dependencies

- Required: `pandas`, `pyyaml`, `numpy`
- `DNAflexpy[plot]`: `matplotlib`
- `DNAflexpy[bed]`: `pyfaidx`
- `DNAflexpy[stats]`: `scipy` (per-position testing and FDR)
- `DNAflexpy[bigwig]`: `pyBigWig`
- `DNAflexpy[all]`: everything above
- Dev: `pytest`

The container installs `[all]`; the shuffle is implemented in numpy and needs no
extra.

## Container

A single full-analysis image containing every extra, so all subcommands work
without a "which feature is missing" question. The Dockerfile lives in the repo
and is built locally; **nothing is published to any registry**. Publishing can be
added later without changing the image.

```
Docker/
  Dockerfile
  README.md      build + run + Apptainer instructions
  .dockerignore
```

Multi-stage: a builder stage compiles wheels, the runtime stage copies only the
installed packages, keeping the image near ~400 MB rather than carrying build
toolchains.

`DNAflexpy` is **not** published to any index, so the image installs it from the
build context — never by name from PyPI:

```dockerfile
FROM python:3.12-slim AS builder
COPY . /src
RUN pip wheel --no-cache-dir --wheel-dir=/wheels "/src[plot,bed]"

FROM python:3.12-slim
COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links=/wheels \
      DNAflexpy[plot,bed] && rm -rf /wheels
ENV MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mpl \
    PYTHONDONTWRITEBYTECODE=1
WORKDIR /data
ENTRYPOINT ["DNAflexpy"]
```

### Apptainer/Singularity compatibility

The image must run unprivileged on an HPC cluster via
`apptainer pull docker://...`. This is cheap to honour from the start and
painful to retrofit, so it constrains the build:

- **No `USER` directive and no assumptions about UID or `$HOME`.** Apptainer runs
  the container as the invoking user, who may have no passwd entry and no
  writable home inside the image.
- **All writes go to the bind mount** (`/data`). Nothing is written to the image
  filesystem at runtime, which is read-only under Apptainer.
- `MPLCONFIGDIR=/tmp/mpl` — matplotlib otherwise tries to write a config
  directory under `$HOME` and warns or fails when that is unwritable.
- `PYTHONDONTWRITEBYTECODE=1` — prevents `.pyc` writes into a read-only layer.
- **No runtime `chown`/`chmod`**; packages install to system site-packages, not
  user site, so they remain importable as an arbitrary UID.

```bash
# Docker, locally
docker build -t dnaflexpy:dev -f Docker/Dockerfile .
docker run --rm -v "$PWD":/data dnaflexpy:dev profile in.fa --feature DNaseI --outfile out.tsv

# Apptainer, on a cluster, no root
apptainer build dnaflexpy.sif docker-daemon://dnaflexpy:dev
apptainer exec --bind "$PWD":/data dnaflexpy.sif DNAflexpy profile /data/in.fa --feature DNaseI
```

A container smoke test runs the CLI end-to-end (profile → encode → plot) on a
packaged test FASTA, so a broken image is caught at build time rather than on the
cluster.

## Implementation phases

Ordered so the performance fix lands early and does not wait on plotting.

0. **Archive** — `git mv DNAflexpy rxv/DNAflexpy`, add `rxv/__init__.py`, fix the
   forced `files("rxv.DNAflexpy.data")` import path, migrate packaging to
   `pyproject.toml`, update `scripts/plot_profiles.py`. Gate: the archived
   package still runs and reproduces the checked-in TSVs. **This must land before
   any new code, since every later phase is verified against it.**
1. **Scaffold + FeatureTable** — package skeleton, YAML loading and validation,
   k inference, unit tests.
2. **Profiling core** — `FlexProfiler`, `FlexProfile`, TSV output, Pool
   initializer, `NaN` masking for ambiguous nucleotides. Gate: differential suite
   byte-matches the archived package on all allowed cases.
3. **Inputs** — `ProfileSet` for multi-feature runs, user-supplied lookup tables,
   bare-string and list/dict entry points, the memoised module-level one-liner,
   labelled TSV input (`profile_table`, `FlexProfile.y`), and BED/genome input.
4. **Encoding** — `encode.py` and normalisation.
5. **Plotting** — heatmap, metaprofile with background, trackplot.
6. **Backgrounds and statistics** — dinucleotide shuffle, per-position testing
   with BH correction, significance shading on the metaprofile.
7. **Scale and provenance** — streaming input, provenance sidecar,
   `ProfileSet.correlate()`, bigWig/BED export.
8. **CLI + docs** — `DNAflexpy` command, populate the empty `docs/api_reference.md`,
   add it to the mkdocs nav, update `CLAUDE.md` to describe both packages.
9. **Container** — Dockerfile, Apptainer constraints, smoke test, `Docker/README.md`.

## Out of scope

- Modifying, deleting, or archiving `DNAflexpy/` — the user will archive it.
- Updating `Notebooks/*.ipynb` — they continue to work against the old package.
- Publishing the image to any container registry, and publishing the package to
  PyPI. Both are outward-facing and deferred to an explicit user decision.
- Methylation support (`convertMethFile`). The lookup tables contain only ACGT
  k-mers with no methylated states, so there is nothing to look up. Adding it
  would require new experimental parameter tables, not a code change.
