# dnaflex — Greenfield Rewrite Design

**Date:** 2026-08-12
**Status:** Awaiting user review
**Supersedes:** nothing — `DNAflexpy/` remains untouched and installable

## Goal

Rewrite DNAflexpy as a new package, `dnaflex`, built from scratch. The existing
`DNAflexpy/` package is left byte-for-byte untouched and will be archived by the
user later. The rewrite fixes a measured performance defect, replaces the
function-based API with classes, and adds ML feature encoding, plotting, and
genomic input modelled on the R/Bioconductor package DNAshapeR.

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

Therefore `dnaflex` infers k from the table and has no hardcoded map. This single
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

```
DNAflexpy/            <- UNTOUCHED. Still installs, still exposes `DNAflexpy` CLI.
dnaflex/              <- NEW
  __init__.py           FlexProfiler, FlexProfile, ProfileSet, __version__, __all__
  lookup.py             FeatureTable: load + validate + infer k
  io.py                 read_fasta, from_bed
  core.py               FlexProfiler: profiling engine
  profile.py            FlexProfile / ProfileSet: results, TSV, DataFrame
  encode.py             ML feature matrices
  plotting.py           heatmap / metaprofile / trackplot (lazy matplotlib)
  cli.py                `dnaflex` command
  data/lookupNEW.yaml   own copy — no shared state with the old package
```

Both packages install side by side; import names never collide. When the user
archives the old package, `dnaflex` is renamed to `DNAflexpy` in one commit.

Version lives only in `dnaflex/__init__.py.__version__`, read by packaging
metadata — avoiding the setup.py/requirements.txt drift in the old package.

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
- **warn**, not fail, when a table has fewer than `4**k` entries, because lookup
  zero-fills missing k-mers via `.get(subseq, 0)` and would otherwise silently
  bias results toward zero

### FlexProfiler (`core.py`)

```python
FlexProfiler(feature, window_size=10, lookup=None)
    feature: str | list[str]
    lookup:  None | path | dict | FeatureTable
.profile(seqid, seq) -> list                      # fast path, ~0.1 ms
.profile_fasta(path, threads=None) -> FlexProfile | ProfileSet
.from_bed(bed, genome, width=None) -> FlexProfile | ProfileSet
```

Features are validated against the table in `__init__`, so a typo fails
immediately rather than after a full FASTA pass.

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
.to_tsv(path)                   # identical format to the old package
.to_frame(tidy=False)
.zscale(axis="column")          # "column" | "row" | "global"
.normalize(method="minmax")
.encode(...)  .heatmap(...)  .metaprofile(...)  .trackplot(...)
```

`ProfileSet` is a dict-like `{feature: FlexProfile}` with `.to_tidy()`.

### Numerics are frozen

The window loop, `sum(w) / len(w)`, `round(v, 3)`, and `.get(subseq, 0)`
zero-filling are reproduced exactly. In particular the mean is **not** replaced
with a numpy reduction: last-bit float differences would flip values at rounding
boundaries and break byte-comparison against the old output.

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
without a plotting stack. matplotlib is an optional extra, `dnaflex[plot]`.

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

Resolution, without abandoning the choice:

- **With a background set** (the DNAshapeR `plotShape(background=)` case):
  z-score the **pooled** foreground+background matrix column-wise, then plot the
  two group means separately. Column-wise scaling is exactly right here, and the
  divergence between curves is the finding.
- **Without a background set**: default to `zscale="global"`, which preserves
  cross-position structure. `zscale="column"` raises with an explanation rather
  than silently drawing a flat line.

### Trackplot

Single-sequence view with all requested features stacked on a shared x-axis
(`trackShape`).

## Genomic input (`io.py`)

```python
FlexProfiler("DNaseI", window_size=10).from_bed("peaks.bed", genome="TAIR10.fa", width=200)
```

Extracts fixed-width, centre-anchored sequences from a reference FASTA via
pyfaidx (optional extra, `dnaflex[bed]`), mirroring `getFasta`. Fixed width means
output is aligned by construction, satisfying the equal-length requirement of the
plotting and encoding paths. Strand-aware: `-` strand intervals are
reverse-complemented.

## CLI

New `dnaflex` command with subcommands. The old `DNAflexpy` command is untouched,
so `scripts/plot_profiles.py`, which shells out to `python -m DNAflexpy.cli`,
keeps working unchanged.

```bash
dnaflex profile in.fa --feature DNaseI trx --window-size 10 --outfile out.tsv
dnaflex profile peaks.bed --genome TAIR10.fa --width 200 --feature DNaseI
dnaflex encode in.fa --features 1-mer 1-flex --out X.npz
dnaflex plot heatmap out.tsv --nbins 20 --order-rows cv --out hm.png
dnaflex plot meta out.tsv --background bg.tsv --out meta.png
```

## Dependencies

- Required: `pandas`, `pyyaml`, `numpy`
- `dnaflex[plot]`: `matplotlib`
- `dnaflex[bed]`: `pyfaidx`
- Dev: `pytest`

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

```dockerfile
FROM python:3.12-slim AS builder
# build wheels for pandas/numpy/matplotlib/pyfaidx + dnaflex

FROM python:3.12-slim
COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links=/wheels \
      dnaflex[plot,bed] && rm -rf /wheels
ENV MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mpl \
    PYTHONDONTWRITEBYTECODE=1
WORKDIR /data
ENTRYPOINT ["dnaflex"]
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
docker build -t dnaflex:dev -f Docker/Dockerfile .
docker run --rm -v "$PWD":/data dnaflex:dev profile in.fa --feature DNaseI --outfile out.tsv

# Apptainer, on a cluster, no root
apptainer build dnaflex.sif docker-daemon://dnaflex:dev
apptainer exec --bind "$PWD":/data dnaflex.sif dnaflex profile /data/in.fa --feature DNaseI
```

A container smoke test runs the CLI end-to-end (profile → encode → plot) on a
packaged test FASTA, so a broken image is caught at build time rather than on the
cluster.

## Implementation phases

Ordered so the performance fix lands early and does not wait on plotting.

1. **Scaffold + FeatureTable** — package skeleton, YAML loading and validation,
   k inference, unit tests.
2. **Profiling core** — `FlexProfiler`, `FlexProfile`, TSV output, Pool
   initializer. Gate: differential suite byte-matches the old package on all
   allowed cases.
3. **Multi-feature, custom tables, `from_bed`** — `ProfileSet`, user-supplied
   lookups, BED/genome input.
4. **Encoding** — `encode.py` and normalisation.
5. **Plotting** — heatmap, metaprofile with background, trackplot.
6. **CLI + docs** — `dnaflex` command, populate the empty `docs/api_reference.md`,
   add it to the mkdocs nav, update `CLAUDE.md` to describe both packages.
7. **Container** — Dockerfile, Apptainer constraints, smoke test, `Docker/README.md`.

## Out of scope

- Modifying, deleting, or archiving `DNAflexpy/` — the user will archive it.
- Updating `Notebooks/*.ipynb` — they continue to work against the old package.
- Publishing the image to any container registry, and publishing the package to
  PyPI. Both are outward-facing and deferred to an explicit user decision.
- Methylation support (`convertMethFile`). The lookup tables contain only ACGT
  k-mers with no methylated states, so there is nothing to look up. Adding it
  would require new experimental parameter tables, not a code change.
