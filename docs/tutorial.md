# Tutorial

A complete analysis, start to finish, using data that ships with the package.
Nothing to download. Every block on this page is executed as part of the test
suite, so it runs as written.

The question we will answer: **do these promoters have a flexible region that
control sequences do not?**

## The example data

Five files come with the package:

```python
from DNAflexpy import describe_examples, example_path

for name, what in sorted(describe_examples().items()):
    print(f"{name:16} {what}")
```

```
affinity.tsv     the same promoters with a binding score, for model fitting (has a header)
control.fa       12 matched controls: same length and base composition, no positioned element
genome.fa        a two-contig toy genome, chr1 (2000 bp) and chr2 (1200 bp)
peaks.bed        5 intervals into that genome, named and stranded
promoters.fa     12 synthetic promoters, 120 bp each, with an AT-rich element at 40-60
```

`example_path(name)` gives you the full path, wherever the package is
installed. It is the equivalent of DNAshapeR's `system.file("extdata", ...)`.

The sequences are synthetic, but built to behave like real ones. Each promoter
carries an AT-rich stretch at bases 40–60 — what a TATA-box region looks like
to a flexibility profile — and `control.fa` holds sequences of the same length
and the same base composition with no positioned element. So any difference we
find is about **where** the bases are, not how many of each there are.

## Step 1 — profile the sequences

```python
from DNAflexpy import FlexProfiler, example_path

profiler = FlexProfiler(feature="DNaseI", window_size=10)
promoters = profiler.profile_fasta(example_path("promoters.fa"))

print(promoters.frame.shape)
print(promoters.seqids[:3])
```

```
(12, 111)
['promoter_1', 'promoter_2', 'promoter_3']
```

Twelve sequences, 111 positions each. The sequences are 120 bases and the
window is 10, so there are `120 - 10 + 1 = 111` window positions.

`DNaseI` is one of thirteen tables — see [Feature tables](features.md). Build
the profiler once and reuse it: it reads the table at construction, which is
why reusing it is about 130× faster than a fresh one per sequence.

### What a profile is

Each row is one sequence; each column is one position. The value is the mean
of the k-mer parameter values inside the window starting there.

```python
print(promoters.frame.iloc[0, :6].round(3).to_dict())
```

```
{1: -0.03, 2: 0.005, 3: 0.026, 4: 0.008, 5: -0.004, 6: 0.011}
```

Set `window_size=0` to skip the averaging and get one value per k-mer instead.

## Step 2 — look at it

```python
import matplotlib
matplotlib.use("Agg")

ax = promoters.heatmap(nbins=30)
ax.figure.savefig("tutorial_heatmap.png", dpi=150, bbox_inches="tight")
```

Rows are sequences, columns are positions. `nbins=30` averages the 111
positions into 30 bins so the picture stays readable; without it you get one
column per position.

Rows come back sorted with the most variable sequence first. Pass
`order_rows="input"` to keep file order.

One feature per figure, always. `stiffness` runs to 5500 and `DNaseI` to
0.194 — they cannot share a colour scale, so the function will not try.

## Step 3 — compare against a control

A single heatmap tells you what the sequences look like. It does not tell you
whether that is unusual. For that you need something to compare against.

```python
controls = profiler.profile_fasta(example_path("control.fa"))

ax = promoters.metaprofile(background=controls)
ax.figure.savefig("tutorial_meta.png", dpi=150, bbox_inches="tight")
```

The y-axis is now **background standard deviations**. Each position is
standardised against the control set's own mean and spread at that position,
so the grey control line sits flat at zero and the promoter line shows how far
it departs.

Where does it depart most?

```python
import numpy as np
from DNAflexpy.plotting import _matrix

fg, bg = _matrix(promoters), _matrix(controls)
z = (fg.mean(axis=0) - bg.mean(axis=0)) / bg.std(axis=0)

peak = int(np.argmax(np.abs(z)))
print(f"largest departure at position {peak + 1}: {z[peak]:+.2f} SDs")
```

```
largest departure at position 49: -1.88 SDs
```

Position 49, which is inside the AT-rich element we planted at 40–60. The
window is 10 bases wide, so an element at 40–60 shows up between positions 31
and 60. The sign is negative because AT-rich DNA has lower `DNaseI` values in
this table.

That is the analysis working: the difference is positional, since the two sets
have matched base composition.

!!! note "One option is refused here"
    `metaprofile(zscale="column")` **without** a background raises rather than
    drawing. Standardising each position across sequences and then averaging
    down that same position cancels to exactly zero — the line would be flat by
    construction, at every position. Rather than hand you a meaningless chart,
    it explains the problem. With a background, as above, the comparison is
    well defined.

## Step 4 — several features at once

```python
several = FlexProfiler(["DNaseI", "gc", "stiffness"], window_size=10)
result = several.profile_fasta(example_path("promoters.fa"))

print(sorted(result))
fig = result.trackplot(seqid="promoter_1")
fig.savefig("tutorial_tracks.png", dpi=150, bbox_inches="tight")
```

```
['DNaseI', 'gc', 'stiffness']
```

You get a `ProfileSet` — a dict of `{feature: profile}`. They are kept apart
because different tables use different k-mer lengths, so their vectors differ
in length for the same sequence.

`trackplot` stacks one sequence's features on a shared x-axis, each with its
own y-axis.

## Step 5 — fit a model

`affinity.tsv` holds the same promoters with a binding score. One file gives
you both the predictors and the target.

```python
labelled = profiler.profile_table(
    example_path("affinity.tsv"),
    seq_col="sequence", value_col="binding_score", id_col="name",
    header=True,
)
print(len(labelled.y), labelled.y[:3])
```

```
12 [7.801, 7.945, 7.987]
```

Turn it into a design matrix:

```python
fm = labelled.encode(["1-mer", "1-DNaseI"])
print(fm.shape, len(fm.y))
print(fm.columns[0], "...", fm.columns[-1])
```

```
(12, 591) 12
seq.1mer.p1.A ... DNaseI.w10.o1.p111
```

`fm.X` is a NumPy array and `fm.y` is the target list, so they go straight
into scikit-learn:

```python
# from sklearn.linear_model import Ridge
# from sklearn.model_selection import cross_val_score
#
# model = Ridge()
# cross_val_score(model, fm.X, fm.y, cv=3)
```

`1-mer` gives four binary columns per base. `1-DNaseI` gives the profile
values. Add `2-DNaseI` for products of adjacent positions — the interaction
terms a linear model would otherwise miss.

!!! warning "Normalisation and train/test splits"
    `normalize=True` (the default) min-max scales each feature block using
    **this dataset's** range. If you are going to split into train and test
    sets, that leaks. Pass `normalize=False` and scale inside a scikit-learn
    pipeline instead.

## Step 6 — start from genomic intervals

Real peak sets come as BED files. `from_bed` reads them against a reference
genome. Needs the `bed` extra: `pip install 'DNAflexpy[bed]'`.

First copy the example data somewhere you can write:

```python
import pathlib
import shutil

from DNAflexpy import example_files, example_path

work = pathlib.Path("dnaflexpy_examples")
work.mkdir(exist_ok=True)
for name in example_files():
    shutil.copy(example_path(name), work / name)
```

!!! warning "Why copy, rather than read the packaged genome directly"
    `pyfaidx` writes an index file (`genome.fa.fai`) **next to the genome** the
    first time it reads it. Point it at the installed package and it will try
    to write into your `site-packages`, which fails on a read-only install —
    a container, a system-wide install, or an HPC module. This applies to any
    genome, not just this example: make sure the directory holding your FASTA
    is writable, or pre-build the index with `samtools faidx`.

```python
peaks = profiler.from_bed(
    work / "peaks.bed",
    genome=work / "genome.fa",
    width=120,
)
print(peaks.seqids)
print(peaks.frame.shape)
```

```
['peak_1', 'peak_2', 'peak_3', 'peak_4', 'peak_5']
(5, 111)
```

`width=120` re-centres every interval on its midpoint and cuts it to exactly
120 bases, so all rows are the same length. That is what makes position-wise
comparison meaningful — and it is required by both `encode` and the plots.

Intervals on the minus strand are reverse-complemented automatically. An
interval whose centred window runs off the end of a contig is dropped with a
warning; `on_edge="pad"` or `on_edge="error"` change that.

## The same thing from the command line

Every step above has a command-line equivalent:

```bash
# profile
DNAflexpy profile promoters.fa --feature DNaseI --window-size 10

# several features -> one file each
DNAflexpy profile promoters.fa --feature DNaseI gc stiffness

# plots
DNAflexpy plot heatmap promoters.fa --nbins 30 --out heatmap.png
DNAflexpy plot meta promoters.fa --background control.fa --out meta.png
DNAflexpy plot track promoters.fa --feature DNaseI gc stiffness --out tracks.png

# design matrix
DNAflexpy encode affinity.tsv --header --seq-col sequence \
  --value-col binding_score --features 1-mer 1-DNaseI --out X.npz

# from BED
DNAflexpy profile peaks.bed --genome genome.fa --width 120
```

Run them from the `dnaflexpy_examples` directory created in step 6:

```python
print(sorted(p.name for p in work.iterdir() if p.suffix != ".fai"))
```

```
['affinity.tsv', 'control.fa', 'genome.fa', 'peaks.bed', 'promoters.fa']
```

## Where to go next

- [Usage](usage.md) — every option, every input format, in detail
- [Feature tables](features.md) — what each of the thirteen measures
- [API reference](api_reference.md) — every class and function
- [FAQ](faq.md) — masked values, window sizes, ambiguous bases
