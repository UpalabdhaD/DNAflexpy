# FAQ

## Why do I have to say `header=True` or `header=False`?

Because guessing is not safe here. DNA is written in letters that also spell
words: `dna`, `rna`, `tag`, `cat` and `data` are all made only of nucleotide
letters. So a header row can look exactly like a sequence, and any rule that
guesses will sometimes read your column names as a data row and quietly add a
made-up example to your training set.

Two guessing rules were tried during development. Each fixed one case and broke
another. Being asked once is cheaper than finding a fabricated row later.

## Why do sequences with `N` give blank values instead of zero?

Because zero is a real measurement in these tables. `gc` scores `AA` as `0.0`,
and `trx` scores `AA` as `0`. If a gap in your sequence were also scored zero,
you could not tell a missing value from a measured one.

So any k-mer that cannot be looked up is masked instead. The original package
scored them zero, which is the one behaviour this rewrite deliberately changed.

## `n_masked` says 0, but I got a warning about odd letters. Why?

`n_masked` counts a position as masked only when *none* of the k-mers in that
window could be looked up. If a window has some usable k-mers, it is simply
averaged over fewer of them, and does not show up as masked.

At `window_size=0` you see the true count. Above that, `n_masked` under-reports.
The warning is the reliable signal.

## Where is the command-line tool?

The rewritten package does not have one yet. The original still works:

```bash
python -m rxv.DNAflexpy.cli sequences.fa --window-size 10 --feature DNaseI --outfile out.tsv
```

## What is the `rxv/` folder?

The original package, kept exactly as it was. It is frozen: its logic is never
changed.

It is there because the rewrite is checked against it. 230 tests run both
versions on the same input and compare the output **byte for byte**. If a
change alters any number, those tests fail.

That is stricter than comparing values. Six features hold whole numbers, and
writing `18.0` where the original wrote `18` is a real difference that
`==` would miss, because `18 == 18.0` is true in Python.

## Why is my profile shorter than my sequence?

A window of `N` bases starting at every position gives `len(sequence) - N + 1`
windows, not `len(sequence)`. With `window_size=0` you get one value per k-mer,
so `len(sequence) - k + 1`.

If the window is longer than the sequence, none fit, and you get the sequence
name with no values.

## Can I use my own feature values?

Yes:

```python
from DNAflexpy import FlexProfiler, FeatureTable

p = FlexProfiler(feature="mine", lookup={"mine": {"AA": 1.0, "AT": 2.0,
                                                  "TA": 3.0, "TT": 4.0}})
table = FeatureTable.from_yaml("my_features.yaml")
```

The k-mer length is read from your keys, so there is no second place to keep in
step. Mixed lengths inside one feature are rejected.

## Why is BED input a separate install?

It needs `pyfaidx` to read sequences out of a reference genome by position.
Most users never do that, so it is optional: `pip install -e '.[bed]'`. The
rest of the package works without it.
