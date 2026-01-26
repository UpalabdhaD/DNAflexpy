# Window Size Checks (30 bp sequences)

This note documents manual checks for window sizing on two 30 bp sequences.

## Test input

```
>long_seq_1
ATCGATCGATCGATCGATCGATCGATCGAT
>long_seq_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
```

## Commands

```bash
python -m DNAflexpy.cli /tmp/window_size_long_check.fa --window-size 0 --feature DNaseI --outfile /tmp/window_long_w0.tsv
python -m DNAflexpy.cli /tmp/window_size_long_check.fa --window-size 10 --feature DNaseI --outfile /tmp/window_long_w10.tsv
python -m DNAflexpy.cli /tmp/window_size_long_check.fa --window-size 15 --feature DNaseI --outfile /tmp/window_long_w15.tsv
```

## Observed output summary (DNaseI)

The number of values per sequence matches the expected window counts:

- `window_size 0`: per-kmer values for full sequence (`30 - 3 + 1 = 28` values)
- `window_size 10`: overlapping windows (`30 - 10 + 1 = 21` values)
- `window_size 15`: overlapping windows (`30 - 15 + 1 = 16` values)

Sample values (first 5 per sequence):

```
== /tmp/window_long_w0.tsv ==
long_seq_1 values: 28 first5: -0.11, -0.003, -0.003, -0.11, -0.11
long_seq_2 values: 28 first5: 0.017, 0.09, 0.09, 0.017, 0.017

== /tmp/window_long_w10.tsv ==
long_seq_1 values: 21 first5: -0.057, -0.057, -0.057, -0.057, -0.057
long_seq_2 values: 21 first5: 0.053, 0.053, 0.053, 0.053, 0.053

== /tmp/window_long_w15.tsv ==
long_seq_1 values: 16 first5: -0.061, -0.052, -0.052, -0.061, -0.061
long_seq_2 values: 16 first5: 0.051, 0.056, 0.056, 0.051, 0.051
```

## Expected behavior confirmed

- Window size alters the number of windows and output length.
- Window size changes the averaged values (`w10` and `w15` outputs differ).
