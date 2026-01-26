# Window Size Behavior Checks

This note documents manual checks for window sizing after fixing the window-average calculation.

## Test input

```
>short_seq
ATCGAT
>longer_seq
ATCGATCGATCG
```

## Commands

```bash
python -m DNAflexpy.cli /tmp/window_size_check.fa --window-size 0 --feature DNaseI --outfile /tmp/window_w0.tsv
python -m DNAflexpy.cli /tmp/window_size_check.fa --window-size 3 --feature DNaseI --outfile /tmp/window_w3.tsv
python -m DNAflexpy.cli /tmp/window_size_check.fa --window-size 10 --feature DNaseI --outfile /tmp/window_w10.tsv
python -m DNAflexpy.cli /tmp/window_size_check.fa --window-size 2 --feature DNaseI --outfile /tmp/window_w2.tsv
```

## Observed output (DNaseI)

```
== /tmp/window_w0.tsv ==
short_seq    -0.11    -0.003   -0.003   -0.11
longer_seq   -0.11    -0.003   -0.003   -0.11   -0.11   -0.003   -0.003   -0.11   -0.11   -0.003

== /tmp/window_w3.tsv ==
short_seq    -0.11    -0.003   -0.003   -0.11
longer_seq   -0.11    -0.003   -0.003   -0.11   -0.11   -0.003   -0.003   -0.11   -0.11   -0.003

== /tmp/window_w10.tsv ==
short_seq
longer_seq   -0.057   -0.057   -0.057

== /tmp/window_w2.tsv ==
short_seq    0.0      0.0      0.0      0.0      0.0
longer_seq   0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0
```

## Expected behavior confirmed

- `--window-size 0` returns per-kmer values for the full sequence.
- `--window-size 3` matches the per-kmer values for `DNaseI` because the window contains exactly one 3-mer.
- `--window-size 10` generates windows for `longer_seq` and none for `short_seq` (window size > sequence length).
- `--window-size 2` yields 0.0 averages because the window is smaller than the 3-mer feature size.
