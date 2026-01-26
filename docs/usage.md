# Usage
---

## CLI
This section includes examples of how to use DNAflexpy both as a command-line tool and as a library in Python.

Once DNAflexpy is installed, you can run it from the command line:

```bash

DNAflexpy "<path/to/nt/fasta/file>" \
    --window-size 10 \
    --feature "<DNaseI/NPP/twistDisp/trx/stiffness>" \
    --outfile "<path/to/output.tsv>"

```

To print the BibTeX citation:

```bash
DNAflexpy --citation
```

## Run by importing as library

```py

from DNAflexpy.core import DNAflexpy

# Example usage
df = DNAflexpy(
    input_file="input.fasta",
    window_size=10,
    feature="DNaseI",
    threads=4
)

# Display the results
print(df.head())

```

## Feature options

The `--feature` option (or `feature` in the API) accepts the following values:

- DNaseI
- NPP
- twistDisp
- trx
- stiffness

## Window size behavior

- `--window-size 0` disables windowing and returns per-kmer feature values for the full sequence.
- `--window-size N` generates overlapping windows of length N (step size 1), computes per-kmer
  feature values within each window, and reports the mean value per window.
- If `N` equals the feature k-mer length (e.g., `N=3` for `DNaseI`), each window contains a single
  k-mer value, so outputs will match the per-kmer values.
- If `N` is larger than the sequence length, no windows are generated and only the sequence ID
  is returned for that record.
- If `N` is smaller than the feature k-mer length, the window averages are 0.0 because no k-mers
  fit inside the window.

## Repository

Git URL: https://github.com/upalabdhaD/DNAflexpy.git
