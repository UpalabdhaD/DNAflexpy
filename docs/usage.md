# Uage
---

## CLI
This section includes examples of how to use DNAflexpy both as a command-line tool and as a library in Python.

Once DNAflexpy is installed, you can run it from the command line:

```bash
feature-profiler input.fasta --window-size 10 --kmer-len 3 --feature DNaseI --outfile output.tsv
```

## Run by importing as library

```py

from DNAflexpy.core import calculate_features

# Example usage
df = calculate_bendability(
    input_file="input.fasta",
    window_size=10,
    kmer_len=3,
    feature="DNaseI",
    feature_file="data/feature_lookup.yaml",
    threads=4
)

# Display the results
print(df.head())

```