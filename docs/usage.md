# Uage
---

## CLI
This section includes examples of how to use DNAflexpy both as a command-line tool and as a library in Python.

Once DNAflexpy is installed, you can run it from the command line:

```bash

DNAflexpy "<path/to/nt/fasta/file>" \
    --window-size 10 \
    --feature "<DNaseI>/<NPP>" \
    --outfile "<path/to/output.tsv>"

```

## Run by importing as library

```py

from DNAflexpy.core import DNAflexpy

# Example usage
df = DNAflexpy(
    input_file="input.fasta",
    window_size=10,
    kmer_len=3,
    feature="DNaseI",
    threads=4
)

# Display the results
print(df.head())

```