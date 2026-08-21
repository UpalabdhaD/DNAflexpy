# DNAflexpy

**DNAflexpy** is a python package designed to calculate bendability of DNA sequences using trinucleotide feature values, helping researchers analyze and understand sequence properties.

This documentation provides an overview of installation, usage, and API references to help you get started quickly.

Repository: https://github.com/upalabdhaD/DNAflexpy.git


## Features

- **Bendability estimation**: Calculate bendability of DNA sequences based on trinucleotide or dinucleotide features.

- **Customizable**: Specify feature types and window sizes for in-depth analysis.

- **Five ways in**: a single sequence, a list, a FASTA file, a labelled table of
  sequences and values, or a BED file read against a reference genome.

- **Built for model fitting**: a labelled table gives you the features and the
  targets from one file.

- **Python Library**: Import DNAflexpy as a library in Jupyter notebooks or other scripts.

> `pip install -e .` gives you a `DNAflexpy` command:
> `DNAflexpy profile sequences.fa --feature DNaseI --window-size 10`.
> The original one still works too, and produces byte-identical output:
> `python -m rxv.DNAflexpy.cli`.

---

## Table of Contents

- [Installation](installation.md)
- [Usage](usage.md)
- [API Reference](api_reference.md)
- [FAQ](faq.md)
