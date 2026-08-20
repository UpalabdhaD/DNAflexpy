# Welcome to DNAflexpy

Git URL: https://github.com/upalabdhaD/DNAflexpy.git

DNAflexpy is a Python package for computing DNA flexibility profiles from sequence
data using published di- and tri-nucleotide feature lookup tables. It provides
both a command-line interface and a Python API, supports sliding-window
averaging, and outputs results in a tabular format that is easy to integrate
with downstream analyses.

## Key features

- Fast profile generation from multi-FASTA inputs.
- CLI and Python API for scripts, notebooks, or pipelines.
- Sliding-window averaging for local flexibility estimates.
- Multiple supported features from lookup tables.

## Installation

### With pip (GitHub)

```bash
pip install git+https://github.com/upalabdhaD/DNAflexpy.git
```

### From source

```bash
git clone https://github.com/upalabdhaD/DNAflexpy.git
cd DNAflexpy
pip install .
```

## Quick start (CLI)

```bash
DNAflexpy profile sequences.fa --feature DNaseI --window-size 10
```

That writes `sequences_w10nt_DNaseI.tsv`.

```bash
DNAflexpy profile sequences.fa --feature DNaseI gc    # one file per feature
DNAflexpy profile peaks.bed --genome TAIR10.fa --width 200
DNAflexpy profile affinity.tsv --no-header
DNAflexpy profile --seq ATGCGTACGT                    # straight to stdout
DNAflexpy encode sequences.fa --features 1-mer 1-DNaseI --out X.npz
```

### CLI options (summary)

- `input`: a FASTA, BED or table file. The format comes from the extension;
  override with `--format`.
- `--feature`: one or more feature names. Several means one file each.
- `--window-size`: window length; `0` gives per-k-mer values.
- `--threads`: worker processes. Left out, the tool decides.
- `--outfile`: output path, for a single feature only.
- `--genome`, `--width`, `--on-edge`: for BED input.
- `--header` / `--no-header`, `--seq-col`, `--value-col`: for table input.

The original command is untouched and still works. Its output is byte-for-byte
identical to the new one:

```bash
python -m rxv.DNAflexpy.cli "<path/to/fasta>" \
  --window-size 10 \
  --feature "DNaseI" \
  --outfile "<path/to/output.tsv>"
```

## Python API

```py
from DNAflexpy import FlexProfiler

p = FlexProfiler(feature="DNaseI", window_size=10)

p.profile("ATGCGTACGTAGCTAGCGTAGCTAGT")   # one sequence -> array of values
prof = p.profile_fasta("input.fasta")      # a whole file
prof.to_tsv("out.tsv")
```

Or for a single sequence, without building anything:

```py
import DNAflexpy

DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT", feature="DNaseI", window_size=10)
```

Sequences can also come from a labelled table or from a BED file plus a
reference genome. See [docs/usage.md](docs/usage.md).

## Supported features

Thirteen features. Four use 3-mers:

- DNaseI, NPP, bendabilityDNase, bendabilityConcensus

Nine use 2-mers:

- wedge, prop, freeen, gc, twistDisp, stiffness, bendingStiffness, mechen, trx

## Window-size behavior

- `window_size == 0` disables windowing and returns per-kmer feature values for
  the full sequence.
- `window_size == N` generates overlapping windows of length `N` (step size 1),
  computes per-kmer feature values within each window, and reports the mean
  value per window.
- Note: `window_size > 0` shortens the feature vector compared to per-kmer output
  because windows collapse multiple k-mer values into one mean per window.
- If `N` equals the feature k-mer length (e.g., `N=3` for `DNaseI`), each window
  contains a single k-mer value, so outputs match the per-kmer values.
- If `N` is larger than the sequence length, no windows are generated and only
  the sequence ID is returned for that record.
- If `N` is smaller than the feature k-mer length, window averages are 0.0
  because no k-mers fit inside the window.

## Input and output format

Input must be a FASTA or multi-FASTA file. Output is a tab-separated file with
one row per input record. The first column is the sequence identifier followed
by either:

- per-kmer feature values (`window_size == 0`), or
- per-window mean feature values (`window_size > 0`).

## Documentation

See the documentation in `docs/` or build the site with MkDocs.

## Citation

Please refer to `docs/citation.md` for the current citation information, or run
`DNAflexpy --citation` to print the BibTeX placeholder.


Print the citation:

```bash
DNAflexpy --citation
```

## License

See `LICENSE`.
