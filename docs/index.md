# DNAflexpy

**DNA flexibility profiling from sequence.**

DNA is not uniformly stiff. Some sequences bend easily, others resist — and it
varies along a sequence in ways that matter for nucleosome positioning,
promoter architecture and protein binding.

DNAflexpy turns a sequence into a **flexibility profile**: one number per
position, computed from published di- and trinucleotide parameter tables.
Thirteen of them ship with the package. From there you can plot it, compare it
against a control, or turn it into a design matrix and fit a model.

```python
import DNAflexpy

print(DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT"))
```

```
[ 0.011 -0.003 -0.001  0.01   0.017  0.025  0.033  0.039  0.034  0.026
  0.018  0.027  0.027  0.018  0.018  0.027  0.014]
```

## Start here

<div class="grid cards" markdown>

- **[Tutorial](tutorial.md)**

    A complete worked analysis using data that ships with the package. Nothing
    to download. Start here if you are new.

- **[Installation](installation.md)**

    pip, optional extras, and the container.

- **[Usage](usage.md)**

    Every input format, option and output, in detail.

- **[Feature tables](features.md)**

    What each of the thirteen tables measures, and how to add your own.

- **[API reference](api_reference.md)**

    Every class, method and function.

- **[FAQ](faq.md)**

    Masked values, window sizes, ambiguous bases.

</div>

## What it can do

| | |
|---|---|
| **Five ways in** | a bare string, a list or dict, a FASTA file, a labelled table of sequences and values, or BED intervals against a reference genome |
| **Thirteen parameter tables** | bendability, stiffness, propeller twist, free energy, GC content and more — or supply your own |
| **Machine-learning encoding** | one-hot sequence blocks and nth-order interaction terms, straight into scikit-learn |
| **Three figures** | heatmap, metaprofile with a background, stacked trackplot |
| **Command line and API** | the same five inputs from either |
| **Containerised** | Docker and Apptainer, unprivileged, for HPC |

## Two ways to run it

Python:

```python
from DNAflexpy import FlexProfiler, example_path

p = FlexProfiler(feature="DNaseI", window_size=10)
prof = p.profile_fasta(example_path("promoters.fa"))
print(prof.frame.shape)
```

```
(12, 111)
```

Or the command line:

```bash
DNAflexpy profile promoters.fa --feature DNaseI --window-size 10
DNAflexpy plot heatmap promoters.fa --nbins 30 --out heat.png
DNAflexpy encode promoters.fa --features 1-mer 1-DNaseI --out X.npz
```

## Verified against the original

DNAflexpy is a rewrite. The original implementation is preserved, frozen, at
`rxv/DNAflexpy/`, and **230 tests compare the new output against it byte for
byte** on every push — three FASTA files, ten features, seven window sizes,
plus the parallel code path.

The rewrite is also about **133× faster per sequence**: the original re-parsed
its 400-entry parameter table on every call.

The original command still works, and produces identical bytes:

```bash
python -m rxv.DNAflexpy.cli sequences.fa --window-size 10 --feature DNaseI
```

## Citation

See [Citation](citation.md), or run `DNAflexpy --citation`.

Repository: [github.com/UpalabdhaD/DNAflexpy](https://github.com/UpalabdhaD/DNAflexpy)
