<h1 align="center">DNAflexpy</h1>

<p align="center">
  <strong>DNA flexibility profiling from sequence</strong><br>
  Sliding-window profiles, machine-learning features and publication figures,
  from di- and trinucleotide parameter tables.
</p>

<p align="center">
  <a href="https://github.com/UpalabdhaD/DNAflexpy/actions/workflows/ci.yml">
    <img alt="CI" src="https://github.com/UpalabdhaD/DNAflexpy/actions/workflows/ci.yml/badge.svg?branch=main"></a>
  <a href="https://upalabdhad.github.io/DNAflexpy/">
    <img alt="Documentation" src="https://img.shields.io/badge/docs-online-14b8a6"></a>
  <img alt="Python" src="https://img.shields.io/badge/python-3.12%20%7C%203.13-3776ab?logo=python&logoColor=white">
  <a href="LICENSE">
    <img alt="License" src="https://img.shields.io/badge/license-BSD--3--Clause-blue"></a>
</p>

<p align="center">
  <img alt="Tests" src="https://img.shields.io/badge/tests-552%20passing-brightgreen">
  <img alt="Byte-equality gate" src="https://img.shields.io/badge/output-byte--identical%20to%20v1-8b5cf6">
  <img alt="Container" src="https://img.shields.io/badge/container-Docker%20%7C%20Apptainer-2496ed?logo=docker&logoColor=white">
  <img alt="Platform" src="https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey">
</p>

---

## What it does

DNA is not uniformly stiff. Some sequences bend easily, others resist — and
that varies along a sequence in ways that matter for nucleosome positioning,
promoter architecture and protein binding.

DNAflexpy turns a sequence into a **flexibility profile**: one number per
position, computed from published experimental parameter tables. Thirteen of
them ship with the package, covering bendability, propeller twist, free
energy, stiffness and more.

From there you can plot it, compare it against a control, or turn it into a
design matrix and fit a model.

```python
import DNAflexpy

DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT")
# array([ 0.011, -0.003, -0.001,  0.01 ,  0.017, ... ])
```

## Install

```bash
pip install git+https://github.com/UpalabdhaD/DNAflexpy.git
```

For plotting and BED input:

```bash
pip install "DNAflexpy[plot,bed] @ git+https://github.com/UpalabdhaD/DNAflexpy.git"
```

Python 3.12 or newer. From a clone: `pip install -e '.[plot,bed]'`.

There is also a [container](Docker/README.md) with everything already
installed, which runs under Apptainer on a cluster without root.

## Sixty-second tour

Example data ships with the package, so this runs as written:

```python
from DNAflexpy import FlexProfiler, example_path

p = FlexProfiler(feature="DNaseI", window_size=10)
prof = p.profile_fasta(example_path("promoters.fa"))

prof.frame          # 12 sequences x 111 positions, as a DataFrame
prof.to_tsv("profiles.tsv")
```

Draw it:

```python
ax = prof.heatmap(nbins=30)
ax.figure.savefig("heatmap.png", dpi=150, bbox_inches="tight")
```

Compare it against a matched control:

```python
control = p.profile_fasta(example_path("control.fa"))
ax = prof.metaprofile(background=control)
```

Build a design matrix for scikit-learn:

```python
labelled = p.profile_table(example_path("affinity.tsv"),
                           seq_col="sequence", value_col="binding_score",
                           header=True)
fm = labelled.encode(["1-mer", "1-DNaseI"])
fm.X, fm.y          # features and targets, from one file
```

Or from the shell:

```bash
DNAflexpy profile promoters.fa --feature DNaseI --window-size 10
DNAflexpy plot heatmap promoters.fa --nbins 30 --out heat.png
DNAflexpy encode promoters.fa --features 1-mer 1-DNaseI --out X.npz
```

## Features

| | |
|---|---|
| **Five ways in** | a bare string, a list or dict, a FASTA file, a labelled table, or BED intervals against a reference genome |
| **Thirteen parameter tables** | bendability, stiffness, propeller twist, free energy, GC content and more — or supply your own |
| **Machine-learning encoding** | one-hot sequence blocks and nth-order interaction terms, straight into scikit-learn |
| **Three figures** | heatmap, metaprofile with a background, stacked trackplot |
| **Command line and API** | the same five inputs from either |
| **Parallel** | multiprocessing for large FASTA files, with a threshold so small inputs stay fast |
| **Containerised** | Docker and Apptainer, unprivileged, for HPC |

### The parameter tables

Four use trinucleotides:

`DNaseI` · `NPP` · `bendabilityDNase` · `bendabilityConcensus`

Nine use dinucleotides:

`wedge` · `prop` · `freeen` · `gc` · `twistDisp` · `stiffness` ·
`bendingStiffness` · `mechen` · `trx`

k-mer length is read from each table's own keys, so adding a table to
`lookupNEW.yaml` is the only step needed to use it.

## Two things it deliberately refuses

Both are cases where it would be easy to return a plausible-looking answer.

**A column-scaled metaprofile with no background.** Standardising each position
across sequences and then averaging down that same position cancels to exactly
zero, at every position, by construction. The line would be flat and carry no
information, so it raises and explains why.

**Sorting heatmap rows by coefficient of variation on signed data.** CV is
`std/mean`; several tables take negative values, which inverts the sort. It
raises rather than showing you a confidently wrong order.

## Output is byte-identical to the original

This package is a rewrite. The original is preserved, frozen, at
`rxv/DNAflexpy/`, and **230 tests compare the new output against it byte for
byte** on every push — three FASTA files, ten features, seven window sizes,
plus the parallel code path.

Byte comparison, not value comparison. That distinction caught a real bug:
six tables hold integers, and coercing them to float wrote `18.0` where the
original writes `18`. `assert a == b` passes on `18 == 18.0`; only
serialisation exposes it.

The rewrite is also **133× faster per sequence** — the original re-parsed its
400-entry YAML table on every call.

## Documentation

**[upalabdhad.github.io/DNAflexpy](https://upalabdhad.github.io/DNAflexpy/)**

- [Tutorial](https://upalabdhad.github.io/DNAflexpy/tutorial/) — a full worked
  analysis using the packaged example data
- [Usage](https://upalabdhad.github.io/DNAflexpy/usage/) — every input, option
  and output format
- [Feature tables](https://upalabdhad.github.io/DNAflexpy/features/) — what
  each of the thirteen measures
- [API reference](https://upalabdhad.github.io/DNAflexpy/api_reference/)
- [FAQ](https://upalabdhad.github.io/DNAflexpy/faq/)

## Citation

If DNAflexpy is useful in your work, please cite it. See
[docs/citation.md](docs/citation.md), or run:

```bash
DNAflexpy --citation
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).
