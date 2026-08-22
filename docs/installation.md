# Installation

Follow these steps to install DNAflexpy and its dependencies.

## Prerequisites

Python 3.12 or later. Check yours with:

!!! warning "Why 3.12 and not earlier"

    Python 3.12 changed how the builtin `sum()` adds floating-point numbers.
    The same sequence therefore gives slightly different values on Python 3.11:
    a window that reads `0.011` on 3.12 reads `0.01` on 3.11. Every number in
    this documentation, and every expected result in the test suite, is a 3.12
    value. Requiring 3.12 keeps results reproducible.

## Dependencies

Installed for you:

- pandas
- pyyaml
- numpy

Optional extras:

- `pip install -e '.[bed]'` adds **pyfaidx**, needed only for BED input
- `pip install -e '.[dev]'` adds **pytest** and the BED extra, for running the tests

```bash
python --version
```

## Installation with pip

- Directly use pip to install from github

```bash
pip install git+https://github.com/upalabdhaD/DNAflexpy.git

```
OR

- Use conda to create a environment and clone the repo and then use pip to install the package

```bash
# create env with pip installation
conda create -n DNAflexpy_env pip -c conda-forge
conda activate DNAflexpy_env

# git clone 
git clone https://github.com/upalabdhaD/DNAflexpy.git

# move to the directory root
cd DNAflexpy

# then install with pip 
pip install .
```





## Container

If you would rather not install anything, or you need this on a cluster, build
the image:

```bash
docker build -t dnaflexpy:dev -f Docker/Dockerfile .
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data dnaflexpy:dev \
  profile in.fa --feature DNaseI --window-size 10
```

It has every optional extra already installed, so plotting and BED input work
out of the box.

On an HPC cluster, Apptainer runs it without root:

```bash
apptainer build dnaflexpy.sif docker-daemon://dnaflexpy:dev
apptainer exec --bind "$PWD":/data dnaflexpy.sif \
  DNAflexpy profile /data/in.fa --feature DNaseI --outfile /data/out.tsv
```

Nothing is published to a registry, so there is no image to pull — you build
it from the repository. Full instructions, including why it works unprivileged,
are in [`Docker/README.md`](https://github.com/UpalabdhaD/DNAflexpy/blob/main/Docker/README.md).
