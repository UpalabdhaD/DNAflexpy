# Installation

Follow these steps to install DNAflexpy and its dependencies.

## Prerequisites

Python 3.11 or later. Check yours with:

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




