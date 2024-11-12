# Installation

Follow these steps to install DNAflexpy and its dependencies.

## Prerequisites

Ensure that Python 3.8 or later is installed on your system. You can check your Python version by running:

## Dependencies
- pandas
- pyyaml

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




