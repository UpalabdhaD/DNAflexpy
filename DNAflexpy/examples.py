"""Small example datasets that ship with the package.

The equivalent of DNAshapeR's `system.file("extdata", ..., package=)`: every
example in the documentation runs against these, so a reader can follow along
without downloading anything.

    from DNAflexpy import example_path

    prof = FlexProfiler("DNaseI").profile_fasta(example_path("promoters.fa"))

The sequences are synthetic, but not noise. Each promoter carries an AT-rich
stretch around position 40 -- what a TATA-box region looks like to a
flexibility profile -- and `control.fa` holds matched sequences with the same
length and base composition but no positioned element. That is what makes the
worked examples show a real feature rather than a flat line.

Regenerate them with `python scripts/make_example_data.py`.
"""
from __future__ import annotations

from importlib.resources import files
from pathlib import Path

_DESCRIPTIONS = {
    "promoters.fa": "12 synthetic promoters, 120 bp each, with an AT-rich element at 40-60",
    "control.fa": "12 matched controls: same length and base composition, no positioned element",
    "affinity.tsv": "the same promoters with a binding score, for model fitting (has a header)",
    "genome.fa": "a two-contig toy genome, chr1 (2000 bp) and chr2 (1200 bp)",
    "peaks.bed": "5 intervals into that genome, named and stranded",
}


def example_path(name: str) -> Path:
    """The full path to a packaged example file.

    Raises with the available names if `name` is not one of them.
    """
    if name not in _DESCRIPTIONS:
        raise ValueError(
            f"no example file named {name!r}. Available: "
            f"{', '.join(sorted(_DESCRIPTIONS))}"
        )
    path = Path(str(files("DNAflexpy").joinpath("data", "examples", name)))
    if not path.exists():
        raise FileNotFoundError(
            f"example file {name!r} is missing from the installed package at "
            f"{path}. Regenerate it with scripts/make_example_data.py, or "
            "reinstall."
        )
    return path


def example_files() -> list[str]:
    """The names of every packaged example file, sorted."""
    return sorted(_DESCRIPTIONS)


def describe_examples() -> dict[str, str]:
    """`{filename: one-line description}` for every example file."""
    return dict(_DESCRIPTIONS)
