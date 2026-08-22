"""Citation metadata, in one place.

The CLI's `--citation` output and `docs/citation.md` both read from here, so
they cannot drift apart. `version` is not stored: it comes from
`DNAflexpy.__version__`, which is itself the single source of the version.
"""
from __future__ import annotations

AUTHORS = ["Upalabdha Dey", "Rajesh Yella", "Aditya Kumar"]
TITLE = "DNAflexpy: DNA flexibility profiling from sequence"
YEAR = "2026"
URL = "https://github.com/UpalabdhaD/DNAflexpy"


def bibtex(version: str) -> str:
    """A BibTeX entry for this installed version."""
    authors = " and ".join(AUTHORS)
    return (
        "@software{DNAflexpy,\n"
        f"  author    = {{{authors}}},\n"
        f"  title     = {{{TITLE}}},\n"
        f"  year      = {{{YEAR}}},\n"
        "  publisher = {GitHub},\n"
        f"  version   = {{{version}}},\n"
        f"  url       = {{{URL}}}\n"
        "}"
    )


def plain(version: str) -> str:
    """A one-line citation for a reference list."""
    authors = ", ".join(AUTHORS)
    return f"{authors} ({YEAR}). {TITLE}. Version {version}. {URL}"
