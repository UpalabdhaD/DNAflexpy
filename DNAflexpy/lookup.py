"""Feature lookup tables: load once, validate, infer k-mer length.

The archive hardcoded a feature -> k-mer-length map that had to be kept in
sync with the YAML by hand, and silently produced None rows for any feature
missing from it. Here k is inferred from the table's own keys, which is
safe because validation rejects any table with mixed key lengths.
"""
from __future__ import annotations

import functools
import warnings
from collections.abc import Mapping
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType

import yaml

_ALPHABET = frozenset("ACGT")


class FeatureTable:
    """A validated collection of k-mer -> value tables, keyed by feature."""

    def __init__(self, tables: dict[str, dict[str, float]], kmer_lens: dict[str, int]):
        self._tables = tables
        self._kmer_lens = kmer_lens

    @property
    def features(self) -> list[str]:
        return list(self._tables)

    def kmer_len(self, feature: str) -> int:
        self._require(feature)
        return self._kmer_lens[feature]

    def table(self, feature: str) -> Mapping[str, float]:
        self._require(feature)
        return MappingProxyType(self._tables[feature])

    def _require(self, feature: str) -> None:
        if feature not in self._tables:
            raise ValueError(
                f"unknown feature {feature!r}. Available: {', '.join(self.features)}"
            )

    @classmethod
    def from_dict(cls, mapping: dict) -> "FeatureTable":
        tables: dict[str, dict[str, float]] = {}
        kmer_lens: dict[str, int] = {}
        for feature, raw in mapping.items():
            tables[feature], kmer_lens[feature] = _validate(feature, raw)
        return cls(tables, kmer_lens)

    @classmethod
    def from_yaml(cls, path=None) -> "FeatureTable":
        if path is None:
            source = files("DNAflexpy.data").joinpath("lookupNEW.yaml")
            text = source.read_text()
        else:
            text = Path(path).read_text()
        try:
            loaded = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise ValueError(f"could not parse feature table: {exc}") from exc
        if not isinstance(loaded, dict):
            raise ValueError("feature table must map feature names to k-mer tables")
        return cls.from_dict(loaded)


def _validate(feature: str, raw) -> tuple[dict[str, float], int]:
    """Return an uppercased table and its inferred k-mer length."""
    if not isinstance(raw, dict) or not raw:
        raise ValueError(f"feature {feature!r} has an empty or malformed table")

    table: dict[str, float] = {}
    for key, value in raw.items():
        kmer = str(key).upper()
        if set(kmer) - _ALPHABET:
            raise ValueError(
                f"feature {feature!r} has non-ACGT k-mer {kmer!r}; "
                "tables must contain unambiguous k-mers only"
            )
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(
                f"feature {feature!r} has non-numeric value {value!r} for {kmer!r}"
            )
        table[kmer] = value

    lengths = {len(k) for k in table}
    if len(lengths) > 1:
        raise ValueError(
            f"feature {feature!r} has mixed k-mer lengths {sorted(lengths)}; "
            "k-mer length is inferred from the keys and must be uniform"
        )
    kmer_len = lengths.pop()

    expected = 4 ** kmer_len
    if len(table) < expected:
        warnings.warn(
            f"feature {feature!r} table is incomplete: {len(table)} of {expected} "
            f"{kmer_len}-mers. Missing k-mers resolve to NaN and are masked.",
            UserWarning,
            stacklevel=3,
        )
    return table, kmer_len


@functools.lru_cache(maxsize=1)
def default_table() -> FeatureTable:
    """The packaged lookup table, parsed at most once per process."""
    return FeatureTable.from_yaml()
