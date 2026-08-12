"""The profiling engine."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from DNAflexpy.kernel import profile_values
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.profile import FlexProfile


class ProfileSet(dict):
    """`{feature: FlexProfile}`.

    Multi-feature runs cannot share one matrix: features have different k,
    so their vectors differ in length for the same sequence.
    """

    def to_tidy(self) -> pd.DataFrame:
        return pd.concat([p.to_frame(tidy=True) for p in self.values()],
                         ignore_index=True)


class FlexProfiler:
    """Profiles sequences against one or more flexibility features.

    The lookup table is loaded once, at construction. The archive reloaded
    it on every per-sequence call, which cost 8 ms against 0.1 ms of actual
    computation.
    """

    def __init__(self, feature, window_size: int = 10, lookup=None):
        self._table = _resolve_lookup(lookup)
        self._features = [feature] if isinstance(feature, str) else list(feature)
        if not self._features:
            raise ValueError("at least one feature is required")
        # Validate now so a typo fails immediately, not after a full pass.
        for name in self._features:
            self._table.kmer_len(name)
        if window_size < 0:
            raise ValueError(f"window_size must be >= 0, got {window_size}")
        self.window_size = window_size

    @property
    def features(self) -> list[str]:
        return list(self._features)

    def kmer_len(self, feature: str | None = None) -> int:
        return self._table.kmer_len(feature or self._features[0])

    def profile(self, seq: str) -> np.ndarray:
        """Values for one sequence. Always a 1-D array, never a result object."""
        if len(self._features) != 1:
            raise ValueError(
                "profile() requires a single feature; use profile_seqs() for "
                f"multi-feature profilers (has {len(self._features)})"
            )
        return np.asarray(self._values(self._features[0], seq), dtype=float)

    def profile_seqs(self, seqs):
        """Profile a list of sequences or an {id: sequence} mapping."""
        if isinstance(seqs, dict):
            pairs = list(seqs.items())
        else:
            pairs = [(f"seq_{i}", s) for i, s in enumerate(seqs)]
        return self._build(pairs)

    def _values(self, feature: str, seq: str) -> list[float]:
        return profile_values(
            seq, self._table.kmer_len(feature), self._table.table(feature),
            self.window_size,
        )

    def _build(self, pairs):
        """Turn (seqid, sequence) pairs into a FlexProfile or ProfileSet."""
        rows_by_feature = {
            feature: [[sid, *self._values(feature, seq)] for sid, seq in pairs]
            for feature in self._features
        }
        return self._assemble(rows_by_feature)

    def _assemble(self, rows_by_feature: dict[str, list[list]]):
        """Turn per-feature row lists into a FlexProfile or ProfileSet.

        Factored out of `_build` so Task 7's `profile_fasta` can reuse this
        same final step when starting from worker-pool results instead of
        sequences.
        """
        built = {
            feature: FlexProfile(
                rows,
                feature=feature,
                window_size=self.window_size,
                kmer_len=self._table.kmer_len(feature),
            )
            for feature, rows in rows_by_feature.items()
        }
        if len(self._features) == 1:
            return built[self._features[0]]
        return ProfileSet(built)


def _resolve_lookup(lookup) -> FeatureTable:
    if lookup is None:
        return default_table()
    if isinstance(lookup, FeatureTable):
        return lookup
    if isinstance(lookup, dict):
        return FeatureTable.from_dict(lookup)
    if isinstance(lookup, (str, Path)):
        return FeatureTable.from_yaml(lookup)
    raise TypeError(
        f"lookup must be None, a path, a dict or a FeatureTable, got {type(lookup).__name__}"
    )
