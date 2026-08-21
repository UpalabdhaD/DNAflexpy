"""Profiling results and their serialisation."""
from __future__ import annotations

import math

import pandas as pd


class FlexProfile:
    """Profile values for a set of sequences, one row per sequence.

    Rows arrive as `[seqid, *values]` and may be ragged: a sequence shorter
    than the window contributes only its ID.
    """

    def __init__(self, rows, feature: str, window_size: int, kmer_len: int,
                 y=None, seqs=None):
        self._rows = [list(r) for r in rows]
        self.feature = feature
        self.window_size = window_size
        self.kmer_len = kmer_len
        self.y = y
        self._seqs = list(seqs) if seqs is not None else None

    @property
    def seqs(self) -> list[str] | None:
        """The input sequences, aligned to `.seqids`, or None.

        Only sequence one-hot encoding needs these; every other operation
        works from the profile values alone.
        """
        return list(self._seqs) if self._seqs is not None else None

    @property
    def seqids(self) -> list[str]:
        return [r[0] for r in self._rows]

    @property
    def n_masked(self) -> dict[str, int]:
        """Masked measurements per sequence, excluding ragged padding.

        Padding is NaN in the frame too, which is why this counts from the
        original rows rather than from the padded frame.
        """
        return {
            row[0]: sum(1 for v in row[1:] if isinstance(v, float) and math.isnan(v))
            for row in self._rows
        }

    @property
    def frame(self) -> pd.DataFrame:
        """Values as a DataFrame indexed by sequence ID, NaN-padded."""
        widest = max((len(r) - 1 for r in self._rows), default=0)
        data = [r[1:] + [float("nan")] * (widest - (len(r) - 1)) for r in self._rows]
        return pd.DataFrame(data, index=self.seqids, columns=range(1, widest + 1))

    def to_frame(self, tidy: bool = False) -> pd.DataFrame:
        if not tidy:
            return self.frame
        long = (
            self.frame.stack(future_stack=True)
            .rename("value")
            .reset_index()
            .rename(columns={"level_0": "seqid", "level_1": "position"})
        )
        long["feature"] = self.feature
        return long[["seqid", "position", "value", "feature"]]

    def encode(self, feature_names, normalize: bool = True):
        """Build a design matrix from this profile. See `DNAflexpy.encode`."""
        from DNAflexpy.encode import encode

        return encode(self, feature_names, normalize=normalize)

    def heatmap(self, **kwargs):
        """Draw this profile as a heatmap. See `DNAflexpy.plotting.heatmap`."""
        from DNAflexpy.plotting import heatmap

        return heatmap(self, **kwargs)

    def metaprofile(self, **kwargs):
        """Draw the position-wise average. See `DNAflexpy.plotting.metaprofile`."""
        from DNAflexpy.plotting import metaprofile

        return metaprofile(self, **kwargs)

    def to_tsv(self, path) -> None:
        """Write the archive's exact format.

        Built from ragged lists so pandas pads with NaN and renders each as
        an empty field, reproducing the archive's trailing tabs. Masked
        values are also NaN and therefore also empty: no `na_rep` could
        separate the two without relabelling padding and breaking byte
        equality.
        """
        pd.DataFrame(self._rows).to_csv(path, index=False, header=False, sep="\t")
