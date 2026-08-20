"""The profiling engine."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from DNAflexpy.kernel import profile_values
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.results import FlexProfile

# Below this many records, process startup dominates: spawning a Pool (macOS
# uses `spawn`, which forks a fresh interpreter and re-imports the main
# module) costs far more than just running the records serially. Measured at
# ~390 ms (Pool) vs ~0.08 ms (serial) for a 2-record FASTA. `threads=None`
# means "decide automatically", not "always spawn a Pool" -- this constant is
# how that decision is made regardless of `threads`.
_MIN_RECORDS_FOR_POOL = 64


class ProfileSet(dict):
    """`{feature: FlexProfile}`.

    Multi-feature runs cannot share one matrix: features have different k,
    so their vectors differ in length for the same sequence.
    """

    def to_tidy(self) -> pd.DataFrame:
        return pd.concat([p.to_frame(tidy=True) for p in self.values()],
                         ignore_index=True)

    def encode(self, feature_names, normalize: bool = True):
        """Build a design matrix across these features. See `DNAflexpy.encode`."""
        from DNAflexpy.encode import encode

        return encode(self, feature_names, normalize=normalize)


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
        from DNAflexpy.io import warn_if_ambiguous

        warn_if_ambiguous([("sequence", seq)], "the given sequence")
        return np.asarray(self._values(self._features[0], seq), dtype=float)

    def profile_seqs(self, seqs):
        """Profile a list of sequences or an {id: sequence} mapping."""
        if isinstance(seqs, dict):
            pairs = list(seqs.items())
        else:
            pairs = [(f"seq_{i}", s) for i, s in enumerate(seqs)]
        from DNAflexpy.io import warn_if_ambiguous

        warn_if_ambiguous(pairs, "the given sequences")
        return self._build(pairs)

    def profile_fasta(self, path, threads: int | None = None):
        """Profile every record in a FASTA file.

        One read and one pool spawn regardless of how many features are
        requested. With differing k the features cannot share a matrix, so
        the saving is I/O and process startup, not a shared k-mer scan.

        `threads=None` ("decide automatically") does NOT mean "always spawn
        a Pool". Below `_MIN_RECORDS_FOR_POOL` records, this always runs
        serially, regardless of `threads` -- process startup dominates at
        that size, so a Pool would be strictly slower (measured ~4900x on a
        2-record file). Pass `threads=1` to force serial explicitly on any
        input size; pass an explicit `threads > 1` to force a Pool once the
        input is large enough to make it worthwhile.
        """
        from DNAflexpy.io import read_fasta, warn_if_ambiguous

        records = list(read_fasta(path))
        warn_if_ambiguous(records, path)
        if (threads is not None and threads <= 1) or len(records) < _MIN_RECORDS_FOR_POOL:
            return self._build(records)

        import multiprocessing

        with multiprocessing.Pool(
            processes=threads,
            initializer=_init_worker,
            initargs=(self._table, self._features, self.window_size),
        ) as pool:
            results = pool.map(_profile_record, records)

        # pool.map preserves input order, so output stays deterministic.
        rows_by_feature = {
            feature: [[seqid, *values[feature]] for seqid, values in results]
            for feature in self._features
        }
        # The workers return values only, never the sequences, so this is the
        # one path where `seqs` has to be re-attached from the parent's copy
        # of `records`. `to_tsv` ignores `seqs`, so the byte-equality gate
        # cannot catch it going missing here.
        return self._assemble(rows_by_feature, seqs=[s for _, s in records])

    def profile_table(self, path, seq_col=0, value_col=1, id_col=None,
                      header=None, sep="\t"):
        """Profile every sequence in a labelled table.

        The table's numeric column becomes `FlexProfile.y`, aligned to
        `.seqids`, so one file supplies both the features and the targets.
        `header` must be passed explicitly as `True` or `False`; omitting it
        raises, because row 1 cannot be reliably classified as data or
        column names from its content alone.
        """
        from DNAflexpy.io import read_table

        records = read_table(
            path, seq_col=seq_col, value_col=value_col, id_col=id_col,
            header=header, sep=sep,
        )
        pairs = [(seqid, seq) for seqid, seq, _ in records]
        y = [value for _, _, value in records]
        return self._build(pairs, y=y)

    def from_bed(self, bed, genome, width=None, on_edge="drop"):
        """Profile intervals from a BED file against a reference genome.

        With `width`, intervals are re-centred and cut to that many bases,
        so every row is the same length. `on_edge` decides what happens to
        an interval whose centred window runs past a contig boundary.
        """
        from DNAflexpy.io import extract_intervals, read_bed

        pairs = extract_intervals(
            read_bed(bed), genome, width=width, on_edge=on_edge
        )
        return self._build(pairs)

    def _values(self, feature: str, seq: str) -> list[float]:
        return _feature_values(self._table, feature, seq, self.window_size)

    def _build(self, pairs, y=None):
        """Turn (seqid, sequence) pairs into a FlexProfile or ProfileSet."""
        rows_by_feature = {
            feature: [[sid, *self._values(feature, seq)] for sid, seq in pairs]
            for feature in self._features
        }
        return self._assemble(rows_by_feature, y=y, seqs=[s for _, s in pairs])

    def _assemble(self, rows_by_feature: dict[str, list[list]], y=None, seqs=None):
        """Turn per-feature row lists into a FlexProfile or ProfileSet.

        Factored out of `_build` so Task 7's `profile_fasta` can reuse this
        same final step when starting from worker-pool results instead of
        sequences.

        `y` is the aligned label vector from labelled input, or None. Every
        profile in a multi-feature set shares the same labels.

        `seqs` is the aligned input sequences, kept so `encode` can build
        one-hot sequence features. It is a parameter here rather than only
        on `_build` because the pooled `profile_fasta` branch calls this
        directly and would otherwise leave it None on that path alone.
        """
        built = {
            feature: FlexProfile(
                rows,
                feature=feature,
                window_size=self.window_size,
                kmer_len=self._table.kmer_len(feature),
                y=y,
                seqs=seqs,
            )
            for feature, rows in rows_by_feature.items()
        }
        if len(self._features) == 1:
            return built[self._features[0]]
        return ProfileSet(built)


def _feature_values(table: FeatureTable, feature: str, seq: str, window_size: int) -> list[float]:
    """The one call site both the serial and pooled paths use.

    `FlexProfiler._values` (serial) and `_profile_record` (pooled worker)
    used to invoke `profile_values` independently; two copies of the same
    call is exactly the kind of drift that could quietly break byte
    equality on one path while the other still passes.
    """
    return profile_values(seq, table.kmer_len(feature), table.table(feature), window_size)


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


_WORKER_STATE: dict = {}


def _init_worker(table: FeatureTable, features: list[str], window_size: int) -> None:
    """Seed each worker once with the already-parsed table.

    Deliberately receives the loaded table, not a path: reloading per worker
    would reintroduce the 8 ms x N cost this design removes.
    """
    _WORKER_STATE["table"] = table
    _WORKER_STATE["features"] = features
    _WORKER_STATE["window_size"] = window_size


def _profile_record(record: tuple[str, str]) -> tuple[str, dict[str, list]]:
    seqid, seq = record
    table = _WORKER_STATE["table"]
    window_size = _WORKER_STATE["window_size"]
    return seqid, {
        feature: _feature_values(table, feature, seq, window_size)
        for feature in _WORKER_STATE["features"]
    }
