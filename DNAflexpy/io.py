"""Input readers."""
from __future__ import annotations

import math
import warnings
from pathlib import Path
from typing import Iterator

import pandas as pd


def read_fasta(path) -> Iterator[tuple[str, str]]:
    """Yield `(record_id, sequence)` for each record in a FASTA file.

    Handles wrapped sequence lines and a missing trailing newline.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def read_table(path, seq_col=0, value_col=1, id_col=None, header=None, sep="\t"):
    """Read a labelled table of sequences and numeric values.

    Returns `(seqid, sequence, value)` triples in file order. This is the
    natural input for model fitting, since it carries X and y in one file.

    A row whose value is missing or non-numeric raises rather than being
    dropped: silently discarding labelled rows would corrupt a training set.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"table not found: {path}")

    if header is None:
        raise ValueError(
            f"pass header=True or header=False for {path}: whether row 1 holds "
            "column names or data cannot be guessed reliably, because words "
            "like 'dna', 'rna' and 'tag' are made only of nucleotide letters"
        )

    try:
        frame = pd.read_csv(path, sep=sep, header=0 if header else None, dtype=str)
    except pd.errors.EmptyDataError:
        raise ValueError(f"table has no data rows: {path}") from None
    if frame.empty:
        raise ValueError(f"table has no data rows: {path}")

    seqs = _pick_column(frame, seq_col, "seq_col")
    values = _pick_column(frame, value_col, "value_col")
    ids = _pick_column(frame, id_col, "id_col") if id_col is not None else None

    out = []
    for i in range(len(frame)):
        line_no = i + 1 + (1 if header else 0)
        raw = values.iloc[i]
        if raw is None or (isinstance(raw, float) and math.isnan(raw)) or str(raw).strip() == "":
            raise ValueError(f"line {line_no} of {path} has a missing value in value_col")
        try:
            value = float(raw)
        except (TypeError, ValueError):
            raise ValueError(
                f"line {line_no} of {path} has non-numeric value {raw!r} in value_col"
            ) from None
        seqid = str(ids.iloc[i]) if ids is not None else f"seq_{i}"
        out.append((seqid, str(seqs.iloc[i]), value))
    warn_if_ambiguous([(sid, seq) for sid, seq, _ in out], path)
    return out


def _pick_column(frame, col, argname):
    """Select a column by position or by name, with a useful error."""
    if isinstance(col, int):
        if col >= len(frame.columns):
            raise ValueError(
                f"{argname}={col} is out of range; the table has "
                f"{len(frame.columns)} column(s)"
            )
        return frame.iloc[:, col]
    if col not in frame.columns:
        raise ValueError(
            f"{argname}={col!r} not found; columns are {list(frame.columns)}"
        )
    return frame[col]


_ACGTN = frozenset("ACGTN")


def warn_if_ambiguous(records, source):
    """Warn when sequences carry IUPAC codes beyond ACGTN.

    The feature tables hold only ACGT k-mers, so any k-mer covering one of
    these resolves to NaN and is masked. N is an ordinary placeholder and
    is not worth warning about; the rarer codes usually are.
    """
    found = {}
    for seqid, sequence in records:
        odd = sorted(set(sequence.upper()) - _ACGTN)
        if odd:
            found[seqid] = odd
    if not found:
        return
    letters = sorted({c for codes in found.values() for c in codes})
    examples = ", ".join(list(found)[:3])
    warnings.warn(
        f"{len(found)} sequence(s) in {source} contain non-ACGTN letters "
        f"({', '.join(letters)}), for example {examples}. The feature tables "
        "hold only ACGT k-mers, so every k-mer covering one of these resolves "
        "to NaN and is masked. See FlexProfile.n_masked.",
        UserWarning,
        stacklevel=3,
    )
