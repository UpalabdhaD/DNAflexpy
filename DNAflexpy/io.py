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


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_bed(path):
    """Read a BED file into `(chrom, start, end, name, strand)` tuples.

    Coordinates stay exactly as BED defines them: 0-based, half-open.
    `track`, `browser` and `#` lines are skipped, as are blank lines.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"BED file not found: {path}")

    out = []
    with path.open() as handle:
        for number, line in enumerate(handle):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"line {number} of {path} has {len(fields)} column(s); "
                    "BED needs at least 3 (chrom, start, end)"
                )
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 and fields[3] != "." else None
            strand = fields[5] if len(fields) > 5 and fields[5] in ("+", "-") else "+"
            out.append((chrom, start, end, name, strand))
    return out


def extract_intervals(intervals, genome, width=None, on_edge="drop"):
    """Pull sequences for BED intervals out of a reference genome FASTA.

    With `width`, every interval is re-centred on its midpoint and cut to
    exactly that many bases, so all outputs are the same length. That is
    what makes position-wise comparison downstream meaningful.

    A centred window can run off the end of a contig. `on_edge` decides:
    `"drop"` skips it and warns with a count, `"error"` raises naming the
    intervals, `"pad"` pads with `N`. Padded positions do not resolve to a
    k-mer, so they are masked as NaN and counted in `FlexProfile.n_masked`
    - they are not silently scored as zero.
    """
    if on_edge not in ("drop", "error", "pad"):
        raise ValueError(
            f"on_edge must be 'drop', 'error' or 'pad', got {on_edge!r}"
        )
    try:
        from pyfaidx import Fasta
    except ImportError:
        raise ImportError(
            "reading BED input needs pyfaidx, which is an optional extra. "
            "Install it with: pip install 'DNAflexpy[bed]'"
        ) from None

    fasta = Fasta(str(genome))
    out, dropped, padded = [], [], []

    for chrom, start, end, name, strand in intervals:
        if chrom not in fasta:
            raise ValueError(
                f"contig {chrom!r} is not in {genome}; "
                f"it has {len(fasta.keys())} contig(s)"
            )
        length = len(fasta[chrom])

        if width is not None:
            centre = (start + end) // 2
            start = centre - width // 2
            end = start + width

        seqid = name if name is not None else f"{chrom}:{start}-{end}"

        if start < 0 or end > length:
            if on_edge == "error":
                raise ValueError(
                    f"interval {seqid} runs past the bounds of {chrom} "
                    f"(length {length}); pass on_edge='drop' or 'pad'"
                )
            if on_edge == "drop":
                dropped.append(seqid)
                continue
            left = "N" * max(0, -start)
            right = "N" * max(0, end - length)
            body = str(fasta[chrom][max(0, start):min(end, length)])
            sequence = left + body + right
            padded.append(seqid)
        else:
            sequence = str(fasta[chrom][start:end])

        if strand == "-":
            sequence = sequence.translate(_COMPLEMENT)[::-1]
        out.append((seqid, sequence))

    if dropped:
        warnings.warn(
            f"{len(dropped)} interval(s) dropped for running past a contig "
            f"boundary: {', '.join(dropped[:5])}"
            + (" ..." if len(dropped) > 5 else ""),
            UserWarning,
            stacklevel=2,
        )
    if padded:
        warnings.warn(
            f"{len(padded)} interval(s) padded with N at a contig boundary. "
            "Padded positions do not resolve to a k-mer and are masked as "
            "NaN, so those windows average fewer values.",
            UserWarning,
            stacklevel=2,
        )
    warn_if_ambiguous(out, genome)
    return out
