"""Input readers."""
from __future__ import annotations

from pathlib import Path
from typing import Iterator


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
