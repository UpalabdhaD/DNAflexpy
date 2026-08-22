"""Generate the example data shipped in `DNAflexpy/data/examples/`.

Run this only to regenerate the files; they are committed, so a normal
install does not need it.

    python scripts/make_example_data.py

The sequences are synthetic but not random noise. Each "promoter" carries an
AT-rich stretch around position 40, which is what a real TATA-box region looks
like to a flexibility profile: AT-rich DNA is more bendable, so the example
metaprofile shows an actual feature instead of a flat line. A tutorial whose
figures show nothing teaches nothing.

The matched control set has the same length and the same overall GC content
but no positioned element, so `metaprofile(background=...)` has something
honest to compare against.
"""
import pathlib
import random

HERE = pathlib.Path(__file__).resolve().parent.parent
OUT = HERE / "DNAflexpy" / "data" / "examples"

SEED = 20260822
LENGTH = 120
N_PROMOTERS = 12
ELEMENT_AT = 40          # where the AT-rich stretch starts
ELEMENT_LEN = 20


def _background(rng, n, gc=0.45):
    """Random sequence at a given GC fraction."""
    return "".join(
        rng.choice("GC") if rng.random() < gc else rng.choice("AT")
        for _ in range(n)
    )


def _promoter(rng):
    """Background sequence with an AT-rich element dropped into it."""
    seq = list(_background(rng, LENGTH))
    for i in range(ELEMENT_AT, ELEMENT_AT + ELEMENT_LEN):
        seq[i] = rng.choice("AT")
    return "".join(seq)


def _control(rng):
    """Same length and GC, no positioned element.

    The AT-rich bases are still there, just scattered, so a difference in the
    metaprofile reflects *position* rather than composition.
    """
    body = _background(rng, LENGTH - ELEMENT_LEN)
    scattered = "".join(rng.choice("AT") for _ in range(ELEMENT_LEN))
    mixed = list(body + scattered)
    rng.shuffle(mixed)
    return "".join(mixed)


def _fasta(records):
    return "".join(f">{name}\n{seq}\n" for name, seq in records)


def main() -> None:
    rng = random.Random(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    promoters = [(f"promoter_{i + 1}", _promoter(rng)) for i in range(N_PROMOTERS)]
    controls = [(f"control_{i + 1}", _control(rng)) for i in range(N_PROMOTERS)]

    (OUT / "promoters.fa").write_text(_fasta(promoters))
    (OUT / "control.fa").write_text(_fasta(controls))

    # A labelled table: sequence plus a made-up binding score. The score is
    # built to correlate with AT content in the element, so a model fitted on
    # it in the tutorial actually has signal to find.
    rows = []
    for name, seq in promoters:
        element = seq[ELEMENT_AT:ELEMENT_AT + ELEMENT_LEN]
        at = sum(1 for b in element if b in "AT") / len(element)
        score = round(2.0 + 6.0 * at + rng.uniform(-0.3, 0.3), 3)
        rows.append(f"{name}\t{seq}\t{score}")
    (OUT / "affinity.tsv").write_text(
        "name\tsequence\tbinding_score\n" + "\n".join(rows) + "\n")

    # A tiny two-contig genome and a BED file pointing into it, so the BED
    # example needs no external download.
    chr1 = _background(rng, 2000)
    chr2 = _background(rng, 1200)
    (OUT / "genome.fa").write_text(_fasta([("chr1", chr1), ("chr2", chr2)]))
    (OUT / "peaks.bed").write_text(
        "chr1\t300\t420\tpeak_1\t0\t+\n"
        "chr1\t900\t1020\tpeak_2\t0\t-\n"
        "chr1\t1500\t1620\tpeak_3\t0\t+\n"
        "chr2\t200\t320\tpeak_4\t0\t+\n"
        "chr2\t700\t820\tpeak_5\t0\t-\n"
    )

    for path in sorted(OUT.iterdir()):
        print(f"  {path.name:16} {path.stat().st_size:6} bytes")
    print(f"wrote example data to {OUT}")


if __name__ == "__main__":
    main()
