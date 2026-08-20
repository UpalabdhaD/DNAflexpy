"""Write the input files the documentation examples refer to.

`scripts/check_doc_examples.py` runs every python block in the docs against a
working directory. Several blocks open files by name -- `affinity.tsv`,
`peaks.bed` and so on -- which a reader would supply themselves. This writes a
set with those names so the checker can run:

    python scripts/make_doc_fixtures.py /tmp/docfix
    python scripts/check_doc_examples.py /tmp/docfix

Keep this in step with the docs: a new example that opens a new filename needs
its fixture added here, or the checker reports a FileNotFoundError that looks
like a broken example.
"""
import pathlib
import sys

DINUCS = [a + b for a in "ACGT" for b in "ACGT"]

FIXTURES = {
    # docs/usage.md -- labelled tables
    "affinity.tsv": "ATGCGTACGT\t1.5\nCGTAGCTAGT\t2.5\n",
    "affinity.csv": "sequence,affinity\nATGCGTACGT,1.5\nCGTAGCTAGT,2.5\n",
    "labelled.tsv": "sequence\taffinity\nATGCGTACGT\t1.5\nCGTAGCTAGT\t2.5\n",
    "named.tsv": (
        "name\tsequence\taffinity\n"
        "p1\tATGCGTACGT\t1.5\n"
        "p2\tCGTAGCTAGT\t2.5\n"
    ),
    # FASTA -- usage.md says "sequences.fa", README.md says "input.fasta"
    "sequences.fa": (
        ">a\nATGCGTACGTAGCTAGCGTAGCTAGT\n"
        ">b\nCGTAGCTAGTATGCGTACGTAGCTA\n"
    ),
    "input.fasta": ">a\nATGCGTACGTAGCTAGCGTAGCTAGT\n",
    # BED against a reference genome
    "TAIR10.fa": ">chr1\n" + "ACGT" * 200 + "\n",
    "peaks.bed": "chr1\t100\t300\tpeakA\nchr1\t400\t600\tpeakB\n",
    # A user-supplied lookup table. Two features, because usage.md names
    # "my_feature" and faq.md names "mine" -- both read this one file.
    "my_features.yaml": "".join(
        f"{feature}:\n" + "".join(f"  {kmer}: {i}.0\n" for i, kmer in enumerate(DINUCS))
        for feature in ("my_feature", "mine")
    ),
}


def main() -> None:
    target = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    target.mkdir(parents=True, exist_ok=True)
    for name, text in FIXTURES.items():
        (target / name).write_text(text)
    print(f"wrote {len(FIXTURES)} fixture(s) to {target}")


if __name__ == "__main__":
    main()
