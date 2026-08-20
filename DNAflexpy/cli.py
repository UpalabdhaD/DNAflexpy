"""The `DNAflexpy` command.

This is a front end and nothing else. It parses arguments, builds a
`FlexProfiler`, and writes results through `FlexProfile.to_tsv()`. It never
formats profile values itself: `to_tsv` is what the 230-case byte-equality
suite covers, so a path that serialised values here would be a path nothing
verifies.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from DNAflexpy import __version__

# Author and title are still placeholders. The rest is filled in from the
# package itself, so at least the version and URL are never stale.
CITATION_TEMPLATE = """@software{{DNAflexpy,
  author    = {{ADD_AUTHOR_LIST}},
  title     = {{ADD_TITLE}},
  year      = {{ADD_YEAR}},
  publisher = {{GitHub}},
  version   = {{{version}}},
  url       = {{https://github.com/UpalabdhaD/DNAflexpy}}
}}"""

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}
_TABLE_SUFFIXES = {".tsv", ".csv", ".txt", ".tab"}
_BED_SUFFIXES = {".bed"}


def _infer_format(path: Path) -> str:
    """Guess the input format from the file extension.

    Deliberately not content sniffing: a wrong guess here would silently
    profile the wrong thing, and `--format` is one word to type.
    """
    suffix = path.suffix.lower()
    if suffix == ".gz":
        suffix = Path(path.stem).suffix.lower()
    if suffix in _FASTA_SUFFIXES:
        return "fasta"
    if suffix in _BED_SUFFIXES:
        return "bed"
    if suffix in _TABLE_SUFFIXES:
        return "table"
    raise SystemExit(
        f"cannot tell what kind of file {path.name} is from its extension. "
        "Pass --format fasta, --format bed or --format table."
    )


def _default_outfile(path: Path, window_size: int, feature: str) -> str:
    """The archive's own naming pattern, kept so habits carry over."""
    return f"{path.stem}_w{window_size}nt_{feature}.tsv"


def _profile_input(args, features):
    """Run the right entry point for the input format. Returns a mapping."""
    from DNAflexpy import FlexProfiler

    profiler = FlexProfiler(features, window_size=args.window_size,
                            lookup=args.lookup)

    if args.seq:
        result = profiler.profile_seqs({"sequence": args.seq})
        return profiler, _as_set(result, features)

    path = Path(args.input)
    if not path.exists():
        raise SystemExit(f"input file not found: {path}")
    fmt = args.format or _infer_format(path)

    if fmt == "fasta":
        result = profiler.profile_fasta(path, threads=args.threads)
    elif fmt == "bed":
        if not args.genome:
            raise SystemExit(
                "BED input needs a reference genome to read sequences from. "
                "Pass --genome <genome.fa>."
            )
        result = profiler.from_bed(path, genome=args.genome, width=args.width,
                                   on_edge=args.on_edge)
    elif fmt == "table":
        if args.header is None:
            raise SystemExit(
                f"pass --header or --no-header for {path.name}: whether row 1 "
                "holds column names or data cannot be guessed reliably, "
                "because words like 'dna', 'rna' and 'tag' are made only of "
                "nucleotide letters."
            )
        result = profiler.profile_table(
            path, seq_col=args.seq_col, value_col=args.value_col,
            id_col=args.id_col, header=args.header, sep=args.sep,
        )
    else:
        raise SystemExit(f"unknown --format {fmt!r}")

    return profiler, _as_set(result, features)


def _as_set(result, features):
    """Always hand back a ProfileSet, so callers have one shape to handle."""
    from DNAflexpy.core import ProfileSet

    if isinstance(result, ProfileSet):
        return result
    return ProfileSet({features[0]: result})


def _cmd_profile(args) -> int:
    features = args.feature
    try:
        _, profiles = _profile_input(args, features)
    except ValueError as exc:
        raise SystemExit(str(exc)) from None

    if args.seq:
        # One bare sequence goes to stdout: there is no file to name, and
        # piping it somewhere is the point. Still written through to_tsv,
        # not formatted here -- str(float("nan")) gives "nan" where to_tsv
        # writes an empty field, so a hand-rolled print would be a second,
        # untested serialisation of the same values.
        for feature, profile in profiles.items():
            if len(profiles) > 1:
                print(f"# {feature}")
            profile.to_tsv(sys.stdout)
        return 0

    source = Path(args.input)
    if args.outfile and len(profiles) > 1:
        raise SystemExit(
            f"--outfile names one file, but {len(profiles)} features were "
            "requested. Features have different k-mer lengths, so they cannot "
            "share a table. Drop --outfile to get one file per feature."
        )

    for feature, profile in profiles.items():
        out = args.outfile or _default_outfile(source, args.window_size, feature)
        profile.to_tsv(out)
        print(f"wrote {out}", file=sys.stderr)
    return 0


def _features_in(requests) -> list[str]:
    """Which lookup features a set of encode requests actually needs.

    `1-DNaseI` implies profiling `DNaseI`; `1-mer` implies nothing. Deriving
    this from the request means there is no second flag to keep in sync.
    """
    found = []
    for request in requests:
        head, _, what = request.partition("-")
        if head.isdigit() and what and what != "mer" and what not in found:
            found.append(what)
    return found


def _cmd_encode(args) -> int:
    import numpy as np

    features = _features_in(args.features)
    if not features:
        # A sequence-only request still needs a profile, because the profile
        # is what carries the sequences. Any feature will do; this one is
        # cheap and always present.
        features = ["gc"]

    try:
        _, profiles = _profile_input(args, features)
        matrix = profiles.encode(args.features, normalize=not args.no_normalize)
    except (ValueError, TypeError) as exc:
        raise SystemExit(str(exc)) from None

    payload = {
        "X": matrix.X,
        "columns": np.array(matrix.columns),
        "seqids": np.array(matrix.seqids),
        "window_size": np.array(matrix.window_size),
        "feature_names": np.array(matrix.feature_names),
    }
    if matrix.y is not None:
        payload["y"] = np.array(matrix.y, dtype=float)

    np.savez(args.out, **payload)
    rows, cols = matrix.shape
    labels = "with labels" if matrix.y is not None else "no labels"
    print(f"wrote {args.out}: {rows} x {cols}, {labels}", file=sys.stderr)
    return 0


def _add_input_arguments(parser) -> None:
    """Arguments shared by every subcommand that reads sequences."""
    parser.add_argument("input", nargs="?", help="FASTA, BED or table file")
    parser.add_argument("--seq", help="profile one sequence given inline, instead of a file")
    parser.add_argument("--format", choices=("fasta", "bed", "table"),
                        help="override the format guessed from the extension")
    parser.add_argument("--window-size", type=int, default=10,
                        help="sliding window in bases; 0 gives per-k-mer values [default: 10]")
    parser.add_argument("--threads", type=int, default=None,
                        help="worker processes; omit to decide automatically")
    parser.add_argument("--lookup", default=None,
                        help="a YAML feature table of your own; omit for the packaged one")

    bed = parser.add_argument_group("BED input")
    bed.add_argument("--genome", help="reference genome FASTA, required for BED")
    bed.add_argument("--width", type=int, default=None,
                     help="re-centre every interval and cut it to this many bases")
    bed.add_argument("--on-edge", choices=("drop", "error", "pad"), default="drop",
                     help="what to do when a window runs off a contig [default: drop]")

    table = parser.add_argument_group("table input")
    header = table.add_mutually_exclusive_group()
    header.add_argument("--header", dest="header", action="store_true", default=None,
                        help="row 1 holds column names")
    header.add_argument("--no-header", dest="header", action="store_false",
                        help="row 1 is data")
    table.add_argument("--seq-col", default=0, help="sequence column, by position or name")
    table.add_argument("--value-col", default=1, help="label column, by position or name")
    table.add_argument("--id-col", default=None, help="name column, by position or name")
    table.add_argument("--sep", default="\t", help="column separator [default: tab]")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="DNAflexpy",
        description="DNA flexibility profiling from sequence.",
    )
    parser.add_argument("--version", action="version",
                        version=f"DNAflexpy {__version__}")
    parser.add_argument("--citation", action="store_true",
                        help="print a BibTeX entry and exit")
    subparsers = parser.add_subparsers(dest="command")

    profile = subparsers.add_parser(
        "profile", help="profile sequences and write a TSV per feature")
    _add_input_arguments(profile)
    profile.add_argument("--feature", nargs="+", default=["DNaseI"],
                         help="one or more features [default: DNaseI]")
    profile.add_argument("--outfile",
                         help="output file; only valid with a single feature")
    profile.set_defaults(func=_cmd_profile)

    encode = subparsers.add_parser(
        "encode", help="build a machine-learning design matrix")
    _add_input_arguments(encode)
    encode.add_argument("--features", nargs="+", required=True,
                        help="blocks to build, e.g. 1-mer 1-DNaseI 2-stiffness")
    encode.add_argument("--out", default="X.npz", help="output .npz [default: X.npz]")
    encode.add_argument("--no-normalize", action="store_true",
                        help="skip min-max scaling of feature blocks")
    encode.set_defaults(func=_cmd_encode)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.citation:
        print(CITATION_TEMPLATE.format(version=__version__))
        print(
            "\nNote: author and title are placeholders. Fill them in before "
            "citing this.", file=sys.stderr,
        )
        return 0

    if not args.command:
        parser.print_help()
        return 1
    if not args.input and not args.seq:
        raise SystemExit(
            f"{args.command} needs an input file, or --seq with a sequence."
        )
    if args.input and args.seq:
        raise SystemExit("pass an input file or --seq, not both.")

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
