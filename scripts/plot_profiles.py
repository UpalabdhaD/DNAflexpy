#!/usr/bin/env python3
import argparse
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile


def read_profile(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if not parts:
            continue
        seqid = parts[0]
        values = [p for p in parts[1:] if p != ""]
        rows[seqid] = [float(v) for v in values]
    return rows


def aggregate_mean(profile_dict):
    lengths = [len(v) for v in profile_dict.values() if v]
    if not lengths:
        return []
    min_len = min(lengths)
    if len(set(lengths)) > 1:
        print("Warning: sequence profiles have different lengths; truncating to min length")
    means = []
    for idx in range(min_len):
        values = [v[idx] for v in profile_dict.values() if len(v) > idx]
        means.append(sum(values) / len(values))
    return means


def zscale(values):
    if not values:
        return values
    mean = statistics.mean(values)
    std = statistics.pstdev(values)
    if std == 0:
        return [0.0 for _ in values]
    return [(v - mean) / std for v in values]


def generate_random_fasta(path, count, length, seed):
    rng = random.Random(seed)
    alphabet = ["A", "C", "G", "T"]
    lines = []
    for i in range(1, count + 1):
        seq = "".join(rng.choice(alphabet) for _ in range(length))
        lines.append(f">random_seq_{i}")
        lines.append(seq)
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Plot DNAflexpy profiles from TSV outputs."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        help="List of DNAflexpy TSV outputs to plot",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        help="Optional labels matching --inputs order",
    )
    parser.add_argument(
        "--seqid",
        help="Sequence ID to plot (defaults to first row in first input)",
    )
    parser.add_argument(
        "--feature",
        help="Feature name for labeling the y-axis",
    )
    parser.add_argument(
        "--zscale",
        action="store_true",
        help="Plot z-scored values for each profile",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output image path (e.g., /tmp/plot.png)",
    )
    parser.add_argument(
        "--generate-random-fasta",
        help="Generate a random multi-FASTA file at the given path",
    )
    parser.add_argument(
        "--random-count",
        type=int,
        default=5,
        help="Number of random sequences to generate [default: 5]",
    )
    parser.add_argument(
        "--random-length",
        type=int,
        default=30,
        help="Length of each random sequence [default: 30]",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=1,
        help="Random seed for reproducibility [default: 1]",
    )
    parser.add_argument(
        "--window-sizes",
        nargs="+",
        type=int,
        help="Window sizes to run DNAflexpy on generated FASTA",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads to pass to DNAflexpy",
    )
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Error: matplotlib is required to plot ({exc})")
        sys.exit(1)

    inputs = args.inputs or []

    if args.generate_random_fasta:
        generate_random_fasta(
            args.generate_random_fasta,
            args.random_count,
            args.random_length,
            args.random_seed,
        )

        if not args.window_sizes or not args.feature:
            print("Error: --window-sizes and --feature are required with --generate-random-fasta")
            sys.exit(1)

        outdir = Path(tempfile.mkdtemp(prefix="dnaflexpy_plot_"))
        inputs = []
        for w in args.window_sizes:
            outpath = outdir / f"random_w{w}.tsv"
            cmd = [
                sys.executable,
                "-m",
                "DNAflexpy.cli",
                args.generate_random_fasta,
                "--window-size",
                str(w),
                "--feature",
                args.feature,
                "--outfile",
                str(outpath),
            ]
            if args.threads:
                cmd += ["--threads", str(args.threads)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(result.stderr.strip() or result.stdout.strip())
                print("Error: DNAflexpy CLI failed")
                sys.exit(result.returncode)
            inputs.append(str(outpath))

    if not inputs:
        print("Error: --inputs is required unless --generate-random-fasta is used")
        sys.exit(1)

    labels = args.labels or [
        f"w{w}" for w in (args.window_sizes or [])
    ] or [Path(p).stem for p in inputs]

    if len(labels) != len(inputs):
        print("Error: --labels must match the number of --inputs")
        sys.exit(1)

    profiles = [read_profile(p) for p in inputs]

    plt.figure(figsize=(10, 4))
    for label, profile in zip(labels, profiles):
        if args.seqid:
            if args.seqid not in profile:
                print(f"Warning: {args.seqid} not found in {label}, skipping")
                continue
            y = profile[args.seqid]
        else:
            y = aggregate_mean(profile)
        if args.zscale:
            y = zscale(y)
        x = list(range(1, len(y) + 1))
        plt.plot(x, y, label=label, linewidth=1.8)

    title_seq = args.seqid or "mean across sequences"
    plt.title(f"DNAflexpy profile: {title_seq}")
    plt.xlabel("Window index")
    if args.feature:
        ylabel = f"{args.feature} mean value"
    else:
        ylabel = "Mean feature value"
    if args.zscale:
        ylabel = f"Z-scored {ylabel}"
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)


if __name__ == "__main__":
    main()
