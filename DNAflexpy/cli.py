import argparse
import multiprocessing
import os
from .core import DNAflexpyMP



def main():
    """
    Command-line interface for the feature profiling tool.
    """
    parser = argparse.ArgumentParser(description="Process a multifasta file and calculate bendability of DNA sequence")

    parser.add_argument("input_file", help="Path to the input multifasta file")
    parser.add_argument("--window-size", type=int, default=0, help="Size of the processing window [default: 10]")
    parser.add_argument("--feature", default="DNaseI", help="Feature(s) to calculate ('DNaseI', 'NPP') [default: DNaseI]")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads")
    parser.add_argument("--outfile", help="Output file name [optional]")
    parser.add_argument("--feature-file", default="data/lookupNEW.yaml", help="Path to the YAML file containing feature data")

    args = parser.parse_args()

    # Determine output file name if not specified
    out_base_name = os.path.splitext(os.path.basename(args.input_file))[0]
    output_filename = args.outfile or f"{out_base_name}_w{args.window_size}nt_{args.feature}.tsv"

    # Run the calculation
    DNAflexpyMP(
        input_file=args.input_file,
        window_size=args.window_size,
        feature=args.feature,
        feature_file=args.feature_file,
        threads=args.threads,
        outfile=output_filename
    )

if __name__ == "__main__":
    main()
