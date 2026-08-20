import argparse
import multiprocessing
import os
from .core import DNAflexpyMP


CITATION_BIBTEX = """@software{DNAflexpy,
  author = {ADD_AUTHOR_LIST},
  title = {ADD_TITLE},
  year = {ADD_YEAR},
  publisher = {GitHub},
  version = {ADD_VERSION},
  url = {https://github.com/upalabdhaD/DNAflexpy}
}"""



def main():
    """
    Command-line interface for the feature profiling tool.
    """
    parser = argparse.ArgumentParser(description="Process a multifasta file and calculate bendability of DNA sequence")

    parser.add_argument("input_file", nargs="?", help="Path to the input multifasta file")
    parser.add_argument("--window-size", type=int, default=10, help="Size of the processing window [default: 10]")
    parser.add_argument(
        "--feature",
        default="DNaseI",
        help=(
            "Feature to calculate ('DNaseI', 'NPP', 'twistDisp', 'trx', 'stiffness') "
            "[default: DNaseI]"
        ),
    )
    parser.add_argument("--threads", type=int, default= multiprocessing.cpu_count(), help="Number of threads")
    parser.add_argument("--outfile", help="Output file name [optional]")
    parser.add_argument(
        "--citation",
        action="store_true",
        help="Print BibTeX citation and exit",
    )
    args = parser.parse_args()

    if args.citation:
        print(CITATION_BIBTEX)
        return

    if not args.input_file:
        parser.error("the following arguments are required: input_file")

    # Determine output file name if not specified
    out_base_name = os.path.splitext(os.path.basename(args.input_file))[0]
    output_filename = args.outfile or f"{out_base_name}_w{args.window_size}nt_{args.feature}.tsv"

    # Run the calculation
    DNAflexpyMP(
        input_file=args.input_file,
        window_size=args.window_size,
        feature=args.feature,
        threads=args.threads,
        outfile=output_filename
    )

if __name__ == "__main__":
    main()
