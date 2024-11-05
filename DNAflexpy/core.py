import pandas as pd

import multiprocessing
import logging
from typing import Optional
from .utils import read_fasta, load_feature_data, seq_to_numeric_profile

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_sequence(record, 
                    window_size, 
                    kmer_len, 
                    feature, 
                    feature_lookup):
    """
        Takes a sequence, 
        window size,
        kmer length,
        feature,
        feature_lookup
        
    """

    feature_value = seq_to_numeric_profile(record, window_size, kmer_len, feature, feature_lookup)
    return feature_value


def calculate_bendability(input_file: str, window_size: int, prop: str, kmer_len: int, feature: str, feature_file: str, threads: int, outfile: Optional[str] = None) -> Optional[pd.DataFrame]:
    """
    Main function to calculate feature profiles from a FASTA file.

    Args:
        input_file (str): Path to the input multifasta file.
        window_size (int): Size of the processing window.
        prop (int): Dinucleotide/Trinucleotide.
        feature (str): Feature to calculate.
        feature_file (str): Path to the YAML file containing feature data.
        threads (int): Number of threads to use.
        outfile (str, optional): Output file name.

    Returns:
        Optional[pd.DataFrame]: Feature profiles if outfile is not specified, otherwise None.
    """
    try:
        # Load feature data
        feature_lookup = load_feature_data(prop, feature_file)
        
        # Parse the fasta file and process sequences
#         records = list(read_fasta(input_file))
        with multiprocessing.Pool(processes=threads) as pool:
            results = pool.starmap(process_sequence, [(seq, window_size, kmer_len, feature, feature_lookup) for _, seq in read_fasta(input_file)])

        # Create DataFrame from results
        df = pd.DataFrame(results)
        
        if outfile:
            df.to_csv(outfile, index=False, header=False, sep="\t")
            logging.info("Process completed! Output written to %s", outfile)
        else:
            return df  # Return DataFrame if no output file is specified
    except Exception as e:
        logging.error("An error occurred: %s", e)
        raise e
    

# def calculate_features(input_file: str, window_size: int, kmer_len: int, feature: str, feature_file: str, threads: int, outfile: Optional[str] = None) -> Optional[pd.DataFrame]:
#     """
#     Main function to calculate feature profiles from a FASTA file.

#     Args:
#         input_file (str): Path to the input multifasta file.
#         window_size (int): Size of the processing window.
#         kmer_len (int): Length of k-mers.
#         feature (str): Feature to calculate.
#         feature_file (str): Path to the YAML file containing feature data.
#         threads (int): Number of threads to use.
#         outfile (str, optional): Output file name.

#     Returns:
#         Optional[pd.DataFrame]: Feature profiles if outfile is not specified, otherwise None.
#     """
#     try:
#         # Load feature data
#         feature_lookup = load_feature_data(feature_file, kmer_len)
        
#         # Parse the fasta file and process sequences
#         records = list(read_fasta(input_file))
#         with multiprocessing.Pool(processes=threads) as pool:
#             results = pool.starmap(process_sequence, [(seq, window_size, kmer_len, feature, feature_lookup) for _, seq in records])

#         # Create DataFrame from results
#         df = pd.DataFrame(results)
#         if outfile:
#             df.to_csv(outfile, index=False, header=False, sep="\t")
#             logging.info("Process completed! Output written to %s", outfile)
#         else:
#             return df  # Return DataFrame if no output file is specified
#     except Exception as e:
#         logging.error("An error occurred: %s", e)
#         raise e
    
    
# import pandas as pd
# import multiprocessing
# from .utils import read_fasta, load_feature_data, process_sequence


# def calculate_features(input_file, window_size, kmer_len, feature, feature_file, threads, outfile=None):
#     """
#     Main function to calculate feature profiles from a FASTA file.

#     Args:
#         input_file (str): Path to the input multifasta file.
#         window_size (int): Size of the processing window.
#         kmer_len (int): Length of k-mers.
#         feature (str): Feature to calculate.
#         feature_file (str): Path to the YAML file containing feature data.
#         threads (int): Number of threads to use.
#         outfile (str, optional): Output file name.

#     Returns:
#         DataFrame: Feature profiles if outfile is not specified, otherwise None.
#     """
#     # Load feature data
#     feature_lookup = load_feature_data(feature_file, kmer_len)
    
#     # Parse the fasta file and process sequences
#     records = list(read_fasta(input_file))
#     with multiprocessing.Pool(processes=threads) as pool:
#         results = pool.starmap(process_sequence, [(seq, window_size, kmer_len, feature, feature_lookup) for _, seq in records])

#     # Create DataFrame from results
#     df = pd.DataFrame(results)
#     if outfile:
#         df.to_csv(outfile, index=False, header=False, sep="\t")
#         print("Process completed! Output written to", outfile)
#     else:
#         return df  # Return DataFrame if no output file is specified
