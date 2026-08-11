import pandas as pd
import multiprocessing
import logging
from typing import Optional
from .utils import read_fasta, get_kmer_len, load_feature_data, seq_to_numeric_profile

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def DNAflexpy(seqid:str, record:str, window_size:int, feature:str, feature_lookup: str):
    """
        Takes a sequence, 
        window size,
        kmer length,
        feature,
        feature_lookup
    """
    try:
        feature_lookup = load_feature_data()
        kmer_len = get_kmer_len(feature)
        feature_value = seq_to_numeric_profile(seqid, record, kmer_len, window_size, feature, feature_lookup)
        return feature_value
    
    except Exception as e:
        print(f"Error occured while processing sequence: {e}")


def DNAflexpy_for_CLI(seqid:str, record:str, kmer_len:int, window_size:int, feature:str, feature_lookup: str):
    """Takes a sequence, 
        window size,
        kmer length,
        feature,
        feature_lookup
    """
    try:
        feature_value = seq_to_numeric_profile(seqid, record, kmer_len, window_size, feature, feature_lookup)
        return feature_value
    
    except Exception as e:
        print(f"Error occured while processing sequence: {e}")


def DNAflexpyMP(input_file: str, 
                          window_size: int, 
                          feature: str,
                          threads: int,
                          outfile: Optional[str] = None) -> Optional[pd.DataFrame]:
    """
    Main function to calculate feature profilesß from a FASTA file.

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
        feature_lookup = load_feature_data()
        kmer_len = get_kmer_len(feature)

        sequence_data = list(read_fasta(input_file))  # List of tuples (id, seq)
        
        with multiprocessing.Pool(processes=threads) as pool:
            processed_results = pool.starmap(DNAflexpy_for_CLI, [(seqid, seq, kmer_len, window_size, feature, feature_lookup) for seqid, seq in sequence_data])

        df = pd.DataFrame(processed_results)
        
        if outfile:
            df.to_csv(outfile, index=False, header=False, sep="\t")
            logging.info("Process completed! Output written to %s", outfile)
        else:
            return df  # Return DataFrame if no output file is specified
    except Exception as e:
        logging.error("An error occurred: %s", e)
        raise e