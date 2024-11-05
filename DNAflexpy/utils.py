import yaml
from typing import Generator, Tuple, Dict

def read_fasta(filepath: str) -> Generator[Tuple[str, str], None, None]:
    """
    Parses a FASTA file and yields record names and sequences.
    Args:
        filepath (str): Path to the FASTA file.
    Yields:
        tuple: A tuple containing the record name and sequence.
    """
    try:
        with open(filepath, 'r') as file:
            name, sequence = None, []
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if name:
                        yield name, ''.join(sequence)
                    name, sequence = line[1:], []
                else:
                    sequence.append(line)
            if name:
                yield name, ''.join(sequence)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")
    except Exception as e:
        raise RuntimeError(f"An error occurred while reading the FASTA file: {e}")

def load_feature_data(prop: str="trinucleotide", feature_file: str="data/lookup.yaml") -> Dict[str, float]:
    """
    Loads feature data from a YAML file for a given k-mer length.

    Args:
        feature_file (str): Path to the YAML file containing feature data.
        kmer_len (int): Length of the k-mer.

    Returns:
        dict: A dictionary of feature data.
    """
    try:
        with open(feature_file, 'r') as f:
            data = yaml.safe_load(f)
            if str(prop) in data:
                return data[str(prop)]
            else:
                raise ValueError(f"k-mer length {kmer_len} not found in feature data.")
    except FileNotFoundError:
        raise FileNotFoundError(f"Feature file not found: {feature_file}")
    except yaml.YAMLError as e:
        raise RuntimeError(f"An error occurred while parsing the YAML file: {e}")
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading feature data: {e}")

def transform_seq_to_feat(sequence, kmer_len, feature, feature_lookup):
    """
    Calculates average structural profile values for a given sequence window.

    Args:
        sequence (str): DNA sequence to process.
        kmer_len (int): Length of k-mers.
        feature (str): Feature to calculate.
        feature_lookup (dict): Dictionary containing feature data.

    Returns:
        float: Average feature value for the given sequence window.
    """
    sequence = sequence.upper()
    feature_data = feature_lookup.get(feature, {})

    if not feature_data:
        print(f"**FATAL** ==>{feature}<== not in lookup data (.yaml file). Check spelling or add it.")
        return 0

    ls_values_w = []
    for i in range(len(sequence) - kmer_len + 1):
        subseq = sequence[i: i + kmer_len]
        value = feature_data.get(subseq)
        if value is not None:
            ls_values_w.append(value)
        else:
            print(f"**Warning** Subsequence {subseq} not found in the dictionary for {feature}")

    avg_w = sum(ls_values_w) / (len(sequence) - 1) if ls_values_w else 0
    
    return avg_w

def seq_to_numeric_profile(sequence, 
                        window_size, 
                        kmer_len, 
                        feature, 
                        feature_lookup):
    """
        Desc:
            - Operates on a sequence to aggregate average feature value of all the overlapping window

        arguments:
            - sequence
            - window size to make overlapping window
            - k-mer length
            - feature: Feature to calculate
            - feature_lookup: Dictionary containing feature data
            
        returns:
            - a list of values for a sequence

    """
    ls_window_avg_for_seq = []
    for w_start in range(len(sequence) - window_size + 1):
        current_w_seq = sequence[w_start : w_start + window_size]
        avg_w = transform_seq_to_feat(current_w_seq, kmer_len, feature, feature_lookup)
        ls_window_avg_for_seq.append(round(avg_w, 3))
    
    return ls_window_avg_for_seq

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
