import yaml
from typing import Generator, Tuple, Dict
from collections import defaultdict
from statistics import mean
from importlib.resources import files as ir_files

def seq_to_numeric_profile(seqid: str, sequence: str, kmer_len: int, window_size: int,
                           feature: str, feature_lookup: dict):
    """
    Operates on a sequence to aggregate average feature values over overlapping windows,
    or returns the per-position values when window_size == 0.
    """
    try:
        if window_size > 0:
            ls_window_avg_for_seq = [seqid]

            for w_start in range(len(sequence) - window_size + 1):
                current_w_seq = sequence[w_start: w_start + window_size]
                w_vals = transform_seq_to_feat(current_w_seq, kmer_len, feature, feature_lookup)
                avg_w = sum(w_vals) / len(w_vals) if w_vals else 0.0   # <-- fixed denominator
                ls_window_avg_for_seq.append(round(avg_w, 3))

            return ls_window_avg_for_seq

        elif window_size == 0:
            ls_for_seq = [seqid]
            per_pos_vals = transform_seq_to_feat(sequence, kmer_len, feature, feature_lookup)
            ls_for_seq.extend(per_pos_vals)
            return ls_for_seq

    except Exception as e:
        print(f"Error occurd while generating profile from sequence: {e}")


def transform_seq_to_feat(sequence: str, kmer_len: int, feature: str, feature_lookup: dict):
    """
    Returns a list of numeric values for each overlapping k-mer in `sequence`.
    Unknown/underscored k-mers (e.g., 'A_G') and any k-mer missing from the table
    receive the MEAN value of the selected feature’s lookup table.
    """
    try:
        sequence = sequence.upper()
        feature_data = feature_lookup.get(feature, {})
        if not feature_data:
            print(f"{feature} feature not in lookup data (.yaml file). Check the spelling or add it.")
            return []

        # default value = mean of the entire lookup table for this feature
        vals = list(feature_data.values())
        default_val = float(mean(vals)) if vals else 0.0

        # defaultdict returns default_val for any unknown k-mer (incl. ones with '_')
        lookup = defaultdict(lambda: default_val, feature_data)

        ls_values_w = []
        for i in range(len(sequence) - kmer_len + 1):
            subseq = sequence[i: i + kmer_len]
            # Any subseq not present (including containing '_') -> default mean
            value = lookup[subseq]
            ls_values_w.append(value)

        return ls_values_w

    except Exception as e:
        print(f"Error occured while transforming sequence to feature space: {e}")
        return []



# def process_sequence(seqid:str, record:str, window_size:int, feature:str, feature_lookup: str):
#     """
#         Takes a sequence, 
#         window size,
#         kmer length,
#         feature,
#         feature_lookup
#     """
#     try:
#         kmer_len = get_kmer_len(feature)
#         feature_value = seq_to_numeric_profile(seqid, record, kmer_len, window_size, feature, feature_lookup)
#         return feature_value
#     except Exception as e:
#         print(f"Error occured while processing sequence: {e}")



# def seq_to_numeric_profile(seqid:str, sequence:str, kmer_len:int, window_size:int, feature:str, feature_lookup:dict):
#     """
#         Desc:
#             - Operates on a sequence to aggregate average feature value of all the overlapping window

#         arguments:
#             - sequence
#             - window size to make overlapping window
#             - k-mer length
#             - feature: Feature to calculate
#             - feature_lookup: Dictionary containing feature data
            
#         returns:
#             - a list of values for a sequence

#     """
    
#     try:
#         if window_size > 0:
#             ls_window_avg_for_seq = [seqid]
        
#             for w_start in range(len(sequence) - window_size + 1):
#                 current_w_seq = sequence[w_start : w_start + window_size]
#                 w = transform_seq_to_feat(current_w_seq, kmer_len, feature, feature_lookup)
#                 avg_w = sum(w) / (len(sequence) - 1) if w else 0
#                 ls_window_avg_for_seq.append(round(avg_w, 3))

#             return ls_window_avg_for_seq

#         elif window_size == 0:
#             ls_for_seq = [seqid]
#             feature = transform_seq_to_feat(sequence, kmer_len, feature, feature_lookup)
#             for i in feature:
#                 ls_for_seq.append(i)
            
#             return ls_for_seq
        
#     except Exception as e:
#         print(f"Error occurd while generating profile from sequence: {e}")
 

# def transform_seq_to_feat(sequence, kmer_len, feature, feature_lookup):
#     """
#     Calculates average structural profile values for a given sequence window.

#     Args:
#         sequence (str): DNA sequence to process.
#         kmer_len (int): Length of k-mers.
#         feature (str): Feature to calculate.
#         feature_lookup (dict): Dictionary containing feature data.

#     Returns:
#         float: Average feature value for the given sequence window.
#     """
    
#     try:
#         sequence = sequence.upper()
#         feature_data = feature_lookup.get(feature, {})
    
#     except Exception as e:    
#         print(f"{feature} feature not in lookup data (.yaml file). Check the spelling or add it.")
#         return 0

    
    
#     try:
#         ls_values_w = []

#         for i in range(len(sequence) - kmer_len + 1):
#             subseq = sequence[i: i + kmer_len]
#             value = feature_data.get(subseq, 0)
#             # value = feature_data.get(subseq, 0)  # Default to 0 if missing
#             # ls_values_w.append(value)
#             if value is not None:
#                 ls_values_w.append(value)
#             else:
#                 print(f"**Warning** Subsequence {subseq} not found in the dictionary for {feature}")
    
#         return ls_values_w
    
#     except Exception as e:
#         print(f"Error occured while transforming sequence to feature space: {e}")



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


def get_kmer_len(feature:str):
    
    kmerlen_to_feature = {3:['DNaseI', 'NPP', 'bendabilityDNase', 'bendabilityConcensus'],
                          2: ['wedge', 'prop', 'twistDisp', 'stiffness', 'bendingStiffness', 'trx']}
        
    try:
        for kmer_len, features in kmerlen_to_feature.items():
            if feature in features:
                return int(kmer_len)
            
    except Exception as e:
        print(f"Error occured while determining appropriate kmer length for {feature}: {e}")



# def load_feature_data() -> Dict[str, float]:
#     """
#     Loads feature data from a YAML file for a given k-mer length.

#     Args:
#         feature_file (str): Path to the YAML file containing feature data.
#         kmer_len (int): Length of the k-mer.

#     Returns:
#         dict: A dictionary of feature data.
#     """
#     try:
#         yamlfilepath = files("DNAflexpy.data").joinpath("lookupNEW.yaml")

#         with open(yamlfilepath, 'r') as f:
#             return yaml.safe_load(f)
        
#     except FileNotFoundError:
#         raise FileNotFoundError(f"Feature file not found: {feature_file}")
        
#     except yaml.YAMLError as e:
#         raise RuntimeError(f"An error occurred while parsing the YAML file: {e}")

#     except Exception as e:
#         raise RuntimeError(f"An error occurred while loading feature data: {e}")


def load_feature_data() -> Dict[str, float]:
    """
    Loads feature data from the packaged YAML file.

    Returns:
        dict: Feature lookup dictionary.
    """
    yaml_resource = ir_files("DNAflexpy.data").joinpath("lookupNEW.yaml")

    try:
        with yaml_resource.open('r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
            return data or {}  # handle empty YAML gracefully
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Feature file not found: {yaml_resource}") from e
    except yaml.YAMLError as e:
        raise RuntimeError(f"Error parsing YAML from {yaml_resource}: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Error loading feature data from {yaml_resource}: {e}") from e