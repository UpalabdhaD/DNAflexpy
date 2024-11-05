import unittest
from utils import calculate_window_averages




class TestCalculateWindowAverages(unittest.TestCase):

    def setUp(self):
        self.feature_lookup = {
            "example_feature": {
                "AT": 1.0,
                "TC": 2.0,
                "CG": 3.0,
                "GA": 4.0
            }
        }

    def test_simple_sequence(self):
        sequence = "ATCGATCG"
        window_size = 4
        kmer_len = 2
        feature = "example_feature"
        expected_output = [2.5, 2.5, 2.5, 2.5, 2.5]
        result = calculate_window_averages(sequence, window_size, kmer_len, feature, self.feature_lookup)
        self.assertEqual(result, expected_output)

    def test_no_matching_kmers(self):
        sequence = "ATCGATCG"
        window_size = 4
        kmer_len = 2
        feature = "example_feature"
        feature_lookup = {"example_feature": {"GG": 1.0, "CC": 2.0}}
        expected_output = [0, 0, 0, 0, 0]
        result = calculate_window_averages(sequence, window_size, kmer_len, feature, feature_lookup)
        self.assertEqual(result, expected_output)

    def test_empty_sequence(self):
        sequence = ""
        window_size = 4
        kmer_len = 2
        feature = "example_feature"
        expected_output = []
        result = calculate_window_averages(sequence, window_size, kmer_len, feature, self.feature_lookup)
        self.assertEqual(result, expected_output)

    def test_sequence_shorter_than_window(self):
        sequence = "ATC"
        window_size = 4
        kmer_len = 2
        feature = "example_feature"
        expected_output = []
        result = calculate_window_averages(sequence, window_size, kmer_len, feature, self.feature_lookup)
        self.assertEqual(result, expected_output)

if __name__ == '__main__':
    unittest.main()