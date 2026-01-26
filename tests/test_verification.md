# DNAflexpy Window-Size Functionality Verification

**Date:** 2026-01-18
**Version Tested:** 0.2.0
**Test File:** DNAflexpy/data/test_fasta.fa

## Test Overview

Manual verification was performed to confirm that DNAflexpy produces the intended results when using different window sizes on the same FASTA file. This test validates that the window-size bug fix is working correctly.

## Test Sequences

The test file contains two sequences:
- **sequence1:** ATGCGTACGTAGCTAGCGTAGCTAGT (26 nucleotides)
- **sequence2:** CGTAGCTAGTACGATCGTACGTAGCT (26 nucleotides)

## Test Parameters

- **Feature:** DNaseI (k-mer length = 3)
- **Window sizes tested:** 0, 3, 10, 15
- **Threads:** Default (all cores)

## Commands Executed

```bash
DNAflexpy DNAflexpy/data/test_fasta.fa --window-size 0 --feature DNaseI --outfile test_outputs/output_w0.tsv
DNAflexpy DNAflexpy/data/test_fasta.fa --window-size 3 --feature DNaseI --outfile test_outputs/output_w3.tsv
DNAflexpy DNAflexpy/data/test_fasta.fa --window-size 10 --feature DNaseI --outfile test_outputs/output_w10.tsv
DNAflexpy DNAflexpy/data/test_fasta.fa --window-size 15 --feature DNaseI --outfile test_outputs/output_w15.tsv
```

## Results Summary

### Output Dimensions

| Window Size | Columns | Data Points per Sequence | Expected | Status |
|-------------|---------|-------------------------|----------|---------|
| 0           | 25      | 24                      | 26-3+1=24 | ✅ Correct |
| 3           | 25      | 24                      | 26-3+1=24 | ✅ Correct |
| 10          | 18      | 17                      | 26-10+1=17 | ✅ Correct |
| 15          | 13      | 12                      | 26-15+1=12 | ✅ Correct |

### Sample Values Comparison (sequence1, first 5 data points)

| Window Size | Value 1 | Value 2 | Value 3 | Value 4 | Value 5 | Interpretation |
|-------------|---------|---------|---------|---------|---------|----------------|
| 0           | 0.134   | 0.076   | -0.077  | -0.033  | 0.025   | Per k-mer values (no averaging) |
| 3           | 0.134   | 0.076   | -0.077  | -0.033  | 0.025   | Same as w0 (1 k-mer per window) |
| 10          | 0.011   | -0.003  | -0.001  | 0.010   | 0.017   | Averaged over 8 k-mers/window |
| 15          | 0.025   | 0.021   | 0.017   | 0.017   | 0.017   | Averaged over 13 k-mers/window |

## Detailed Analysis

### Window-Size 0 (No Windowing)
```
sequence1: 0.134  0.076  -0.077  -0.033  0.025  0.025  -0.033  -0.033  0.025  0.09  ...
sequence2: -0.033  0.025  0.09  0.017  0.017  0.09  0.09  -0.183  0.025  0.025  ...
```

- Returns raw per-k-mer feature values
- 24 values per sequence (26 nucleotides - 3 k-mer length + 1)
- Values range from -0.183 to 0.134
- **Behavior:** ✅ Expected - returns individual k-mer values without averaging

### Window-Size 3 (Window = K-mer Length)
```
sequence1: 0.134  0.076  -0.077  -0.033  0.025  0.025  -0.033  -0.033  0.025  0.09  ...
sequence2: -0.033  0.025  0.09  0.017  0.017  0.09  0.09  -0.183  0.025  0.025  ...
```

- **CRITICAL TEST:** Window size equals k-mer length
- Each window contains exactly 1 k-mer (3-3+1=1)
- Average of 1 value = the value itself
- Results are **identical** to window-size 0
- **Behavior:** ✅ Expected - confirms averaging formula is correct

### Window-Size 10 (Small Window)
```
sequence1: 0.011  -0.003  -0.001  0.010  0.017  0.025  0.033  0.039  0.034  0.026  ...
sequence2: 0.014  0.021  0.021  0.006  0.004  -0.012  -0.037  -0.049  -0.03  -0.03  ...
```

- Creates 17 windows per sequence (26-10+1=17)
- Each window contains 8 k-mers (10-3+1=8)
- Values are averaged over 8 k-mers
- Results show moderate smoothing
- Value ranges are narrower than window-size 0/3
- **Behavior:** ✅ Expected - shows proper averaging over multiple k-mers

### Window-Size 15 (Large Window)
```
sequence1: 0.025  0.021  0.017  0.017  0.017  0.017  0.022  0.026  0.030  0.035  ...
sequence2: 0.001  -0.005  -0.007  -0.016  -0.016  -0.015  -0.024  -0.034  -0.018  -0.013  ...
```

- Creates 12 windows per sequence (26-15+1=12)
- Each window contains 13 k-mers (15-3+1=13)
- Values are averaged over 13 k-mers
- Results show stronger smoothing than window-size 10
- Value ranges are even narrower
- **Behavior:** ✅ Expected - shows increased smoothing with larger windows

## Mathematical Verification

### Formula Used (Correct)
```python
avg_w = (sum(w) / len(w)) if w else 0.0
```

Where:
- `w` is the list of k-mer feature values within the current window
- `len(w)` = window_size - kmer_len + 1

### Example Calculation for Window-Size 10

For the first window of sequence1 (positions 0-9: "ATGCGTACGT"):

**K-mers in window:**
1. ATG → 0.134
2. TGC → 0.076
3. GCG → -0.077
4. CGT → -0.033
5. GTA → 0.025
6. TAC → 0.025
7. ACG → -0.033
8. CGT → -0.033

**Calculation:**
```
sum = 0.134 + 0.076 + (-0.077) + (-0.033) + 0.025 + 0.025 + (-0.033) + (-0.033) = 0.084
len(w) = 8
avg = 0.084 / 8 = 0.0105 ≈ 0.011 (rounded to 3 decimals)
```

**Result:** 0.011 ✅ Matches output!

## Comparison to Original Bug

### Original (Incorrect) Formula
```python
avg_w = sum(w) / (len(sequence) - 1) if w else 0
```

This would give:
```
avg = 0.084 / (26 - 1) = 0.084 / 25 = 0.00336
```

This is **very different** from the correct value of 0.011, and would make all window sizes produce nearly identical (and incorrect) results.

### Bug Fix Validation

The test results confirm:
1. ✅ Window-size 0 and 3 are identical (correct, since each w3 window has 1 k-mer)
2. ✅ Window-size 10 produces averaged values (not scaled-down versions of w0)
3. ✅ Window-size 15 produces more smoothed values than w10
4. ✅ Output dimensions match mathematical expectations
5. ✅ Manual calculation matches actual output

## Smoothing Effect Visualization

The smoothing effect can be clearly seen:

**Sequence1 Value Ranges:**
- Window-size 0/3: -0.183 to 0.134 (range: 0.317)
- Window-size 10: -0.003 to 0.039 (range: 0.042)
- Window-size 15: 0.017 to 0.035 (range: 0.018)

As window size increases:
- Range decreases (more smoothing)
- Extreme values are averaged out
- Profile becomes smoother and less variable

## Edge Case Handling

### Window Size Smaller Than K-mer Length
Not tested in this verification, but documented behavior:
- Window size < 3 would result in 0.0 averages (no k-mers fit in window)

### Window Size Larger Than Sequence
Not applicable for these test sequences, but documented behavior:
- Returns only sequence ID with warning

## Conclusion

✅ **ALL TESTS PASSED**

The DNAflexpy package produces the intended and mathematically correct results with different window sizes. The window-size bug fix is working correctly:

1. **Window size 0:** Returns per-k-mer values (no averaging)
2. **Window size = k-mer length:** Returns same values as window size 0 (correct)
3. **Window size > k-mer length:** Returns properly averaged values with increasing smoothing as window size increases

The fix ensures that:
- Each window is averaged correctly over the k-mers it contains
- Different window sizes produce appropriately different results
- The mathematical formula `avg = sum(values) / count(values)` is correctly implemented
- Output dimensions match the expected formula: `num_windows = len(sequence) - window_size + 1`

## Recommendations

1. ✅ The package is ready for use with confidence in window-size functionality
2. ✅ Users can now vary window size to achieve desired smoothing levels
3. ✅ The bug that prevented window-size changes from having effect has been eliminated
4. Consider adding automated tests to prevent regression of this bug
5. Consider adding visual output examples in documentation showing smoothing effects

## Test Environment

- **OS:** Linux 6.8.0-79-generic
- **Python:** 3.12
- **DNAflexpy Version:** 0.2.0
- **Installation Method:** Development mode (`pip install -e .`)
