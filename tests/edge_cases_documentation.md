# DNAflexpy Edge Cases Documentation

**Date:** 2026-01-18
**Version:** 0.2.0

## Overview

This document comprehensively tests and documents the behavior of DNAflexpy under various edge case conditions. Edge cases are critical to understand for proper usage and to ensure the package handles unusual inputs gracefully.

---

## Edge Case 1: Window Size Smaller Than K-mer Length

### Test Setup
- **Sequences:** Various lengths (2-26 nucleotides)
- **Feature:** DNaseI (k-mer length = 3)
- **Window sizes tested:** 1, 2

### Behavior

When `window_size < kmer_len`, no complete k-mers can fit within any window, so the package:

1. **Issues a warning** for each affected sequence:
   ```
   **Warning** window_size (1) < kmer_len (3) for sequence_name; windows will yield 0.0 averages.
   ```

2. **Returns 0.0 values** for all windows:
   ```
   short_seq    0.0    0.0    0.0    0.0
   normal_seq   0.0    0.0    0.0    0.0    0.0    0.0    ...
   ```

3. **Generates the correct number of windows** based on the formula: `num_windows = len(sequence) - window_size + 1`

### Example Output (window-size 1, DNaseI)

| Sequence | Length | Windows Created | Output |
|----------|--------|----------------|---------|
| very_short | 2 | 2 | `very_short  0.0  0.0` |
| short_seq | 4 | 4 | `short_seq  0.0  0.0  0.0  0.0` |
| medium_seq | 12 | 12 | `medium_seq  0.0  0.0  0.0  ...` (12 zeros) |
| normal_seq | 26 | 26 | `normal_seq  0.0  0.0  0.0  ...` (26 zeros) |

### Rationale

This behavior is mathematically correct because:
- A window of size 1 with k-mer length 3 means we need 3 nucleotides but only have 1
- No k-mers can be extracted → empty list → average of empty list = 0.0

### User Guidance

✅ **Expected behavior** - Package warns users about this condition
⚠️ **User action:** Choose `window_size >= kmer_len` for meaningful results

---

## Edge Case 2: Window Size Larger Than Sequence Length

### Test Setup
- **Sequences:** Various lengths (2-26 nucleotides)
- **Feature:** DNaseI (k-mer length = 3)
- **Window size:** 30 (larger than all test sequences)

### Behavior

When `window_size > len(sequence)`, no windows can be created, so the package:

1. **Issues a warning** for each affected sequence:
   ```
   **Warning** window_size (30) > sequence length (4) for short_seq; no windows generated.
   ```

2. **Returns only the sequence ID** with no data values:
   ```
   short_seq
   normal_seq
   very_short
   medium_seq
   ```

3. **Does not crash or produce errors** - handles gracefully

### Example Output (window-size 30, DNaseI)

```
short_seq
normal_seq
very_short
medium_seq
```

Each line contains only the sequence identifier, no numeric values.

### Mixed Scenario (window-size 10, mixed sequence lengths)

When processing multiple sequences with varying lengths:

```
short_seq          (length 4 < 10) → only ID
normal_seq  0.011  (length 26 > 10) → 17 values
very_short         (length 2 < 10) → only ID
medium_seq  0.011  (length 12 > 10) → 3 values
```

### Rationale

Cannot create a window larger than the sequence itself. The package appropriately:
- Warns the user
- Returns the sequence ID for record keeping
- Continues processing other sequences

### User Guidance

✅ **Expected behavior** - Package handles gracefully with warnings
⚠️ **User action:** Ensure `window_size <= sequence_length` for all sequences

---

## Edge Case 3: Window Size Equal to Sequence Length

### Test Setup
- **Sequences:**
  - `exactly_3nt`: ATG (length 3)
  - `sequence1`: 26 nucleotides
- **Features:** DNaseI (k-mer length = 3)
- **Window sizes:** 3, 26 (equal to sequence lengths)

### Behavior

When `window_size == len(sequence)`, exactly **one window** is created that spans the entire sequence:

1. **Number of windows:** 1 (from formula: `26 - 26 + 1 = 1`)

2. **Window content:** All k-mers from the entire sequence

3. **Output:** Single averaged value representing the entire sequence

### Examples

#### Small Sequence (length 3, window 3, k-mer 3)
```
exactly_3nt    0.134
```
- 1 window containing 1 k-mer (3-3+1=1)
- Average of 1 value = the value itself

#### Normal Sequence (length 26, window 26, k-mer 3)
```
sequence1    0.019
sequence2    0.0
```
- 1 window containing 24 k-mers (26-3+1=24)
- Single value is the average of all 24 k-mers

### Comparison to Window-Size 0

| Setting | Output for 26nt sequence |
|---------|-------------------------|
| window-size 0 | 24 individual k-mer values |
| window-size 26 | 1 value (average of all 24 k-mers) |

### Rationale

This is the **extreme case of window averaging** where the window encompasses the entire sequence, yielding a single global average.

### User Guidance

✅ **Expected behavior** - Produces one value summarizing the entire sequence
💡 **Use case:** When you want a single metric for each sequence

---

## Edge Case 4: Sequences Shorter Than K-mer Length

### Test Setup
- **Sequences:**
  - `very_short`: AT (length 2)
  - `exactly_2nt`: AT (length 2)
- **Features:**
  - DNaseI (k-mer length = 3)
  - trx (k-mer length = 2)

### Behavior: With Tri-nucleotide Feature (DNaseI, k-mer=3)

When `len(sequence) < kmer_len`, no k-mers can be extracted:

```
very_short
exactly_2nt
```

**Output:** Only sequence ID, no values

### Behavior: With Di-nucleotide Feature (trx, k-mer=2)

When `len(sequence) == kmer_len`:

```
very_short    0.0
exactly_2nt   0.0
```

**Output:** Sequence ID + exactly 1 value (the single k-mer that can be extracted)

### Comparison Table

| Sequence Length | K-mer Length | Can Extract K-mers? | Output |
|----------------|--------------|---------------------|---------|
| 2 | 3 | No | ID only |
| 2 | 2 | Yes (1 k-mer) | ID + 1 value |
| 3 | 3 | Yes (1 k-mer) | ID + 1 value |
| 4 | 3 | Yes (2 k-mers) | ID + 2 values |

### Rationale

Number of k-mers = `len(sequence) - kmer_len + 1`
- If result ≤ 0, no k-mers can be extracted
- If result > 0, that many k-mers can be extracted

### User Guidance

✅ **Expected behavior** - Package handles gracefully
⚠️ **User action:** Ensure sequences are at least as long as the k-mer length of your chosen feature
💡 **Tip:** Use di-nucleotide features (k-mer=2) for very short sequences

---

## Edge Case 5: Different Features (Di-nucleotide vs Tri-nucleotide)

### Test Setup
- **Sequences:** Various lengths (2-26 nucleotides)
- **Features:**
  - DNaseI (tri-nucleotide, k-mer = 3)
  - trx (di-nucleotide, k-mer = 2)
- **Window sizes:** 0, 1, 2

### Behavior Differences

#### Tri-nucleotide Feature (DNaseI, k-mer=3)

| Window Size | Behavior |
|------------|----------|
| 0 | Per k-mer values (seq_len - 3 + 1 values) |
| 1 | ⚠️ Warning: 0.0 averages (window < k-mer) |
| 2 | ⚠️ Warning: 0.0 averages (window < k-mer) |
| 3 | ✅ Same as window 0 (window = k-mer) |

#### Di-nucleotide Feature (trx, k-mer=2)

| Window Size | Behavior |
|------------|----------|
| 0 | Per k-mer values (seq_len - 2 + 1 values) |
| 1 | ⚠️ Warning: 0.0 averages (window < k-mer) |
| 2 | ✅ Same as window 0 (window = k-mer) |

### Example: 4-nucleotide Sequence (ATGC)

#### With DNaseI (k-mer=3)
```
window 0: ATGC → [ATG, TGC] → 2 values
window 1: ATGC → 0.0, 0.0, 0.0, 0.0 (4 windows, all zeros)
window 2: ATGC → 0.0, 0.0, 0.0 (3 windows, all zeros)
window 3: ATGC → [ATG, TGC] → 2 values (average) → 1 value
window 4: ATGC → average of 2 k-mers → 1 value
```

#### With trx (k-mer=2)
```
window 0: ATGC → [AT, TG, GC] → 3 values (0.0, 42.0, 25.0)
window 1: ATGC → 0.0, 0.0, 0.0, 0.0 (4 windows, all zeros)
window 2: ATGC → [AT, TG, GC] → 3 values (0.0, 42.0, 25.0)
```

### Very Short Sequence (AT, length 2)

| Feature | K-mer | Output |
|---------|-------|--------|
| DNaseI | 3 | `very_short` (ID only) |
| trx | 2 | `very_short  0.0` (ID + 1 value) |

### User Guidance

✅ **Key insight:** Di-nucleotide features work with shorter sequences
💡 **Recommendation:**
- Use tri-nucleotide features (DNaseI, NPP) for sequences ≥3 nt
- Use di-nucleotide features (trx, twistDisp, stiffness) for sequences ≥2 nt

---

## Edge Case 6: Window Size = K-mer Length (Critical Test)

### Test Setup
- **Features:** Both DNaseI (k-mer=3) and trx (k-mer=2)
- **Window sizes:** Set equal to k-mer length

### Behavior

When `window_size == kmer_len`, each window contains **exactly 1 k-mer**:

#### DNaseI (k-mer=3, window=3)
```
window 0: 0.134  0.076  -0.077  -0.033  0.025  ...
window 3: 0.134  0.076  -0.077  -0.033  0.025  ...
```
**Result:** Identical outputs ✅

#### trx (k-mer=2, window=2)
```
window 0: 0.0  42.0  25.0  ...
window 2: 0.0  42.0  25.0  ...
```
**Result:** Identical outputs ✅

### Mathematical Proof

For a window of size W containing k-mers of length K:
- Number of k-mers in window = `W - K + 1`
- When `W = K`: `K - K + 1 = 1`
- Average of 1 value = the value itself

Therefore: `window(K) output = window(0) output` ✅

### Significance

This is a **critical validation** of the averaging formula fix:
- Confirms `avg = sum(w) / len(w)` is correct
- Proves the old bug `avg = sum(w) / (len(sequence) - 1)` would fail this test

### User Guidance

✅ **Expected behavior** - This confirms the bug fix is working correctly
💡 **Testing tip:** This is the easiest way to verify the package is working correctly

---

## Summary Table: All Edge Cases

| Edge Case | Condition | Behavior | Output | Warning |
|-----------|-----------|----------|--------|---------|
| Window < K-mer | `W < K` | Windows created but all 0.0 | ID + zeros | ✅ Yes |
| Window > Sequence | `W > L` | No windows created | ID only | ✅ Yes |
| Window = Sequence | `W = L` | One window (entire seq) | ID + 1 value | No |
| Sequence < K-mer | `L < K` | No k-mers extractable | ID only | No |
| Sequence = K-mer | `L = K` | Exactly 1 k-mer | ID + 1 value | No |
| Window = K-mer | `W = K` | 1 k-mer per window | Same as W=0 | No |
| Very short + di-nuc | `L=2, K=2` | 1 k-mer extractable | ID + 1 value | No |
| Very short + tri-nuc | `L=2, K=3` | No k-mers | ID only | No |

**Legend:**
- W = window size
- K = k-mer length
- L = sequence length

---

## Recommendations for Users

### ✅ Good Practices

1. **Match window size to analysis goals:**
   - `window = 0`: Per-k-mer resolution (maximum detail)
   - `window = kmer_len`: Same as window 0 but validates implementation
   - `window > kmer_len`: Smoothed profile (less noise)
   - `window = sequence_len`: Single summary value per sequence

2. **Choose appropriate features for sequence length:**
   - Sequences ≥3 nt: Any feature
   - Sequences 2 nt: Di-nucleotide features only (trx, twistDisp, stiffness)
   - Sequences <2 nt: No features will work

3. **Pre-filter sequences:**
   - Remove sequences shorter than your feature's k-mer length
   - Or use `--window-size 0` to identify which sequences are too short

4. **Understand warnings:**
   - Warnings are informational, not errors
   - Package continues processing other sequences
   - Check output to see which sequences were affected

### ⚠️ Avoid These Pitfalls

1. ❌ **Don't use** `window_size < kmer_len` expecting meaningful results
   - You'll get all zeros

2. ❌ **Don't ignore warnings** about window size vs sequence length
   - Those sequences won't have useful data

3. ❌ **Don't use tri-nucleotide features** on very short sequences (<3 nt)
   - Use di-nucleotide features instead

4. ❌ **Don't expect output** from sequences shorter than k-mer length
   - These will only have sequence IDs

### 💡 Advanced Use Cases

1. **Quality control:** Use `window_size = sequence_length` to get a single quality metric per sequence

2. **Multi-scale analysis:** Run multiple window sizes (0, 5, 10, 20) to analyze features at different scales

3. **Short sequence handling:** Use di-nucleotide features (k-mer=2) for maximum compatibility with short sequences

---

## Testing Validation

All edge cases were tested and documented. The package handles edge cases gracefully:

✅ Appropriate warnings are issued
✅ No crashes or errors
✅ Mathematically correct behavior
✅ Continues processing despite edge cases
✅ Output format remains consistent

## Files Generated During Testing

- `test_edge_cases.fasta` - Test sequences of various lengths
- `test_exact_kmer.fasta` - Sequences exactly k-mer length
- `edge_case_outputs/` - Directory containing all test outputs

---

**Conclusion:** DNAflexpy handles edge cases robustly and predictably. Users should be aware of these behaviors to interpret results correctly and choose appropriate parameters for their analysis.
