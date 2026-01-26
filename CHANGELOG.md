# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project adheres to
Semantic Versioning.

## [Unreleased]

## [0.2.0] - 2026-01-17

### Added
- Detailed README for publication use, including CLI/API usage and windowing behavior.
- Window-size behavior documentation and manual verification notes.
- BSD 3-Clause license text.
- Docs navigation entry for window size checks.

### Changed
- Window averaging now uses per-window k-mer counts (window size affects results).
- CLI default window size set to 10 and feature list clarified.

### Fixed
- Window-size edge cases now warn on `window_size > sequence length` and
  `window_size < kmer_len`.
