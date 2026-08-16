# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **This repo contains a new, class-based `DNAflexpy` package plus a frozen archive of the original.**
>
> `DNAflexpy/` is the new package: rewritten from scratch, class-based, and the one to import in new code. `rxv/DNAflexpy/` is the original package, moved there verbatim and **frozen — never change its logic**. It must keep running because `tests/test_differential.py` byte-compares the new package's output against it; if the archive stops working, or stops reading its own data, that gate stops meaning anything. The only sanctioned edit to the archive is the forced import-path fix at `rxv/DNAflexpy/utils.py:175` (`files("DNAflexpy.data")` → `files("rxv.DNAflexpy.data")`), required by the move. See [the design spec](docs/superpowers/specs/2026-08-12-dnaflex-rewrite-design.md) and [the phase plan](docs/superpowers/plans/2026-08-12-dnaflexpy-phases-0-2.md) for the full rationale.
>
> Active work happens on the `dev` branch.

## Commands

```bash
pip install -e .                    # editable install; see "No console script yet" below
python -m pytest -q                 # 319 passed (see "Testing" below)
python -m pytest tests/test_differential.py -k byte_identical -q   # just the byte-equality gate

mkdocs serve                        # docs preview; mkdocs build to render
python scripts/plot_profiles.py --generate-random-fasta /tmp/r.fa \
  --window-sizes 0 5 10 --feature DNaseI --out /tmp/plot.png   # needs matplotlib (not a declared dep)
```

**No console script yet.** There is no `DNAflexpy` CLI entry point today — `DNAflexpy.cli` does not exist. Phase 8 of the rewrite restores it. Until then the archive's CLI still works: `python -m rxv.DNAflexpy.cli <fasta> --window-size 10 --feature DNaseI --outfile out.tsv`.

## Testing

Tests run and pass: `python -m pytest -q` reports **319 passed**. This includes two byte-equality matrices in `tests/test_differential.py`:

- 210 cases (3 FASTAs x 10 features x 7 window sizes), run with `threads=1`, comparing `FlexProfiler(...).profile_fasta(...).to_tsv(...)` byte-for-byte against the archive's `seq_to_numeric_profile` output.
- 20 cases forcing the pooled `multiprocessing.Pool` code path (10 features x 2 window sizes on the `edge` FASTA), added because the 210-case matrix alone never exercises `_init_worker`/`_profile_record` in `DNAflexpy/core.py`.

**The byte-equality contract is load-bearing and fragile.** `DNAflexpy/kernel.py` deliberately uses builtin `sum()` over plain Python lists, left to right — never a numpy or `math.fsum` reduction. Those change the last bits of the float and can flip a value at the `round(v, 3)` boundary the archive also rounds to, silently breaking byte equality. If you touch `kernel.py`'s arithmetic, re-run the differential gate, not just `==` comparisons (`test_kernel.py::test_window_zero_bytes_match_the_archive_for_integer_features` explains why `==` can hide a `18 == 18.0` divergence that only serialisation exposes).

`tests/test_archive.py::test_archive_reads_its_own_lookup_table` guards the one thing that would let the differential gate pass while testing nothing: if the archive's import-path fix were ever reverted, it would read the *new* package's lookup table instead of its own. Both YAMLs currently parse identically, so no behavioural test could catch that reversion — only inspecting the archive's source for the exact `files("rxv.DNAflexpy.data")` call site can. Don't weaken that assertion to a loose substring check.

## Architecture

The new package, `DNAflexpy/`:

- `DNAflexpy/kernel.py` — pure numeric core, no I/O or pandas: `kmer_values`, `window_means`, `profile_values`. This is where the archive-compatible arithmetic lives.
- `DNAflexpy/lookup.py` — `FeatureTable`: loads and validates a feature -> k-mer -> value table, infers k-mer length from the keys themselves (no hand-maintained map to keep in sync, unlike the archive's `get_kmer_len`), and rejects mixed-length or non-ACGT keys. `default_table()` is `functools.lru_cache`d so the packaged YAML parses at most once per process.
- `DNAflexpy/core.py` — `FlexProfiler`: the class-based engine. `.profile(seq)` for one sequence, `.profile_seqs(...)` for a list/dict of sequences, `.profile_fasta(path, threads=...)` for a FASTA file, `.profile_table(path, seq_col=..., value_col=...)` for a labelled table (fills `FlexProfile.y`), `.from_bed(bed, genome, width=..., on_edge=...)` for BED intervals against a reference genome. `_MIN_RECORDS_FOR_POOL` (currently 64) governs when `profile_fasta` actually spawns a `multiprocessing.Pool`: below that record count it always runs serially regardless of `threads`, because process startup dominates at small sizes (~4900x slower measured on a 2-record file). `threads=None` means "decide automatically", not "always spawn a Pool". `threads=1` forces serial on any input size; an explicit `threads > 1` on a large-enough input still forces a Pool.
- `DNAflexpy/results.py` — `FlexProfile` (and `ProfileSet`, a `{feature: FlexProfile}` dict, defined in `core.py`): holds per-sequence rows and serialises them with `.to_tsv()`, reproducing the archive's exact ragged/NaN-padded TSV format.
- `DNAflexpy/io.py` — `read_fasta`, `read_table`, `read_bed`, `extract_intervals`. `extract_intervals` imports `pyfaidx` lazily, since it is only needed by this one entry point and ships as the optional `DNAflexpy[bed]` extra.
- `DNAflexpy/__init__.py` — re-exports the above and defines the module-level convenience function `profile(seq, feature="DNaseI", window_size=10)`. **`DNAflexpy.profile` is that function, not the `results` submodule** — the submodule is named `results.py`, not `profile.py`, specifically to avoid the function shadowing it as a package attribute. Import the class directly (`from DNAflexpy.results import FlexProfile`) rather than via `DNAflexpy.profile`.

The archive, `rxv/DNAflexpy/` (frozen, do not modify logic):

- `rxv/DNAflexpy/utils.py` — `seq_to_numeric_profile` (windowing), `transform_seq_to_feat` (k-mer -> float lookup), `get_kmer_len` (hardcoded feature -> k map), `load_feature_data`, `read_fasta`.
- `rxv/DNAflexpy/core.py` — `DNAflexpyMP` (file-level, pool-based) and `DNAflexpy` (per-sequence).
- `rxv/DNAflexpy/cli.py` — the only working CLI right now (see "No console script yet" above).
- `rxv/DNAflexpy/data/lookupNEW.yaml` — the archive's own copy of the feature table; byte-identical today to `DNAflexpy/data/lookupNEW.yaml`, but the two are loaded through separate `importlib.resources` calls and must stay that way.

## Conventions

- Version is declared once, in `DNAflexpy/__init__.py` (`__version__`); `pyproject.toml` resolves it dynamically via `[tool.setuptools.dynamic] version = { attr = "DNAflexpy.__version__" }`. There is no `setup.py` — it was deleted; `setup.py` beside a `[project]` table in `pyproject.toml` makes setuptools error, which `tests/test_packaging.py::test_setup_py_is_gone` guards.
- `pyproject.toml` explicitly lists `packages = ["DNAflexpy", "rxv", "rxv.DNAflexpy"]` rather than using `find_packages()`, and ships both YAMLs via `package-data`. Keep `requirements.txt` and `pyproject.toml`'s `dependencies` in sync when either changes (`requirements.txt` currently has no floor tighter than `pyproject.toml`'s).
- New feature tables belong in `lookupNEW.yaml` only — `DNAflexpy/lookup.py` infers k-mer length from the table's own keys, so (unlike the archive's `get_kmer_len`) there is no second place to update by hand.
- `--threads` (via `profile_fasta(threads=...)`) spawns processes, not threads, via `multiprocessing.Pool` with `spawn` on macOS. Never start a `multiprocessing.Pool` from a `python - <<EOF` stdin heredoc when testing this by hand — `spawn` re-imports the main module, and a heredoc can't be re-imported, causing a runaway respawn loop. Use a real `.py` file with an `if __name__ == "__main__":` guard, or pass `threads=1`.
- `docs/superpowers/specs/` and `docs/superpowers/plans/` are internal planning docs, not published documentation. `mkdocs.yml` excludes the whole tree via an `exclude_docs` block listing `superpowers/`, so `mkdocs gh-deploy` (run by `.github/workflows/ci.yml` on push to `main`) never publishes them — mkdocs renders every `.md` under `docs/` regardless of what's listed in `nav`, so omitting them from `nav` alone would not have been enough.
- No Cursor (`.cursor/`, `.cursorrules`) or Copilot instruction files exist in this repo.
