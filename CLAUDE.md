# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **This repo contains a new, class-based `DNAflexpy` package plus a frozen archive of the original.**
>
> `DNAflexpy/` is the new package: rewritten from scratch, class-based, and the one to import in new code. `rxv/DNAflexpy/` is the original package, moved there verbatim and **frozen — never change its logic**. It must keep running because `tests/test_differential.py` byte-compares the new package's output against it; if the archive stops working, or stops reading its own data, that gate stops meaning anything. The only sanctioned edit to the archive is the forced import-path fix at `rxv/DNAflexpy/utils.py:175` (`files("DNAflexpy.data")` → `files("rxv.DNAflexpy.data")`), required by the move. See [the design spec](docs/superpowers/specs/2026-08-12-dnaflex-rewrite-design.md) and [the phase plan](docs/superpowers/plans/2026-08-12-dnaflexpy-phases-0-2.md) for the full rationale.
>
> Active work happens on the `dev` branch.

## Commands

```bash
pip install -e .                    # editable install; also puts the DNAflexpy command on PATH
pip install -e '.[dev]'             # test deps; BED tests need pyfaidx, included here
python -m pytest -q                 # 448 passed (see "Testing" below)
python -m pytest tests/test_differential.py -k byte_identical -q   # just the byte-equality gate

mkdocs serve                        # docs preview; mkdocs build to render
python scripts/plot_profiles.py --generate-random-fasta /tmp/r.fa \
  --window-sizes 0 5 10 --feature DNaseI --out /tmp/plot.png   # needs matplotlib (not a declared dep)
```

**The `DNAflexpy` command.** `pip install -e .` puts it on PATH via `[project.scripts]`. Subcommands only — there is no flat form:

```bash
DNAflexpy profile in.fa --feature DNaseI gc --window-size 10
DNAflexpy profile peaks.bed --genome TAIR10.fa --width 200
DNAflexpy profile affinity.tsv --no-header --value-col 1
DNAflexpy profile --seq ATGCGTACGT            # to stdout, no file needed
DNAflexpy encode in.fa --features 1-mer 1-DNaseI --out X.npz
DNAflexpy plot heatmap|meta|track in.fa --feature DNaseI --out fig.png
```

`DNAflexpy/cli.py` is a front end and nothing more: **every profile byte it writes goes through `FlexProfile.to_tsv()`**, because that method is what the 230-case gate covers. A CLI path that formatted values itself would be a path nothing verifies. Multi-feature runs therefore write one file per feature (`{base}_w{window}nt_{feature}.tsv`, the archive's own pattern) rather than one combined table, since `ProfileSet` has no `to_tsv` and adding one would mean new uncovered serialisation. `--threads` defaults to `None` ("decide automatically"), not `cpu_count()` as the archive did — passing a count explicitly bypasses `_MIN_RECORDS_FOR_POOL`. Table input requires `--header`/`--no-header`, mirroring `read_table`'s deliberate refusal to guess.

The archive's CLI is untouched and still runs: `python -m rxv.DNAflexpy.cli <fasta> --window-size 10 --feature DNaseI --outfile out.tsv`. Verified byte-identical to `DNAflexpy profile` on the same input.

## Testing

Tests run and pass: `python -m pytest -q` reports **448 passed**. This includes two byte-equality matrices in `tests/test_differential.py`:

- 210 cases (3 FASTAs x 10 features x 7 window sizes), run with `threads=1`, comparing `FlexProfiler(...).profile_fasta(...).to_tsv(...)` byte-for-byte against the archive's `seq_to_numeric_profile` output.
- 20 cases forcing the pooled `multiprocessing.Pool` code path (10 features x 2 window sizes on the `edge` FASTA), added because the 210-case matrix alone never exercises `_init_worker`/`_profile_record` in `DNAflexpy/core.py`.

**The byte-equality contract is load-bearing and fragile.** `DNAflexpy/kernel.py` deliberately uses builtin `sum()` over plain Python lists — the same call the archive makes, so both shift together on any interpreter. Never substitute a numpy reduction. Note that on Python 3.12+ builtin `sum()` is itself compensated (Neumaier) and equals `math.fsum` for floats; before 3.12 it added naively and gave different last bits, which is why `requires-python` is `>=3.12` and why the checked-in expected outputs only reproduce there (`tests/test_kernel.py::test_builtin_sum_is_compensated_on_this_interpreter` pins this). Those change the last bits of the float and can flip a value at the `round(v, 3)` boundary the archive also rounds to, silently breaking byte equality. If you touch `kernel.py`'s arithmetic, re-run the differential gate, not just `==` comparisons (`test_kernel.py::test_window_zero_bytes_match_the_archive_for_integer_features` explains why `==` can hide a `18 == 18.0` divergence that only serialisation exposes).

`tests/test_archive.py::test_archive_reads_its_own_lookup_table` guards the one thing that would let the differential gate pass while testing nothing: if the archive's import-path fix were ever reverted, it would read the *new* package's lookup table instead of its own. Both YAMLs currently parse identically, so no behavioural test could catch that reversion — only inspecting the archive's source for the exact `files("rxv.DNAflexpy.data")` call site can. Don't weaken that assertion to a loose substring check.

## Architecture

The new package, `DNAflexpy/`:

- `DNAflexpy/kernel.py` — pure numeric core, no I/O or pandas: `kmer_values`, `window_means`, `profile_values`. This is where the archive-compatible arithmetic lives.
- `DNAflexpy/lookup.py` — `FeatureTable`: loads and validates a feature -> k-mer -> value table, infers k-mer length from the keys themselves (no hand-maintained map to keep in sync, unlike the archive's `get_kmer_len`), and rejects mixed-length or non-ACGT keys. `default_table()` is `functools.lru_cache`d so the packaged YAML parses at most once per process.
- `DNAflexpy/core.py` — `FlexProfiler`: the class-based engine. `.profile(seq)` for one sequence, `.profile_seqs(...)` for a list/dict of sequences, `.profile_fasta(path, threads=...)` for a FASTA file, `.profile_table(path, seq_col=..., value_col=..., header=...)` for a labelled table (`header` is required, not guessed; fills `FlexProfile.y`), `.from_bed(bed, genome, width=..., on_edge=...)` for BED intervals against a reference genome. `_MIN_RECORDS_FOR_POOL` (currently 64) governs when `profile_fasta` actually spawns a `multiprocessing.Pool`: below that record count it always runs serially regardless of `threads`, because process startup dominates at small sizes (~4900x slower measured on a 2-record file). `threads=None` means "decide automatically", not "always spawn a Pool". `threads=1` forces serial on any input size; an explicit `threads > 1` on a large-enough input still forces a Pool.
- `DNAflexpy/results.py` — `FlexProfile` (and `ProfileSet`, a `{feature: FlexProfile}` dict, defined in `core.py`): holds per-sequence rows and serialises them with `.to_tsv()`, reproducing the archive's exact ragged/NaN-padded TSV format. It also carries `.seqs`, the input sequences, purely so one-hot encoding has something to encode. **`seqs` is a parameter of `_assemble`, not just `_build`**, because the pooled `profile_fasta` branch calls `_assemble` directly — and since `to_tsv` never reads `seqs`, the byte-equality gate cannot catch that path leaving it `None`.
- `DNAflexpy/plotting.py` — `heatmap`, `metaprofile`, `trackplot`, reached as `FlexProfile.heatmap/.metaprofile` and `ProfileSet.trackplot`. **matplotlib is imported inside each function, never at module level**, and `__init__.py` must never import this module — doing so would make the optional `plot` extra a hard dependency of `import DNAflexpy` (`tests/test_plotting.py::test_importing_dnaflexpy_does_not_import_matplotlib` guards this in a subprocess). The numeric preparation (`_matrix`, `_bin`, `_zscale`, `_order`) is deliberately split from the drawing so tests assert on the plotted numbers rather than on pixels. Two refusals are load-bearing, not fussiness: `metaprofile(zscale="column")` without a background raises because column-scaling then averaging down the same column cancels to exactly zero — the line would be flat by construction; and `order_rows="cv"` raises on non-positive data because `std/mean` inverts on negative means, which `DNaseI`, `prop`, `freeen` and `mechen` all have. Binning happens **before** scaling, and `tests/test_plotting.py::test_bin_then_scale_is_not_the_same_as_scale_then_bin` pins that order with numbers a shape-only test would miss.
- `DNAflexpy/cli.py` — the `DNAflexpy` command. `plot` takes the same inputs as `profile` rather than reading a TSV back: `to_tsv` writes masked values and ragged padding both as empty fields, so a TSV reader could not tell them apart — exactly the distinction a heatmap must render.
- `DNAflexpy/encode.py` — `encode(profiles, feature_names, normalize=True)` and `FeatureMatrix`: builds a scikit-learn design matrix from a `FlexProfile` or `ProfileSet`. A request is a list of names — `"{k}-mer"` for one-hot sequence blocks, `"{n}-{feature}"` for nth-order product terms — concatenated in the order given. Min-max normalisation runs per feature block (features span 0–1 to 0–5500; one shared scale would flatten the small ones) and never touches one-hot blocks. NaN is preserved, never filled: Phase 2 made NaN mean "unresolvable" as distinct from a real zero, and a `fill=` argument would reinstate exactly that ambiguity one layer up. Equal sequence length is checked on `.seqs`, **not** on row widths — two sequences shorter than the window both produce empty rows, so a row-width check passes on input that one-hot encoding cannot handle.
- `DNAflexpy/io.py` — `read_fasta`, `read_table`, `read_bed`, `extract_intervals`. `extract_intervals` imports `pyfaidx` lazily, since it is only needed by this one entry point and ships as the optional `DNAflexpy[bed]` extra.
- `DNAflexpy/__init__.py` — re-exports the above and defines the module-level convenience function `profile(seq, feature="DNaseI", window_size=10)`. **`DNAflexpy.profile` is that function, not the `results` submodule** — the submodule is named `results.py`, not `profile.py`, specifically to avoid the function shadowing it as a package attribute. Import the class directly (`from DNAflexpy.results import FlexProfile`) rather than via `DNAflexpy.profile`. The same trap is avoided the other way round for encoding: `DNAflexpy.encode` stays the **submodule**, and the `encode` function is deliberately *not* re-exported at package level, so `from DNAflexpy.encode import encode` always resolves. The normal way to reach it is the method, `prof.encode([...])`.

The archive, `rxv/DNAflexpy/` (frozen, do not modify logic):

- `rxv/DNAflexpy/utils.py` — `seq_to_numeric_profile` (windowing), `transform_seq_to_feat` (k-mer -> float lookup), `get_kmer_len` (hardcoded feature -> k map), `load_feature_data`, `read_fasta`.
- `rxv/DNAflexpy/core.py` — `DNAflexpyMP` (file-level, pool-based) and `DNAflexpy` (per-sequence).
- `rxv/DNAflexpy/cli.py` — the archive's CLI, still runnable as `python -m rxv.DNAflexpy.cli`. Superseded by `DNAflexpy/cli.py`, and verified to produce byte-identical output.
- `rxv/DNAflexpy/data/lookupNEW.yaml` — the archive's own copy of the feature table; byte-identical today to `DNAflexpy/data/lookupNEW.yaml`, but the two are loaded through separate `importlib.resources` calls and must stay that way.

## Container

`Docker/Dockerfile` builds one full-analysis image with the `plot` and `bed` extras installed. **Nothing is published to any registry** — it is built from the repository, because DNAflexpy is not on PyPI and the builder installs it from the build context, never by name from an index.

**There is no Docker on the maintainer's machine.** The `container` CI job is therefore the only thing that verifies the Dockerfile, and it runs on pull requests as well as pushes. Never claim the image builds without pointing at a green run of that job.

Four Dockerfile decisions exist for Apptainer, where the container runs as the invoking user (who may have no passwd entry and no writable home) against a read-only filesystem: no `USER` directive, packages in system site-packages rather than user site, `MPLCONFIGDIR=/tmp/mpl` (matplotlib otherwise writes under `$HOME`), and `PYTHONDONTWRITEBYTECODE=1`. `Docker/smoke_test.sh` pins this by running one case as UID 12345 with `HOME=/nonexistent`.

`.dockerignore` excludes `tests/`, `docs/` and `Notebooks/`, but **must never exclude `rxv/`, `pyproject.toml` or `README.md`** — `pyproject.toml` lists `rxv` and `rxv.DNAflexpy` as shipped packages and reads `README.md` at build time.

## Conventions

- Version is declared once, in `DNAflexpy/__init__.py` (`__version__`); `pyproject.toml` resolves it dynamically via `[tool.setuptools.dynamic] version = { attr = "DNAflexpy.__version__" }`. There is no `setup.py` — it was deleted; `setup.py` beside a `[project]` table in `pyproject.toml` makes setuptools error, which `tests/test_packaging.py::test_setup_py_is_gone` guards.
- `pyproject.toml` explicitly lists `packages = ["DNAflexpy", "rxv", "rxv.DNAflexpy"]` rather than using `find_packages()`, and ships both YAMLs via `package-data`. Keep `requirements.txt` and `pyproject.toml`'s `dependencies` in sync when either changes (`requirements.txt` currently has no floor tighter than `pyproject.toml`'s).
- New feature tables belong in `lookupNEW.yaml` only — `DNAflexpy/lookup.py` infers k-mer length from the table's own keys, so (unlike the archive's `get_kmer_len`) there is no second place to update by hand.
- `--threads` (via `profile_fasta(threads=...)`) spawns processes, not threads, via `multiprocessing.Pool` with `spawn` on macOS. Never start a `multiprocessing.Pool` from a `python - <<EOF` stdin heredoc when testing this by hand — `spawn` re-imports the main module, and a heredoc can't be re-imported, causing a runaway respawn loop. Use a real `.py` file with an `if __name__ == "__main__":` guard, or pass `threads=1`.
- `docs/superpowers/specs/` and `docs/superpowers/plans/` are internal planning docs, not published documentation. `mkdocs.yml` excludes the whole tree via an `exclude_docs` block listing `superpowers/`, so `mkdocs gh-deploy` (run by `.github/workflows/ci.yml` on push to `main`) never publishes them — mkdocs renders every `.md` under `docs/` regardless of what's listed in `nav`, so omitting them from `nav` alone would not have been enough.
- No Cursor (`.cursor/`, `.cursorrules`) or Copilot instruction files exist in this repo.
