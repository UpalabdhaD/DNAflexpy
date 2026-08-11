# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **`DNAflexpy/` is frozen — do not modify it.** It is being replaced by a new `dnaflex/` package written from scratch; see [the design spec](docs/superpowers/specs/2026-08-12-dnaflex-rewrite-design.md). The old package must keep working, because the new one is verified by byte-comparing its output against it. The user will archive `DNAflexpy/` themselves. The rest of this file documents that frozen package. Active work happens on the `dev` branch.

## Commands

```bash
pip install -e .                    # installs the `DNAflexpy` console script
python -m DNAflexpy.cli <fasta> --window-size 10 --feature DNaseI --outfile out.tsv
DNAflexpy <fasta> --window-size 10 --feature DNaseI    # same thing, needs the install above
DNAflexpy --citation

pytest tests/                                                   # see caveat below
pytest tests/test_utils.py::TestCalculateWindowAverages::test_simple_sequence   # single test

mkdocs serve                        # docs preview; mkdocs build to render
python scripts/plot_profiles.py --generate-random-fasta /tmp/r.fa \
  --window-sizes 0 5 10 --feature DNaseI --out /tmp/plot.png   # needs matplotlib (not a declared dep)
```

**Testing is currently manual, not automated.** `pytest` fails at collection: [tests/test_utils.py:2](tests/test_utils.py#L2) does `from utils import calculate_window_averages`, a module path and a function name that no longer exist (the function was split into `seq_to_numeric_profile` / `transform_seq_to_feat`). [tests/test_core.py](tests/test_core.py) is empty. CI ([.github/workflows/ci.yml](.github/workflows/ci.yml)) only runs `mkdocs gh-deploy` — there is no test job. Verification to date lives in prose plus checked-in expected outputs: [tests/test_verification.md](tests/test_verification.md), [tests/edge_cases_documentation.md](tests/edge_cases_documentation.md), [docs/window_size_checks.md](docs/window_size_checks.md), and the TSVs under [tests/test_outputs/](tests/test_outputs/) and [tests/edge_case_outputs/](tests/edge_case_outputs/). When changing numeric behavior, re-run the commands in those markdown files and diff against the recorded values.

## Architecture

Call path, CLI to values:

[cli.py](DNAflexpy/cli.py) → `DNAflexpyMP` fans records out over a `multiprocessing.Pool` → `DNAflexpy_for_CLI` per record → `seq_to_numeric_profile` (windowing) → `transform_seq_to_feat` (k-mer → float lookup).

**Two public entry points with confusingly similar names.** `DNAflexpyMP(input_file, window_size, feature, threads, outfile)` in [core.py:42](DNAflexpy/core.py#L42) is the file-level API and the one the CLI uses. `DNAflexpy(seqid, record, window_size, feature, feature_lookup)` at [core.py:9](DNAflexpy/core.py#L9) is **per-sequence**, reloads the entire YAML on every call, and overwrites the `feature_lookup` argument it was handed. The README and [docs/usage.md](docs/usage.md) show `DNAflexpy(input_file=..., threads=...)`, which does not match either signature and raises `TypeError` — use `DNAflexpyMP` for file-level work.

**Feature → k-mer length is a hand-maintained map that must stay in sync with the YAML.** `get_kmer_len` ([utils.py:148](DNAflexpy/utils.py#L148)) hardcodes which features are trinucleotide (k=3: DNaseI, NPP, bendabilityDNase, bendabilityConcensus) vs dinucleotide (k=2: wedge, prop, twistDisp, stiffness, bendingStiffness, trx). Adding a feature to the YAML without adding it here yields `None`, which then fails silently. `lookupNEW.yaml` currently holds 13 top-level features; `freeen`, `gc`, and `mechen` are unmapped and produce all-`None` rows today. The CLI help and docs advertise only 5 of the 10 that work.

**Two YAML lookup tables, only one live.** `load_feature_data` ([utils.py:163](DNAflexpy/utils.py#L163)) reads `DNAflexpy/data/lookupNEW.yaml` via `importlib.resources`, and that is the only file shipped by [MANIFEST.in](MANIFEST.in). Its schema is `feature → kmer → value`. `DNAflexpy/data/lookup.yaml` is legacy, unused, and uses an incompatible schema grouped by `trinucleotide` / `dinucleotide`.

**Error model: exceptions are swallowed, failures surface as data.** Nearly every function wraps its body in a broad `try/except` that prints and returns `None` or `0`. Bad features, unknown k-mers, and malformed input therefore show up as `None`/NaN rows in the output TSV rather than as a raised error. When debugging a wrong-looking profile, suspect a swallowed exception before suspecting the math; the only place that re-raises is `DNAflexpyMP`'s outer handler.

**Row shape drives the output format.** Each record becomes a flat list `[seqid, *values]`. Lists are ragged across records (a sequence shorter than the window yields `[seqid]` alone), so `pd.DataFrame(rows)` NaN-pads to the widest row, then writes with `header=False, index=False, sep="\t"`. This is why short sequences appear as ID-only lines and why trailing empty fields show up in the TSV.

**Windowing semantics** (`seq_to_numeric_profile`, [utils.py:22](DNAflexpy/utils.py#L22)): `window_size == 0` returns per-k-mer values for the whole sequence; `window_size == N` slides overlapping length-N windows with step 1 and emits the mean of the k-mer values inside each, rounded to 3 decimals, giving `len(seq) - N + 1` values. `N == kmer_len` reproduces the `N == 0` output exactly (one k-mer per window) — that identity is the cheapest regression check for the averaging formula, and is what the `a9bd674` window-size fix restored.

## Conventions

- `DNAflexpy/__init__.py` is empty and exports nothing — always import from the submodule: `from DNAflexpy.core import DNAflexpyMP`.
- Version is declared only in [setup.py](DNAflexpy/../setup.py) (`0.2.0`); [pyproject.toml](pyproject.toml) carries build-system config only. `setup.py`'s `install_requires` pins (`pandas>=1.0`, `pyyaml>=5.4`) are looser than [requirements.txt](requirements.txt) (`pandas>=2.2.3`, `pyyaml>=6.0.2`) — update both together.
- `--threads` spawns processes, not threads; the full lookup dict is pickled to every worker.
- New feature tables belong in `lookupNEW.yaml` **and** `get_kmer_len` **and** the CLI `--feature` help text.
- No Cursor (`.cursor/`, `.cursorrules`) or Copilot instruction files exist in this repo.
