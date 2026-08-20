# DNAflexpy — Restoring the CLI on the new package

**Goal:** A working `DNAflexpy` command backed by the rewritten package. The archive's CLI (`python -m rxv.DNAflexpy.cli`) is the only one that works today; this replaces it without touching it.

**Architecture:** `DNAflexpy/cli.py` holds an `argparse` front end with subcommands. It does no computation and no serialisation of its own — every path constructs a `FlexProfiler` and writes through `FlexProfile.to_tsv()`.

**Ordering note:** brought forward from Phase 8 at the user's request, ahead of Phase 5 (plotting).

## Global Constraints

- **Branch:** `dev`. Commits stay **local** until PR #2 (Phase 4) is merged on GitHub by the user — pushing early would sweep unfinished CLI work into that PR.
- **Commit messages:** no `Co-Authored-By` trailers, no AI attribution.
- **Never modify `rxv/DNAflexpy/`.** The archive's CLI keeps working exactly as it does now.
- **Every byte of profile output must come from `FlexProfile.to_tsv()`.** The 230-case gate covers that method and nothing else. A CLI path that formats values itself is a path nothing verifies. After implementation, grep `cli.py`: if any line writes profile values other than via `.to_tsv()`, that path is untested.
- `python -m pytest -q` currently reports **450 passed**, including 230 byte-equality cases. Nothing here should change any numeric output.

## Decisions

1. **Subcommands only** (user's choice). `DNAflexpy profile ...`, `DNAflexpy encode ...`. The old flat form (`DNAflexpy in.fa --feature X`) is *not* accepted. There is no compatibility debt: `pyproject.toml` has no `[project.scripts]` today, so no `DNAflexpy` command has ever been installed. The archive stays reachable as `python -m rxv.DNAflexpy.cli`.
2. **No `plot` subcommand yet.** A subcommand that exists and errors is worse than one that does not. It arrives with Phase 5.
3. **Multi-feature writes one file per feature**, named with the archive's own pattern `{base}_w{window}nt_{feature}.tsv`. A single combined file would need new serialisation code that the byte gate does not cover. `ProfileSet` has no `to_tsv`, and this is why.
4. **`--header` / `--no-header` is required for table input.** Phase 3 decided deliberately that header-ness cannot be guessed (`dna`, `rna`, `tag` are spelled with nucleotide letters). Extension-based format inference cannot recover it. Defaulting either way at the command line would quietly undo that decision.
5. **`--threads` defaults to `None`, not `cpu_count()`.** The archive defaulted to `cpu_count()`; passing that explicitly forces a Pool on any input and bypasses `_MIN_RECORDS_FOR_POOL`. `None` means "decide automatically", which is the new package's own logic.
6. **`--citation` keeps the template but says it is a template.** The archive prints literal `ADD_AUTHOR_LIST` / `ADD_TITLE` / `ADD_YEAR` / `ADD_VERSION`. Version and URL are derivable and get filled in; author and title stay as placeholders with a clear note that they need completing. The user will supply them later.
7. **`encode` derives which features to profile from `--features`.** `1-DNaseI` implies profiling `DNaseI`. No second flag to keep in sync. A request of only `k-mer` terms still needs a profile, because that is what carries the sequences — it falls back to one cheap feature and says so in the help text.

## File structure

| File | Responsibility |
|---|---|
| `DNAflexpy/cli.py` (create) | argparse front end, `profile` and `encode` subcommands. |
| `pyproject.toml` (modify) | `[project.scripts] DNAflexpy = "DNAflexpy.cli:main"`. |
| `tests/test_cli.py` (create) | Calls `main(argv)` directly; asserts files and exit codes. |
| `CLAUDE.md`, `docs/usage.md`, `README.md` (modify) | All three currently say there is no CLI. All three must change. |

**Verification gap this closes:** `scripts/check_doc_examples.py` runs ```` ```python ```` blocks only, so the bash examples in the docs are never executed. `tests/test_cli.py` covers them instead, by calling `main()` with an argv list rather than shelling out.

---

### Task 1: `profile` subcommand

- [x] `DNAflexpy profile <input>` with `--feature` (one or more), `--window-size`, `--threads`, `--outfile`, `--format`.
- [x] Format inferred from extension: `.fa`/`.fasta`/`.fna`/`.fa.gz` → FASTA, `.bed` → BED, `.tsv`/`.csv`/`.txt` → table. `--format` overrides.
- [x] BED needs `--genome`; also `--width`, `--on-edge`.
- [x] Table needs `--header`/`--no-header`; also `--seq-col`, `--value-col`, `--id-col`, `--sep`.
- [x] `--seq ATGC...` profiles one bare string to stdout, no file needed.
- [x] Default output name matches the archive: `{base}_w{window}nt_{feature}.tsv`.
- [x] Every write goes through `FlexProfile.to_tsv()`.

**Tests:** FASTA round trip; the written bytes equal `to_tsv` on the same input; multi-feature writes two correctly-named files; table without `--header` exits non-zero with a message naming the flag; BED without `--genome` exits non-zero; `--seq` prints values to stdout; unknown feature exits non-zero listing the valid ones.

### Task 2: `encode` subcommand

- [x] `DNAflexpy encode <input> --features 1-mer 1-DNaseI --out X.npz`.
- [x] Features to profile derived from the `--features` request.
- [x] `--no-normalize` maps to `normalize=False`.
- [x] Writes `.npz` holding `X`, `columns`, `seqids`, `y` (when present), `window_size`, `feature_names`.
- [x] Table input fills `y`, so one command turns a labelled file into a training set.

**Tests:** round trip through `np.load` gives back the same `X` and columns; `y` present for table input and absent otherwise; a malformed feature name exits non-zero.

### Task 3: Packaging, docs, and the three stale claims

- [x] `[project.scripts]` added; `pip install -e .` puts `DNAflexpy` on PATH.
- [x] `CLAUDE.md` — replace the "No console script yet" block.
- [x] `docs/usage.md` — rewrite `## Command line`.
- [x] `README.md` lines ~36–40 — replace "There is no command-line tool for the rewritten package yet."
- [x] `--version` and `--citation` behave as decided above.

## Definition of done

- [x] `DNAflexpy profile` and `DNAflexpy encode` work on all input types.
- [x] `python -m pytest -q` green; 230 byte-equality cases unchanged.
- [x] No profile bytes written outside `FlexProfile.to_tsv()`.
- [x] No document still claims the rewritten package has no CLI.
- [x] `python scripts/check_doc_examples.py` still clean; `mkdocs build --strict` passes.

## Outcome

`python -m pytest -q` reports **486 passed** (450 before, 36 new). The 230
byte-equality cases are unchanged. `DNAflexpy profile` output was additionally
compared against `python -m rxv.DNAflexpy.cli` on the same input and is
**byte-identical**.

### Caught during implementation

- **`--seq` had a second serialiser.** It printed values with
  `"\t".join(str(v))`, which writes `nan` where `to_tsv` writes an empty
  field — two formats for the same data, only one of them covered by the gate.
  Now routed through `to_tsv(sys.stdout)`, with a test comparing stdout against
  the file byte for byte on an `N`-containing sequence.
- **Three more documents claimed there was no CLI** beyond the three the plan
  listed: two further spots in `CLAUDE.md` and one in `docs/index.md`. Found by
  grepping for the claim rather than trusting the list.
- **The doc fixture FASTA had unequal-length records**, so the documented
  `DNAflexpy encode sequences.fa` example could not actually run. Fixed in
  `scripts/make_doc_fixtures.py`, which now also builds the `X.npz` that
  `usage.md` loads back.

### Left for later

- `--citation` still prints `ADD_AUTHOR_LIST` and `ADD_TITLE`, and says so on
  stderr. Version and URL are filled in from the package. The user will supply
  author and title.
- No `plot` subcommand. It arrives with Phase 5.
