# DNAflexpy Phase 5 — Plotting

**Goal:** Three figures — heatmap, metaprofile, trackplot — modelled on DNAshapeR's `heatShape`, `plotShape` and `trackShape`.

**Architecture:** `DNAflexpy/plotting.py`, with the numeric preparation split from the drawing. `_prepare_*` functions return the matrix that is actually plotted; the draw functions take that and hand it to matplotlib. Tests assert on the prepared numbers, not on pixels. This is the same discipline `kernel.py` follows, and it is what makes plotting testable at all.

**Tech Stack:** matplotlib as an optional extra, `DNAflexpy[plot]`. No new hard dependency.

## Global Constraints

- **Branch:** `dev`. Commits stay **local** until the user says to push, so any open PR keeps its reviewed scope.
- **Commit messages:** no `Co-Authored-By` trailers, no AI attribution.
- **Never modify `rxv/DNAflexpy/`.**
- **The byte-equality gate is not in this path.** Nothing here writes profile values, so the `to_tsv` audit discipline the CLI needed does not apply. The 230 cases must still pass, but no plotting code should be able to touch them.
- `python -m pytest -q` currently reports **486 passed**.
- **Phase 5 stops at a user-supplied background.** The dinucleotide shuffle, `prof.test()` and `mark_significant=` are Phase 6. A background here is a second `FlexProfile` built from a control set — the DNAshapeR `plotShape(background=)` case, which stands on its own.

## Decisions

1. **`DNAflexpy/__init__.py` must not import `plotting`.** It imports `encode` at module level, which is safe because numpy is a hard dependency. The same line for plotting would make matplotlib mandatory for `import DNAflexpy` and break every install that did not ask for the extra. Import inside the method, the way `io.py` defers `pyfaidx`. **Verify in a clean venv without matplotlib that `import DNAflexpy` still works.**
2. **`DNAflexpy plot` takes the same inputs as `profile`, not a TSV.** The spec sketched `DNAflexpy plot heatmap out.tsv`. That cannot work: `to_tsv` writes masked values and ragged padding both as empty fields, and its own docstring records that no `na_rep` can separate them without breaking byte equality. A TSV reader could not tell a padded tail from a masked position — exactly the distinction a heatmap has to render. So `plot` profiles and draws in one pass. No new reader, no lossy round trip.
3. **`order_rows="cv"` raises unless the values are strictly positive.** CV is `std/mean`, meaningful only for a positive quantity. `DNaseI` spans −0.28 to 0.194, so a row mean near zero sends CV to infinity; `prop`, `freeen` and `mechen` are entirely negative, so CV comes back negative and the sort silently inverts. This follows the precedent already set for `zscale="column"` on a single-set metaprofile: raise with an explanation rather than draw something wrong. Options are `"input"`, `"std"` (default) and `"cv"`.
4. **Bin first, then z-scale.** The two orders give different figures and the spec fixes this one, so that statistics are computed on the values actually drawn. A test hand-computes a 2x6 matrix binned to 3 and scaled; a shape-only test would pass with the order reversed.
5. **Equal row widths is the right check here** — the opposite of `encode`, where sequence length was correct. The columns *are* profile positions. Also require width > 0: two sequences shorter than the window both give empty rows, and there is nothing to draw.
6. **The background must match the foreground** on feature, window size and row width. A DNaseI foreground against a gc background is silently meaningless. Dividing by the background's per-column standard deviation needs a zero guard, the same shape as `_minmax`'s constant-block branch.
7. **`ax=None` and return the `Axes`** for heatmap and metaprofile. Trackplot builds its own stacked figure, so it returns the `Figure`. The caller does the `savefig`; the CLI is what writes files.

## File structure

| File | Responsibility |
|---|---|
| `DNAflexpy/plotting.py` (create) | `_prepare_*` and the three draw functions. |
| `DNAflexpy/results.py` (modify) | `FlexProfile.heatmap`, `.metaprofile` — lazy import inside the method. |
| `DNAflexpy/core.py` (modify) | `ProfileSet.trackplot`. |
| `DNAflexpy/cli.py` (modify) | `plot` subcommand: `heatmap`, `meta`, `track`. |
| `pyproject.toml` (modify) | `plot = ["matplotlib>=3.8"]`, and `DNAflexpy[plot]` added to `dev` — CI installs `.[dev]`, so without that the plotting tests never run. |
| `tests/conftest.py` (modify) | `matplotlib.use("Agg")` before anything imports pyplot. |
| `tests/test_plotting.py` (create) | All Phase 5 tests. |
| `docs/usage.md`, `docs/api_reference.md` (modify) | Examples must end in `savefig`, never `show()` — `check_doc_examples.py` really executes them and `show()` would try to open a window. |

---

### Task 1: Shared preparation

- [x] `_matrix(profile)` — values as a 2-D array; raises on ragged rows naming the widths, and on zero width pointing at the window size.
- [x] `_bin(matrix, nbins)` — mean within equal-width position bins; raises if `nbins` exceeds the number of positions.
- [x] `_zscale(matrix, mode, background=None)` — `"column"`, `"global"`, `None`, or signal-vs-control against a background's column mean and standard deviation, with a zero-std guard.
- [x] `_order(matrix, mode)` — `"input"`, `"std"`, `"cv"`, with `"cv"` raising on non-positive data.
- [x] `_check_background(profile, background)` — feature, window size and width must match.

### Task 2: Heatmap

- [x] `heatmap(profile, nbins=None, order_rows="std", zscale="column", cmap=None, ax=None)`.
- [x] Rows are sequences, columns are positions. One feature per figure; `stiffness` (1230–5500) and `DNaseI` (−0.28–0.194) cannot share a colour scale.
- [x] A diverging colormap centred at 0 when the data is signed; sequential otherwise.
- [x] Returns the `Axes`.

### Task 3: Metaprofile

- [x] `metaprofile(profile, background=None, nbins=None, zscale=None, ax=None)`.
- [x] **Without a background, `zscale="column"` raises.** Column-wise scaling then averaging down the column gives exactly 0 at every position, by construction — a flat line at zero. Default is `"global"`.
- [x] With a background, default to the signal-vs-control z-score: standardise each column against the background's own mean and standard deviation, so the background sits flat at ~0 and the foreground's departure is readable in units of background standard deviations.
- [x] Both curves drawn when a background is given, with a legend.
- [x] Returns the `Axes`.

### Task 4: Trackplot, CLI, packaging and docs

- [x] `trackplot(profiles, seqid=None, nbins=None)` — one sequence, every feature stacked on a shared x-axis. Returns the `Figure`.
- [x] `DNAflexpy plot heatmap|meta|track <input> --feature ... --out fig.png`, taking the same input arguments as `profile`.
- [x] `pyproject.toml` extras; `tests/conftest.py` backend.
- [x] Grep for the claim that there is no plot subcommand — `CLAUDE.md` and the CLI plan both say it. The no-CLI grep found three places a hand-written list had missed.

## Definition of done

- [x] `import DNAflexpy` works with matplotlib absent, verified in a clean venv.
- [x] All three figures render from every input type.
- [x] `python -m pytest -q` green; 230 byte-equality cases unchanged.
- [x] `check_doc_examples.py` clean with `MPLBACKEND=Agg`; `mkdocs build --strict` passes.
- [x] No document still says there is no plotting.

## Outcome

`python -m pytest -q` reports **537 passed** (486 before, 51 new). The 230
byte-equality cases are unchanged. `mkdocs build --strict` passes and 34/34
documentation examples run clean under `MPLBACKEND=Agg`.

The optional-extra boundary was verified the hard way, in a fresh venv with no
matplotlib: `import DNAflexpy` succeeds, and `prof.heatmap()` raises the
install hint rather than an ImportError traceback.

### Caught during implementation

- **A documentation example was wrong**, and the doc runner caught it:
  `prof.heatmap(nbins=20)` on the fixture FASTA raised, because 26-base
  records at `window_size=10` give only 17 positions. The fixtures are now
  200 bases, and `shuffled.fa` was added for the background example.
- **`tests/test_cli.py::test_there_is_no_plot_subcommand_yet` kept passing
  after `plot` was added** - `main(["plot", "out.tsv"])` still raised, but as
  an argparse "invalid choice" rather than an unknown subcommand. A test that
  passes for the wrong reason is worse than no test; it was replaced with ten
  real ones.
- **One stale claim**, found by grepping for `plot subcommand` rather than
  trusting the plan's list. Same method that found three misses last time.

### Known issues left open at the end of Phase 5

- **The heatmap x-axis is 0-based** (matplotlib `imshow` index), while
  positions are 1-based everywhere else in the package. Cosmetic, but
  inconsistent.
- **`trackplot` reads `profiles[feature]._rows` per feature** and assumes the
  same row order across the set. That holds for anything built by a single
  `FlexProfiler`, which is the only way to get a `ProfileSet` today, but a
  hand-assembled set could break it silently. `encode` validates this;
  `trackplot` does not.
- **No `metaprofile` significance shading.** `mark_significant=` needs the
  per-position testing from Phase 6.
