"""Run every python example in the docs, in order, the way a reader would.

Blocks share one namespace per file, because a reader following along has
whatever the earlier blocks defined.

Two checks, not one:

1. Every ```python block runs without raising.
2. When a python block is immediately followed by a plain ``` block, that
   block is treated as the expected stdout and compared against what the code
   actually printed.

The second check exists because the first is not enough. A block can run
perfectly while the output printed underneath it in the documentation is
something the author made up -- which is exactly the mistake this was added
to catch. Code that runs and a claim that is false is worse than code that
does not run, because nothing flags it.

To show output that is genuinely not reproducible (a timing, a path), put it
in a fenced block with a language tag such as ```text so it is not compared.

    python scripts/make_doc_fixtures.py /tmp/docfix
    MPLBACKEND=Agg python scripts/check_doc_examples.py /tmp/docfix
"""
import contextlib
import io
import os
import pathlib
import re
import sys
import traceback
import warnings

REPO = pathlib.Path(__file__).resolve().parent.parent
os.chdir(sys.argv[1] if len(sys.argv) > 1 else ".")

FILES = sorted((REPO / "docs").glob("*.md")) + [REPO / "README.md"]

# A python block, optionally followed by an untagged block holding its output.
# The output block must start on the very next line to count.
BLOCK = re.compile(
    r"```(?:python|py)\n(?P<code>.*?)```"          # the code
    r"(?:\n\n?(?P<expected>```\n.*?```))?",        # its output, if given
    re.S,
)
FENCED = re.compile(r"^```\n(.*?)```$", re.S)


def _normalise(text):
    """Compare on visible content, not on trailing whitespace."""
    return [line.rstrip() for line in text.strip().splitlines()]


fails = mismatches = total = checked_output = 0

for path in FILES:
    if not path.exists():
        continue
    namespace = {}
    for i, match in enumerate(BLOCK.finditer(path.read_text()), 1):
        code = match.group("code")
        total += 1
        captured = io.StringIO()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with contextlib.redirect_stdout(captured):
                    exec(compile(code, f"{path.name}:block{i}", "exec"), namespace)
        except Exception:
            fails += 1
            print(f"\nFAIL  {path.name} block {i} raised")
            print("  " + "\n  ".join(code.strip().splitlines()[:6]))
            print("  ->", traceback.format_exc().strip().splitlines()[-1])
            continue

        expected_block = match.group("expected")
        if not expected_block:
            continue
        expected = FENCED.match(expected_block.strip())
        if not expected:
            continue

        checked_output += 1
        want = _normalise(expected.group(1))
        got = _normalise(captured.getvalue())
        if want != got:
            mismatches += 1
            print(f"\nMISMATCH  {path.name} block {i} printed something else")
            print("  documented:")
            for line in want[:6]:
                print(f"    {line}")
            print("  actual:")
            for line in got[:6] or ["(nothing)"]:
                print(f"    {line}")

print(f"\n{total - fails}/{total} python examples run clean")
print(f"{checked_output - mismatches}/{checked_output} documented outputs match")
sys.exit(1 if fails or mismatches else 0)
