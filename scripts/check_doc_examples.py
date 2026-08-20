"""Run every python example in the docs, in order, the way a reader would.

Blocks share one namespace per file, because a reader following along has
whatever the earlier blocks defined.
"""
import io as _io, os, pathlib, re, sys, traceback, warnings

REPO = pathlib.Path(__file__).resolve().parent.parent
os.chdir(sys.argv[1] if len(sys.argv) > 1 else ".")

FILES = sorted((REPO / "docs").glob("*.md")) + [REPO / "README.md"]
BLOCK = re.compile(r"```(?:python|py)\n(.*?)```", re.S)

# Examples that name a file the reader supplies; we made fixtures with these names.
fails = total = 0
for path in FILES:
    if not path.exists():
        continue
    ns = {}
    for i, code in enumerate(BLOCK.findall(path.read_text()), 1):
        # skip pure import-listing blocks that are illustrative only
        total += 1
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                exec(compile(code, f"{path.name}:block{i}", "exec"), ns)
        except Exception:
            fails += 1
            print(f"\nFAIL  {path.name} block {i}")
            print("  " + "\n  ".join(code.strip().splitlines()[:6]))
            print("  ->", traceback.format_exc().strip().splitlines()[-1])

print(f"\n{total - fails}/{total} python examples run clean")
sys.exit(1 if fails else 0)
