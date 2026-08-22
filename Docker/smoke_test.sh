#!/usr/bin/env bash
# End-to-end check that a built image actually works.
#
# Runs every subcommand against data the script generates itself, so it needs
# no bind mount of repository files and no packaged fixture. Run it against an
# image:
#
#   ./Docker/smoke_test.sh dnaflexpy:dev
#
# A broken image should fail here, at build time, rather than on a cluster.
set -euo pipefail

IMAGE="${1:-dnaflexpy:dev}"

# Two directories, deliberately. Docker runs as root by default, so everything
# written into WORK ends up root-owned. The unprivileged case at the end needs
# a directory it can actually read and write, so it gets a fresh one that this
# script owns -- which is also what the real Apptainer case looks like: a user
# running against their own data directory.
WORK="$(mktemp -d)"
UWORK="$(mktemp -d)"
trap 'rm -rf "$WORK" "$UWORK" 2>/dev/null || true' EXIT

# Equal-length records: encode and all three plots need equal lengths.
{
  echo ">a"
  echo "ATGCGTACGTAGCTAGCGTAGCTAGTATGCGTACGTAGCTAGCGTAGCTAG"
  echo ">b"
  echo "CGTAGCTAGTATGCGTACGTAGCTAGCGTAGCTAGTATGCGTACGTAGCTA"
} > "$WORK/in.fa"
printf 'ATGCGTACGT\t1.5\nCGTAGCTAGT\t2.5\n' > "$WORK/affinity.tsv"

run() { docker run --rm -v "$WORK":/data "$IMAGE" "$@"; }

check() {  # check <full-path> <description> [magic-string]
  local path="$1"
  [ -s "$path" ] || { echo "FAIL: $2 produced nothing at $path"; exit 1; }
  if [ $# -ge 3 ]; then
    # A fixed 16 bytes, never ${#3}: a PNG header is \x89PNG\r\n\x1a\n, so the
    # literal "PNG" starts at byte 2. Reading only as many bytes as the pattern
    # is long would cut it off and fail on a valid file. -a because it is
    # binary input.
    head -c 16 "$path" | grep -qa "$3" \
      || { echo "FAIL: $path is not a $3 file"; exit 1; }
  fi
  echo "  ok: $2"
}

echo "Smoke-testing $IMAGE"

echo "--- version and help ---"
run --version
run --help > /dev/null
echo "  ok: the entrypoint runs"

echo "--- profile ---"
run profile in.fa --feature DNaseI --window-size 10 --outfile out.tsv
check "$WORK/out.tsv" "profile from FASTA"

run profile in.fa --feature DNaseI gc --window-size 5
check "$WORK/in_w5nt_DNaseI.tsv" "multi-feature profile (DNaseI)"
check "$WORK/in_w5nt_gc.tsv" "multi-feature profile (gc)"

run profile --seq ATGCGTACGTAGCTAGCGTAGCTAGT | grep -q "sequence" \
  || { echo "FAIL: --seq wrote nothing to stdout"; exit 1; }
echo "  ok: profile --seq to stdout"

run profile affinity.tsv --no-header --outfile labelled.tsv
check "$WORK/labelled.tsv" "profile from a labelled table"

echo "--- encode ---"
run encode in.fa --features 1-mer 1-DNaseI --window-size 0 --out X.npz
check "$WORK/X.npz" "encode to npz"

echo "--- plot (needs the plot extra) ---"
run plot heatmap in.fa --feature DNaseI --nbins 10 --out heat.png
check "$WORK/heat.png" "heatmap" "PNG"

run plot meta in.fa --feature DNaseI --out meta.png
check "$WORK/meta.png" "metaprofile" "PNG"

run plot track in.fa --feature DNaseI gc stiffness --out track.png
check "$WORK/track.png" "trackplot" "PNG"

echo "--- the packaged example data ships ---"
docker run --rm --entrypoint python "$IMAGE" -c "
from DNAflexpy import example_path, example_files
names = example_files()
assert 'promoters.fa' in names, names
assert example_path('promoters.fa').read_text().startswith('>'), 'promoters.fa is empty'
print(f'  ok: {len(names)} example file(s) ship with the package')
"

echo "--- optional extras are really installed ---"
docker run --rm --entrypoint python "$IMAGE" -c \
  "import pyfaidx, matplotlib; print('  ok: pyfaidx and matplotlib import')"

echo "--- the frozen archive ships and still reads its own table ---"
docker run --rm --entrypoint python "$IMAGE" -c "
from rxv.DNAflexpy.utils import load_feature_data
assert load_feature_data(), 'the archive could not read its lookup table'
print('  ok: rxv.DNAflexpy loads its own lookup table')
"

echo "--- unprivileged, no writable home: the Apptainer case ---"
# A clean directory this script owns, holding only the input. Reusing WORK
# would not work: the runs above wrote into it as root, and a non-root caller
# cannot chmod root-owned files.
cp "$WORK/in.fa" "$UWORK/in.fa"
chmod 0777 "$UWORK"
chmod 0644 "$UWORK/in.fa"
docker run --rm --user 12345:12345 -e HOME=/nonexistent -v "$UWORK":/data "$IMAGE" \
  plot heatmap in.fa --feature DNaseI --out asuser.png
check "$UWORK/asuser.png" "heatmap as an arbitrary UID with no home" "PNG"

echo
echo "All smoke tests passed for $IMAGE"
