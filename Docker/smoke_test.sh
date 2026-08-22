#!/usr/bin/env bash
# End-to-end check that a built image actually works.
#
# Runs every subcommand against a FASTA the script generates itself, so it
# needs no bind mount and no packaged fixture. Run it against an image:
#
#   ./Docker/smoke_test.sh dnaflexpy:dev
#
# A broken image should fail here, at build time, rather than on a cluster.
set -euo pipefail

IMAGE="${1:-dnaflexpy:dev}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Two equal-length records: encode and the plots both need equal lengths.
{
  echo ">a"
  echo "ATGCGTACGTAGCTAGCGTAGCTAGTATGCGTACGTAGCTAGCGTAGCTAG"
  echo ">b"
  echo "CGTAGCTAGTATGCGTACGTAGCTAGCGTAGCTAGTATGCGTACGTAGCTA"
} > "$WORK/in.fa"
printf 'ATGCGTACGT\t1.5\nCGTAGCTAGT\t2.5\n' > "$WORK/affinity.tsv"

run() { docker run --rm -v "$WORK":/data "$IMAGE" "$@"; }

check() {  # check <file> <description> [magic-bytes]
  local path="$WORK/$1"
  [ -s "$path" ] || { echo "FAIL: $2 produced nothing at $1"; exit 1; }
  if [ $# -ge 3 ]; then
    # Read a fixed 16 bytes, not ${#3}. A PNG header is \x89PNG\r\n\x1a\n, so
    # the literal "PNG" starts at byte 2 -- reading only as many bytes as the
    # pattern is long would cut it off and fail on a perfectly good file.
    # grep -a because the input is binary.
    head -c 16 "$path" | grep -qa "$3" || { echo "FAIL: $1 is not a $3 file"; exit 1; }
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
check out.tsv "profile from FASTA"

run profile in.fa --feature DNaseI gc --window-size 5
check in_w5nt_DNaseI.tsv "multi-feature profile (DNaseI)"
check in_w5nt_gc.tsv "multi-feature profile (gc)"

run profile --seq ATGCGTACGTAGCTAGCGTAGCTAGT | grep -q "sequence" \
  || { echo "FAIL: --seq wrote nothing to stdout"; exit 1; }
echo "  ok: profile --seq to stdout"

run profile affinity.tsv --no-header --outfile labelled.tsv
check labelled.tsv "profile from a labelled table"

echo "--- encode ---"
run encode in.fa --features 1-mer 1-DNaseI --window-size 0 --out X.npz
check X.npz "encode to npz"

echo "--- plot (needs the plot extra) ---"
run plot heatmap in.fa --feature DNaseI --nbins 10 --out heat.png
check heat.png "heatmap" "PNG"

run plot meta in.fa --feature DNaseI --out meta.png
check meta.png "metaprofile" "PNG"

run plot track in.fa --feature DNaseI gc stiffness --out track.png
check track.png "trackplot" "PNG"

echo "--- optional extras are really installed ---"
docker run --rm --entrypoint python "$IMAGE" -c "import pyfaidx, matplotlib; print('  ok: pyfaidx and matplotlib import')"

echo "--- the frozen archive ships and still reads its own table ---"
docker run --rm --entrypoint python "$IMAGE" -c "
from rxv.DNAflexpy.utils import load_feature_data
assert load_feature_data(), 'the archive could not read its lookup table'
print('  ok: rxv.DNAflexpy loads its own lookup table')
"

echo "--- unprivileged, no writable home: the Apptainer case ---"
# mktemp -d creates the directory 0700, owned by whoever runs this script, and
# the runs above wrote into it as root. An arbitrary UID could not even stat
# in.fa. Open it up first: that is a property of the host bind mount, not of
# the image, and under Apptainer the user owns their own data directory
# anyway. Without this the check fails on a perfectly good image.
chmod -R a+rwX "$WORK"
docker run --rm --user 12345:12345 -e HOME=/nonexistent -v "$WORK":/data "$IMAGE" \
  plot heatmap in.fa --feature DNaseI --out asuser.png
check asuser.png "heatmap as an arbitrary UID with no home" "PNG"

echo
echo "All smoke tests passed for $IMAGE"
