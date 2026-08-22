# Running DNAflexpy in a container

One image with every optional extra installed, so all subcommands work and
there is no "which feature is missing" question.

The image is built from this repository. **Nothing is published to any
registry**, so there is no `docker pull` — you build it yourself. Publishing can
be added later without changing the image.

## Build

From the repository root, not from this directory:

```bash
docker build -t dnaflexpy:dev -f Docker/Dockerfile .
```

DNAflexpy is not on PyPI, so the build installs it from the files you just
checked out. That is why the build context is the repository root.

## Run

Mount the directory holding your data at `/data`. That is the working
directory inside the container, so filenames are relative to it:

```bash
docker run --rm -v "$PWD":/data dnaflexpy:dev \
  profile in.fa --feature DNaseI --window-size 10 --outfile out.tsv
```

The entrypoint is `DNAflexpy` itself, so everything after the image name is
just the normal command line:

```bash
docker run --rm -v "$PWD":/data dnaflexpy:dev --version
docker run --rm -v "$PWD":/data dnaflexpy:dev profile in.fa --feature DNaseI gc
docker run --rm -v "$PWD":/data dnaflexpy:dev encode in.fa --features 1-mer 1-DNaseI --out X.npz
docker run --rm -v "$PWD":/data dnaflexpy:dev plot heatmap in.fa --nbins 20 --out heat.png
```

Run it with no arguments and it prints the help.

### Files come out owned by root

Docker runs the container as root by default, so anything written to the bind
mount belongs to root on your machine. Pass your own IDs to avoid that:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data dnaflexpy:dev \
  profile in.fa --feature DNaseI
```

This works because nothing in the image needs root, and nothing is written
outside `/data` and `/tmp`.

## Apptainer / Singularity, on a cluster

The image is built to run unprivileged. You do not need root on the cluster,
only somewhere to build or fetch the `.sif`.

```bash
# From a local Docker image
apptainer build dnaflexpy.sif docker-daemon://dnaflexpy:dev

# Then, on the cluster
apptainer exec --bind "$PWD":/data dnaflexpy.sif \
  DNAflexpy profile /data/in.fa --feature DNaseI --outfile /data/out.tsv
```

Note the two differences from Docker:

- `apptainer exec` ignores the image entrypoint, so name `DNAflexpy`
  explicitly.
- Apptainer does not apply `WORKDIR`, so give absolute paths under `/data`.

Everything else works the same.

### Why it survives an unprivileged, read-only environment

Apptainer runs the container as **you** — a UID that may have no entry in the
image's `/etc/passwd` and no writable home. The image filesystem is read-only.
Four decisions in the Dockerfile exist for that:

| Decision | Why |
|---|---|
| No `USER` directive | Apptainer overrides it anyway; a `USER` line would only mislead |
| Packages in system site-packages | They stay importable as an arbitrary UID |
| `MPLCONFIGDIR=/tmp/mpl` | matplotlib writes a config directory under `$HOME` otherwise, and warns or fails when that is unwritable |
| `PYTHONDONTWRITEBYTECODE=1` | Stops `.pyc` writes into a read-only layer |

There is also no `chown` or `chmod` at runtime, and every output goes to the
bind mount rather than into the image.

`Docker/smoke_test.sh` checks this directly: one case runs the container as
UID 12345 with `HOME=/nonexistent` and asserts a plot still renders.

## Verifying an image

```bash
./Docker/smoke_test.sh dnaflexpy:dev
```

That runs `profile`, `encode` and all three plots end to end against a FASTA it
generates itself, checks the outputs are real files (PNGs are checked for the
PNG magic bytes, not just for existing), confirms `pyfaidx` and `matplotlib`
are genuinely installed, confirms the frozen `rxv.DNAflexpy` archive still
reads its own lookup table, and runs the unprivileged case described above.

CI runs this on every push and every pull request, so a broken image is caught
there rather than on a cluster.

## What is inside

- Python 3.12 slim
- DNAflexpy with the `plot` and `bed` extras — matplotlib and pyfaidx
- The frozen `rxv.DNAflexpy` archive, which the differential test suite
  compares against

Two stages: the first compiles wheels, the second installs them and keeps no
build toolchain, so the image stays small.
