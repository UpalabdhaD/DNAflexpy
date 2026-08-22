# A single full-analysis image: every optional extra is installed, so every
# subcommand works and there is no "which feature is missing" question.
#
# Built locally from this repository. Nothing is published to any registry.
#
# Two stages. The builder compiles wheels; the runtime copies only the
# installed packages, so the build toolchain never ships.

FROM python:3.12-slim AS builder

# DNAflexpy is not on PyPI, so it is installed from the build context and
# never by name from an index. `pip wheel` resolves the dependencies too, so
# the runtime stage can install with --no-index.
COPY . /src
RUN pip wheel --no-cache-dir --wheel-dir=/wheels "/src[plot,bed]"


FROM python:3.12-slim

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links=/wheels \
        "DNAflexpy[plot,bed]" \
    && rm -rf /wheels

# These three exist for Apptainer, where the image filesystem is read-only and
# the container runs as the invoking user -- who may have no passwd entry and
# no writable home.
#
#   MPLBACKEND=Agg          there is no display; the default backend would fail
#   MPLCONFIGDIR=/tmp/mpl   matplotlib otherwise writes under $HOME and warns
#                           or fails when that is unwritable
#   PYTHONDONTWRITEBYTECODE stops .pyc writes into a read-only layer
ENV MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/mpl \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# There is deliberately no USER directive. Apptainer runs the container as the
# invoking user regardless, and a USER line here would only mislead. Packages
# install to system site-packages, not user site, so they stay importable as
# an arbitrary UID. Nothing is chowned or chmodded at runtime.

# Everything the user reads or writes lives on the bind mount.
WORKDIR /data

ENTRYPOINT ["DNAflexpy"]
CMD ["--help"]
