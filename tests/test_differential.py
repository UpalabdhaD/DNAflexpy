"""Phase 2 gate: new output must be byte-identical to the archive."""
import contextlib
import io
import math

import pytest

from DNAflexpy.core import FlexProfiler
from rxv.DNAflexpy.utils import (
    get_kmer_len,
    load_feature_data,
    read_fasta as archive_read_fasta,
    seq_to_numeric_profile,
)

from tests.conftest import ROOT, FASTAS, archive_bytes

# The 10 features the archive's hardcoded map can reach. gc, freeen and
# mechen are excluded on purpose: they are allow-listed divergences below.
ARCHIVE_FEATURES = [
    "DNaseI", "NPP", "bendabilityDNase", "bendabilityConcensus",
    "wedge", "prop", "twistDisp", "stiffness", "bendingStiffness", "trx",
]
WINDOWS = [0, 2, 3, 10, 15, 26, 30]


def new_bytes(fasta, feature, window, tmp_path):
    out = tmp_path / "new.tsv"
    FlexProfiler(feature, window_size=window).profile_fasta(fasta, threads=1).to_tsv(out)
    return out.read_bytes()


@pytest.mark.parametrize("fasta_key", sorted(FASTAS))
@pytest.mark.parametrize("feature", ARCHIVE_FEATURES)
@pytest.mark.parametrize("window", WINDOWS)
def test_byte_identical_to_archive(fasta_key, feature, window, tmp_path):
    fasta = FASTAS[fasta_key]
    assert new_bytes(fasta, feature, window, tmp_path) == archive_bytes(
        fasta_key, feature, window
    )


@pytest.mark.parametrize("feature", ["gc", "freeen", "mechen"])
def test_allowlisted_unlocked_features(feature, tmp_path):
    """Archive: None rows, because get_kmer_len has no entry. New: numbers."""
    lookup = load_feature_data()
    assert get_kmer_len(feature) is None
    with contextlib.redirect_stdout(io.StringIO()):
        archived = seq_to_numeric_profile(
            "s", "ATGCGTACGT", get_kmer_len(feature), 10, feature, lookup
        )
    assert archived is None

    values = FlexProfiler(feature, window_size=0).profile("ATGCGTACGT")
    assert len(values) == 9
    assert not any(math.isnan(v) for v in values)

    # And the whole-file path produces a real TSV, not a column of "None".
    text = new_bytes(FASTAS["test_fasta"], feature, 10, tmp_path).decode()
    assert "None" not in text
    assert len(text.splitlines()) == 2


def test_allowlisted_unknown_feature_now_raises():
    """Archive: get_kmer_len misses, the TypeError is swallowed, None row.
    New: raises at construction."""
    assert get_kmer_len("not_a_feature") is None
    with contextlib.redirect_stdout(io.StringIO()):
        archived = seq_to_numeric_profile(
            "s", "ATGCGTACGT", get_kmer_len("not_a_feature"), 10,
            "not_a_feature", load_feature_data(),
        )
    assert archived is None

    with pytest.raises(ValueError, match="unknown feature"):
        FlexProfiler("not_a_feature")


def test_allowlisted_negative_window_now_raises():
    """The archive silently returns None; the new package refuses."""
    lookup = load_feature_data()
    with contextlib.redirect_stdout(io.StringIO()):
        assert seq_to_numeric_profile("s", "ATGCGT", 3, -1, "DNaseI", lookup) is None
    with pytest.raises(ValueError, match="window_size must be >= 0"):
        FlexProfiler("DNaseI", window_size=-1)


def test_allowlisted_ambiguity_masking(tmp_path):
    """Archive scores N as 0; the new package masks it. Off the matrix above.

    `tests/test_ambiguous.fasta` is deliberately not one of the three FASTAs
    in `FASTAS` (see `test_no_test_fasta_contains_ambiguity` below), so the
    shared `archive_bytes(fasta_key, ...)` helper -- keyed only on
    "test_fasta"/"edge"/"exact" -- doesn't apply here. Rather than widen
    `FASTAS` (which would silently pull this fixture into the 3x10x7 matrix
    above and its no-ambiguity guard test) or redefine a second general
    archive-reproduction helper (the drift risk the brief's own
    `archive_bytes` override exists to avoid), this calls the same archive
    primitives directly -- exactly what the brief's own
    `test_allowlisted_unlocked_features` and
    `test_allowlisted_negative_window_now_raises` already do for cases
    outside the matrix.
    """
    fasta = ROOT / "tests/test_ambiguous.fasta"
    prof = FlexProfiler("DNaseI", window_size=0).profile_fasta(fasta, threads=1)

    assert prof.n_masked["clean"] == 0
    assert prof.n_masked["has_n"] == 3
    assert prof.n_masked["has_iupac"] == 3
    assert prof.n_masked["all_n"] == 8

    lookup = load_feature_data()
    kmer_len = get_kmer_len("DNaseI")
    seqs = dict(archive_read_fasta(str(fasta)))
    profiler = FlexProfiler("DNaseI", window_size=0)

    # Position-exact divergence check per sequence: same length throughout,
    # and the archive has a literal fabricated 0 at exactly the positions
    # where the new package masks to NaN -- everywhere else the two agree.
    with contextlib.redirect_stdout(io.StringIO()):
        for seqid, seq in seqs.items():
            archived_row = seq_to_numeric_profile(
                seqid, seq, kmer_len, 0, "DNaseI", lookup
            )[1:]
            new_row = profiler.profile(seq)
            assert len(archived_row) == len(new_row), seqid

            masked_positions = {i for i, v in enumerate(new_row) if math.isnan(v)}
            assert len(masked_positions) == prof.n_masked[seqid]

            for i, (archived_v, new_v) in enumerate(zip(archived_row, new_row)):
                if i in masked_positions:
                    assert archived_v == 0, f"{seqid} pos {i}: expected fabricated 0"
                else:
                    assert archived_v == pytest.approx(new_v), f"{seqid} pos {i}"

    # And the brief's own framing, verified concretely: the archive's row
    # for has_n literally carries three consecutive fabricated zeros.
    with contextlib.redirect_stdout(io.StringIO()):
        archived_has_n = seq_to_numeric_profile(
            "has_n", seqs["has_n"], kmer_len, 0, "DNaseI", lookup
        )
    n_row = "\t".join(str(v) for v in archived_has_n)
    assert "\t0\t0\t0\t" in n_row


def test_no_test_fasta_contains_ambiguity():
    """Guards the claim that masking cannot affect the differential matrix."""
    for fasta in FASTAS.values():
        seqs = "".join(s for _, s in archive_read_fasta(str(fasta))).upper()
        assert not set(seqs) - set("ACGT"), fasta
