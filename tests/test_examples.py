"""The example data shipped with the package."""
import pytest

from DNAflexpy import FlexProfiler, describe_examples, example_files, example_path
from DNAflexpy.io import read_bed, read_fasta, read_table


def test_every_named_example_actually_ships():
    for name in example_files():
        assert example_path(name).stat().st_size > 0, name


def test_an_unknown_name_lists_what_exists():
    with pytest.raises(ValueError, match="promoters.fa"):
        example_path("nosuch.fa")


def test_descriptions_cover_every_file():
    assert set(describe_examples()) == set(example_files())


def test_the_promoters_are_equal_length_so_every_feature_works():
    """encode and all three plots need equal lengths. If this file ever
    stops being uniform, half the documentation stops running."""
    records = list(read_fasta(example_path("promoters.fa")))
    assert len(records) == 12
    assert len({len(seq) for _, seq in records}) == 1


def test_the_control_set_matches_the_promoters():
    """A background is only a control if it lines up with the foreground."""
    fg = list(read_fasta(example_path("promoters.fa")))
    bg = list(read_fasta(example_path("control.fa")))
    assert {len(s) for _, s in fg} == {len(s) for _, s in bg}


def test_the_control_set_has_the_same_base_composition():
    """The difference the examples show must be about position, not
    composition -- otherwise the tutorial teaches the wrong lesson."""
    def gc(records):
        seq = "".join(s for _, s in records)
        return sum(1 for b in seq if b in "GC") / len(seq)

    fg = gc(list(read_fasta(example_path("promoters.fa"))))
    bg = gc(list(read_fasta(example_path("control.fa"))))
    assert abs(fg - bg) < 0.02, f"GC differs too much: {fg:.3f} vs {bg:.3f}"


def test_the_promoters_carry_a_real_signal():
    """The examples exist to show something. If the AT-rich element stops
    being detectable, the documented figures become meaningless."""
    import numpy as np

    from DNAflexpy.plotting import _matrix

    p = FlexProfiler("DNaseI", window_size=10)
    fg = _matrix(p.profile_fasta(example_path("promoters.fa")))
    bg = _matrix(p.profile_fasta(example_path("control.fa")))
    z = (fg.mean(axis=0) - bg.mean(axis=0)) / bg.std(axis=0)
    strongest = int(np.argmax(np.abs(z)))
    # The element sits at bases 40-60; with a 10-base window that shows up
    # between positions 31 and 60.
    assert 31 <= strongest + 1 <= 60, f"strongest signal at {strongest + 1}"
    assert abs(z[strongest]) > 1.0, f"signal is only {z[strongest]:.2f} SDs"


def test_the_labelled_table_loads_and_has_a_header():
    rows = read_table(example_path("affinity.tsv"), seq_col="sequence",
                      value_col="binding_score", id_col="name", header=True)
    assert len(rows) == 12
    assert all(v > 0 for _, _, v in rows)


def test_the_bed_intervals_are_inside_the_toy_genome():
    genome = dict(read_fasta(example_path("genome.fa")))
    for chrom, start, end, name, strand in read_bed(example_path("peaks.bed")):
        assert chrom in genome, chrom
        assert end <= len(genome[chrom]), f"{name} runs past {chrom}"


def test_bed_input_works_end_to_end_from_the_examples():
    pytest.importorskip("pyfaidx")
    prof = FlexProfiler("DNaseI", window_size=10).from_bed(
        example_path("peaks.bed"), genome=example_path("genome.fa"), width=120)
    assert prof.seqids == [f"peak_{i}" for i in range(1, 6)]
    assert len({len(r) - 1 for r in prof._rows}) == 1
