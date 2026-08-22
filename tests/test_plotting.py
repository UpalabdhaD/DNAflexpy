"""Phase 5: heatmap, metaprofile and trackplot.

The numeric preparation is tested directly; the drawing gets one smoke
assertion each. No image comparison - it would tell you a figure changed
without telling you whether it was ever right.
"""
import numpy as np
import pytest

from DNAflexpy import FlexProfiler
from DNAflexpy.plotting import (
    _bin,
    _check_background,
    _matrix,
    _order,
    _zscale,
    heatmap,
    metaprofile,
    trackplot,
)
from DNAflexpy.results import FlexProfile

SEQS = {
    "s1": "ATGCGTACGTAGCTAGCGTAGCTAGT",
    "s2": "CGTAGCTAGTATGCGTACGTAGCTAG",
    "s3": "TTTTTTTTTTTTTTTTTTTTTTTTTT",
}
GRID = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                 [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]])


def _profile(feature="DNaseI", window_size=10, seqs=SEQS):
    return FlexProfiler(feature, window_size=window_size).profile_seqs(seqs)


# --- _matrix ----------------------------------------------------------------


def test_matrix_is_one_row_per_sequence():
    matrix = _matrix(_profile())
    assert matrix.shape == (3, 26 - 10 + 1)


def test_ragged_rows_raise_and_point_at_from_bed():
    prof = _profile(seqs=["ATGCGTACGTAGCTAGCGTAGCTAGT", "ATGCGTACGTAGC"])
    with pytest.raises(ValueError, match="from_bed"):
        _matrix(prof)


def test_all_rows_empty_blames_the_window_size():
    """Two short sequences give equal row widths, both zero. There is
    nothing to draw, and the window size is why."""
    prof = _profile(window_size=20, seqs=["ACGTA", "ACGTC"])
    assert [len(r) - 1 for r in prof._rows] == [0, 0]
    with pytest.raises(ValueError, match="window_size=20"):
        _matrix(prof)


def test_an_empty_profile_says_so():
    with pytest.raises(ValueError, match="nothing to plot"):
        _matrix(FlexProfile([], feature="gc", window_size=0, kmer_len=2))


# --- _bin -------------------------------------------------------------------


def test_binning_averages_within_equal_width_bins():
    np.testing.assert_allclose(_bin(GRID, 3), [[1.5, 3.5, 5.5], [3.0, 7.0, 11.0]])


def test_no_binning_leaves_the_matrix_alone():
    assert _bin(GRID, None) is GRID


def test_more_bins_than_positions_is_refused():
    with pytest.raises(ValueError, match="exceeds"):
        _bin(GRID, 7)


def test_zero_bins_is_refused():
    with pytest.raises(ValueError, match="nbins must be >= 1"):
        _bin(GRID, 0)


def test_binning_ignores_masked_values():
    grid = np.array([[1.0, np.nan, 3.0, 5.0]])
    np.testing.assert_allclose(_bin(grid, 2), [[1.0, 4.0]])


# --- ordering ---------------------------------------------------------------


def test_bin_then_scale_is_not_the_same_as_scale_then_bin():
    """The spec fixes bin-first so the colours describe what is drawn. A
    shape-only test would pass with the order reversed; these numbers do
    not."""
    correct = _zscale(_bin(GRID, 3), "global")
    reversed_order = _bin(_zscale(GRID, "global"), 3)
    np.testing.assert_allclose(
        correct, [[-1.202246, -0.561048, 0.080150],
                  [-0.721348, 0.561048, 1.843444]], atol=1e-6)
    assert not np.allclose(correct, reversed_order)


# --- _zscale ----------------------------------------------------------------


def test_column_scaling_centres_each_position():
    scaled = _zscale(GRID, "column")
    np.testing.assert_allclose(np.nanmean(scaled, axis=0), np.zeros(6), atol=1e-12)


def test_global_scaling_keeps_positions_different():
    scaled = _zscale(GRID, "global")
    assert not np.allclose(np.nanmean(scaled, axis=0), 0.0)
    np.testing.assert_allclose(np.nanmean(scaled), 0.0, atol=1e-12)


def test_no_scaling_returns_the_values():
    np.testing.assert_allclose(_zscale(GRID, None), GRID)


def test_a_constant_column_does_not_divide_by_zero():
    grid = np.array([[3.0, 1.0], [3.0, 5.0]])
    scaled = _zscale(grid, "column")
    assert np.isfinite(scaled).all()
    np.testing.assert_allclose(scaled[:, 0], [0.0, 0.0])


def test_an_unknown_scale_mode_raises():
    with pytest.raises(ValueError, match="zscale must be"):
        _zscale(GRID, "sideways")


def test_a_background_puts_itself_flat_at_zero():
    """Signal-vs-control: the background is its own reference, so it
    standardises to a mean of 0 at every position."""
    control = np.array([[1.0, 2.0], [3.0, 4.0]])
    scaled_control = _zscale(control, None, background=control)
    np.testing.assert_allclose(np.nanmean(scaled_control, axis=0), [0.0, 0.0])


def test_the_foreground_is_measured_in_background_sds():
    control = np.array([[0.0, 0.0], [2.0, 2.0]])       # mean 1, sd 1
    foreground = np.array([[3.0, 5.0]])
    np.testing.assert_allclose(
        _zscale(foreground, None, background=control), [[2.0, 4.0]])


def test_a_flat_background_column_does_not_divide_by_zero():
    control = np.array([[2.0, 1.0], [2.0, 5.0]])       # column 0 has no spread
    scaled = _zscale(np.array([[4.0, 3.0]]), None, background=control)
    assert np.isfinite(scaled).all()
    np.testing.assert_allclose(scaled[0][0], 2.0)      # difference kept, unscaled


# --- _order -----------------------------------------------------------------


def test_input_order_is_the_identity():
    np.testing.assert_array_equal(_order(GRID, "input"), [0, 1])


def test_std_puts_the_most_variable_row_first():
    grid = np.array([[1.0, 1.0, 1.0], [0.0, 5.0, 10.0]])
    np.testing.assert_array_equal(_order(grid, "std"), [1, 0])


def test_cv_works_on_strictly_positive_data():
    grid = np.array([[10.0, 10.0, 10.0], [1.0, 5.0, 9.0]])
    np.testing.assert_array_equal(_order(grid, "cv"), [1, 0])


def test_cv_refuses_data_that_can_be_negative():
    """std/mean inverts on negative means and blows up near zero. Several
    real features (DNaseI, prop, freeen, mechen) take negative values."""
    with pytest.raises(ValueError, match="strictly positive"):
        _order(np.array([[-1.0, 2.0], [3.0, 4.0]]), "cv")


def test_cv_refuses_a_real_signed_feature():
    with pytest.raises(ValueError, match="order_rows='std'"):
        _profile("DNaseI").heatmap(order_rows="cv")


def test_an_unknown_order_mode_raises():
    with pytest.raises(ValueError, match="order_rows must be"):
        _order(GRID, "alphabetical")


# --- background validation --------------------------------------------------


def test_a_background_of_a_different_feature_is_refused():
    with pytest.raises(ValueError, match="meaningless"):
        _check_background(_profile("DNaseI"), _profile("gc"))


def test_a_background_with_a_different_window_is_refused():
    with pytest.raises(ValueError, match="window_size"):
        _check_background(_profile(window_size=10), _profile(window_size=5))


# --- heatmap ----------------------------------------------------------------


def test_heatmap_returns_an_axes_with_the_right_shape():
    ax = _profile().heatmap()
    image = ax.images[0].get_array()
    assert image.shape == (3, 17)
    assert ax.get_ylabel() == "sequence"
    assert "DNaseI" in ax.get_title()


def test_heatmap_binning_narrows_the_image():
    ax = _profile().heatmap(nbins=5)
    assert ax.images[0].get_array().shape == (3, 5)
    assert "5 bins" in ax.get_xlabel()


def test_heatmap_rows_follow_the_requested_order():
    ax = _profile().heatmap(order_rows="input", zscale=None)
    assert [t.get_text() for t in ax.get_yticklabels()] == ["s1", "s2", "s3"]


def test_a_signed_feature_gets_a_symmetric_diverging_scale():
    ax = _profile("DNaseI").heatmap(zscale=None)
    low, high = ax.images[0].get_clim()
    assert ax.images[0].get_cmap().name == "RdBu_r"
    np.testing.assert_allclose(low, -high)


def test_a_positive_feature_gets_a_sequential_scale():
    ax = _profile("stiffness").heatmap(zscale=None)
    assert ax.images[0].get_cmap().name == "viridis"


# --- metaprofile ------------------------------------------------------------


def test_column_scaling_without_a_background_is_refused():
    with pytest.raises(ValueError, match="flat line at zero"):
        _profile().metaprofile(zscale="column")


def test_that_refusal_is_not_superstition():
    """Column-scale then average down the same column and it cancels to
    zero exactly. This is why the option raises instead of drawing."""
    scaled = _zscale(_matrix(_profile()), "column")
    np.testing.assert_allclose(np.nanmean(scaled, axis=0),
                               np.zeros(scaled.shape[1]), atol=1e-12)


def test_metaprofile_defaults_to_global_scaling():
    ax = _profile().metaprofile()
    assert "whole matrix" in ax.get_ylabel()
    assert len(ax.lines) == 1


def test_metaprofile_draws_both_curves_with_a_background():
    fg = _profile(seqs={"a": SEQS["s1"], "b": SEQS["s2"]})
    bg = _profile(seqs={"c": SEQS["s3"], "d": SEQS["s3"]})
    ax = fg.metaprofile(background=bg)
    assert len(ax.lines) == 3           # foreground, background, zero rule
    assert ax.get_ylabel() == "background SDs"
    assert [t.get_text() for t in ax.get_legend().get_texts()] == [
        "foreground", "background"]


def test_the_background_curve_sits_at_zero():
    fg = _profile(seqs={"a": SEQS["s1"]})
    bg = _profile(seqs={"c": SEQS["s2"], "d": SEQS["s3"]})
    ax = fg.metaprofile(background=bg)
    background_curve = ax.lines[1].get_ydata()
    np.testing.assert_allclose(background_curve, np.zeros_like(background_curve),
                               atol=1e-12)


def test_a_background_of_a_different_length_is_refused():
    fg = _profile(seqs={"a": SEQS["s1"]})
    bg = _profile(seqs={"c": "ATGCGTACGTAGCTAGCG"})
    with pytest.raises(ValueError, match="line up"):
        fg.metaprofile(background=bg)


# --- trackplot --------------------------------------------------------------


def test_trackplot_stacks_one_axis_per_feature():
    profiles = FlexProfiler(["DNaseI", "gc", "stiffness"],
                            window_size=10).profile_seqs(SEQS)
    fig = profiles.trackplot()
    assert len(fig.axes) == 3
    assert [ax.get_ylabel() for ax in fig.axes] == ["DNaseI", "gc", "stiffness"]
    assert "s1" in fig.axes[0].get_title()


def test_trackplot_can_pick_the_sequence():
    profiles = FlexProfiler(["DNaseI", "gc"], window_size=10).profile_seqs(SEQS)
    assert "s2" in profiles.trackplot(seqid="s2").axes[0].get_title()


def test_trackplot_names_the_sequences_it_has():
    profiles = FlexProfiler(["DNaseI", "gc"], window_size=10).profile_seqs(SEQS)
    with pytest.raises(ValueError, match="s1"):
        profiles.trackplot(seqid="nosuch")


def test_trackplot_of_nothing_is_refused():
    from DNAflexpy.core import ProfileSet

    with pytest.raises(ValueError, match="nothing to plot"):
        trackplot(ProfileSet())


# --- the optional extra -----------------------------------------------------


def test_importing_dnaflexpy_does_not_import_matplotlib():
    """matplotlib is an optional extra. If `import DNAflexpy` pulled it in,
    every install without the extra would break."""
    import subprocess
    import sys

    code = (
        "import sys; import DNAflexpy; "
        "sys.exit(1 if 'matplotlib' in sys.modules else 0)"
    )
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0


# --- the two refusals, measured rather than asserted ------------------------


def test_cv_would_invert_the_order_on_negative_means():
    """Why `cv` raises instead of sorting. Row 0 varies far more than row 1,
    but both have negative means, so std/mean comes back negative and the
    sort puts the *least* variable row first."""
    grid = np.array([[-10.0, -1.0], [-5.0, -4.0]])
    spread = np.nanstd(grid, axis=1)
    assert spread[0] > spread[1]                       # row 0 is more variable
    cv = spread / np.nanmean(grid, axis=1)
    assert cv[0] < cv[1]                               # ... but ranks lower
    assert list(np.argsort(-cv)) == [1, 0]             # the inverted order
    with pytest.raises(ValueError, match="strictly positive"):
        _order(grid, "cv")


def test_binning_always_returns_exactly_nbins_columns():
    """An empty bin would give fewer columns than the axis label claims."""
    grid = np.zeros((2, 17))
    for nbins in range(1, 18):
        assert _bin(grid, nbins).shape[1] == nbins


def test_trackplot_refuses_a_set_with_mismatched_sequences():
    """Stacking these would draw rows from different sequences as if they
    were one. `encode` already refuses the same input."""
    from DNAflexpy.core import ProfileSet

    mixed = ProfileSet({
        "gc": FlexProfile([["a", 1.0, 2.0]], feature="gc",
                          window_size=0, kmer_len=2),
        "wedge": FlexProfile([["b", 1.0, 2.0]], feature="wedge",
                             window_size=0, kmer_len=2),
    })
    with pytest.raises(ValueError, match="same sequences"):
        trackplot(mixed)
