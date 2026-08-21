"""Figures: heatmap, metaprofile, trackplot.

matplotlib is imported inside the drawing functions, never at module level,
so `import DNAflexpy` still works without the optional `plot` extra.

The numeric preparation is deliberately split from the drawing. `_matrix`,
`_bin`, `_zscale` and `_order` return the numbers that are actually plotted,
and the tests assert on those. Pixel comparison would tell you a figure
changed but not whether it was ever right.
"""
from __future__ import annotations

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError(
            "plotting needs matplotlib, which is an optional extra. "
            "Install it with: pip install 'DNAflexpy[plot]'"
        ) from None
    return plt


def _matrix(profile) -> np.ndarray:
    """Profile values as a 2-D array, one row per sequence.

    Equal row widths is the right check here, unlike `encode`, which checks
    sequence length: the columns of these figures *are* profile positions.
    """
    widths = {len(row) - 1 for row in profile._rows}
    if not widths:
        raise ValueError("there is nothing to plot: the profile has no sequences")
    if len(widths) > 1:
        raise ValueError(
            f"plotting needs one value per position for every sequence, but "
            f"the rows have widths {sorted(widths)}. Position-wise comparison "
            "is only meaningful for aligned sequences: use "
            "from_bed(..., width=N) to cut fixed-width windows, or trim the "
            "sequences yourself."
        )
    if widths == {0}:
        raise ValueError(
            f"every sequence is shorter than window_size={profile.window_size}, "
            "so there are no values to plot. Lower the window size."
        )
    return np.array([row[1:] for row in profile._rows], dtype=float)


def _bin(matrix: np.ndarray, nbins: int | None) -> np.ndarray:
    """Average positions into `nbins` equal-width bins.

    Applied *before* z-scaling, so the statistics describe the values that
    are actually drawn. The other order gives a different figure.
    """
    if nbins is None:
        return matrix
    positions = matrix.shape[1]
    if nbins < 1:
        raise ValueError(f"nbins must be >= 1, got {nbins}")
    if nbins > positions:
        raise ValueError(
            f"nbins={nbins} exceeds the {positions} position(s) available; "
            "binning can only make the matrix narrower"
        )
    edges = np.linspace(0, positions, nbins + 1).astype(int)
    with np.errstate(invalid="ignore"):
        return np.column_stack([
            np.nanmean(matrix[:, lo:hi], axis=1) if hi > lo else matrix[:, lo:lo]
            for lo, hi in zip(edges[:-1], edges[1:])
        ])


def _zscale(matrix: np.ndarray, mode, background: np.ndarray | None = None) -> np.ndarray:
    """Standardise the matrix.

    `"column"` scales each position across sequences. `"global"` uses one
    mean and standard deviation for the whole matrix, preserving the
    differences between positions. `None` leaves the values alone.

    With a `background`, this becomes a signal-vs-control z-score: every
    column is standardised against the *background's* mean and standard
    deviation, so the background sits flat at ~0 and the foreground's
    departure from it is readable in units of background standard
    deviations.
    """
    if background is not None:
        centre = np.nanmean(background, axis=0)
        spread = np.nanstd(background, axis=0)
        # A background column with no variation cannot scale anything. Leave
        # the difference unscaled rather than dividing by zero.
        spread = np.where(spread > 0, spread, 1.0)
        return (matrix - centre) / spread
    if mode is None:
        return matrix
    if mode == "column":
        centre = np.nanmean(matrix, axis=0)
        spread = np.nanstd(matrix, axis=0)
        spread = np.where(spread > 0, spread, 1.0)
        return (matrix - centre) / spread
    if mode == "global":
        centre = np.nanmean(matrix)
        spread = np.nanstd(matrix)
        return (matrix - centre) / (spread if spread > 0 else 1.0)
    raise ValueError(
        f"zscale must be 'column', 'global' or None, got {mode!r}"
    )


def _order(matrix: np.ndarray, mode: str) -> np.ndarray:
    """Row order: as given, by spread, or by coefficient of variation.

    `"cv"` is `std/mean`, which only means anything for a strictly positive
    quantity. `DNaseI` runs from -0.28 to 0.194, so a row mean near zero
    sends CV to infinity, and `prop`, `freeen` and `mechen` are entirely
    negative, which flips the sort. Rather than silently produce a wrong
    order it raises, and `"std"` is the default.
    """
    if mode == "input":
        return np.arange(matrix.shape[0])
    if mode == "std":
        return np.argsort(-np.nan_to_num(np.nanstd(matrix, axis=1)))
    if mode == "cv":
        if not np.all(np.nan_to_num(matrix, nan=1.0) > 0):
            raise ValueError(
                "order_rows='cv' needs strictly positive values: the "
                "coefficient of variation is std/mean, which is meaningless "
                "when the mean can be zero or negative. Several features "
                "(DNaseI, prop, freeen, mechen) take negative values. Use "
                "order_rows='std' instead."
            )
        means = np.nanmean(matrix, axis=1)
        return np.argsort(-(np.nanstd(matrix, axis=1) / means))
    raise ValueError(
        f"order_rows must be 'input', 'std' or 'cv', got {mode!r}"
    )


def _check_background(profile, background) -> None:
    """A background is only a control if it measures the same thing."""
    if background.feature != profile.feature:
        raise ValueError(
            f"background measures {background.feature!r} but the foreground "
            f"measures {profile.feature!r}; comparing them is meaningless"
        )
    if background.window_size != profile.window_size:
        raise ValueError(
            f"background used window_size={background.window_size} but the "
            f"foreground used {profile.window_size}; a position means a "
            "different span in each"
        )


def _diverging(matrix: np.ndarray) -> bool:
    """Signed data wants a colormap centred at zero."""
    return bool(np.nanmin(matrix) < 0 < np.nanmax(matrix))


def heatmap(profile, nbins=None, order_rows="std", zscale="column",
            cmap=None, ax=None):
    """Rows are sequences, columns are positions.

    Strictly one feature per figure: `stiffness` runs to 5500 and `DNaseI`
    to 0.194, so they cannot share a colour scale.

    `nbins` averages positions into equal-width bins, which keeps a large
    matrix renderable. Binning happens before scaling, so the colours
    describe the values actually drawn.

    `order_rows="std"` puts the most variable sequences first. `"cv"` is
    available but raises on features that take negative values; `"input"`
    keeps file order.

    Returns the `Axes`, so the caller decides whether to save or show.
    """
    plt = _require_matplotlib()

    matrix = _bin(_matrix(profile), nbins)
    matrix = _zscale(matrix, zscale)
    order = _order(matrix, order_rows)
    matrix = matrix[order]
    labels = [profile.seqids[i] for i in order]

    if ax is None:
        _, ax = plt.subplots(figsize=(10, max(2.0, 0.25 * len(labels))))
    if cmap is None:
        cmap = "RdBu_r" if _diverging(matrix) else "viridis"

    limit = np.nanmax(np.abs(matrix)) if _diverging(matrix) else None
    image = ax.imshow(
        matrix, aspect="auto", cmap=cmap, interpolation="nearest",
        vmin=-limit if limit else None, vmax=limit if limit else None,
    )
    ax.set_xlabel("position" if nbins is None else f"position ({nbins} bins)")
    ax.set_ylabel("sequence")
    ax.set_title(f"{profile.feature} (window {profile.window_size})")
    if len(labels) <= 40:
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
    else:
        ax.set_yticks([])
    ax.figure.colorbar(image, ax=ax, label=_scale_label(zscale))
    return ax


def _scale_label(zscale) -> str:
    return {
        "column": "z-score (per position)",
        "global": "z-score (whole matrix)",
        None: "value",
    }.get(zscale, "value")


def metaprofile(profile, background=None, nbins=None, zscale="auto", ax=None):
    """The position-wise average across sequences, as a line.

    **Without a background, `zscale="column"` raises.** Scaling each column
    across sequences and then averaging down that same column gives exactly
    0 at every position, by construction - the line would be flat at zero
    and tell you nothing. The default is `"global"`, which keeps the
    differences between positions.

    **With a background** - a second profile built from a control set - the
    default becomes a signal-vs-control z-score: each column is
    standardised against the background's own mean and standard deviation.
    The background then sits flat at ~0 and the foreground's departure from
    it is a real quantity, in units of background standard deviations.

    Returns the `Axes`.
    """
    plt = _require_matplotlib()

    matrix = _bin(_matrix(profile), nbins)
    control = None
    if background is not None:
        _check_background(profile, background)
        control = _bin(_matrix(background), nbins)
        if control.shape[1] != matrix.shape[1]:
            raise ValueError(
                f"background has {control.shape[1]} position(s) but the "
                f"foreground has {matrix.shape[1]}; they must line up"
            )

    if zscale == "auto":
        zscale = None if background is not None else "global"
    elif zscale == "column" and background is None:
        raise ValueError(
            "zscale='column' with no background would draw a flat line at "
            "zero. Standardising each position across sequences and then "
            "averaging down that same position cancels exactly, by "
            "construction. Use zscale='global', or supply background=."
        )

    if control is not None:
        scaled = _zscale(matrix, None, background=control)
        scaled_bg = _zscale(control, None, background=control)
    else:
        scaled = _zscale(matrix, zscale)
        scaled_bg = None

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 3.5))

    positions = np.arange(1, scaled.shape[1] + 1)
    ax.plot(positions, np.nanmean(scaled, axis=0), label="foreground", lw=1.8)
    if scaled_bg is not None:
        ax.plot(positions, np.nanmean(scaled_bg, axis=0), label="background",
                lw=1.2, color="grey")
        ax.axhline(0, color="black", lw=0.6, ls=":")
        ax.legend(frameon=False)
        ax.set_ylabel("background SDs")
    else:
        ax.set_ylabel(_scale_label(zscale))

    ax.set_xlabel("position" if nbins is None else f"position ({nbins} bins)")
    ax.set_title(f"{profile.feature} (window {profile.window_size})")
    return ax


def trackplot(profiles, seqid=None, nbins=None, figsize=None):
    """One sequence, every feature stacked on a shared x-axis.

    Features have different units, so each gets its own y-axis rather than
    a shared one. Returns the `Figure`.
    """
    plt = _require_matplotlib()

    if not profiles:
        raise ValueError("there is nothing to plot: no features were given")
    features = list(profiles)
    first = profiles[features[0]]
    if seqid is None:
        seqid = first.seqids[0]
    if seqid not in first.seqids:
        raise ValueError(
            f"sequence {seqid!r} is not in this profile. Available: "
            f"{', '.join(first.seqids[:5])}"
            + (" ..." if len(first.seqids) > 5 else "")
        )
    row = first.seqids.index(seqid)

    if figsize is None:
        figsize = (10, 1.6 * len(features))
    fig, axes = plt.subplots(len(features), 1, sharex=True, figsize=figsize,
                             squeeze=False)
    for ax, feature in zip(axes[:, 0], features):
        values = _bin(_matrix(profiles[feature]), nbins)[row]
        ax.plot(np.arange(1, len(values) + 1), values, lw=1.4)
        ax.set_ylabel(feature, fontsize=8)
        ax.margins(x=0)
    axes[-1, 0].set_xlabel("position" if nbins is None else f"position ({nbins} bins)")
    axes[0, 0].set_title(f"{seqid} (window {first.window_size})")
    fig.tight_layout()
    return fig
