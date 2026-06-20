"""
Plotting utilities — enforce CLAUDE.md graph standards:
  - 300+ DPI, perceptually uniform colormaps, physical units on all axes
  - Dual-stack rule: linear scale first, log/dB scale second
  - Figure numbering and export
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from .config import DPI, COLORMAP, FIGURE_FORMATS, GRAPH_DIR, FORBIDDEN_COLORMAPS


def _apply_defaults(fig):
    """Apply project-wide figure defaults."""
    fig.set_dpi(DPI)  # enforce minimum DPI
    fig.tight_layout()  # prevent label clipping


def _validate_colormap(cmap_name):
    """Reject forbidden colormaps (jet, rainbow, hsv)."""
    assert (
        cmap_name not in FORBIDDEN_COLORMAPS
    ), f"Colormap '{cmap_name}' is forbidden — use viridis, inferno, plasma, or magma"


def save_figure(fig, lab_number, fig_id, raster_only=False):
    """
    Save a figure to results/graphs/labN/ in PNG (always) and SVG (unless raster_only).
    Use raster_only=True for spectrograms and heatmaps — SVG of dense 2D data
    produces massive files that crash renderers.

    Parameters
    ----------
    fig         : matplotlib Figure
    lab_number  : int — lab number (1–8)
    fig_id      : str — figure identifier (e.g. "01", "02a")
    raster_only : bool — if True, skip SVG export (for spectrograms/heatmaps)
    """
    lab_dir = os.path.join(GRAPH_DIR, f"lab{lab_number}")  # target directory
    os.makedirs(lab_dir, exist_ok=True)  # ensure it exists
    formats = ["png"] if raster_only else FIGURE_FORMATS  # skip SVG for heatmaps
    paths = []  # collect saved paths
    for fmt in formats:  # save in each format
        filename = f"figure_B_{fig_id}.{fmt}"  # construct filename
        filepath = os.path.join(lab_dir, filename)  # full path
        fig.savefig(filepath, dpi=DPI, bbox_inches="tight")  # save
        print(f"Saved: {filepath}")  # confirm to user
        paths.append(filepath)
    return paths


# --- Time-domain plot ---
def plot_time_domain(
    t, x, xlabel="Time (s)", ylabel="Amplitude", title=None, fig_id=None, t_range=None
):
    """
    Plot signal in the time domain. Always shown before any spectral analysis.

    Parameters
    ----------
    t       : ndarray — time axis in seconds
    x       : ndarray — signal samples
    xlabel  : str — x-axis label with units
    ylabel  : str — y-axis label with units
    title   : str — figure title (optional)
    fig_id  : str — figure number for labeling (e.g. "Figure B.1")
    t_range : tuple — (t_start, t_end) to zoom into a time range (optional)

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 3))  # wide, short aspect
    if t_range is not None:  # zoom to range if specified
        mask = (t >= t_range[0]) & (t <= t_range[1])  # boolean mask
        ax.plot(t[mask], x[mask], linewidth=0.5)  # plot zoomed segment
    else:
        ax.plot(t, x, linewidth=0.5)  # plot full signal
    ax.set_xlabel(xlabel)  # x label with units
    ax.set_ylabel(ylabel)  # y label with units
    if title:
        ax.set_title(title)  # optional title
    if fig_id:
        ax.set_title(
            fig_id, loc="left", fontsize=9, fontstyle="italic"  # figure number top-left
        )
    ax.grid(True, alpha=0.3)  # light grid
    _apply_defaults(fig)  # DPI and layout
    return fig, ax


# --- Dual-stack spectrum plot (linear + log/dB) ---
def plot_dual_stack_spectrum(
    freqs,
    power,
    xlabel="Frequency (Hz)",
    ylabel_linear="Power",
    ylabel_db="Power (dB)",
    title=None,
    fig_id=None,
    f_range=None,
):
    """
    Dual-stack plot: linear scale on top, log/dB scale on bottom.
    Linear is primary (physical units); dB is secondary (dynamic range).

    Parameters
    ----------
    freqs        : ndarray — frequency axis in Hz
    power        : ndarray — power values (linear scale, physical units)
    xlabel       : str — x-axis label with units
    ylabel_linear: str — y-axis label for linear plot
    ylabel_db    : str — y-axis label for dB plot
    title        : str — overall title (optional)
    fig_id       : str — figure number for labeling (e.g. "Figure B.2")
    f_range      : tuple — (f_min, f_max) to zoom frequency axis (optional)

    Returns
    -------
    fig, (ax_lin, ax_db) : matplotlib Figure and two Axes
    """
    fig, (ax_lin, ax_db) = plt.subplots(
        2, 1, figsize=(12, 6), sharex=True  # two vertically stacked
    )
    # --- Top: linear scale (primary) ---
    ax_lin.plot(freqs, power, linewidth=0.5)  # linear power
    ax_lin.set_ylabel(ylabel_linear)  # label with physical units
    ax_lin.grid(True, alpha=0.3)  # light grid
    if title:
        ax_lin.set_title(title)  # overall title on top panel
    if fig_id:
        ax_lin.set_title(
            fig_id, loc="left", fontsize=9, fontstyle="italic"  # figure number
        )

    # --- Bottom: log/dB scale (secondary) ---
    power_db = 10 * np.log10(np.maximum(power, 1e-20))  # dB, floor at -200 dB
    ax_db.plot(freqs, power_db, linewidth=0.5)  # dB power
    ax_db.set_ylabel(ylabel_db)  # label with dB units
    ax_db.set_xlabel(xlabel)  # frequency label on bottom
    ax_db.grid(True, alpha=0.3)  # light grid

    if f_range is not None:  # zoom frequency axis
        ax_lin.set_xlim(f_range)  # apply to both (sharex)

    _apply_defaults(fig)  # DPI and layout
    return fig, (ax_lin, ax_db)


# --- Spectrogram plot (dual-stack: linear + dB) ---
def plot_spectrogram(
    t_stft,
    f_stft,
    Sxx,
    xlabel="Time (s)",
    ylabel="Frequency (Hz)",
    title=None,
    fig_id=None,
    cmap=COLORMAP,
    f_range=None,
    t_range=None,
):
    """
    Dual-stack spectrogram: linear power on top, dB on bottom.

    Parameters
    ----------
    t_stft  : ndarray — time axis of the STFT (s)
    f_stft  : ndarray — frequency axis of the STFT (Hz)
    Sxx     : 2D ndarray — spectrogram power (shape: freq × time)
    xlabel  : str — x-axis label
    ylabel  : str — y-axis label
    title   : str — overall title (optional)
    fig_id  : str — figure number (optional)
    cmap    : str — colormap name (must be perceptually uniform)
    f_range : tuple — (f_min, f_max) frequency zoom (optional)
    t_range : tuple — (t_min, t_max) time zoom (optional)

    Returns
    -------
    fig, (ax_lin, ax_db) : matplotlib Figure and two Axes
    """
    _validate_colormap(cmap)  # reject jet/rainbow

    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)

    # --- Top: linear scale ---
    im_lin = ax_lin.pcolormesh(
        t_stft, f_stft, Sxx, shading="gouraud", cmap=cmap  # linear power
    )
    ax_lin.set_ylabel(ylabel)  # frequency label
    fig.colorbar(im_lin, ax=ax_lin, label="Power (linear)")  # colorbar with units
    if title:
        ax_lin.set_title(title)
    if fig_id:
        ax_lin.set_title(fig_id, loc="left", fontsize=9, fontstyle="italic")

    # --- Bottom: dB scale ---
    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))  # dB conversion
    im_db = ax_db.pcolormesh(
        t_stft, f_stft, Sxx_db, shading="gouraud", cmap=cmap  # dB power
    )
    ax_db.set_ylabel(ylabel)  # frequency label
    ax_db.set_xlabel(xlabel)  # time label on bottom
    fig.colorbar(im_db, ax=ax_db, label="Power (dB)")  # colorbar with units

    if f_range is not None:  # zoom frequency
        ax_lin.set_ylim(f_range)
    if t_range is not None:  # zoom time
        ax_lin.set_xlim(t_range)

    _apply_defaults(fig)  # DPI and layout
    return fig, (ax_lin, ax_db)
