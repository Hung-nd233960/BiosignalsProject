"""
Lab 5: Two-Tone Resolution on the Spectrogram (↔ A.3, A.5, Lab 3)
-------------------------------------------------------------------
Confirms the resolution limit Δf_min ≈ β·fs/M visually:
two stationary tones on a spectrogram — too close = one line,
far enough = two lines.

Figures B.26–B.38 in the report.
"""

import numpy as np  # numerical computing
import matplotlib.pyplot as plt  # plotting
from scipy import signal as sp_signal  # STFT
import sys
import os  # path manipulation

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (  # project utilities
    FS,
    DURATION,
    DPI,
    COLORMAP,
    make_mixed_tones,
    plot_time_domain,
    save_figure,
)
from src.common.windows import rectangular, hann, blackman

# ============================================================
# Parameters (Table B.8 in the report)
# ============================================================
F1 = 10.0  # first tone frequency (Hz) — alpha band
M = 256  # window length (samples) — typical EEG epoch
DURATION_LAB = DURATION  # 1200 s signal duration
LAB_NUMBER = 5  # lab number for figure saving
FIG_START = 26  # first figure number (B.26)

# Resolution limits per window at M=256, fs=250
# Δf_min = β · fs / M
WINDOWS = [
    ("Rectangular", "hann", rectangular, 1, 1 * FS / M),  # 0.977 Hz
    ("Hann", "hann", hann, 2, 2 * FS / M),  # 1.953 Hz
    ("Blackman", "hann", blackman, 3, 3 * FS / M),  # 2.930 Hz
]

# Tone separations to demonstrate (Hz)
SEPARATIONS = [0.5, 2.0, 3.0, 5.0]


# ============================================================
# 1. Helper: compute spectrogram with a specific window array
# ============================================================
def compute_stft_custom_window(x, window_array, fs=FS):
    """
    Compute spectrogram using a custom window array.

    Parameters
    ----------
    x            : ndarray — input signal
    window_array : ndarray — window samples
    fs           : int — sampling frequency

    Returns
    -------
    t_stft : ndarray — time axis (s)
    f_stft : ndarray — frequency axis (Hz)
    Sxx    : 2D ndarray — power spectrogram
    """
    nperseg = len(window_array)  # window length
    noverlap = nperseg // 2  # 50% overlap
    f_stft, t_stft, Sxx = sp_signal.spectrogram(
        x,
        fs=fs,
        window=window_array,  # pass actual window samples
        nperseg=nperseg,
        noverlap=noverlap,
    )
    return t_stft, f_stft, Sxx


# ============================================================
# 2. Plot one spectrogram panel (single, not dual-stack)
# ============================================================
def plot_single_spectrogram(ax, t_stft, f_stft, Sxx, title="", f_range=(0, 20)):
    """
    Plot a single spectrogram on a given axes (for grid layout).

    Parameters
    ----------
    ax      : matplotlib Axes
    t_stft  : ndarray — time axis
    f_stft  : ndarray — frequency axis
    Sxx     : 2D ndarray — power spectrogram
    title   : str — subplot title
    f_range : tuple — frequency zoom
    """
    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))  # dB scale
    ax.pcolormesh(t_stft, f_stft, Sxx_db, shading="gouraud", cmap=COLORMAP)
    ax.set_ylim(f_range)
    ax.set_ylabel("Freq (Hz)")
    ax.set_title(title, fontsize=9)


# ============================================================
# 3. Main experiment: separation × window grid
# ============================================================
def experiment_resolution_grid():
    """
    For each separation, generate two-tone signal and compute
    spectrogram with each window. Arrange as a grid:
    rows = separations, columns = windows.
    """
    print("\nTwo-tone resolution grid")
    print("-" * 40)

    n_seps = len(SEPARATIONS)  # number of separations
    n_wins = len(WINDOWS)  # number of windows

    # --- Print resolution limits ---
    for name, _, _, beta, df_min in WINDOWS:
        print(f"  {name}: β={beta}, Δf_min = {df_min:.3f} Hz")

    # --- Time-domain figure of the beat pattern at each separation ---
    fig_td, axes_td = plt.subplots(n_seps, 1, figsize=(12, 2.5 * n_seps), sharex=True)
    for j, sep in enumerate(SEPARATIONS):
        f2 = F1 + sep  # second tone frequency
        x, _, t = make_mixed_tones([F1, f2], amplitudes=[1.0, 1.0], duration=10.0)
        axes_td[j].plot(t, x, linewidth=0.5)  # plot 10 s of beat pattern
        axes_td[j].set_ylabel("Amplitude")
        axes_td[j].set_title(
            f"Δ = {sep:.1f} Hz  (f₁={F1:.0f}, f₂={f2:.1f} Hz)",
            fontsize=9,
        )
        axes_td[j].grid(True, alpha=0.3)
    axes_td[-1].set_xlabel("Time (s)")
    fig_td.suptitle(
        f"Figure B.{FIG_START} — Two-tone time domain at each separation",
        fontsize=12,
        fontstyle="italic",
    )
    fig_td.set_dpi(DPI)
    fig_td.tight_layout()

    # --- Spectrogram grid: one figure per separation ---
    figs_grid = []
    for j, sep in enumerate(SEPARATIONS):
        f2 = F1 + sep  # second tone
        x, _, t = make_mixed_tones(
            [F1, f2], amplitudes=[1.0, 1.0], duration=DURATION_LAB
        )

        fig, axes = plt.subplots(1, n_wins, figsize=(5 * n_wins, 4), sharey=True)

        for i, (name, _, win_func, beta, df_min) in enumerate(WINDOWS):
            w = win_func(M)  # generate window
            t_stft, f_stft, Sxx = compute_stft_custom_window(x, w)

            resolved = "YES" if sep >= df_min else "NO"  # prediction
            plot_single_spectrogram(
                axes[i],
                t_stft,
                f_stft,
                Sxx,
                title=f"{name} (β={beta})\nΔf_min={df_min:.2f} Hz → {resolved}",
                f_range=(F1 - 3, F1 + sep + 3),  # zoom around the tones
            )
            if i == 0:
                axes[i].set_ylabel("Frequency (Hz)")
            axes[i].set_xlabel("Time (s)")

            # --- Mark the true tone frequencies ---
            axes[i].axhline(F1, color="white", ls="--", lw=0.5, alpha=0.6)
            axes[i].axhline(f2, color="white", ls="--", lw=0.5, alpha=0.6)

        fig_num = FIG_START + 1 + j
        fig.suptitle(
            f"Figure B.{fig_num} — Δ = {sep:.1f} Hz "
            f"(f₁={F1:.0f}, f₂={f2:.1f} Hz), M={M}",
            fontsize=12,
            fontstyle="italic",
        )
        fig.set_dpi(DPI)
        fig.tight_layout()
        figs_grid.append(fig)

        print(f"\n  Δ = {sep:.1f} Hz:")
        for name, _, _, beta, df_min in WINDOWS:
            resolved = "resolved" if sep >= df_min else "merged"
            print(f"    {name} (β={beta}): Δf_min={df_min:.3f} Hz → {resolved}")

    return fig_td, figs_grid


# ============================================================
# 4. Run all experiments
# ============================================================
def run_lab5(save=False):
    """Run all Lab 5 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 5: Two-Tone Resolution on the Spectrogram")
    print("=" * 60)

    fig_td, figs_grid = experiment_resolution_grid()

    # --- Tag each figure: (fig, raster_only?) ---
    all_figs = [
        (fig_td, False),  # B.26: time domain (line plot)
        *[(f, True) for f in figs_grid],  # B.27–B.30: spectrogram grids
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(
                fig,
                lab_number=LAB_NUMBER,
                fig_id=fig_id,
                raster_only=raster_only,
            )

    plt.show()
    return [fig for fig, _ in all_figs]


if __name__ == "__main__":
    run_lab5(save=True)
