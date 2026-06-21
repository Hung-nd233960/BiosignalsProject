"""
Lab 8: The Smoothed Pseudo Wigner-Ville Distribution (SPWVD) and its Tradeoffs (↔ A.8)
--------------------------------------------------------------------------------------
Tests: PWVD vs SPWVD ghost suppression, independent control of time/frequency
       smoothing, and the resolution-vs-ghost tradeoff.

Figures B.51–B.54 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.signal.windows import hann, gaussian

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (
    FS,
    SEED,
    DPI,
    COLORMAP,
    make_tone,
    make_chirp,
    make_time_axis,
    save_figure,
    wigner_ville,
    smoothed_pseudo_wigner_ville,
)

# ============================================================
# Parameters (Table B.11 in the report)
# ============================================================
DURATION_LAB = 2.0  # 2 seconds
N_FFT = 512  # FFT resolution
LAB_NUMBER = 8
FIG_START = 51


# ============================================================
# 1. Experiment A: WVD -> PWVD -> SPWVD Progression
# ============================================================
def experiment_progression():
    """
    Compare raw WVD, PWVD (lag-smoothed only), and SPWVD (time and lag smoothed)
    on the multi-component signal (chirp + tone) from Lab 7.
    """
    print("\nExperiment A: WVD -> PWVD -> SPWVD Progression")
    print("-" * 40)

    # Component 1: Linear Chirp (10 Hz to 90 Hz)
    f0 = 10.0
    f_end = 90.0
    mu = (f_end - f0) / DURATION_LAB  # 40 Hz/s
    x_chirp, _, t = make_chirp(f0, mu, A=1.0, duration=DURATION_LAB, fs=FS)

    # Component 2: Stationary Tone at 40 Hz
    f_tone = 40.0
    x_tone, _, _ = make_tone(f_tone, A=1.0, duration=DURATION_LAB, fs=FS)

    # Combined signal
    x = x_chirp + x_tone

    # --- 1. Compute Raw WVD ---
    wvd, _, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- 2. Compute PWVD (Frequency/Lag smoothed only) ---
    # PWVD is SPWVD with time smoothing window g of length 1 (no time smoothing)
    h_lag = hann(51, sym=True)  # lag window
    g_time_none = np.array([1.0])  # no time smoothing
    pwvd, _, _ = smoothed_pseudo_wigner_ville(x, FS, h=h_lag, g=g_time_none, n_fft=N_FFT)

    # --- 3. Compute SPWVD (Both time and lag smoothed) ---
    g_time_hann = hann(15, sym=True)  # time smoothing window
    spwvd, t_spwvd, f_spwvd = smoothed_pseudo_wigner_ville(x, FS, h=h_lag, g=g_time_hann, n_fft=N_FFT)

    # --- Plot the progression ---
    fig, (ax_wvd, ax_pwvd, ax_spwvd) = plt.subplots(3, 1, figsize=(12, 10), sharex=True, sharey=True)

    # 1. Raw WVD
    im_wvd = ax_wvd.pcolormesh(t, f_wvd, wvd, shading="gouraud", cmap=COLORMAP)
    ax_wvd.set_ylabel("Frequency (Hz)")
    ax_wvd.set_title("1. Raw Wigner-Ville Distribution (sharp lines, heavy midpoint cross-terms)")
    ax_wvd.set_title(f"Figure B.{FIG_START}", loc="left", fontsize=9, fontstyle="italic")
    fig.colorbar(im_wvd, ax=ax_wvd, label="Power")

    # 2. PWVD
    im_pwvd = ax_pwvd.pcolormesh(t, f_wvd, pwvd, shading="gouraud", cmap=COLORMAP)
    ax_pwvd.set_ylabel("Frequency (Hz)")
    ax_pwvd.set_title("2. Pseudo Wigner-Ville Distribution (lag window h=51, frequency-oscillating ghosts reduced)")
    fig.colorbar(im_pwvd, ax=ax_pwvd, label="Power")

    # 3. SPWVD
    im_spwvd = ax_spwvd.pcolormesh(t_spwvd, f_spwvd, spwvd, shading="gouraud", cmap=COLORMAP)
    ax_spwvd.set_ylabel("Frequency (Hz)")
    ax_spwvd.set_xlabel("Time (s)")
    ax_spwvd.set_title("3. Smoothed Pseudo WVD (g=15, h=51, both time and frequency ghosts suppressed)")
    fig.colorbar(im_spwvd, ax=ax_spwvd, label="Power")

    # Set limits
    ax_wvd.set_ylim(0, 100)
    ax_wvd.set_xlim(0, DURATION_LAB)

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# 2. Experiment B: Time vs. Frequency Smoothing Window Sweep
# ============================================================
def experiment_smoothing_sweep():
    """
    Sweep different window lengths for h (lag window) and g (time window)
    to demonstrate independent control over the time-frequency smoothing tradeoff.
    """
    print("\nExperiment B: Time vs. Frequency Smoothing Window Sweep")
    print("-" * 40)

    # Component 1: Linear Chirp (10 Hz to 90 Hz)
    f0 = 10.0
    f_end = 90.0
    mu = (f_end - f0) / DURATION_LAB  # 40 Hz/s
    x_chirp, _, t = make_chirp(f0, mu, A=1.0, duration=DURATION_LAB, fs=FS)

    # Component 2: Stationary Tone at 40 Hz
    f_tone = 40.0
    x_tone, _, _ = make_tone(f_tone, A=1.0, duration=DURATION_LAB, fs=FS)
    x = x_chirp + x_tone

    # We will test two extreme cases of smoothing to compare:
    # Case 1: Strong frequency smoothing (h=101), weak time smoothing (g=5)
    # Case 2: Weak frequency smoothing (h=25), strong time smoothing (g=31)
    
    h1 = hann(101, sym=True)
    g1 = hann(5, sym=True)
    spwvd1, t1, f1 = smoothed_pseudo_wigner_ville(x, FS, h=h1, g=g1, n_fft=N_FFT)

    h2 = hann(25, sym=True)
    g2 = hann(31, sym=True)
    spwvd2, t2, f2 = smoothed_pseudo_wigner_ville(x, FS, h=h2, g=g2, n_fft=N_FFT)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)

    # Plot Case 1
    im1 = ax1.pcolormesh(t1, f1, spwvd1, shading="gouraud", cmap=COLORMAP)
    ax1.set_ylabel("Frequency (Hz)")
    ax1.set_title("Case 1: Strong Frequency Smoothing, Weak Time Smoothing (h=101, g=5)")
    ax1.set_title(f"Figure B.{FIG_START + 1}", loc="left", fontsize=9, fontstyle="italic")
    fig.colorbar(im1, ax=ax1, label="Power")

    # Plot Case 2
    im2 = ax2.pcolormesh(t2, f2, spwvd2, shading="gouraud", cmap=COLORMAP)
    ax2.set_ylabel("Frequency (Hz)")
    ax2.set_xlabel("Time (s)")
    ax2.set_title("Case 2: Weak Frequency Smoothing, Strong Time Smoothing (h=25, g=31)")
    fig.colorbar(im2, ax=ax2, label="Power")

    ax1.set_ylim(0, 100)
    ax1.set_xlim(0, DURATION_LAB)

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# 3. Run Lab 8
# ============================================================
def run_lab8(save=False):
    """Run all Lab 8 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 8: The Smoothed Pseudo Wigner-Ville Distribution (SPWVD) and its Tradeoffs")
    print("=" * 60)

    fig_comp_a = experiment_progression()  # Experiment A
    fig_comp_b = experiment_smoothing_sweep()  # Experiment B

    all_figs = [
        (fig_comp_a, True),   # B.51: WVD -> PWVD -> SPWVD progression
        (fig_comp_b, True),   # B.52: Time vs Frequency smoothing sweep
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
    run_lab8(save=True)
