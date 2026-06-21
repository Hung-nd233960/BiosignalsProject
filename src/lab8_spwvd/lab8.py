"""
Lab 8: The SPWVD and its Tradeoffs (↔ A.8)
--------------------------------------------
Tests: WVD→PWVD→SPWVD progression, independent two-knob control,
       and the duality: PWVD kills frequency-oscillating ghosts (impulses),
       SPWVD kills both (time + frequency oscillating).

Figures B.45-B.49 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (                                   # project utilities
    FS,
    DPI,
    COLORMAP,
    make_tone,
    make_chirp,
    make_transient,
    plot_time_domain,
    save_figure,
    wigner_ville,
    smoothed_pseudo_wigner_ville,
)
from src.common.windows import hann                        # our window, not scipy's

# ============================================================
# Parameters
# ============================================================
# Duration is 2.0 s — same justification as Lab 7 (WVD is O(N²))
DURATION_LAB = 2.0        # signal duration (s)
LAB_NUMBER = 8            # lab number for figure saving
FIG_START = 45            # first figure number (B.45), continuing from Lab 7

# Signal: same chirp + tone from Lab 7 Experiment B
F0_CHIRP = 10.0           # chirp start (Hz)
F_END_CHIRP = 90.0        # chirp end (Hz)
MU_CHIRP = (F_END_CHIRP - F0_CHIRP) / DURATION_LAB  # 40 Hz/s
F_TONE = 40.0             # stationary tone (Hz)

# Window parameters
N_FFT = 512               # frequency bins


# ============================================================
# Experiment A: WVD → PWVD → SPWVD Progression
# ============================================================
def experiment_progression():
    """
    Same chirp+tone signal, processed with WVD, PWVD, and SPWVD.
    Shows step-by-step ghost suppression.
    """
    print("\nExperiment A: WVD → PWVD → SPWVD Progression")
    print("-" * 40)

    # --- Generate signal (same as Lab 7 Exp B) ---
    x_chirp, _, t = make_chirp(F0_CHIRP, MU_CHIRP,         # 10-90 Hz chirp
                                A=1.0, duration=DURATION_LAB)
    x_tone, _, _ = make_tone(F_TONE, A=1.0,                # 40 Hz tone
                              duration=DURATION_LAB)
    x = x_chirp + x_tone                                   # combined signal

    # --- Time-domain (always first) ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Chirp ({F0_CHIRP}-{F_END_CHIRP} Hz) + Tone ({F_TONE} Hz) - same as Lab 7",
        fig_id=f"Figure B.{FIG_START}",
    )

    # --- 1. Raw WVD ---
    wvd_raw, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- 2. PWVD (lag window only, no time smoothing) ---
    h_lag = hann(51)                                       # lag window: Hann 51 (0.204 s)
    g_none = np.array([1.0])                               # no time smoothing (length 1)
    pwvd, _, _ = smoothed_pseudo_wigner_ville(             # PWVD = SPWVD with g=1
        x, FS, h=h_lag, g=g_none, n_fft=N_FFT,
    )

    # --- 3. SPWVD (both windows) ---
    g_time = hann(15)                                      # time window: Hann 15 (0.060 s)
    spwvd, t_spwvd, f_spwvd = smoothed_pseudo_wigner_ville(
        x, FS, h=h_lag, g=g_time, n_fft=N_FFT,
    )

    print(f"  Signal: chirp {F0_CHIRP}-{F_END_CHIRP} Hz + tone {F_TONE} Hz")
    print(f"  Cross-term type: TIME-oscillating (components separated in frequency)")
    print(f"  PWVD lag window h: Hann 51 ({51/FS:.3f} s) - smooths FREQUENCY axis")
    print(f"  SPWVD time window g: Hann 15 ({15/FS:.3f} s) - smooths TIME axis")

    # --- Plot progression: 3 rows, dual-stack (linear + dB) ---
    fig, axes = plt.subplots(3, 2, figsize=(14, 12),
                              sharex=True, sharey="col")

    distributions = [
        ("Raw WVD", wvd_raw, t_wvd, f_wvd),
        ("PWVD (h=51, no g)", pwvd, t_wvd, f_wvd),
        ("SPWVD (h=51, g=15)", spwvd, t_spwvd, f_spwvd),
    ]

    for i, (name, dist, t_d, f_d) in enumerate(distributions):
        # Linear (left column)
        im_lin = axes[i, 0].pcolormesh(t_d, f_d, dist,
                                        shading="gouraud", cmap=COLORMAP)
        axes[i, 0].set_ylabel(f"{name}\nFreq (Hz)")
        axes[i, 0].set_ylim(0, 100)
        fig.colorbar(im_lin, ax=axes[i, 0], label="WVD value")

        # dB (right column)
        dist_db = 10 * np.log10(np.maximum(np.abs(dist), 1e-20))
        im_db = axes[i, 1].pcolormesh(t_d, f_d, dist_db,
                                       shading="gouraud", cmap=COLORMAP)
        axes[i, 1].set_ylabel("Freq (Hz)")
        axes[i, 1].set_ylim(0, 100)
        fig.colorbar(im_db, ax=axes[i, 1], label="Power (dB)")

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    axes[0, 0].set_title("Linear scale")
    axes[0, 1].set_title("dB scale")

    fig.suptitle(f"Figure B.{FIG_START + 1} - WVD → PWVD → SPWVD progression",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig_td, fig


# ============================================================
# Experiment B: Two-Knob Sweep
# ============================================================
def experiment_smoothing_sweep():
    """
    Two extreme SPWVD tunings on the same chirp+tone signal.
    Shows independent control of time and frequency resolution.
    """
    print("\nExperiment B: Two-Knob Sweep")
    print("-" * 40)

    # --- Same signal ---
    x_chirp, _, t = make_chirp(F0_CHIRP, MU_CHIRP,
                                A=1.0, duration=DURATION_LAB)
    x_tone, _, _ = make_tone(F_TONE, A=1.0, duration=DURATION_LAB)
    x = x_chirp + x_tone

    # --- Case 1: strong frequency smoothing, weak time ---
    h1 = hann(101)                                         # long lag window (0.404 s)
    g1 = hann(5)                                           # short time window (0.020 s)
    spwvd1, t1, f1 = smoothed_pseudo_wigner_ville(
        x, FS, h=h1, g=g1, n_fft=N_FFT,
    )

    # --- Case 2: weak frequency smoothing, strong time ---
    h2 = hann(25)                                          # short lag window (0.100 s)
    g2 = hann(31)                                          # long time window (0.124 s)
    spwvd2, t2, f2 = smoothed_pseudo_wigner_ville(
        x, FS, h=h2, g=g2, n_fft=N_FFT,
    )

    print(f"  Case 1: h=101 ({101/FS:.3f} s), g=5 ({5/FS:.3f} s) - strong freq, weak time")
    print(f"  Case 2: h=25 ({25/FS:.3f} s), g=31 ({31/FS:.3f} s) - weak freq, strong time")

    # --- Dual-stack: 2 rows × 2 cols ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 8),
                              sharex=True, sharey="col")

    cases = [
        (f"Case 1: h=101 ({101/FS:.3f} s), g=5 ({5/FS:.3f} s)", spwvd1, t1, f1),
        (f"Case 2: h=25 ({25/FS:.3f} s), g=31 ({31/FS:.3f} s)", spwvd2, t2, f2),
    ]

    for i, (name, dist, t_d, f_d) in enumerate(cases):
        # Linear
        im = axes[i, 0].pcolormesh(t_d, f_d, dist,
                                    shading="gouraud", cmap=COLORMAP)
        axes[i, 0].set_ylabel(f"{name}\nFreq (Hz)")
        axes[i, 0].set_ylim(0, 100)
        fig.colorbar(im, ax=axes[i, 0], label="WVD value")

        # dB
        dist_db = 10 * np.log10(np.maximum(np.abs(dist), 1e-20))
        im_db = axes[i, 1].pcolormesh(t_d, f_d, dist_db,
                                       shading="gouraud", cmap=COLORMAP)
        axes[i, 1].set_ylabel("Freq (Hz)")
        axes[i, 1].set_ylim(0, 100)
        fig.colorbar(im_db, ax=axes[i, 1], label="Power (dB)")

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    axes[0, 0].set_title("Linear scale")
    axes[0, 1].set_title("dB scale")

    fig.suptitle(f"Figure B.{FIG_START + 2} - SPWVD two-knob sweep",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# Experiment C: Frequency-Oscillating Ghosts (Two Impulses)
# ============================================================
def experiment_freq_oscillating():
    """
    Two impulses separated in TIME produce FREQUENCY-oscillating ghosts.
    The PWVD (lag window) suppresses these — demonstrating the other half
    of the duality from Table A.12.
    """
    print("\nExperiment C: Frequency-Oscillating Ghosts (Two Impulses)")
    print("-" * 40)

    # --- Two impulses separated in time ---
    n0_1 = int(0.5 * FS)                                  # first impulse at t = 0.5 s
    n0_2 = int(1.5 * FS)                                  # second impulse at t = 1.5 s
    sigma_t = int(0.02 * FS)                               # very short (0.02 s = 5 samples)

    x1, _, t = make_transient(n0_1, sigma_t, f0=0.0,       # impulse 1 (baseband)
                               A=1.0, duration=DURATION_LAB)
    x2, _, _ = make_transient(n0_2, sigma_t, f0=0.0,       # impulse 2 (baseband)
                               A=1.0, duration=DURATION_LAB)
    x = x1 + x2                                           # two impulses

    dt_sep = (n0_2 - n0_1) / FS                           # time separation (s)
    print(f"  Impulse 1: t = {n0_1/FS:.1f} s")
    print(f"  Impulse 2: t = {n0_2/FS:.1f} s")
    print(f"  Separation: Δt = {dt_sep:.1f} s")
    print(f"  Cross-term type: FREQUENCY-oscillating (components separated in time)")
    print(f"  Cross-term midpoint: t_c = {(n0_1 + n0_2)/(2*FS):.1f} s")
    print(f"  Oscillation rate in frequency: Δt = {dt_sep:.1f} s")

    # --- Time-domain ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Two impulses at t = {n0_1/FS:.1f} s and t = {n0_2/FS:.1f} s",
        fig_id=f"Figure B.{FIG_START + 3}",
    )

    # --- WVD, PWVD, SPWVD ---
    wvd_raw, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    h_lag = hann(51)                                       # lag window
    g_none = np.array([1.0])                               # no time smoothing
    pwvd, _, _ = smoothed_pseudo_wigner_ville(
        x, FS, h=h_lag, g=g_none, n_fft=N_FFT,
    )

    g_time = hann(15)                                      # time window
    spwvd, t_spwvd, f_spwvd = smoothed_pseudo_wigner_ville(
        x, FS, h=h_lag, g=g_time, n_fft=N_FFT,
    )

    # --- Plot: 3 rows (WVD, PWVD, SPWVD) × 2 cols (linear, dB) ---
    fig, axes = plt.subplots(3, 2, figsize=(14, 12),
                              sharex=True, sharey="col")

    distributions = [
        ("Raw WVD", wvd_raw, t_wvd, f_wvd),
        ("PWVD (h=51)", pwvd, t_wvd, f_wvd),
        ("SPWVD (h=51, g=15)", spwvd, t_spwvd, f_spwvd),
    ]

    for i, (name, dist, t_d, f_d) in enumerate(distributions):
        im_lin = axes[i, 0].pcolormesh(t_d, f_d, dist,
                                        shading="gouraud", cmap=COLORMAP)
        axes[i, 0].set_ylabel(f"{name}\nFreq (Hz)")
        axes[i, 0].set_ylim(0, 100)
        fig.colorbar(im_lin, ax=axes[i, 0], label="WVD value")

        dist_db = 10 * np.log10(np.maximum(np.abs(dist), 1e-20))
        im_db = axes[i, 1].pcolormesh(t_d, f_d, dist_db,
                                       shading="gouraud", cmap=COLORMAP)
        axes[i, 1].set_ylabel("Freq (Hz)")
        axes[i, 1].set_ylim(0, 100)
        fig.colorbar(im_db, ax=axes[i, 1], label="Power (dB)")

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    axes[0, 0].set_title("Linear scale")
    axes[0, 1].set_title("dB scale")

    fig.suptitle(f"Figure B.{FIG_START + 4} - Two impulses: "
                 f"frequency-oscillating ghosts, PWVD suppresses them",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig_td, fig


# ============================================================
# Run all
# ============================================================
def run_lab8(save=False):
    """Run all Lab 8 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 8: The SPWVD and its Tradeoffs")
    print("=" * 60)

    fig_td_a, fig_prog = experiment_progression()           # Experiment A
    fig_sweep = experiment_smoothing_sweep()                # Experiment B
    fig_td_c, fig_impulse = experiment_freq_oscillating()   # Experiment C

    all_figs = [
        (fig_td_a, False),      # B.45: chirp+tone time domain
        (fig_prog, True),       # B.46: WVD→PWVD→SPWVD progression
        (fig_sweep, True),      # B.47: two-knob sweep
        (fig_td_c, False),      # B.48: two impulses time domain
        (fig_impulse, True),    # B.49: freq-oscillating ghost, PWVD works
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return [fig for fig, _ in all_figs]


if __name__ == "__main__":
    run_lab8(save=True)
