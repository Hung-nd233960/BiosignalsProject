"""
Lab 7: The Wigner-Ville Distribution (WVD) and its Tradeoffs (↔ A.7)
--------------------------------------------------------------------
Tests: WVD vs STFT resolution for a chirp, cross-terms in multi-component
signals, and the necessity of the analytic signal.

Figures B.39-B.44 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.signal import stft as scipy_stft             # scipy STFT for comparison

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (                                  # project utilities
    FS,
    DPI,
    COLORMAP,
    make_tone,
    make_chirp,
    plot_time_domain,
    save_figure,
    wigner_ville,                                         # WVD from src/common/wvd.py
)

# ============================================================
# Parameters
# ============================================================
# Duration is 2.0 s, NOT the 1200 s lab default.
# Justification: WVD has O(N^2) complexity. At N=300,000 (1200 s)
# the computation would require ~90 billion operations and exceed
# available memory. 2.0 s (N=500) is standard for WVD demonstrations.
DURATION_LAB = 2.0        # signal duration (s)
LAB_NUMBER = 7            # lab number for figure saving
FIG_START = 39            # first figure number (B.39), continuing from Lab 6

# Experiment A: single chirp
F0_CHIRP_A = 10.0         # start frequency (Hz)
F_END_CHIRP_A = 80.0      # end frequency (Hz)
MU_CHIRP_A = (F_END_CHIRP_A - F0_CHIRP_A) / DURATION_LAB  # 35 Hz/s

# Experiment B: chirp + tone
F0_CHIRP_B = 10.0         # chirp start (Hz)
F_END_CHIRP_B = 90.0      # chirp end (Hz)
MU_CHIRP_B = (F_END_CHIRP_B - F0_CHIRP_B) / DURATION_LAB  # 40 Hz/s
F_TONE_B = 40.0           # stationary tone (Hz)

# Experiment C: analytic signal test
F_TONE_C = 30.0           # test tone (Hz)

# STFT parameters for comparison
STFT_NPERSEG = 64         # 64 samples (0.256 s at 250 Hz)
STFT_NOVERLAP = 32        # 50% overlap
N_FFT = 512               # frequency bins for both STFT and WVD


# ============================================================
# Helper: dual-stack time-frequency comparison (STFT vs WVD)
# ============================================================
def plot_stft_vs_wvd(t_stft, f_stft, Sxx, t_wvd, f_wvd, wvd,
                      title_stft="STFT", title_wvd="WVD",
                      fig_id=None, f_range=(0, 100)):
    """
    Plot STFT and WVD side by side, each with dual-stack (linear + dB).

    Parameters
    ----------
    t_stft, f_stft, Sxx : STFT time, freq, power
    t_wvd, f_wvd, wvd   : WVD time, freq, power
    title_stft, title_wvd: subplot titles
    fig_id               : figure label
    f_range              : frequency zoom

    Returns
    -------
    fig : matplotlib Figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10),       # 2 rows x 2 cols
                              sharex=True)

    # Top-left: STFT linear
    im1 = axes[0, 0].pcolormesh(t_stft, f_stft, Sxx,       # STFT linear power
                                 shading="gouraud", cmap=COLORMAP)
    axes[0, 0].set_ylabel("Frequency (Hz)")
    axes[0, 0].set_title(f"{title_stft} (linear)")
    axes[0, 0].set_ylim(f_range)
    fig.colorbar(im1, ax=axes[0, 0], label="Power (linear)")

    # Bottom-left: STFT dB
    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))         # dB conversion
    im2 = axes[1, 0].pcolormesh(t_stft, f_stft, Sxx_db,
                                 shading="gouraud", cmap=COLORMAP)
    axes[1, 0].set_ylabel("Frequency (Hz)")
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_title(f"{title_stft} (dB)")
    axes[1, 0].set_ylim(f_range)
    fig.colorbar(im2, ax=axes[1, 0], label="Power (dB)")

    # Top-right: WVD linear
    im3 = axes[0, 1].pcolormesh(t_wvd, f_wvd, wvd,         # WVD linear
                                 shading="gouraud", cmap=COLORMAP)
    axes[0, 1].set_ylabel("Frequency (Hz)")
    axes[0, 1].set_title(f"{title_wvd} (linear)")
    axes[0, 1].set_ylim(f_range)
    fig.colorbar(im3, ax=axes[0, 1], label="WVD value")

    # Bottom-right: WVD dB
    wvd_db = 10 * np.log10(np.maximum(np.abs(wvd), 1e-20))  # dB of |WVD| (abs for negative values)
    im4 = axes[1, 1].pcolormesh(t_wvd, f_wvd, wvd_db,
                                 shading="gouraud", cmap=COLORMAP)
    axes[1, 1].set_ylabel("Frequency (Hz)")
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_title(f"{title_wvd} (dB)")
    axes[1, 1].set_ylim(f_range)
    fig.colorbar(im4, ax=axes[1, 1], label="|WVD| (dB)")

    if fig_id:
        axes[0, 0].set_title(fig_id, loc="left", fontsize=9, fontstyle="italic")

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig


# ============================================================
# Experiment A: Single Chirp - WVD vs STFT
# ============================================================
def experiment_single_chirp():
    """
    Compare WVD and STFT on a single linear chirp.
    WVD should show a sharp diagonal; STFT should show a blurred one.
    """
    print("\nExperiment A: Single Chirp - WVD vs STFT")
    print("-" * 40)

    # --- Generate chirp ---
    x, _, t = make_chirp(F0_CHIRP_A, MU_CHIRP_A,           # 10-80 Hz chirp
                          A=1.0, duration=DURATION_LAB)

    print(f"  Chirp: {F0_CHIRP_A}-{F_END_CHIRP_A} Hz, µ = {MU_CHIRP_A:.1f} Hz/s")
    print(f"  N = {len(x)}, duration = {DURATION_LAB} s")

    # --- Time-domain plot (always first) ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Linear chirp {F0_CHIRP_A}-{F_END_CHIRP_A} Hz",
        fig_id=f"Figure B.{FIG_START}",
    )

    # --- STFT ---
    f_stft, t_stft, Zxx = scipy_stft(                      # scipy STFT
        x, FS,
        window="hann",
        nperseg=STFT_NPERSEG,                              # 64 samples (0.256 s)
        noverlap=STFT_NOVERLAP,                            # 50% overlap
        nfft=N_FFT,                                        # 512 frequency bins
    )
    Sxx = np.abs(Zxx)**2                                   # power spectrogram

    # --- WVD (uses analytic signal + interpolation internally) ---
    wvd, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)  # from src/common/wvd.py

    # --- Dual-stack comparison ---
    fig_comp = plot_stft_vs_wvd(
        t_stft, f_stft, Sxx,
        t_wvd, f_wvd, wvd,
        title_stft=f"STFT (M={STFT_NPERSEG}, {STFT_NPERSEG/FS:.3f} s)",
        title_wvd="WVD (analytic signal)",
        fig_id=f"Figure B.{FIG_START + 1}",
    )

    return fig_td, fig_comp


# ============================================================
# Experiment B: Chirp + Tone - Cross-Terms
# ============================================================
def experiment_cross_terms():
    """
    WVD of a two-component signal (chirp + stationary tone).
    Cross-terms appear at the midpoint between the two components.
    """
    print("\nExperiment B: Chirp + Tone - Cross-Terms")
    print("-" * 40)

    # --- Generate signals ---
    x_chirp, _, t = make_chirp(F0_CHIRP_B, MU_CHIRP_B,     # 10-90 Hz chirp
                                A=1.0, duration=DURATION_LAB)
    x_tone, _, _ = make_tone(F_TONE_B,                     # 40 Hz tone
                              A=1.0, duration=DURATION_LAB)
    x = x_chirp + x_tone                                   # combined signal

    print(f"  Chirp: {F0_CHIRP_B}-{F_END_CHIRP_B} Hz, µ = {MU_CHIRP_B:.1f} Hz/s")
    print(f"  Tone: {F_TONE_B} Hz")
    print(f"  Expected cross-term midpoint: {(F0_CHIRP_B + F_TONE_B)/2:.0f} + "
          f"{MU_CHIRP_B/2:.0f}t Hz")

    # --- Time-domain plot (always first) ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Chirp ({F0_CHIRP_B}-{F_END_CHIRP_B} Hz) + Tone ({F_TONE_B} Hz)",
        fig_id=f"Figure B.{FIG_START + 2}",
    )

    # --- STFT ---
    f_stft, t_stft, Zxx = scipy_stft(
        x, FS,
        window="hann",
        nperseg=STFT_NPERSEG,                              # 64 samples (0.256 s)
        noverlap=STFT_NOVERLAP,
        nfft=N_FFT,
    )
    Sxx = np.abs(Zxx)**2

    # --- WVD ---
    wvd, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- Dual-stack comparison ---
    fig_comp = plot_stft_vs_wvd(
        t_stft, f_stft, Sxx,
        t_wvd, f_wvd, wvd,
        title_stft="STFT (clean superposition)",
        title_wvd="WVD (sharp + cross-terms)",
        fig_id=f"Figure B.{FIG_START + 3}",
    )

    return fig_td, fig_comp


# ============================================================
# Experiment C: Real vs Analytic Signal (DC Self-Ghost)
# ============================================================
def experiment_analytic_signal():
    """
    Compare WVD of a real-valued signal vs its analytic form.
    Real signal creates a DC self-ghost; analytic signal removes it.
    """
    print("\nExperiment C: Real vs Analytic Signal")
    print("-" * 40)

    # --- Generate tone ---
    x, _, t = make_tone(F_TONE_C, A=1.0,                   # 30 Hz tone
                         duration=DURATION_LAB)

    print(f"  Tone: {F_TONE_C} Hz")
    print(f"  Expected DC ghost: oscillating at 2×{F_TONE_C} = {2*F_TONE_C} Hz")

    # --- Time-domain ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Tone at {F_TONE_C} Hz (real-valued)",
        fig_id=f"Figure B.{FIG_START + 4}",
    )

    # --- WVD of REAL signal (manual, bypass analytic conversion) ---
    N = len(x)                                             # signal length
    wvd_real = np.zeros((N_FFT, N))                        # allocate WVD matrix
    for n in range(N):                                     # for each time sample
        L = min(n, N - 1 - n, N_FFT // 2 - 1)             # max lag
        r = np.zeros(N_FFT, dtype=float)                   # autocorrelation vector
        r[0] = x[n] * x[n]                                # lag 0
        for m in range(1, L + 1):                          # positive and negative lags
            val = x[n + m] * x[n - m]                      # instantaneous autocorrelation
            r[m] = val                                     # positive lag
            r[N_FFT - m] = val                             # negative lag (symmetric)
        wvd_real[:, n] = 2.0 * np.real(np.fft.fft(r))     # DFT over lag

    f_real = np.linspace(0, FS / 2, N_FFT // 2)           # frequency axis
    wvd_real_half = wvd_real[:N_FFT // 2, :]               # positive frequencies only

    # --- WVD of ANALYTIC signal (standard function) ---
    wvd_analytic, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- Dual-stack: real vs analytic ---
    fig_comp, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)

    # Top-left: real WVD linear
    im1 = axes[0, 0].pcolormesh(t, f_real, wvd_real_half,
                                 shading="gouraud", cmap=COLORMAP)
    axes[0, 0].set_ylabel("Frequency (Hz)")
    axes[0, 0].set_title("Real signal WVD (linear)")
    axes[0, 0].set_ylim(0, 100)
    fig_comp.colorbar(im1, ax=axes[0, 0], label="WVD value")
    axes[0, 0].set_title(f"Figure B.{FIG_START + 5}", loc="left",
                          fontsize=9, fontstyle="italic")

    # Bottom-left: real WVD dB
    wvd_real_db = 10 * np.log10(np.maximum(np.abs(wvd_real_half), 1e-20))
    im2 = axes[1, 0].pcolormesh(t, f_real, wvd_real_db,
                                 shading="gouraud", cmap=COLORMAP)
    axes[1, 0].set_ylabel("Frequency (Hz)")
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_title("Real signal WVD (dB) - DC ghost visible")
    axes[1, 0].set_ylim(0, 100)
    fig_comp.colorbar(im2, ax=axes[1, 0], label="|WVD| (dB)")

    # Top-right: analytic WVD linear
    im3 = axes[0, 1].pcolormesh(t_wvd, f_wvd, wvd_analytic,
                                 shading="gouraud", cmap=COLORMAP)
    axes[0, 1].set_ylabel("Frequency (Hz)")
    axes[0, 1].set_title("Analytic signal WVD (linear)")
    axes[0, 1].set_ylim(0, 100)
    fig_comp.colorbar(im3, ax=axes[0, 1], label="WVD value")

    # Bottom-right: analytic WVD dB
    wvd_an_db = 10 * np.log10(np.maximum(np.abs(wvd_analytic), 1e-20))
    im4 = axes[1, 1].pcolormesh(t_wvd, f_wvd, wvd_an_db,
                                 shading="gouraud", cmap=COLORMAP)
    axes[1, 1].set_ylabel("Frequency (Hz)")
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_title("Analytic signal WVD (dB) - clean")
    axes[1, 1].set_ylim(0, 100)
    fig_comp.colorbar(im4, ax=axes[1, 1], label="|WVD| (dB)")

    fig_comp.set_dpi(DPI)
    fig_comp.tight_layout()

    return fig_td, fig_comp


# ============================================================
# Run all
# ============================================================
def run_lab7(save=False):
    """Run all Lab 7 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 7: The WVD and its Tradeoffs")
    print("=" * 60)

    fig_td_a, fig_comp_a = experiment_single_chirp()        # Experiment A
    fig_td_b, fig_comp_b = experiment_cross_terms()         # Experiment B
    fig_td_c, fig_comp_c = experiment_analytic_signal()     # Experiment C

    all_figs = [
        (fig_td_a, False),     # B.39: chirp time domain
        (fig_comp_a, True),    # B.40: STFT vs WVD (single chirp)
        (fig_td_b, False),     # B.41: chirp+tone time domain
        (fig_comp_b, True),    # B.42: STFT vs WVD (cross-terms)
        (fig_td_c, False),     # B.43: tone time domain
        (fig_comp_c, True),    # B.44: real vs analytic WVD
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return [fig for fig, _ in all_figs]


if __name__ == "__main__":
    run_lab7(save=True)
