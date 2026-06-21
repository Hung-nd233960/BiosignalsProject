"""
Lab 7: The Wigner-Ville Distribution (WVD) and its Tradeoffs (↔ A.7)
--------------------------------------------------------------------
Tests: WVD vs STFT resolution for a chirp, cross-terms in multi-component signals,
       and the necessity of the analytic signal (removing DC/Nyquist ghosts).

Figures B.41–B.44 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.signal import hilbert

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (
    FS,
    SEED,
    DPI,
    COLORMAP,
    make_tone,
    make_chirp,
    make_time_axis,
    plot_time_domain,
    save_figure,
    wigner_ville,
)

# ============================================================
# Parameters (Table B.10 in the report)
# ============================================================
DURATION_LAB = 2.0  # 2 seconds (short duration to prevent memory/CPU overhead)
N_FFT = 512  # FFT resolution for time-frequency analysis
LAB_NUMBER = 7
FIG_START = 41


# ============================================================
# 1. Experiment A: Single Chirp - WVD vs STFT
# ============================================================
def experiment_single_chirp():
    """
    Compare WVD and STFT on a single linear chirp.
    WVD gives a razor-sharp delta trajectory, whereas STFT is blurred.
    """
    print("\nExperiment A: Single Chirp - WVD vs STFT")
    print("-" * 40)

    # --- Generate Chirp (10 Hz to 80 Hz over 2s) ---
    f0 = 10.0
    f_end = 80.0
    mu = (f_end - f0) / DURATION_LAB  # chirp rate: 35 Hz/s
    x, _, t = make_chirp(f0, mu, A=1.0, duration=DURATION_LAB, fs=FS)

    # --- Time-domain plot ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Linear Chirp (f₀={f0} Hz to f_end={f_end} Hz)",
        fig_id=f"Figure B.{FIG_START}",
    )

    # --- STFT Spectrogram ---
    # We use a short window (64 samples = 0.256 s) to show STFT blur
    from scipy.signal import stft
    nperseg = 64
    noverlap = nperseg // 2
    f_stft, t_stft, Zxx = stft(x, FS, window='hann', nperseg=nperseg, noverlap=noverlap, nfft=N_FFT)
    Sxx = np.abs(Zxx)**2

    # --- WVD ---
    wvd, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- Plot comparison ---
    fig_comp, (ax_stft, ax_wvd) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)

    # STFT Plot (linear scale, zoomed 0-100 Hz)
    im_stft = ax_stft.pcolormesh(t_stft, f_stft, Sxx, shading="gouraud", cmap=COLORMAP)
    ax_stft.set_ylabel("Frequency (Hz)")
    ax_stft.set_title("STFT Spectrogram (Hanning window, blurred trajectory)")
    ax_stft.set_title(f"Figure B.{FIG_START + 1}", loc="left", fontsize=9, fontstyle="italic")
    fig_comp.colorbar(im_stft, ax=ax_stft, label="Power (linear)")

    # WVD Plot (linear scale, zoomed 0-100 Hz)
    im_wvd = ax_wvd.pcolormesh(t_wvd, f_wvd, wvd, shading="gouraud", cmap=COLORMAP)
    ax_wvd.set_ylabel("Frequency (Hz)")
    ax_wvd.set_xlabel("Time (s)")
    ax_wvd.set_title("Wigner-Ville Distribution (analytic signal, razor-sharp trajectory)")
    fig_comp.colorbar(im_wvd, ax=ax_wvd, label="Power (linear)")

    # Set limits
    ax_stft.set_ylim(0, 100)
    ax_stft.set_xlim(0, DURATION_LAB)

    fig_comp.set_dpi(DPI)
    fig_comp.tight_layout()

    return fig_td, fig_comp


# ============================================================
# 2. Experiment B: Two Components - Chirp + Tone (Cross-Terms)
# ============================================================
def experiment_cross_terms():
    """
    Apply WVD to a two-component signal (chirp + stationary tone).
    Demonstrates the generation of cross-term artifacts at the midpoint.
    """
    print("\nExperiment B: Two Components - Chirp + Tone")
    print("-" * 40)

    # Component 1: Linear Chirp (10 Hz to 90 Hz over 2s)
    f0 = 10.0
    f_end = 90.0
    mu = (f_end - f0) / DURATION_LAB  # 40 Hz/s
    x_chirp, _, t = make_chirp(f0, mu, A=1.0, duration=DURATION_LAB, fs=FS)

    # Component 2: Stationary Tone at 40 Hz
    f_tone = 40.0
    x_tone, _, _ = make_tone(f_tone, A=1.0, duration=DURATION_LAB, fs=FS)

    # Combined signal
    x = x_chirp + x_tone

    # --- Compute WVD and STFT ---
    from scipy.signal import stft
    nperseg = 64
    noverlap = nperseg // 2
    f_stft, t_stft, Zxx = stft(x, FS, window='hann', nperseg=nperseg, noverlap=noverlap, nfft=N_FFT)
    Sxx = np.abs(Zxx)**2

    wvd, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- Plot comparison ---
    fig, (ax_stft, ax_wvd) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)

    # STFT Plot (linear, clean superposition)
    im_stft = ax_stft.pcolormesh(t_stft, f_stft, Sxx, shading="gouraud", cmap=COLORMAP)
    ax_stft.set_ylabel("Frequency (Hz)")
    ax_stft.set_title("STFT Spectrogram: Chirp + Tone (blurred but free of cross-terms)")
    ax_stft.set_title(f"Figure B.{FIG_START + 2}", loc="left", fontsize=9, fontstyle="italic")
    fig.colorbar(im_stft, ax=ax_stft, label="Power (linear)")

    # WVD Plot (linear, sharp lines but with strong oscillatory cross-terms)
    im_wvd = ax_wvd.pcolormesh(t_wvd, f_wvd, wvd, shading="gouraud", cmap=COLORMAP)
    ax_wvd.set_ylabel("Frequency (Hz)")
    ax_wvd.set_xlabel("Time (s)")
    ax_wvd.set_title("WVD: Chirp + Tone (sharp trajectories, but corrupted by midpoint cross-terms)")
    fig.colorbar(im_wvd, ax=ax_wvd, label="Power (linear)")

    # Set limits
    ax_stft.set_ylim(0, 100)
    ax_stft.set_xlim(0, DURATION_LAB)

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# 3. Experiment C: Real vs. Analytic Signal (Hilbert Transform)
# ============================================================
def experiment_analytic_signal():
    """
    Compare WVD computed using a real-valued signal vs. an analytic signal.
    The real-valued signal creates a strong cross-term at DC and Nyquist frequencies.
    The analytic signal (using the Hilbert transform) removes these artifacts.
    """
    print("\nExperiment C: Real vs. Analytic Signal")
    print("-" * 40)

    # Stationary tone at 30 Hz
    f0 = 30.0
    x, _, t = make_tone(f0, A=1.0, duration=DURATION_LAB, fs=FS)

    # --- WVD using REAL signal ---
    # We bypass hilbert by computing WVD directly on x
    N = len(x)
    wvd_real = np.zeros((N_FFT, N))
    for n in range(N):
        L = min(n, N - 1 - n)
        r = np.zeros(N_FFT, dtype=float)  # real-valued correlation
        r[0] = x[n] * x[n]
        for m in range(1, L + 1):
            val = x[n + m] * x[n - m]
            r[m] = val
            r[N_FFT - m] = val
        wvd_real[:, n] = 2.0 * np.real(np.fft.fft(r))

    f_real = np.linspace(0, FS / 2, N_FFT // 2)
    wvd_real_half = wvd_real[:N_FFT // 2, :]

    # --- WVD using ANALYTIC signal (our standard function) ---
    wvd_analytic, t_wvd, f_wvd = wigner_ville(x, FS, n_fft=N_FFT)

    # --- Plot comparison ---
    fig, (ax_real, ax_analytic) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)

    # Real Signal WVD (contains DC ghost and self-interference)
    im_real = ax_real.pcolormesh(t, f_real, wvd_real_half, shading="gouraud", cmap=COLORMAP)
    ax_real.set_ylabel("Frequency (Hz)")
    ax_real.set_title("WVD of REAL signal (corrupted by DC self-ghost and spectral symmetry)")
    ax_real.set_title(f"Figure B.{FIG_START + 3}", loc="left", fontsize=9, fontstyle="italic")
    fig.colorbar(im_real, ax=ax_real, label="Power (linear)")

    # Analytic Signal WVD (clean)
    im_analytic = ax_analytic.pcolormesh(t_wvd, f_wvd, wvd_analytic, shading="gouraud", cmap=COLORMAP)
    ax_analytic.set_ylabel("Frequency (Hz)")
    ax_analytic.set_xlabel("Time (s)")
    ax_analytic.set_title("WVD of ANALYTIC signal (Hilbert transform, clean tone at 30 Hz)")
    fig.colorbar(im_analytic, ax=ax_analytic, label="Power (linear)")

    # Set limits
    ax_real.set_ylim(0, 100)
    ax_real.set_xlim(0, DURATION_LAB)

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# 4. Run Lab 7
# ============================================================
def run_lab7(save=False):
    """Run all Lab 7 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 7: The Wigner-Ville Distribution (WVD) and its Tradeoffs")
    print("=" * 60)

    fig_td_a, fig_comp_a = experiment_single_chirp()  # Experiment A
    fig_comp_b = experiment_cross_terms()  # Experiment B
    fig_comp_c = experiment_analytic_signal()  # Experiment C

    all_figs = [
        (fig_td_a, False),    # B.41: Single chirp time-domain
        (fig_comp_a, True),   # B.42: WVD vs STFT on single chirp
        (fig_comp_b, True),   # B.43: Cross-terms of chirp + tone
        (fig_comp_c, True),   # B.44: Real vs Analytic WVD
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
    run_lab7(save=True)
