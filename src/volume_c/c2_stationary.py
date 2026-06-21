"""
C.2 Stationary Characterization - DFT and Band Power
------------------------------------------------------
Windowed DFT of CZ, resolution limit applied to real bands,
log-log PSD to check for 1/f vs rhythmic delta.

Figures C.8-C.11 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sp_signal
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (
    DPI,
    COLORMAP,
    save_figure,
    plot_time_domain,
    plot_dual_stack_spectrum,
)
from src.common.config import EEG_BANDS, DATA_DIR, FIGURE_FORMATS
from src.common.eeg import load_eeg, compute_band_power, get_channel_data
from src.common.windows import hann, blackman

# ============================================================
# Parameters
# ============================================================
EDF_FILE = os.path.join(DATA_DIR, "sub-NORB00055_ses-1_task-EEG_eeg.edf")
LAB_NUMBER = "volume_c/c2"
FIG_START = 8

PRIMARY_CHANNEL = "CZ"
WELCH_NPERSEG_S = 5.0
WELCH_OVERLAP = 0.5


# ============================================================
# C.2.1 - Windowed DFT with resolution limit
# ============================================================
def c21_windowed_dft(data, ch_names, fs):
    """Full-length windowed DFT of CZ, with resolution limit marked."""
    print("=" * 60)
    print("C.2.1 - Windowed DFT")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    N = len(x)

    # --- Apply Hann window to full signal ---
    w = hann(N)                                           # Hann window (periodic)
    x_windowed = x * w                                    # windowed signal

    # --- DFT ---
    X = np.fft.fft(x_windowed)                            # full DFT
    freqs = np.fft.fftfreq(N, d=1/fs)                     # frequency axis
    pos = freqs >= 0                                      # positive half
    P = np.abs(X[pos])**2 / N                             # power spectrum (µV²)

    # --- Resolution limit ---
    beta = 2                                              # Hann β
    delta_f = fs / N                                      # bin spacing
    delta_f_min = beta * delta_f                          # resolution limit
    print(f"Channel: {PRIMARY_CHANNEL}")
    print(f"N = {N}, fs = {fs} Hz")
    print(f"Bin spacing: Δf = {delta_f:.6f} Hz")
    print(f"Resolution limit (Hann): Δf_min = {delta_f_min:.6f} Hz")
    print(f"EEG band boundaries vs resolution:")
    for name, (fl, fh) in EEG_BANDS.items():
        print(f"  {name}: {fl}-{fh} Hz, gap = {fh-fl:.1f} Hz >> Δf_min = {delta_f_min:.4f} Hz")

    # --- Dual-stack PSD (0-50 Hz) ---
    fig_dft, (ax_lin, ax_db) = plot_dual_stack_spectrum(
        freqs[pos], P,
        xlabel="Frequency (Hz)",
        ylabel_linear="Power (µV²)",
        ylabel_db="Power (dB)",
        title=f"Channel {PRIMARY_CHANNEL} - Windowed DFT (Hann, N={N})",
        fig_id=f"Figure C.{FIG_START}",
        f_range=(0, 50),
    )

    # --- Mark band boundaries ---
    band_colors = {"delta": "#2563eb", "theta": "#059669", "alpha": "#dc2626",
                   "beta": "#f59e0b", "gamma": "#7c3aed"}
    for name, (fl, fh) in EEG_BANDS.items():
        for ax in (ax_lin, ax_db):
            ax.axvline(fl, color=band_colors[name], ls="--", lw=0.5, alpha=0.5)
            ax.axvline(fh, color=band_colors[name], ls="--", lw=0.5, alpha=0.5)

    return fig_dft


# ============================================================
# C.2.2 - Log-log PSD: 1/f check
# ============================================================
def c22_loglog_psd(data, ch_names, fs):
    """Log-log PSD to check if delta is 1/f noise or rhythmic."""
    print()
    print("=" * 60)
    print("C.2.2 - Log-Log PSD (1/f Check)")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    nperseg = int(WELCH_NPERSEG_S * fs)

    # --- Welch PSD ---
    freqs, psd = sp_signal.welch(x, fs=fs, nperseg=nperseg,
                                  noverlap=int(nperseg * WELCH_OVERLAP),
                                  window="hann")

    # --- Log-log plot ---
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    valid = freqs > 0.3                                   # exclude near-DC
    ax.loglog(freqs[valid], psd[valid], linewidth=0.8, color="#2563eb",
              label="Welch PSD")

    # --- Fit 1/f line to 5-40 Hz range (clearly above any rhythmic peak) ---
    fit_mask = (freqs >= 5) & (freqs <= 40)
    log_f = np.log10(freqs[fit_mask])
    log_p = np.log10(psd[fit_mask])

    from scipy.stats import linregress
    slope, intercept, r_value, _, _ = linregress(log_f, log_p)
    fit_line = 10**(intercept) * freqs[valid]**slope

    ax.loglog(freqs[valid], fit_line, "--", color="#dc2626", linewidth=1.0,
              label=f"1/f fit (5-40 Hz): slope = {slope:.2f}, R² = {r_value**2:.3f}")

    # --- Mark EEG bands ---
    band_colors = {"δ": (0.5, 4, "#2563eb"), "θ": (4, 8, "#059669"),
                   "α": (8, 13, "#dc2626"), "β": (13, 30, "#f59e0b")}
    for name, (fl, fh, color) in band_colors.items():
        ax.axvspan(fl, fh, alpha=0.1, color=color)
        ax.text(np.sqrt(fl * fh), ax.get_ylim()[1] * 0.5, name,
                fontsize=12, ha="center", color=color, alpha=0.7)

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("PSD (µV²/Hz)")
    ax.set_title(
        f"Figure C.{FIG_START + 1} - Log-log PSD, channel {PRIMARY_CHANNEL}",
        fontsize=12, fontstyle="italic", loc="left",
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.15, which="both")
    ax.set_xlim(0.3, 50)
    fig.set_dpi(DPI)
    fig.tight_layout()

    # --- Check: does delta sit above the 1/f trend? ---
    delta_mask = (freqs >= 0.5) & (freqs <= 4)
    delta_psd_actual = np.mean(psd[delta_mask])
    delta_psd_predicted = np.mean(10**(intercept) * freqs[delta_mask]**slope)
    excess_ratio = delta_psd_actual / delta_psd_predicted

    print(f"1/f fit (5-40 Hz): slope = {slope:.2f}, R² = {r_value**2:.3f}")
    print(f"Delta band (0.5-4 Hz):")
    print(f"  Actual mean PSD: {delta_psd_actual:.2f} µV²/Hz")
    print(f"  1/f prediction:  {delta_psd_predicted:.2f} µV²/Hz")
    print(f"  Excess ratio:    {excess_ratio:.1f}x")
    if excess_ratio > 2:
        print(f"  -> Delta power is {excess_ratio:.1f}x above 1/f trend — likely rhythmic, not just noise")
    else:
        print(f"  -> Delta power follows the 1/f trend — may be colored noise, not rhythmic")

    return fig, slope, r_value**2, excess_ratio


# ============================================================
# C.2.3 - Zoomed delta sub-structure
# ============================================================
def c23_delta_zoom(data, ch_names, fs):
    """Zoom into 0-5 Hz to see if there is structure within delta."""
    print()
    print("=" * 60)
    print("C.2.3 - Delta Sub-Structure")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    # --- Use longer segments for finer resolution within delta ---
    nperseg_long = int(20.0 * fs)                         # 20 s segments → Δf = 0.05 Hz
    noverlap_long = nperseg_long // 2

    freqs, psd = sp_signal.welch(x, fs=fs, nperseg=nperseg_long,
                                  noverlap=noverlap_long, window="hann")

    delta_f = freqs[1] - freqs[0]
    print(f"Long segments: M = {nperseg_long} ({nperseg_long/fs:.1f} s), Δf = {delta_f:.4f} Hz")

    # --- Dual-stack zoomed to 0-8 Hz ---
    fig_zoom, _ = plot_dual_stack_spectrum(
        freqs, psd,
        xlabel="Frequency (Hz)",
        ylabel_linear="PSD (µV²/Hz)",
        ylabel_db="PSD (dB re 1 µV²/Hz)",
        title=f"Channel {PRIMARY_CHANNEL} - Delta/theta zoom (20 s segments, Δf = {delta_f:.3f} Hz)",
        fig_id=f"Figure C.{FIG_START + 2}",
        f_range=(0, 8),
    )

    # --- Find peaks in delta range ---
    delta_mask = (freqs >= 0.3) & (freqs <= 4)
    delta_freqs = freqs[delta_mask]
    delta_psd = psd[delta_mask]
    peak_idx = np.argmax(delta_psd)
    peak_freq = delta_freqs[peak_idx]
    peak_power = delta_psd[peak_idx]

    print(f"Peak in delta band: {peak_freq:.2f} Hz ({peak_power:.2f} µV²/Hz)")

    # --- Is there a distinct spectral peak or just monotonic decay? ---
    # Check: does the PSD have a local maximum, or does it just decrease from DC?
    from scipy.signal import argrelextrema
    local_maxima = argrelextrema(delta_psd, np.greater, order=3)[0]
    if len(local_maxima) > 0:
        print(f"Local maxima in delta band: {len(local_maxima)}")
        for idx in local_maxima[:5]:
            print(f"  {delta_freqs[idx]:.2f} Hz: {delta_psd[idx]:.2f} µV²/Hz")
    else:
        print("No local maxima — PSD is monotonically decreasing in delta band")

    return fig_zoom


# ============================================================
# C.2.4 - Window comparison on real EEG
# ============================================================
def c24_window_comparison(data, ch_names, fs):
    """Compare Hann vs Blackman on real EEG — does side-lobe suppression matter?"""
    print()
    print("=" * 60)
    print("C.2.4 - Window Comparison on Real EEG")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    nperseg = int(WELCH_NPERSEG_S * fs)

    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    for name, win, color in [("Hann", "hann", "#2563eb"), ("Blackman", "blackman", "#dc2626")]:
        freqs, psd = sp_signal.welch(x, fs=fs, nperseg=nperseg,
                                      noverlap=int(nperseg * WELCH_OVERLAP),
                                      window=win)
        ax_lin.plot(freqs, psd, linewidth=0.8, color=color, label=name)
        ax_db.plot(freqs, 10*np.log10(np.maximum(psd, 1e-20)),
                   linewidth=0.8, color=color, label=name)

    ax_lin.set_ylabel("PSD (µV²/Hz)")
    ax_lin.set_title(
        f"Figure C.{FIG_START + 3} - Hann vs Blackman on {PRIMARY_CHANNEL}",
        fontsize=12, fontstyle="italic", loc="left",
    )
    ax_lin.legend(fontsize=9)
    ax_lin.grid(True, alpha=0.3)

    ax_db.set_ylabel("PSD (dB re 1 µV²/Hz)")
    ax_db.set_xlabel("Frequency (Hz)")
    ax_db.legend(fontsize=9)
    ax_db.grid(True, alpha=0.3)
    ax_lin.set_xlim(0, 50)

    fig.set_dpi(DPI)
    fig.tight_layout()

    print(f"Comparison: Hann (β=2) vs Blackman (β=3) on real EEG")
    print(f"  If the two curves overlap in the theta/alpha/beta bands,")
    print(f"  delta leakage is not significant and Hann is sufficient.")

    return fig


# ============================================================
# Run all
# ============================================================
def run_c2(save=False):
    """Run C.2 stationary characterization."""
    print("=" * 60)
    print("C.2 - Stationary Characterization")
    print("=" * 60)

    data, ch_names, fs, times = load_eeg(EDF_FILE)

    fig_dft = c21_windowed_dft(data, ch_names, fs)
    fig_loglog, slope, r2, excess = c22_loglog_psd(data, ch_names, fs)
    fig_zoom = c23_delta_zoom(data, ch_names, fs)
    fig_win = c24_window_comparison(data, ch_names, fs)

    all_figs = [
        (fig_dft, False),      # C.8: windowed DFT
        (fig_loglog, False),   # C.9: log-log 1/f check
        (fig_zoom, False),     # C.10: delta sub-structure
        (fig_win, False),      # C.11: window comparison
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return all_figs


if __name__ == "__main__":
    run_c2(save=True)
