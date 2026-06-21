"""
C.3 Time-Varying Characterization - STFT Spectrogram
------------------------------------------------------
Does delta persist continuously or come in bursts?
Heisenberg choice made for this specific signal.

Figures C.12-C.16 in the report.
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
)
from src.common.config import EEG_BANDS, DATA_DIR, FIGURE_FORMATS
from src.common.eeg import load_eeg, get_channel_data

# ============================================================
# Parameters
# ============================================================
EDF_FILE = os.path.join(DATA_DIR, "sub-NORB00055_ses-1_task-EEG_eeg.edf")
LAB_NUMBER = "volume_c/c3"
FIG_START = 12

PRIMARY_CHANNEL = "CZ"


# ============================================================
# C.3.1 - Full-recording spectrogram
# ============================================================
def c31_full_spectrogram(data, ch_names, fs):
    """STFT spectrogram of CZ over the full 19 minutes."""
    print("=" * 60)
    print("C.3.1 - Full-Recording Spectrogram")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    # --- STFT: M = 1000 (5.0 s), Hann, 50% overlap ---
    M = int(5.0 * fs)                                     # 1000 samples (5.0 s)
    noverlap = M // 2                                     # 50% overlap (Hann COLA)
    delta_f = fs / M                                      # 0.2 Hz
    delta_t = M / fs                                      # 5.0 s

    print(f"Channel: {PRIMARY_CHANNEL}")
    print(f"M = {M} ({M/fs:.1f} s), Δf = {delta_f:.2f} Hz, Δt = {delta_t:.1f} s")
    print(f"Δt·Δf = {delta_t * delta_f * 2:.1f} (Hann β=2)")

    f_stft, t_stft, Sxx = sp_signal.spectrogram(
        x, fs=fs, nperseg=M, noverlap=noverlap, window="hann",
    )

    # --- Dual-stack: linear + dB, 0-10 Hz ---
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(16, 8),
                                         sharex=True, sharey=True)

    im_lin = ax_lin.pcolormesh(t_stft, f_stft, Sxx,
                                shading="gouraud", cmap=COLORMAP)
    ax_lin.set_ylabel("Frequency (Hz)")
    ax_lin.set_ylim(0, 10)
    fig.colorbar(im_lin, ax=ax_lin, label="PSD (µV²/Hz)")
    ax_lin.set_title(f"Figure C.{FIG_START}", loc="left",
                      fontsize=9, fontstyle="italic")
    ax_lin.set_title(f"Channel {PRIMARY_CHANNEL} - STFT spectrogram, "
                     f"M={M} ({M/fs:.1f} s), full recording")

    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))
    im_db = ax_db.pcolormesh(t_stft, f_stft, Sxx_db,
                              shading="gouraud", cmap=COLORMAP)
    ax_db.set_ylabel("Frequency (Hz)")
    ax_db.set_xlabel("Time (s)")
    ax_db.set_ylim(0, 10)
    fig.colorbar(im_db, ax=ax_db, label="PSD (dB re 1 µV²/Hz)")

    fig.set_dpi(DPI)
    fig.tight_layout()

    # --- Observations ---
    delta_power = Sxx[(f_stft >= 0.5) & (f_stft <= 4), :].sum(axis=0)
    print(f"\nDelta power time course:")
    print(f"  Mean: {delta_power.mean():.1f}")
    print(f"  Std:  {delta_power.std():.1f}")
    print(f"  Min:  {delta_power.min():.1f}")
    print(f"  Max:  {delta_power.max():.1f}")
    print(f"  Max/Min ratio: {delta_power.max()/delta_power.min():.1f}x")

    return fig, t_stft, f_stft, Sxx


# ============================================================
# C.3.2 - Zoomed 60 s segment
# ============================================================
def c32_zoomed_segment(data, ch_names, fs):
    """Zoomed spectrogram of a 60 s segment showing burst structure."""
    print()
    print("=" * 60)
    print("C.3.2 - Zoomed 60 s Segment")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    M = int(2.0 * fs)                                     # 400 samples (2.0 s) - better time resolution
    noverlap = M // 2
    delta_f = fs / M
    delta_t = M / fs

    print(f"M = {M} ({M/fs:.1f} s), Δf = {delta_f:.2f} Hz, Δt = {delta_t:.1f} s")

    f_stft, t_stft, Sxx = sp_signal.spectrogram(
        x, fs=fs, nperseg=M, noverlap=noverlap, window="hann",
    )

    # --- Dual-stack, zoomed to 0-60 s, 0-10 Hz ---
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(16, 8),
                                         sharex=True, sharey=True)

    im_lin = ax_lin.pcolormesh(t_stft, f_stft, Sxx,
                                shading="gouraud", cmap=COLORMAP)
    ax_lin.set_ylabel("Frequency (Hz)")
    ax_lin.set_ylim(0, 10)
    ax_lin.set_xlim(0, 60)
    fig.colorbar(im_lin, ax=ax_lin, label="PSD (µV²/Hz)")
    ax_lin.set_title(f"Figure C.{FIG_START + 1}", loc="left",
                      fontsize=9, fontstyle="italic")
    ax_lin.set_title(f"Channel {PRIMARY_CHANNEL} - zoomed 0-60 s, "
                     f"M={M} ({M/fs:.1f} s)")

    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))
    im_db = ax_db.pcolormesh(t_stft, f_stft, Sxx_db,
                              shading="gouraud", cmap=COLORMAP)
    ax_db.set_ylabel("Frequency (Hz)")
    ax_db.set_xlabel("Time (s)")
    ax_db.set_ylim(0, 10)
    ax_db.set_xlim(0, 60)
    fig.colorbar(im_db, ax=ax_db, label="PSD (dB re 1 µV²/Hz)")

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# C.3.3 - Heisenberg comparison
# ============================================================
def c33_heisenberg_comparison(data, ch_names, fs):
    """Two window lengths on the same 60 s segment."""
    print()
    print("=" * 60)
    print("C.3.3 - Heisenberg Comparison")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    windows = [
        (int(1.0 * fs), "1.0 s"),                        # M=200, Δf=1.0 Hz, Δt=1.0 s
        (int(5.0 * fs), "5.0 s"),                        # M=1000, Δf=0.2 Hz, Δt=5.0 s
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    for col, (M, label) in enumerate(windows):
        noverlap = M // 2
        delta_f = 2 * fs / M                              # Hann β=2
        delta_t = M / fs

        print(f"M = {M} ({label}): Δf = {delta_f:.2f} Hz, Δt = {delta_t:.1f} s, "
              f"Δt·Δf = {delta_t * delta_f:.1f}")

        f_stft, t_stft, Sxx = sp_signal.spectrogram(
            x, fs=fs, nperseg=M, noverlap=noverlap, window="hann",
        )

        # Linear (top row)
        im = axes[0, col].pcolormesh(t_stft, f_stft, Sxx,
                                      shading="gouraud", cmap=COLORMAP)
        axes[0, col].set_ylim(0, 10)
        axes[0, col].set_xlim(0, 60)
        axes[0, col].set_ylabel("Frequency (Hz)")
        axes[0, col].set_title(f"M={M} ({label}), Δf={delta_f:.1f} Hz (linear)")
        fig.colorbar(im, ax=axes[0, col], label="µV²/Hz")

        # dB (bottom row)
        Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))
        im_db = axes[1, col].pcolormesh(t_stft, f_stft, Sxx_db,
                                         shading="gouraud", cmap=COLORMAP)
        axes[1, col].set_ylim(0, 10)
        axes[1, col].set_xlim(0, 60)
        axes[1, col].set_ylabel("Frequency (Hz)")
        axes[1, col].set_xlabel("Time (s)")
        axes[1, col].set_title(f"M={M} ({label}), Δt={delta_t:.1f} s (dB)")
        fig.colorbar(im_db, ax=axes[1, col], label="dB")

    fig.suptitle(f"Figure C.{FIG_START + 2} - Heisenberg comparison on {PRIMARY_CHANNEL}",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# C.3.4 - Delta power time course
# ============================================================
def c34_delta_time_course(data, ch_names, fs):
    """Extract delta band power over time from the spectrogram."""
    print()
    print("=" * 60)
    print("C.3.4 - Delta Power Time Course")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    M = int(2.0 * fs)                                     # 2 s for time resolution
    noverlap = M // 2

    f_stft, t_stft, Sxx = sp_signal.spectrogram(
        x, fs=fs, nperseg=M, noverlap=noverlap, window="hann",
    )

    # --- Integrate delta and theta power at each time step ---
    df = f_stft[1] - f_stft[0]
    delta_mask = (f_stft >= 0.5) & (f_stft <= 4)
    theta_mask = (f_stft >= 4) & (f_stft <= 8)

    delta_power = np.sum(Sxx[delta_mask, :], axis=0) * df  # µV² per time step
    theta_power = np.sum(Sxx[theta_mask, :], axis=0) * df

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 6), sharex=True)

    # Delta power over time
    ax1.plot(t_stft, delta_power, linewidth=0.5, color="#2563eb", label="Delta (0.5-4 Hz)")
    ax1.set_ylabel("Band power (µV²)")
    ax1.set_title(f"Figure C.{FIG_START + 3}", loc="left",
                   fontsize=9, fontstyle="italic")
    ax1.set_title(f"Delta and theta power over time, channel {PRIMARY_CHANNEL}")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Theta power over time
    ax2.plot(t_stft, theta_power, linewidth=0.5, color="#059669", label="Theta (4-8 Hz)")
    ax2.set_ylabel("Band power (µV²)")
    ax2.set_xlabel("Time (s)")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.set_dpi(DPI)
    fig.tight_layout()

    # --- Burst detection: when does delta power exceed 2x its median? ---
    median_delta = np.median(delta_power)
    burst_threshold = 2.0 * median_delta
    burst_mask = delta_power > burst_threshold
    burst_pct = 100 * np.sum(burst_mask) / len(burst_mask)
    quiet_pct = 100 - burst_pct

    print(f"Delta power statistics:")
    print(f"  Median: {median_delta:.1f} µV²")
    print(f"  Threshold (2x median): {burst_threshold:.1f} µV²")
    print(f"  Time above threshold (burst): {burst_pct:.1f}%")
    print(f"  Time below threshold (quiet): {quiet_pct:.1f}%")
    print(f"  Max/median ratio: {delta_power.max()/median_delta:.1f}x")

    # --- Correlation between delta and theta ---
    corr = np.corrcoef(delta_power, theta_power)[0, 1]
    print(f"\nDelta-theta correlation: {corr:.3f}")
    if corr > 0.7:
        print(f"  -> Strong positive: theta rises and falls with delta (broadband bursts)")
    elif corr < 0.3:
        print(f"  -> Weak: theta is independent of delta")

    return fig, delta_power, theta_power, t_stft


# ============================================================
# C.3.5 - Time-domain overlay with burst markers
# ============================================================
def c35_burst_overlay(data, ch_names, fs, delta_power, t_stft):
    """Overlay delta power on the time-domain signal."""
    print()
    print("=" * 60)
    print("C.3.5 - Burst Overlay")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    times = np.arange(len(x)) / fs

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 6), sharex=True)

    # Time domain (first 120 s)
    mask = times <= 120
    ax1.plot(times[mask], x[mask], linewidth=0.3, color="#2563eb")
    ax1.set_ylabel("Amplitude (µV)")
    ax1.set_title(f"Figure C.{FIG_START + 4}", loc="left",
                   fontsize=9, fontstyle="italic")
    ax1.set_title(f"Channel {PRIMARY_CHANNEL} - time domain with delta power overlay")
    ax1.grid(True, alpha=0.3)

    # Delta power (same time range)
    stft_mask = t_stft <= 120
    ax2.plot(t_stft[stft_mask], delta_power[stft_mask],
             linewidth=1.0, color="#dc2626")
    median_delta = np.median(delta_power)
    ax2.axhline(2 * median_delta, color="black", ls="--", lw=0.8,
                label=f"Burst threshold (2× median = {2*median_delta:.0f} µV²)")
    ax2.set_ylabel("Delta power (µV²)")
    ax2.set_xlabel("Time (s)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# C.3.6 - Window comparison on spectrogram
# ============================================================
def c36_window_comparison(data, ch_names, fs):
    """Hann vs Hamming vs Blackman spectrograms on the same 60 s segment."""
    print()
    print("=" * 60)
    print("C.3.6 - Window Comparison (Spectrogram)")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    M = int(2.0 * fs)                                     # 400 samples (2.0 s)
    noverlap = M // 2

    windows = [
        ("Hann", "hann", "#2563eb"),
        ("Hamming", "hamming", "#059669"),
        ("Blackman", "blackman", "#dc2626"),
    ]

    # --- Compute delta power time course for each window ---
    fig, ax = plt.subplots(1, 1, figsize=(16, 5))

    for name, win, color in windows:
        f_stft, t_stft, Sxx = sp_signal.spectrogram(
            x, fs=fs, nperseg=M, noverlap=noverlap, window=win,
        )
        df = f_stft[1] - f_stft[0]
        delta_mask = (f_stft >= 0.5) & (f_stft <= 4)
        delta_power = np.sum(Sxx[delta_mask, :], axis=0) * df

        t_mask = t_stft <= 60                              # first 60 s
        ax.plot(t_stft[t_mask], delta_power[t_mask],
                linewidth=0.8, color=color, label=name, alpha=0.8)

        # Compare to Hann with fair metrics
        if name != "Hann":
            f_h, t_h, Sxx_h = sp_signal.spectrogram(
                x, fs=fs, nperseg=M, noverlap=noverlap, window="hann",
            )
            delta_h = np.sum(Sxx_h[delta_mask, :], axis=0) * df
            abs_diff = np.abs(delta_power - delta_h)       # absolute difference per time step
            safe = np.maximum(delta_h, 1)                  # avoid division by zero
            rel_local = 100 * abs_diff / safe              # relative to local Hann value
            print(f"  {name} vs Hann:")
            print(f"    Max abs diff:       {np.max(abs_diff):.0f} µV²")
            print(f"    Median rel diff:    {np.median(rel_local):.1f}% (typical time step)")
            print(f"    Mean rel diff:      {np.mean(rel_local):.1f}%")
            print(f"    Max rel diff:       {np.max(rel_local):.1f}% (worst-case, at burst edges)")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Delta band power (µV²)")
    ax.set_title(f"Figure C.{FIG_START + 5}", loc="left",
                  fontsize=9, fontstyle="italic")
    ax.set_title(f"Delta power: Hann vs Hamming vs Blackman, M={M} ({M/fs:.1f} s)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.set_dpi(DPI)
    fig.tight_layout()

    # --- Three-window spectrogram comparison (dual-stack: linear + dB) ---
    betas = {"Hann": 2, "Hamming": 2, "Blackman": 3}
    fig_spec, axes_spec = plt.subplots(2, 3, figsize=(18, 9),
                                        sharex=True, sharey=True,
                                        constrained_layout=True)

    # Compute all spectrograms first for shared color ranges
    spec_data = []
    vmin_lin, vmax_lin = np.inf, -np.inf
    vmin_db, vmax_db = np.inf, -np.inf

    for name, win, color in windows:
        f_s, t_s, Sxx_s = sp_signal.spectrogram(
            x, fs=fs, nperseg=M, noverlap=noverlap, window=win,
        )
        Sxx_db = 10 * np.log10(np.maximum(Sxx_s, 1e-20))
        spec_data.append((name, t_s, f_s, Sxx_s, Sxx_db))
        vmax_lin = max(vmax_lin, np.percentile(Sxx_s, 99))
        vmin_db = min(vmin_db, np.percentile(Sxx_db, 5))
        vmax_db = max(vmax_db, np.percentile(Sxx_db, 99))

    for i, (name, t_s, f_s, Sxx_s, Sxx_db) in enumerate(spec_data):
        t_mask_spec = t_s <= 60

        # Top row: linear scale
        im_lin = axes_spec[0, i].pcolormesh(
            t_s[t_mask_spec], f_s, Sxx_s[:, t_mask_spec],
            shading="gouraud", cmap=COLORMAP,
            vmin=0, vmax=vmax_lin,
        )
        axes_spec[0, i].set_ylim(0, 10)
        axes_spec[0, i].set_title(f"{name} (β={betas[name]})")
        if i == 0:
            axes_spec[0, i].set_ylabel("Frequency (Hz)\n(linear)")

        # Bottom row: dB scale
        im_db = axes_spec[1, i].pcolormesh(
            t_s[t_mask_spec], f_s, Sxx_db[:, t_mask_spec],
            shading="gouraud", cmap=COLORMAP,
            vmin=vmin_db, vmax=vmax_db,
        )
        axes_spec[1, i].set_ylim(0, 10)
        axes_spec[1, i].set_xlabel("Time (s)")
        if i == 0:
            axes_spec[1, i].set_ylabel("Frequency (Hz)\n(dB)")

    # Shared colorbars
    fig_spec.colorbar(im_lin, ax=axes_spec[0, :].tolist(),
                       label="PSD (µV²/Hz)", fraction=0.02, pad=0.04)
    fig_spec.colorbar(im_db, ax=axes_spec[1, :].tolist(),
                       label="PSD (dB re 1 µV²/Hz)", fraction=0.02, pad=0.04)
    fig_spec.suptitle(
        f"Figure C.{FIG_START + 6} - Spectrogram comparison: "
        f"Hann vs Hamming vs Blackman, M={M} ({M/fs:.1f} s), 0-60 s",
        fontsize=12, fontstyle="italic",
    )
    fig_spec.set_dpi(DPI)

    return fig, fig_spec


# ============================================================
# Run all
# ============================================================
def run_c3(save=False):
    """Run C.3 time-varying characterization."""
    print("=" * 60)
    print("C.3 - Time-Varying Characterization")
    print("=" * 60)

    data, ch_names, fs, times = load_eeg(EDF_FILE)

    fig_full, t_stft, f_stft, Sxx = c31_full_spectrogram(data, ch_names, fs)
    fig_zoom = c32_zoomed_segment(data, ch_names, fs)
    fig_heis = c33_heisenberg_comparison(data, ch_names, fs)
    fig_power, delta_power, theta_power, t_stft2 = c34_delta_time_course(data, ch_names, fs)
    fig_overlay = c35_burst_overlay(data, ch_names, fs, delta_power, t_stft2)
    fig_wincomp, fig_winspec = c36_window_comparison(data, ch_names, fs)

    all_figs = [
        (fig_full, True),      # C.12: full spectrogram
        (fig_zoom, True),      # C.13: zoomed 60 s
        (fig_heis, True),      # C.14: Heisenberg comparison
        (fig_power, False),    # C.15: delta/theta power time course
        (fig_overlay, False),  # C.16: time domain with burst overlay
        (fig_wincomp, False),  # C.17: window delta power comparison
        (fig_winspec, True),   # C.18: window spectrogram comparison
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return all_figs


if __name__ == "__main__":
    run_c3(save=True)
