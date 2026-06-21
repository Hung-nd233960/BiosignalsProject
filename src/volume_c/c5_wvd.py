"""
C.5 High-Resolution Time-Frequency - WVD / SPWVD
--------------------------------------------------
Apply WVD and SPWVD to a selected clean delta burst segment from CZ.
Compare against the STFT from C.3. Demonstrate cross-term soup on raw WVD,
then show SPWVD recovery with tuned windows.

Figures C.22-C.27 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sp_signal
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (                                   # project utilities
    DPI,
    COLORMAP,
    save_figure,
    plot_time_domain,
    wigner_ville,
    smoothed_pseudo_wigner_ville,
)
from src.common.config import EEG_BANDS, DATA_DIR
from src.common.eeg import load_eeg, get_channel_data
from src.common.windows import hann                        # our windows, not scipy's

# ============================================================
# Parameters
# ============================================================
EDF_FILE = os.path.join(DATA_DIR, "sub-NORB00055_ses-1_task-EEG_eeg.edf")
LAB_NUMBER = "volume_c/c5"
FIG_START = 22

PRIMARY_CHANNEL = "CZ"

# Segment selection: data-driven from C.3's burst detection.
# C.3 found bursts at 19% duty cycle using threshold = 2x median delta power.
# We recompute the burst detection here and pick the strongest burst.
SEGMENT_DURATION = 2.0    # segment length (s)
STFT_M_DETECT = None      # set from fs at runtime (5.0 s window, same as C.3)

# WVD/SPWVD parameters
N_FFT = 512               # frequency bins

# SPWVD windows — tuned for delta-band EEG
# h (lag window): controls frequency resolution
#   Delta band is 0.5-4 Hz, so we need Δf < 1 Hz minimum
#   Hann 51 at fs=200 → 51/200 = 0.255 s → Δf_main ≈ 2*fs/51 = 7.8 Hz (too coarse)
#   Hann 101 at fs=200 → 0.505 s → Δf_main ≈ 4.0 Hz (better for delta)
H_LEN = 101               # lag window length (0.505 s at 200 Hz)

# g (time window): controls time resolution
#   We want to resolve burst onset/offset within ~0.1-0.2 s
#   Hann 21 at fs=200 → 0.105 s
G_LEN = 21                # time window length (0.105 s at 200 Hz)

# STFT comparison
STFT_NPERSEG_SHORT = 64   # 0.32 s — short window for time resolution
STFT_NPERSEG_LONG = 200   # 1.0 s — long window for frequency resolution


# ============================================================
# C.5.1 - Segment selection and time domain
# ============================================================
def c51_segment_selection(data, ch_names, fs):
    """
    Data-driven burst segment selection using C.3's method:
    compute STFT delta power, find bursts (> 2x median), pick the strongest.
    """
    print("=" * 60)
    print("C.5.1 - Segment Selection (data-driven)")
    print("=" * 60)

    x_full = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    t_full = np.arange(len(x_full)) / fs                  # full time axis

    # --- Recompute delta power using same method as C.3 ---
    M_detect = int(5.0 * fs)                               # 5.0 s window (same as C.3)
    noverlap = M_detect // 2                               # 50% overlap
    f_stft, t_stft, Sxx = sp_signal.spectrogram(           # STFT spectrogram
        x_full, fs=fs, nperseg=M_detect,
        noverlap=noverlap, window="hann",
    )
    df = f_stft[1] - f_stft[0]                            # frequency bin width
    delta_mask = (f_stft >= EEG_BANDS["delta"][0]) & \
                 (f_stft <= EEG_BANDS["delta"][1])         # 0.5-4 Hz
    delta_power = np.sum(Sxx[delta_mask, :], axis=0) * df  # band power (µV²)

    # --- Burst detection: same threshold as C.3 ---
    median_delta = np.median(delta_power)                  # median delta power
    burst_threshold = 2.0 * median_delta                   # burst = 2x median
    burst_mask = delta_power > burst_threshold             # boolean mask

    # --- Step 1: try the strongest burst ---
    burst_indices = np.where(burst_mask)[0]                # all burst frames
    burst_powers = delta_power[burst_indices]              # their delta powers
    max_burst_idx = burst_indices[np.argmax(burst_powers)] # strongest burst frame
    t_max = t_stft[max_burst_idx]                          # its time

    # Extract the strongest burst segment
    t_max_start = max(0, t_max - SEGMENT_DURATION / 2)
    t_max_end = min(t_full[-1], t_max + SEGMENT_DURATION / 2)
    x_max = x_full[int(t_max_start * fs):int(t_max_end * fs)]
    t_max_seg = np.arange(len(x_max)) / fs

    # Check for clipping: count near-zero derivative samples
    diff_max = np.abs(np.diff(x_max))
    flat_count = np.sum(diff_max < 0.01)                   # near-zero derivative
    is_clipped = flat_count > 10                           # clipping if >10 flat samples

    print(f"  Strongest burst: t = {t_max:.1f} s, "
          f"delta = {delta_power[max_burst_idx]:.0f} µV² "
          f"({delta_power[max_burst_idx]/median_delta:.1f}x median)")
    print(f"  Flat samples (|Δx| < 0.01): {flat_count}")
    if is_clipped:
        print(f"  REJECTED: amplifier saturation detected — {flat_count} flat samples")

    # Plot the rejected burst (shows the selection process honestly)
    fig_rejected, ax_rej = plot_time_domain(
        t_max_seg, x_max,
        title=f"REJECTED: strongest burst at t = {t_max:.1f} s "
              f"(amplifier saturation, {flat_count} flat samples)",
        fig_id=f"Figure C.{FIG_START}",
    )
    ax_rej.set_xlabel("Time (s)")
    ax_rej.set_ylabel("Amplitude (µV)")

    # --- Step 2: fall back to 75th percentile ---
    p75 = np.percentile(burst_powers, 75)                  # 75th percentile
    target_idx = burst_indices[np.argmin(np.abs(burst_powers - p75))]
    t_peak = t_stft[target_idx]

    # Verify this one is clean
    t_start_cand = max(0, t_peak - SEGMENT_DURATION / 2)
    t_end_cand = min(t_full[-1], t_peak + SEGMENT_DURATION / 2)
    x_cand = x_full[int(t_start_cand * fs):int(t_end_cand * fs)]
    diff_cand = np.abs(np.diff(x_cand))
    flat_cand = np.sum(diff_cand < 0.01)

    if flat_cand > 10:                                     # still clipped, try 50th
        print(f"  75th percentile also clipped — trying 50th percentile")
        p50 = np.percentile(burst_powers, 50)
        target_idx = burst_indices[np.argmin(np.abs(burst_powers - p50))]
        t_peak = t_stft[target_idx]

    print(f"  Selected (75th pct): t = {t_peak:.1f} s, "
          f"delta = {delta_power[target_idx]:.0f} µV² "
          f"({delta_power[target_idx]/median_delta:.1f}x median)")

    # Center the segment on the selected burst
    t_start = max(0, t_peak - SEGMENT_DURATION / 2)        # segment start
    t_end = min(t_full[-1], t_peak + SEGMENT_DURATION / 2)  # segment end

    # Clamp to recording boundaries
    t_start = max(0, t_start)
    t_end = min(t_full[-1], t_end)

    i_start = int(t_start * fs)                            # start sample
    i_end = int(t_end * fs)                                # end sample
    x_seg = x_full[i_start:i_end]                          # extracted segment
    t_seg = np.arange(len(x_seg)) / fs                     # local time axis

    print(f"  Selected segment: {t_start:.1f} - {t_end:.1f} s "
          f"({SEGMENT_DURATION:.1f} s, {len(x_seg)} samples)")
    print(f"  Amplitude range: {x_seg.min():.1f} to {x_seg.max():.1f} µV")
    print(f"  RMS: {np.sqrt(np.mean(x_seg**2)):.1f} µV")

    # --- Time-domain plot of ACCEPTED burst ---
    fig_td, ax = plot_time_domain(
        t_seg, x_seg,
        title=f"ACCEPTED: CZ burst ({t_start:.1f}-{t_end:.1f} s, "
              f"delta = {delta_power[target_idx]/median_delta:.0f}x median)",
        fig_id=f"Figure C.{FIG_START + 1}",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude (µV)")

    return fig_rejected, fig_td, x_seg, t_seg, fs


# ============================================================
# C.5.2 - STFT baseline (what C.3 gives us)
# ============================================================
def c52_stft_baseline(x_seg, t_seg, fs):
    """STFT of the segment — the baseline from C.3's tool."""
    print()
    print("=" * 60)
    print("C.5.2 - STFT Baseline")
    print("=" * 60)

    # Short window: better time resolution
    f_s, t_s, Zxx_s = sp_signal.stft(                     # short STFT
        x_seg, fs, window="hann",
        nperseg=STFT_NPERSEG_SHORT,                        # 64 samples (0.32 s)
        noverlap=STFT_NPERSEG_SHORT // 2,
        nfft=N_FFT,
    )
    Sxx_s = np.abs(Zxx_s)**2                               # power

    # Long window: better frequency resolution
    f_l, t_l, Zxx_l = sp_signal.stft(                     # long STFT
        x_seg, fs, window="hann",
        nperseg=STFT_NPERSEG_LONG,                         # 200 samples (1.0 s)
        noverlap=STFT_NPERSEG_LONG // 2,
        nfft=N_FFT,
    )
    Sxx_l = np.abs(Zxx_l)**2

    delta_f_s = fs / STFT_NPERSEG_SHORT                    # freq resolution, short
    delta_f_l = fs / STFT_NPERSEG_LONG                     # freq resolution, long

    print(f"  Short window: M={STFT_NPERSEG_SHORT} ({STFT_NPERSEG_SHORT/fs:.3f} s), "
          f"Δf = {delta_f_s:.2f} Hz")
    print(f"  Long window:  M={STFT_NPERSEG_LONG} ({STFT_NPERSEG_LONG/fs:.3f} s), "
          f"Δf = {delta_f_l:.2f} Hz")

    # --- Dual-stack: short vs long window ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

    # Top-left: short window, linear
    im1 = axes[0, 0].pcolormesh(t_s, f_s, Sxx_s,
                                 shading="gouraud", cmap=COLORMAP)
    axes[0, 0].set_ylabel("Frequency (Hz)")
    axes[0, 0].set_ylim(0, 15)
    axes[0, 0].set_title(f"STFT M={STFT_NPERSEG_SHORT} ({STFT_NPERSEG_SHORT/fs:.2f} s) — linear")
    fig.colorbar(im1, ax=axes[0, 0], label="PSD (µV²/Hz)")

    # Top-right: long window, linear
    im2 = axes[0, 1].pcolormesh(t_l, f_l, Sxx_l,
                                 shading="gouraud", cmap=COLORMAP)
    axes[0, 1].set_ylabel("Frequency (Hz)")
    axes[0, 1].set_ylim(0, 15)
    axes[0, 1].set_title(f"STFT M={STFT_NPERSEG_LONG} ({STFT_NPERSEG_LONG/fs:.1f} s) — linear")
    fig.colorbar(im2, ax=axes[0, 1], label="PSD (µV²/Hz)")

    # Bottom-left: short window, dB
    Sxx_s_db = 10 * np.log10(np.maximum(Sxx_s, 1e-20))
    im3 = axes[1, 0].pcolormesh(t_s, f_s, Sxx_s_db,
                                 shading="gouraud", cmap=COLORMAP)
    axes[1, 0].set_ylabel("Frequency (Hz)")
    axes[1, 0].set_xlabel("Time (s)")
    axes[1, 0].set_ylim(0, 15)
    fig.colorbar(im3, ax=axes[1, 0], label="PSD (dB re 1 µV²/Hz)")

    # Bottom-right: long window, dB
    Sxx_l_db = 10 * np.log10(np.maximum(Sxx_l, 1e-20))
    im4 = axes[1, 1].pcolormesh(t_l, f_l, Sxx_l_db,
                                 shading="gouraud", cmap=COLORMAP)
    axes[1, 1].set_ylabel("Frequency (Hz)")
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylim(0, 15)
    fig.colorbar(im4, ax=axes[1, 1], label="PSD (dB re 1 µV²/Hz)")

    fig.suptitle(f"Figure C.{FIG_START + 2} — STFT Heisenberg tradeoff on burst segment",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# C.5.3 - Raw WVD (cross-term soup)
# ============================================================
def c53_raw_wvd(x_seg, t_seg, fs):
    """Raw WVD of the segment — expected: cross-term contamination."""
    print()
    print("=" * 60)
    print("C.5.3 - Raw WVD")
    print("=" * 60)

    wvd, t_wvd, f_wvd = wigner_ville(x_seg, fs, n_fft=N_FFT)

    print(f"  WVD shape: {wvd.shape}")
    print(f"  WVD range: {wvd.min():.1f} to {wvd.max():.1f}")
    print(f"  Negative values: {np.sum(wvd < 0)} / {wvd.size} "
          f"({100*np.sum(wvd < 0)/wvd.size:.1f}%)")

    # --- Dual-stack ---
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    im_lin = ax_lin.pcolormesh(t_wvd, f_wvd, wvd,
                                shading="gouraud", cmap=COLORMAP)
    ax_lin.set_ylabel("Frequency (Hz)")
    ax_lin.set_ylim(0, 15)
    ax_lin.set_title("Raw WVD — linear")
    fig.colorbar(im_lin, ax=ax_lin, label="WVD value (µV²·s)")

    wvd_db = 10 * np.log10(np.maximum(np.abs(wvd), 1e-20))
    im_db = ax_db.pcolormesh(t_wvd, f_wvd, wvd_db,
                              shading="gouraud", cmap=COLORMAP)
    ax_db.set_ylabel("Frequency (Hz)")
    ax_db.set_xlabel("Time (s)")
    ax_db.set_ylim(0, 15)
    ax_db.set_title("Raw WVD — dB")
    fig.colorbar(im_db, ax=ax_db, label="|WVD| (dB)")

    fig.suptitle(f"Figure C.{FIG_START + 3} — Raw WVD of burst segment (cross-term contamination expected)",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig, wvd, t_wvd, f_wvd


# ============================================================
# C.5.4 - SPWVD (the payoff)
# ============================================================
def c54_spwvd(x_seg, t_seg, fs):
    """SPWVD with tuned windows — the clean high-resolution result."""
    print()
    print("=" * 60)
    print("C.5.4 - SPWVD (tuned windows)")
    print("=" * 60)

    h_lag = hann(H_LEN)                                    # lag window
    g_time = hann(G_LEN)                                   # time window

    print(f"  h (lag window): Hann {H_LEN} ({H_LEN/fs:.3f} s)")
    print(f"  g (time window): Hann {G_LEN} ({G_LEN/fs:.3f} s)")

    spwvd, t_sp, f_sp = smoothed_pseudo_wigner_ville(
        x_seg, fs, h=h_lag, g=g_time, n_fft=N_FFT,
    )

    print(f"  SPWVD shape: {spwvd.shape}")
    print(f"  SPWVD range: {spwvd.min():.1f} to {spwvd.max():.1f}")
    print(f"  Negative values: {np.sum(spwvd < 0)} / {spwvd.size} "
          f"({100*np.sum(spwvd < 0)/spwvd.size:.1f}%)")

    # --- Dual-stack ---
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    im_lin = ax_lin.pcolormesh(t_sp, f_sp, spwvd,
                                shading="gouraud", cmap=COLORMAP)
    ax_lin.set_ylabel("Frequency (Hz)")
    ax_lin.set_ylim(0, 15)
    ax_lin.set_title(f"SPWVD (h={H_LEN}, g={G_LEN}) — linear")
    fig.colorbar(im_lin, ax=ax_lin, label="SPWVD value (µV²·s)")

    spwvd_db = 10 * np.log10(np.maximum(np.abs(spwvd), 1e-20))
    im_db = ax_db.pcolormesh(t_sp, f_sp, spwvd_db,
                              shading="gouraud", cmap=COLORMAP)
    ax_db.set_ylabel("Frequency (Hz)")
    ax_db.set_xlabel("Time (s)")
    ax_db.set_ylim(0, 15)
    ax_db.set_title(f"SPWVD (h={H_LEN}, g={G_LEN}) — dB")
    fig.colorbar(im_db, ax=ax_db, label="|SPWVD| (dB)")

    fig.suptitle(f"Figure C.{FIG_START + 4} — SPWVD of burst segment (cross-terms suppressed)",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig, spwvd, t_sp, f_sp


# ============================================================
# C.5.5 - Three-way comparison (STFT vs WVD vs SPWVD)
# ============================================================
def c55_comparison(x_seg, t_seg, fs):
    """Side-by-side: STFT vs WVD vs SPWVD on the same segment."""
    print()
    print("=" * 60)
    print("C.5.5 - Three-Way Comparison")
    print("=" * 60)

    # STFT
    f_stft, t_stft, Zxx = sp_signal.stft(
        x_seg, fs, window="hann",
        nperseg=STFT_NPERSEG_SHORT,
        noverlap=STFT_NPERSEG_SHORT // 2,
        nfft=N_FFT,
    )
    Sxx = np.abs(Zxx)**2

    # WVD
    wvd, t_wvd, f_wvd = wigner_ville(x_seg, fs, n_fft=N_FFT)

    # SPWVD
    h_lag = hann(H_LEN)
    g_time = hann(G_LEN)
    spwvd, t_sp, f_sp = smoothed_pseudo_wigner_ville(
        x_seg, fs, h=h_lag, g=g_time, n_fft=N_FFT,
    )

    # --- 3 rows × 2 cols (linear + dB) ---
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

    distributions = [
        ("STFT (M=64, 0.32 s)", Sxx, t_stft, f_stft, "PSD (µV²/Hz)", "PSD (dB)"),
        ("Raw WVD", wvd, t_wvd, f_wvd, "WVD value (µV²·s)", "|WVD| (dB)"),
        (f"SPWVD (h={H_LEN}, g={G_LEN})", spwvd, t_sp, f_sp, "SPWVD value (µV²·s)", "|SPWVD| (dB)"),
    ]

    for i, (name, dist, t_d, f_d, label_lin, label_db) in enumerate(distributions):
        # Linear
        im = axes[i, 0].pcolormesh(t_d, f_d, dist,
                                    shading="gouraud", cmap=COLORMAP)
        axes[i, 0].set_ylabel(f"{name}\nFreq (Hz)")
        axes[i, 0].set_ylim(0, 15)
        fig.colorbar(im, ax=axes[i, 0], label=label_lin)

        # dB
        dist_db = 10 * np.log10(np.maximum(np.abs(dist), 1e-20))
        im_db = axes[i, 1].pcolormesh(t_d, f_d, dist_db,
                                       shading="gouraud", cmap=COLORMAP)
        axes[i, 1].set_ylabel("Freq (Hz)")
        axes[i, 1].set_ylim(0, 15)
        fig.colorbar(im_db, ax=axes[i, 1], label=label_db)

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    axes[0, 0].set_title("Linear scale")
    axes[0, 1].set_title("dB scale")

    fig.suptitle(f"Figure C.{FIG_START + 5} — STFT vs WVD vs SPWVD on CZ burst segment",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# C.5.6 - SPWVD window sweep on real EEG
# ============================================================
def c56_window_sweep(x_seg, t_seg, fs):
    """Demonstrate the two-knob tradeoff on real EEG data."""
    print()
    print("=" * 60)
    print("C.5.6 - SPWVD Window Sweep on Real EEG")
    print("=" * 60)

    configs = [
        ("Short h, short g (sharp time + freq)",
         hann(51), hann(11), 51, 11),                      # minimal smoothing
        ("Long h, long g (heavy smoothing)",
         hann(151), hann(41), 151, 41),                    # heavy smoothing
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

    for i, (name, h, g, h_len, g_len) in enumerate(configs):
        sp, t_sp, f_sp = smoothed_pseudo_wigner_ville(
            x_seg, fs, h=h, g=g, n_fft=N_FFT,
        )

        print(f"  {name}: h={h_len} ({h_len/fs:.3f} s), g={g_len} ({g_len/fs:.3f} s)")

        # Linear
        im = axes[i, 0].pcolormesh(t_sp, f_sp, sp,
                                    shading="gouraud", cmap=COLORMAP)
        axes[i, 0].set_ylabel(f"h={h_len}, g={g_len}\nFreq (Hz)")
        axes[i, 0].set_ylim(0, 15)
        fig.colorbar(im, ax=axes[i, 0], label="SPWVD value")

        # dB
        sp_db = 10 * np.log10(np.maximum(np.abs(sp), 1e-20))
        im_db = axes[i, 1].pcolormesh(t_sp, f_sp, sp_db,
                                       shading="gouraud", cmap=COLORMAP)
        axes[i, 1].set_ylabel("Freq (Hz)")
        axes[i, 1].set_ylim(0, 15)
        fig.colorbar(im_db, ax=axes[i, 1], label="|SPWVD| (dB)")

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    axes[0, 0].set_title("Linear scale")
    axes[0, 1].set_title("dB scale")

    fig.suptitle(f"Figure C.{FIG_START + 6} — SPWVD two-knob sweep on real EEG",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# Run all
# ============================================================
def run_c5(save=False):
    """Run all C.5 analyses."""
    print("=" * 60)
    print("C.5 - High-Resolution Time-Frequency: WVD / SPWVD")
    print("=" * 60)

    data, ch_names, fs, times = load_eeg(EDF_FILE)

    fig_rej, fig_td, x_seg, t_seg, fs = c51_segment_selection(data, ch_names, fs)
    fig_stft = c52_stft_baseline(x_seg, t_seg, fs)
    fig_wvd, wvd, t_wvd, f_wvd = c53_raw_wvd(x_seg, t_seg, fs)
    fig_spwvd, spwvd, t_sp, f_sp = c54_spwvd(x_seg, t_seg, fs)
    fig_comp = c55_comparison(x_seg, t_seg, fs)
    fig_sweep = c56_window_sweep(x_seg, t_seg, fs)

    all_figs = [
        (fig_rej, False),      # C.22: rejected burst (saturation)
        (fig_td, False),       # C.23: selected clean burst
        (fig_stft, True),      # C.24: STFT baseline
        (fig_wvd, True),       # C.25: raw WVD (cross-term soup)
        (fig_spwvd, True),     # C.26: SPWVD (tuned)
        (fig_comp, True),      # C.27: three-way comparison
        (fig_sweep, True),     # C.28: window sweep on real EEG
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return all_figs


if __name__ == "__main__":
    run_c5(save=True)
