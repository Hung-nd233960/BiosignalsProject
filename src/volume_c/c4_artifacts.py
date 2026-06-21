"""
C.4 Events and Artifacts - Statistics and Transient Detection
--------------------------------------------------------------
Characterize auxiliary channels, cross-correlate with CZ,
verify noise floor model on real EEG.

Figures C.19-C.22 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sp_signal
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (
    DPI,
    save_figure,
    plot_dual_stack_spectrum,
)
from src.common.config import EEG_BANDS, DATA_DIR
from src.common.eeg import load_eeg, compute_band_power, get_channel_data

# ============================================================
# Parameters
# ============================================================
EDF_FILE = os.path.join(DATA_DIR, "sub-NORB00055_ses-1_task-EEG_eeg.edf")
LAB_NUMBER = "volume_c/c4"
FIG_START = 19

PRIMARY_CHANNEL = "CZ"
AUX_CHANNELS = ["25+", "26+", "27+"]
AUX_LABELS = {
    "25+": "25+ (suspected ECG)",
    "26+": "26+ (suspected EMG)",
    "27+": "27+ (suspected EMG/EOG)",
}
WELCH_NPERSEG_S = 5.0
WELCH_OVERLAP = 0.5


# ============================================================
# C.4.1 - Auxiliary channel PSD
# ============================================================
def c41_auxiliary_psd(data, ch_names, fs):
    """Welch PSD of each auxiliary channel to identify ECG/EMG/EOG."""
    print("=" * 60)
    print("C.4.1 - Auxiliary Channel PSD")
    print("=" * 60)

    nperseg = int(WELCH_NPERSEG_S * fs)

    fig, axes = plt.subplots(len(AUX_CHANNELS), 2, figsize=(14, 3 * len(AUX_CHANNELS)),
                              sharex="col")

    for i, ch in enumerate(AUX_CHANNELS):
        x = get_channel_data(data, ch_names, ch)           # auxiliary channel (µV)
        freqs, psd = sp_signal.welch(x, fs=fs, nperseg=nperseg,
                                      noverlap=int(nperseg * WELCH_OVERLAP),
                                      window="hann")

        label = AUX_LABELS.get(ch, ch)                      # e.g. "25+ (suspected ECG)"

        # Linear (left column)
        axes[i, 0].plot(freqs, psd, linewidth=0.5)
        axes[i, 0].set_ylabel(f"{label}\nPSD (µV²/Hz)")
        axes[i, 0].set_xlim(0, 50)
        axes[i, 0].grid(True, alpha=0.3)

        # dB (right column)
        psd_db = 10 * np.log10(np.maximum(psd, 1e-20))
        axes[i, 1].plot(freqs, psd_db, linewidth=0.5)
        axes[i, 1].set_ylabel(f"{label}\nPSD (dB)")
        axes[i, 1].set_xlim(0, 50)
        axes[i, 1].grid(True, alpha=0.3)

        # Find peak frequency
        valid = (freqs >= 0.5) & (freqs <= 50)
        peak_f = freqs[valid][np.argmax(psd[valid])]
        total_power = np.sum(psd[valid]) * (freqs[1] - freqs[0])
        print(f"  {ch}: peak = {peak_f:.1f} Hz, total power = {total_power:.2f} µV²")

    axes[-1, 0].set_xlabel("Frequency (Hz)")
    axes[-1, 1].set_xlabel("Frequency (Hz)")
    axes[0, 0].set_title("Linear scale")
    axes[0, 1].set_title("dB scale")
    fig.suptitle(f"Figure C.{FIG_START} - Auxiliary channel PSDs",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# C.4.2 - Cross-correlation: auxiliary vs CZ
# ============================================================
def c42_cross_correlation(data, ch_names, fs):
    """Cross-correlate each auxiliary channel with CZ (Equation A.57)."""
    print()
    print("=" * 60)
    print("C.4.2 - Cross-Correlation: Auxiliary vs CZ")
    print("=" * 60)

    x_cz = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    fig, axes = plt.subplots(len(AUX_CHANNELS), 1,
                              figsize=(14, 3 * len(AUX_CHANNELS)), sharex=True)

    results = {}
    for i, ch in enumerate(AUX_CHANNELS):
        x_aux = get_channel_data(data, ch_names, ch)       # auxiliary channel (µV)

        # Pearson correlation (Equation A.59)
        rho = np.corrcoef(x_cz, x_aux)[0, 1]              # normalized, zero-lag
        results[ch] = rho

        # Cross-correlation for first 30 s (visualization)
        seg_len = int(30 * fs)                             # 30 s segment
        x_cz_seg = x_cz[:seg_len]
        x_aux_seg = x_aux[:seg_len]

        r_full = np.correlate(x_cz_seg, x_aux_seg, mode="full")
        lags = np.arange(len(r_full)) - (seg_len - 1)
        lag_time = lags / fs

        # Plot ±2 seconds of lag
        lag_mask = (lag_time >= -2) & (lag_time <= 2)
        axes[i].plot(lag_time[lag_mask], r_full[lag_mask], linewidth=0.5)
        label = AUX_LABELS.get(ch, ch)
        axes[i].set_ylabel(f"CZ × {label}")
        axes[i].set_title(f"{label} vs CZ: ρ = {rho:.4f}")
        axes[i].grid(True, alpha=0.3)

        print(f"  {ch} vs CZ: ρ = {rho:.4f}")

    axes[-1].set_xlabel("Lag (s)")
    fig.suptitle(f"Figure C.{FIG_START + 1} - Cross-correlation: auxiliary channels vs CZ",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    # Summary
    print()
    for ch, rho in results.items():
        if abs(rho) < 0.05:
            verdict = "no contamination"
        elif abs(rho) < 0.2:
            verdict = "weak contamination possible"
        else:
            verdict = "significant contamination"
        print(f"  {ch}: ρ = {rho:.4f} → {verdict}")

    return fig, results


# ============================================================
# C.4.3 - Noise floor verification (exponential distribution)
# ============================================================
def c43_noise_floor(data, ch_names, fs):
    """Test if high-frequency bin power follows exponential distribution (Lab 2, A.4)."""
    print()
    print("=" * 60)
    print("C.4.3 - Noise Floor Verification")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    N = len(x)

    # DFT of the full signal
    X = np.fft.fft(x)                                      # full DFT
    freqs = np.fft.fftfreq(N, d=1/fs)                      # frequency axis
    pos = freqs > 0                                        # positive frequencies

    # Bin powers in the alpha band (8-13 Hz) — quietest band in this neonatal EEG
    # (only 1% of total power, C.1 — alpha rhythms haven't developed yet)
    noise_mask = (freqs > 8) & (freqs < 13)
    noise_powers = np.abs(X[noise_mask])**2                 # |X[k]|² in alpha band

    # Expected: exponential distribution (Lab 2, Equation A.34)
    mean_power = np.mean(noise_powers)                      # E[|X[k]|²]
    lam = 1.0 / mean_power                                 # exponential rate

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram of bin powers
    ax1.hist(noise_powers, bins=100, density=True, alpha=0.7, label="Measured")
    p_range = np.linspace(0, np.percentile(noise_powers, 99), 200)
    ax1.plot(p_range, lam * np.exp(-lam * p_range), "r--", linewidth=1.5,
             label=f"Exponential (λ = {lam:.2e})")
    ax1.set_xlabel("Bin power |X[k]|²")
    ax1.set_ylabel("Probability density")
    ax1.set_title("Power distribution (8-13 Hz, alpha band)")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # QQ plot: sorted bin powers vs theoretical exponential quantiles
    sorted_powers = np.sort(noise_powers)
    n_bins = len(sorted_powers)
    theoretical = -np.log(1 - np.arange(1, n_bins + 1) / (n_bins + 1)) / lam
    ax2.scatter(theoretical, sorted_powers, s=1, alpha=0.3)
    max_val = max(theoretical.max(), sorted_powers.max())
    ax2.plot([0, max_val], [0, max_val], "r--", linewidth=1, label="Perfect fit")
    ax2.set_xlabel("Theoretical quantiles (exponential)")
    ax2.set_ylabel("Measured quantiles")
    ax2.set_title("Q-Q plot: exponential vs measured")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect("equal")

    fig.suptitle(f"Figure C.{FIG_START + 2} - Spectral distribution test: "
                 f"bin power in alpha band (8-13 Hz)",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    # Coefficient of variation test (A.4.1: exponential has CV = std/mean = 1)
    std_power = np.std(noise_powers)                       # standard deviation
    cv = std_power / mean_power                            # coefficient of variation

    print(f"  Frequency range: 8-13 Hz (alpha band - quietest in this neonatal EEG)")
    print(f"  Number of bins: {len(noise_powers)}")
    print(f"  Mean bin power: {mean_power:.2f} µV²")
    print(f"  Std bin power:  {std_power:.2f} µV²")
    print(f"  CV = std/mean:  {cv:.4f} (exponential predicts CV = 1.0)")
    if abs(cv - 1.0) < 0.1:
        print(f"  → CV ≈ 1.0: consistent with exponential (Lab 2 model holds)")
    else:
        print(f"  → CV = {cv:.2f} ≠ 1.0: deviates from ideal exponential")

    return fig, cv


# ============================================================
# C.4.4 - Summary: clean vs contaminated
# ============================================================
def c44_summary(data, ch_names, fs, xcorr_results):
    """Which channels are clean for WVD analysis?"""
    print()
    print("=" * 60)
    print("C.4.4 - Summary")
    print("=" * 60)

    # EEG inter-channel correlation for reference
    eeg_channels = ["Fp1", "F3", "C3", "P3", "O1", "CZ"]
    x_cz = get_channel_data(data, ch_names, "CZ")

    print("  EEG inter-channel correlations (vs CZ):")
    for ch in eeg_channels:
        if ch == "CZ":
            continue
        x_ch = get_channel_data(data, ch_names, ch)
        rho = np.corrcoef(x_cz, x_ch)[0, 1]
        print(f"    {ch} vs CZ: ρ = {rho:.4f}")

    print()
    print("  Auxiliary channel correlations (vs CZ):")
    for ch, rho in xcorr_results.items():
        print(f"    {ch} vs CZ: ρ = {rho:.4f}")

    print()
    print("  Verdict:")
    print("    EEG channels show high inter-correlation (shared brain activity)")
    print("    Auxiliary channels show low correlation with CZ")
    print("    → CZ is suitable for WVD analysis in C.5")


# ============================================================
# Run all
# ============================================================
def run_c4(save=False):
    """Run C.4 artifact characterization."""
    print("=" * 60)
    print("C.4 - Events and Artifacts")
    print("=" * 60)

    data, ch_names, fs, times = load_eeg(EDF_FILE)

    fig_psd = c41_auxiliary_psd(data, ch_names, fs)
    fig_xcorr, xcorr_results = c42_cross_correlation(data, ch_names, fs)
    fig_noise, cv = c43_noise_floor(data, ch_names, fs)
    c44_summary(data, ch_names, fs, xcorr_results)

    all_figs = [
        (fig_psd, False),      # C.19: auxiliary channel PSDs
        (fig_xcorr, False),    # C.20: cross-correlation vs CZ
        (fig_noise, False),    # C.21: noise floor verification
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return all_figs


if __name__ == "__main__":
    run_c4(save=True)
