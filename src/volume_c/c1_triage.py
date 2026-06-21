"""
C.1 First Look and Triage
--------------------------
Load the EEG dataset, characterize it, and decide the analysis direction.
Data-first: time domain plot, then broad DFT, then band power.

Figures C.1-C.6 in the report.
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

# ============================================================
# Parameters
# ============================================================
EDF_FILE = os.path.join(DATA_DIR, "sub-NORB00055_ses-1_task-EEG_eeg.edf")
LAB_NUMBER = "volume_c/c1"
FIG_START = 1

# EEG channels (standard 10-20, exclude auxiliary 25+, 26+, 27+)
EEG_CHANNELS = [
    "Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4",
    "O1", "O2", "F7", "F8", "T3", "T4", "T5", "T6",
    "FZ", "CZ", "PZ",
]
AUX_CHANNELS = ["25+", "26+", "27+"]

# Representative channel for detailed analysis
PRIMARY_CHANNEL = "CZ"

# Welch parameters (justified from Lab 2: 5 s segments, Hann, 50% overlap)
WELCH_NPERSEG_S = 5.0
WELCH_OVERLAP = 0.5


# ============================================================
# C.1.1 - Dataset overview
# ============================================================
def c11_dataset_overview():
    """Load the dataset and print its characteristics."""
    print("=" * 60)
    print("C.1.1 - Dataset Overview")
    print("=" * 60)

    data, ch_names, fs, times = load_eeg(EDF_FILE)

    n_channels, n_samples = data.shape
    duration = times[-1] - times[0]

    print(f"File: {os.path.basename(EDF_FILE)}")
    print(f"Subject: NORB00055 (neonatal)")
    print(f"Sampling rate: {fs} Hz")
    print(f"Duration: {duration:.1f} s ({duration/60:.1f} min)")
    print(f"Channels: {n_channels}")
    print(f"EEG channels: {len(EEG_CHANNELS)}")
    print(f"Auxiliary channels: {AUX_CHANNELS}")
    print(f"Samples: {n_samples}")
    print()

    # Per-channel statistics
    print("Channel statistics (µV):")
    print(f"  {'Channel':12s} {'Mean':>8s} {'Std':>8s} {'Min':>10s} {'Max':>10s}")
    for ch in ch_names:
        x = get_channel_data(data, ch_names, ch)
        print(f"  {ch:12s} {x.mean():8.2f} {x.std():8.2f} {x.min():10.2f} {x.max():10.2f}")

    return data, ch_names, fs, times


# ============================================================
# C.1.2 - Time-domain plot (always first)
# ============================================================
def c12_time_domain(data, ch_names, fs, times):
    """Plot the raw EEG in time domain - what is visible by eye."""
    print()
    print("=" * 60)
    print("C.1.2 - Time-Domain Plot")
    print("=" * 60)

    # --- Full recording, primary channel ---
    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)

    fig_full, _ = plot_time_domain(
        times, x,
        xlabel="Time (s)", ylabel="Amplitude (µV)",
        title=f"Channel {PRIMARY_CHANNEL} - full recording",
        fig_id=f"Figure C.{FIG_START}",
    )

    # --- Zoomed: first 30 seconds ---
    fig_zoom1, _ = plot_time_domain(
        times, x,
        xlabel="Time (s)", ylabel="Amplitude (µV)",
        title=f"Channel {PRIMARY_CHANNEL} - first 30 s (zoomed)",
        fig_id=f"Figure C.{FIG_START + 1}",
        t_range=(0, 30),
    )

    # --- Multi-channel overview: 6 representative channels, 30 s ---
    rep_channels = ["Fp1", "F3", "C3", "P3", "O1", "T3"]
    fig_multi, axes = plt.subplots(len(rep_channels), 1,
                                    figsize=(14, 2.5 * len(rep_channels)),
                                    sharex=True)
    for i, ch in enumerate(rep_channels):
        x_ch = get_channel_data(data, ch_names, ch)
        mask = (times >= 0) & (times <= 30)
        axes[i].plot(times[mask], x_ch[mask], linewidth=0.3)
        axes[i].set_ylabel(f"{ch}\n(µV)")
        axes[i].grid(True, alpha=0.3)
        axes[i].set_ylim(-150, 150)
    axes[-1].set_xlabel("Time (s)")
    fig_multi.suptitle(
        f"Figure C.{FIG_START + 2} - Multi-channel EEG, first 30 s",
        fontsize=12, fontstyle="italic",
    )
    fig_multi.set_dpi(DPI)
    fig_multi.tight_layout()

    # --- Auxiliary channels (ECG/EMG/EOG?) ---
    fig_aux, axes_aux = plt.subplots(len(AUX_CHANNELS), 1,
                                      figsize=(14, 2.5 * len(AUX_CHANNELS)),
                                      sharex=True)
    for i, ch in enumerate(AUX_CHANNELS):
        x_ch = get_channel_data(data, ch_names, ch)
        mask = (times >= 0) & (times <= 30)
        axes_aux[i].plot(times[mask], x_ch[mask], linewidth=0.3)
        axes_aux[i].set_ylabel(f"{ch}\n(µV)")
        axes_aux[i].grid(True, alpha=0.3)
    axes_aux[-1].set_xlabel("Time (s)")
    fig_aux.suptitle(
        f"Figure C.{FIG_START + 3} - Auxiliary channels (25+, 26+, 27+), first 30 s",
        fontsize=12, fontstyle="italic",
    )
    fig_aux.set_dpi(DPI)
    fig_aux.tight_layout()

    print(f"Observations from {PRIMARY_CHANNEL}:")
    print(f"  Amplitude range: [{x.min():.1f}, {x.max():.1f}] µV")
    print(f"  Std: {x.std():.1f} µV")
    print(f"  Visual: slow oscillations dominate, possible discontinuities")

    return fig_full, fig_zoom1, fig_multi, fig_aux


# ============================================================
# C.1.3 - Broad DFT and band power
# ============================================================
def c13_spectral_triage(data, ch_names, fs, times):
    """First spectral look: Welch PSD + band power breakdown."""
    print()
    print("=" * 60)
    print("C.1.3 - Spectral Triage")
    print("=" * 60)

    x = get_channel_data(data, ch_names, PRIMARY_CHANNEL)
    nperseg = int(WELCH_NPERSEG_S * fs)
    noverlap = int(nperseg * WELCH_OVERLAP)

    print(f"Channel: {PRIMARY_CHANNEL}")
    print(f"Welch parameters:")
    print(f"  Segment: {nperseg} samples ({WELCH_NPERSEG_S} s)")
    print(f"  Overlap: {WELCH_OVERLAP*100:.0f}%")
    print(f"  Window: Hann")
    print(f"  Justification: 5 s segments give Df = 0.2 Hz (resolves all EEG bands),")
    print(f"    with L ~ {(len(x) - nperseg) // (nperseg - noverlap) + 1} segments for variance reduction (Lab 2)")
    print()

    # --- Welch PSD ---
    band_power, freqs, psd = compute_band_power(
        x, fs=fs, nperseg=nperseg, overlap_frac=WELCH_OVERLAP,
    )

    # --- Dual-stack PSD plot ---
    fig_psd, _ = plot_dual_stack_spectrum(
        freqs, psd,
        xlabel="Frequency (Hz)",
        ylabel_linear="PSD (µV²/Hz)",
        ylabel_db="PSD (dB re 1 µV²/Hz)",
        title=f"Channel {PRIMARY_CHANNEL} - Welch PSD",
        fig_id=f"Figure C.{FIG_START + 4}",
        f_range=(0, 50),
    )

    # --- Band power bar chart ---
    total_power = sum(band_power.values())
    print(f"Band power analysis:")
    print(f"  {'Band':8s} {'Range':>12s} {'Power (µV²)':>14s} {'%':>6s}")
    for name, power in band_power.items():
        fl, fh = EEG_BANDS[name]
        pct = 100 * power / total_power
        print(f"  {name:8s} {fl:4.1f}-{fh:2.0f} Hz {power:14.2f} {pct:5.1f}%")

    fig_bar, ax = plt.subplots(1, 1, figsize=(10, 5))
    band_names = list(band_power.keys())
    band_values = list(band_power.values())
    band_pcts = [100 * v / total_power for v in band_values]
    colors = ["#2563eb", "#059669", "#dc2626", "#f59e0b", "#7c3aed"]

    bars = ax.bar(band_names, band_pcts, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Relative power (%)")
    ax.set_xlabel("EEG band")
    ax.set_title(
        f"Figure C.{FIG_START + 5} - Band power distribution, channel {PRIMARY_CHANNEL}",
        fontsize=12, fontstyle="italic", loc="left",
    )
    for bar, pct in zip(bars, band_pcts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"{pct:.1f}%", ha="center", fontsize=10)
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_ylim(0, 105)
    fig_bar.set_dpi(DPI)
    fig_bar.tight_layout()

    # --- Peak frequency ---
    valid = (freqs >= 0.5) & (freqs <= 50)
    peak_freq = freqs[valid][np.argmax(psd[valid])]
    print(f"\n  Peak frequency: {peak_freq:.2f} Hz")
    print(f"  Dominant band: delta ({band_pcts[0]:.1f}% of total power)")

    return fig_psd, fig_bar, band_power


# ============================================================
# C.1.3b - Multi-channel band power comparison
# ============================================================
def c13b_multichannel_band_power(data, ch_names, fs):
    """Compare band power across all 19 EEG channels."""
    print()
    print("=" * 60)
    print("C.1.3b - Multi-Channel Band Power")
    print("=" * 60)

    nperseg = int(WELCH_NPERSEG_S * fs)                   # 5 s segments

    # --- Compute band power for each EEG channel ---
    results = {}
    for ch in EEG_CHANNELS:
        x = get_channel_data(data, ch_names, ch)          # single channel (µV)
        bp, _, _ = compute_band_power(x, fs=fs, nperseg=nperseg,
                                       overlap_frac=WELCH_OVERLAP)
        total = sum(bp.values())                          # total power
        results[ch] = {band: 100 * power / total for band, power in bp.items()}

    # --- Print table ---
    bands = list(EEG_BANDS.keys())
    print(f"  {'Channel':8s}", end="")
    for b in bands:
        print(f"  {b:>7s}", end="")
    print()
    for ch in EEG_CHANNELS:
        print(f"  {ch:8s}", end="")
        for b in bands:
            print(f"  {results[ch][b]:6.1f}%", end="")
        print()

    # --- Heatmap: channels x bands ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    matrix = np.array([[results[ch][b] for b in bands] for ch in EEG_CHANNELS])
    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd")

    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels(["δ", "θ", "α", "β", "γ"], fontsize=12)
    ax.set_yticks(range(len(EEG_CHANNELS)))
    ax.set_yticklabels(EEG_CHANNELS, fontsize=9)
    ax.set_xlabel("EEG Band")
    ax.set_ylabel("Channel")

    # Annotate each cell with the percentage
    for i in range(len(EEG_CHANNELS)):
        for j in range(len(bands)):
            val = matrix[i, j]
            color = "white" if val > 50 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                    fontsize=8, color=color)

    fig.colorbar(im, ax=ax, label="Relative power (%)", shrink=0.8)
    ax.set_title(
        f"Figure C.{FIG_START + 6} - Band power distribution across all EEG channels",
        fontsize=12, fontstyle="italic", loc="left",
    )
    fig.set_dpi(DPI)
    fig.tight_layout()

    # --- Summary statistics ---
    delta_values = [results[ch]["delta"] for ch in EEG_CHANNELS]
    print(f"\n  Delta power range: {min(delta_values):.1f}% - {max(delta_values):.1f}%")
    print(f"  Delta mean ± std: {np.mean(delta_values):.1f}% ± {np.std(delta_values):.1f}%")

    return fig, results


# ============================================================
# C.1.4 - Triage decision
# ============================================================
def c14_triage_decision(band_power):
    """Map the EEG to the six archetypes and decide analysis direction."""
    print()
    print("=" * 60)
    print("C.1.4 - Triage Decision")
    print("=" * 60)

    total = sum(band_power.values())
    delta_pct = 100 * band_power["delta"] / total
    theta_pct = 100 * band_power["theta"] / total

    print(f"Delta dominance: {delta_pct:.1f}% of total power")
    print(f"Theta secondary: {theta_pct:.1f}%")
    print()
    print("Archetype mapping:")
    print("  - Dominant slow rhythms (delta/theta) -> mixed-tone archetype")
    print("  - Possible discontinuous pattern (neonatal) -> transient archetype")
    print("  - Biological noise floor -> noise archetype")
    print("  - Auxiliary channels (25+, 26+, 27+) -> artifact sources")
    print()
    print("Analysis direction:")
    print("  C.2: Stationary characterization - DFT + band power across channels")
    print("  C.3: Time-varying - STFT spectrogram to check for discontinuity")
    print("  C.4: Artifacts - identify ECG/EMG in auxiliary channels")
    print("  C.5: WVD/SPWVD on a selected clean delta/theta segment")


# ============================================================
# Run all
# ============================================================
def run_c1(save=False):
    """Run the C.1 triage."""
    data, ch_names, fs, times = c11_dataset_overview()
    fig_full, fig_zoom, fig_multi, fig_aux = c12_time_domain(data, ch_names, fs, times)
    fig_psd, fig_bar, band_power = c13_spectral_triage(data, ch_names, fs, times)
    fig_heatmap, multichannel_results = c13b_multichannel_band_power(data, ch_names, fs)
    c14_triage_decision(band_power)

    all_figs = [
        (fig_full, False),    # C.1: full recording
        (fig_zoom, False),    # C.2: zoomed 30 s
        (fig_multi, False),   # C.3: multi-channel
        (fig_aux, False),     # C.4: auxiliary channels
        (fig_psd, False),     # C.5: Welch PSD
        (fig_bar, False),     # C.6: band power bar chart (CZ)
        (fig_heatmap, False), # C.7: multi-channel band power heatmap
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id,
                        raster_only=raster_only)

    plt.show()
    return data, ch_names, fs, times, band_power


if __name__ == "__main__":
    run_c1(save=True)
