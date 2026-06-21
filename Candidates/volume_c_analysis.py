"""
Volume C: High-Resolution Time-Frequency Analysis of Real EEG Data (↔ C.5)
-------------------------------------------------------------------------
Applies WVD and SPWVD to a selected, cleaned segment of real EEG data (from channel O1).
Compares the time-frequency resolution of STFT, raw WVD, and SPWVD.

Figures C.1 and C.2 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.signal import stft
from scipy.signal.windows import hann

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from src.common import (
    DPI,
    COLORMAP,
    save_figure,
    load_preprocessed_eeg,
    get_channel_data,
    wigner_ville,
    smoothed_pseudo_wigner_ville,
)

# ============================================================
# Parameters
# ============================================================
EDF_PATH = r"C:\Users\P1\Downloads\Xlth\EEG.edf"
TARGET_CHANNEL = "O1"
T_START = 100.0  # start time of selected segment (s)
T_DURATION = 2.0  # 2 seconds duration
N_FFT = 512

def analyze_real_eeg():
    print("\nVolume C.5: High-Resolution Time-Frequency Analysis on Real EEG")
    print("-" * 60)

    # 1. Load preprocessed EEG data (notch at 50Hz, bandpass 1-40Hz, average reference)
    if not os.path.exists(EDF_PATH):
        print(f"Error: EDF file not found at {EDF_PATH}")
        return
        
    print(f"Loading data from {EDF_PATH}...")
    data, ch_names, fs, times = load_preprocessed_eeg(
        EDF_PATH, 
        channels=[TARGET_CHANNEL], 
        tmin=T_START, 
        tmax=T_START + T_DURATION
    )
    
    # Get 1D signal for channel O1
    x = get_channel_data(data, ch_names, TARGET_CHANNEL)
    t = times - T_START  # relative time axis starting at 0s
    
    print(f"Segment length: {len(x)} samples, fs = {fs} Hz, duration = {T_DURATION} s")
    
    # 2. Compute STFT Spectrogram
    # Hanning window of length 128 samples = 0.64 s
    nperseg = 128
    noverlap = nperseg // 2
    f_stft, t_stft, Zxx = stft(x, fs, window='hann', nperseg=nperseg, noverlap=noverlap, nfft=N_FFT)
    Sxx = np.abs(Zxx)**2
    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-10))

    # 3. Compute WVD (Raw Wigner-Ville Distribution)
    wvd, t_wvd, f_wvd = wigner_ville(x, fs, n_fft=N_FFT)
    wvd_db = 10 * np.log10(np.maximum(wvd, 1e-10))

    # 4. Compute SPWVD (Smoothed Pseudo Wigner-Ville Distribution)
    # frequency smoothing: lag window h = Hann(51)
    # time smoothing: window g = Hann(15)
    h_lag = hann(51, sym=True)
    g_time = hann(15, sym=True)
    spwvd, t_spwvd, f_spwvd = smoothed_pseudo_wigner_ville(x, fs, h=h_lag, g=g_time, n_fft=N_FFT)
    spwvd_db = 10 * np.log10(np.maximum(spwvd, 1e-10))

    # --- Plot 1: Comparison of STFT vs. Raw WVD vs. SPWVD (Linear Scale) ---
    fig_lin, (ax_stft, ax_wvd, ax_spwvd) = plt.subplots(3, 1, figsize=(12, 10), sharex=True, sharey=True)
    fig_lin.suptitle(f"Real EEG Time-Frequency Analysis — Channel {TARGET_CHANNEL} (Linear Scale)", fontsize=14, fontweight="bold")

    # STFT Plot
    im_stft = ax_stft.pcolormesh(t_stft, f_stft, Sxx, shading="gouraud", cmap=COLORMAP)
    ax_stft.set_ylabel("Frequency (Hz)")
    ax_stft.set_title("1. STFT Spectrogram (Hanning window, blurred due to Heisenberg uncertainty)")
    ax_stft.set_title("Figure C.1", loc="left", fontsize=9, fontstyle="italic")
    fig_lin.colorbar(im_stft, ax=ax_stft, label="Power (linear)")

    # Raw WVD
    im_wvd = ax_wvd.pcolormesh(t_wvd, f_wvd, wvd, shading="gouraud", cmap=COLORMAP)
    ax_wvd.set_ylabel("Frequency (Hz)")
    ax_wvd.set_title("2. Raw Wigner-Ville Distribution (high resolution, but corrupted by cross-term soup)")
    fig_lin.colorbar(im_wvd, ax=ax_wvd, label="Power (linear)")

    # SPWVD
    im_spwvd = ax_spwvd.pcolormesh(t_spwvd, f_spwvd, spwvd, shading="gouraud", cmap=COLORMAP)
    ax_spwvd.set_ylabel("Frequency (Hz)")
    ax_spwvd.set_xlabel("Time (s)")
    ax_spwvd.set_title("3. Smoothed Pseudo WVD (h=51, g=15, cross-terms suppressed, sharp trajectories preserved)")
    fig_lin.colorbar(im_spwvd, ax=ax_spwvd, label="Power (linear)")

    # Zoom in to EEG frequencies (0 to 40 Hz)
    ax_stft.set_ylim(0, 40)
    ax_stft.set_xlim(0, T_DURATION)

    fig_lin.set_dpi(DPI)
    fig_lin.tight_layout()

    # --- Plot 2: Comparison (dB Scale) ---
    fig_db, (ax_stft_db, ax_wvd_db, ax_spwvd_db) = plt.subplots(3, 1, figsize=(12, 10), sharex=True, sharey=True)
    fig_db.suptitle(f"Real EEG Time-Frequency Analysis — Channel {TARGET_CHANNEL} (dB Scale)", fontsize=14, fontweight="bold")

    # STFT dB
    im_stft_db = ax_stft_db.pcolormesh(t_stft, f_stft, Sxx_db, shading="gouraud", cmap=COLORMAP)
    ax_stft_db.set_ylabel("Frequency (Hz)")
    ax_stft_db.set_title("1. STFT Spectrogram (dB Scale)")
    ax_stft_db.set_title("Figure C.2", loc="left", fontsize=9, fontstyle="italic")
    fig_db.colorbar(im_stft_db, ax=ax_stft_db, label="Power (dB)")

    # Raw WVD dB
    im_wvd_db = ax_wvd_db.pcolormesh(t_wvd, f_wvd, wvd_db, shading="gouraud", cmap=COLORMAP)
    ax_wvd_db.set_ylabel("Frequency (Hz)")
    ax_wvd_db.set_title("2. Raw Wigner-Ville Distribution (dB Scale)")
    fig_db.colorbar(im_wvd_db, ax=ax_wvd_db, label="Power (dB)")

    # SPWVD dB
    im_spwvd_db = ax_spwvd_db.pcolormesh(t_spwvd, f_spwvd, spwvd_db, shading="gouraud", cmap=COLORMAP)
    ax_spwvd_db.set_ylabel("Frequency (Hz)")
    ax_spwvd_db.set_xlabel("Time (s)")
    ax_spwvd_db.set_title("3. Smoothed Pseudo WVD (dB Scale)")
    fig_db.colorbar(im_spwvd_db, ax=ax_spwvd_db, label="Power (dB)")

    # Zoom in to EEG frequencies
    ax_stft_db.set_ylim(0, 40)
    ax_stft_db.set_xlim(0, T_DURATION)

    fig_db.set_dpi(DPI)
    fig_db.tight_layout()

    # --- Save figures to results/graphs/volume_c/ ---
    volume_c_dir = os.path.join(os.path.dirname(__file__), "..", "results", "graphs", "volume_c")
    os.makedirs(volume_c_dir, exist_ok=True)
    
    path_lin = os.path.join(volume_c_dir, "figure_C_01.png")
    path_db = os.path.join(volume_c_dir, "figure_C_02.png")
    
    fig_lin.savefig(path_lin, dpi=DPI, bbox_inches="tight")
    fig_db.savefig(path_db, dpi=DPI, bbox_inches="tight")
    
    print(f"Saved Volume C Figures:\n  {path_lin}\n  {path_db}")
    plt.show()

if __name__ == "__main__":
    analyze_real_eeg()
