"""
Appendix B2 - CV of the Six Signal Archetypes
-----------------------------------------------
Compute the coefficient of variation (CV = std/mean) of the DFT bin power
for each of the six signal archetypes from Appendix A.

Figures AB2.1-AB2.2 in the report.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (
    FS,
    DPI,
    save_figure,
    make_tone,
    make_mixed_tones,
    make_chirp,
    make_transient,
    make_noise,
    plot_time_domain,
    plot_dual_stack_spectrum,
)
from src.common.config import FIGURE_FORMATS

# ============================================================
# Parameters
# ============================================================
DURATION = 60.0           # 60 s signals for consistency
LAB_NUMBER = "appendix_b2"
FIG_START = 1


def compute_cv(x, fs=FS):
    """
    Compute the coefficient of variation of |X[k]|² (bin power).
    CV = std(|X[k]|²) / mean(|X[k]|²)
    For exponential distribution (noise): CV = 1.0 exactly.

    Parameters
    ----------
    x  : ndarray - input signal
    fs : float - sampling frequency

    Returns
    -------
    cv    : float - coefficient of variation
    freqs : ndarray - frequency axis (positive, Hz)
    power : ndarray - bin powers |X[k]|² (positive freqs only)
    """
    X = np.fft.fft(x)                                     # full DFT
    freqs = np.fft.fftfreq(len(x), d=1/fs)                # frequency axis
    pos = freqs > 0                                       # positive frequencies (exclude DC)
    power = np.abs(X[pos])**2                              # bin powers
    cv = np.std(power) / np.mean(power)                    # CV = std / mean
    return cv, freqs[pos], power


def run_appendix_b2(save=False):
    """Compute CV for all six archetypes and generate figures."""
    print("=" * 60)
    print("Appendix B2 - CV of the Six Signal Archetypes")
    print("=" * 60)

    # --- Define the six archetypes ---
    archetypes = []

    # 1. Single tone (10 Hz)
    x_tone, _, t = make_tone(10.0, A=1.0, duration=DURATION)
    archetypes.append(("Single tone\n(10 Hz)", x_tone, t))

    # 2. Mixed tones (10 + 20 Hz)
    x_mixed, _, t = make_mixed_tones([10.0, 20.0], amplitudes=[1.0, 1.0], duration=DURATION)
    archetypes.append(("Mixed tones\n(10 + 20 Hz)", x_mixed, t))

    # 3. Chirp (5 -> 25 Hz)
    x_chirp, _, t = make_chirp(5.0, 0.333, A=1.0, duration=DURATION)
    archetypes.append(("Chirp\n(5-25 Hz)", x_chirp, t))

    # 4. Transient (alpha burst at t=30 s)
    n0 = int(30.0 * FS)                                   # center at 30 s
    sigma_t = int(0.5 * FS)                                # 0.5 s width
    x_trans, _, t = make_transient(n0, sigma_t, f0=10.0, A=3.0, duration=DURATION)
    archetypes.append(("Transient\n(0.5 s burst)", x_trans, t))

    # 5. Noise
    x_noise, _, t = make_noise(sigma=1.0, duration=DURATION, seed=42)
    archetypes.append(("Noise\n(σ = 1.0)", x_noise, t))

    # 6. Tone + noise (10 Hz, A=0.5 + noise σ=1.0)
    x_tn, _, t = make_tone(10.0, A=0.5, duration=DURATION)
    x_n, _, _ = make_noise(sigma=1.0, duration=DURATION, seed=42)
    x_tone_noise = x_tn + x_n
    archetypes.append(("Tone + noise\n(A=0.5, σ=1.0)", x_tone_noise, t))

    # --- Compute CV for each ---
    results = []
    print()
    print(f"  {'Archetype':25s}  {'CV':>8s}  {'Interpretation'}")
    print(f"  {'-'*25}  {'-'*8}  {'-'*30}")

    for name, x, t in archetypes:
        cv, freqs, power = compute_cv(x)
        clean_name = name.replace("\n", " ")
        if cv < 1.5:
            interp = "noise-like (energy spread)"
        elif cv < 10:
            interp = "moderate concentration"
        else:
            interp = "highly concentrated (few bins)"
        results.append((name, clean_name, cv, interp))
        print(f"  {clean_name:25s}  {cv:8.2f}  {interp}")

    # --- Figure AB2.1: Time domain + PSD for all six ---
    fig, axes = plt.subplots(6, 2, figsize=(14, 18), constrained_layout=True)

    for i, (name, x, t) in enumerate(archetypes):
        cv = results[i][2]

        # Time domain (left, first 2 seconds)
        mask = t <= 2
        axes[i, 0].plot(t[mask], x[mask], linewidth=0.5)
        axes[i, 0].set_ylabel(name, fontsize=9)
        axes[i, 0].grid(True, alpha=0.3)
        if i == 0:
            axes[i, 0].set_title("Time domain (first 2 s)")
        if i == len(archetypes) - 1:
            axes[i, 0].set_xlabel("Time (s)")

        # PSD (right, 0-50 Hz)
        _, freqs, power = compute_cv(x)
        psd_db = 10 * np.log10(np.maximum(power, 1e-20))
        axes[i, 1].plot(freqs, psd_db, linewidth=0.3)
        axes[i, 1].set_xlim(0, 50)
        axes[i, 1].set_ylabel(f"CV = {cv:.1f}")
        axes[i, 1].grid(True, alpha=0.3)
        if i == 0:
            axes[i, 1].set_title("Power spectrum (dB)")
        if i == len(archetypes) - 1:
            axes[i, 1].set_xlabel("Frequency (Hz)")

    fig.suptitle("Figure AB2.1 - Six signal archetypes: time domain and power spectrum",
                 fontsize=13, fontstyle="italic")
    fig.set_dpi(DPI)

    # --- Figure AB2.2: CV bar chart ---
    fig_bar, ax = plt.subplots(1, 1, figsize=(10, 5))
    names = [r[0] for r in results]
    cvs = [r[2] for r in results]
    colors = ["#2563eb", "#059669", "#f59e0b", "#dc2626", "#94a3b8", "#7c3aed"]

    bars = ax.bar(range(len(names)), cvs, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=8, ha="center")
    ax.set_ylabel("Coefficient of Variation (CV = std / mean)")
    ax.set_title("Figure AB2.2 - CV of bin power for the six signal archetypes",
                 fontsize=12, fontstyle="italic", loc="left")
    ax.axhline(1.0, color="black", ls="--", lw=1, label="CV = 1.0 (exponential / noise)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")

    for bar, cv in zip(bars, cvs):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + max(cvs)*0.02,
                f"{cv:.1f}", ha="center", fontsize=9)

    fig_bar.set_dpi(DPI)
    fig_bar.tight_layout()

    # --- Save ---
    all_figs = [
        (fig, False),       # AB2.1: time domain + PSD
        (fig_bar, False),   # AB2.2: CV bar chart
    ]

    if save:
        for i, (f, raster_only) in enumerate(all_figs):
            save_figure(f, lab_number=LAB_NUMBER, fig_id=f"0{i+1}",
                        raster_only=raster_only)

    plt.show()

    # --- Summary table ---
    print()
    print("Table AB2.1 - CV of the six signal archetypes")
    print(f"  {'Archetype':25s}  {'CV':>8s}")
    for _, clean_name, cv, _ in results:
        print(f"  {clean_name:25s}  {cv:8.2f}")

    return results


if __name__ == "__main__":
    run_appendix_b2(save=True)
