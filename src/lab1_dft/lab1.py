"""
Lab 1: The DFT of Basic Signals (↔ A.1, A.2)
----------------------------------------------
Tests: bin spacing, on-bin vs off-bin placement, resolution vs bin count,
       zero-padding, negative-frequency twin.
"""

import numpy as np  # numerical computing
import matplotlib.pyplot as plt  # plotting
import sys
import os  # path manipulation

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (  # project utilities
    FS,
    DURATION,
    make_tone,
    make_mixed_tones,
    plot_time_domain,
    plot_dual_stack_spectrum,
    save_figure,
)

# ============================================================
# Parameters (Table B.1 in the report)
# ============================================================
F0_ON_BIN = 10.0  # tone frequency that lands exactly on a bin (Hz)
F0_OFF_BIN = 10.5  # maximally off-bin: halfway between bin 10 and bin 11
DURATION_ONOFF = 1.0  # short segment so Δf = 1 Hz and leakage is visible
F1_CHORD = 10.0  # first tone of the dual-tone chord (Hz)
F2_CHORD = 12.0  # second tone of the dual-tone chord (Hz)
A = 1.0  # amplitude for all tones
DURATION_LAB = DURATION  # signal duration (s) — uses project default (1200 s)
ZERO_PAD_FACTOR = 4  # zero-pad to 4× original length


# ============================================================
# 1. Helper: compute DFT and frequency axis
# ============================================================
def compute_dft(x, fs=FS, n_fft=None):
    """
    Compute the DFT of x and return the frequency axis and complex spectrum.

    Parameters
    ----------
    x     : ndarray — input signal
    fs    : int — sampling frequency (Hz)
    n_fft : int — DFT length (None = len(x), i.e. no zero-padding)

    Returns
    -------
    freqs : ndarray — frequency axis in Hz (positive half only)
    X     : ndarray — complex DFT coefficients (positive half only)
    """
    if n_fft is None:
        n_fft = len(x)  # DFT length = signal length
    X_full = np.fft.fft(x, n=n_fft)  # full DFT (Equation A.3)
    freqs_full = np.fft.fftfreq(n_fft, d=1 / fs)  # frequency axis in Hz
    pos = freqs_full >= 0  # positive frequencies only
    return freqs_full[pos], X_full[pos]  # return positive half


def compute_power(X, N):
    """
    Compute power from DFT coefficients: P[k] = |X[k]|² / N

    Parameters
    ----------
    X : ndarray — complex DFT coefficients
    N : int — original signal length (for normalization)

    Returns
    -------
    P : ndarray — power at each bin
    """
    return np.abs(X) ** 2 / N  # Equation (A.29)


# ============================================================
# 2. Experiment A: Single tone — on-bin vs. off-bin
# ============================================================
def experiment_on_vs_off_bin():
    """
    Compare DFT of a tone at an on-bin frequency vs. an off-bin frequency.
    Demonstrates leakage when the tone does not land on a bin.
    """
    # --- Generate signals (short segment: Δf = 1 Hz, leakage visible) ---
    x_on, n_on, t_on = make_tone(F0_ON_BIN, A=A, duration=DURATION_ONOFF)
    x_off, n_off, t_off = make_tone(F0_OFF_BIN, A=A, duration=DURATION_ONOFF)

    # --- Verify bin spacing ---
    N = len(x_on)  # number of samples
    delta_f = FS / N  # bin spacing (Equation A.6)
    print(f"N = {N}, fs = {FS} Hz, duration = {DURATION_ONOFF} s")
    print(f"Bin spacing Δf = fs/N = {FS}/{N} = {delta_f:.4f} Hz")
    print(
        f"On-bin tone:  f0 = {F0_ON_BIN} Hz, f0/Δf = {F0_ON_BIN/delta_f:.1f} (integer = on-bin)"
    )
    print(
        f"Off-bin tone: f0 = {F0_OFF_BIN} Hz, f0/Δf = {F0_OFF_BIN/delta_f:.1f} (half-integer = maximally off-bin)"
    )

    # --- Compute DFTs ---
    freqs_on, X_on = compute_dft(x_on)  # DFT of on-bin tone
    freqs_off, X_off = compute_dft(x_off)  # DFT of off-bin tone
    P_on = compute_power(X_on, N)  # power spectrum (on-bin)
    P_off = compute_power(X_off, N)  # power spectrum (off-bin)

    # --- Time-domain plots (always first) ---
    t_show = (0, 0.5)  # show first 0.5 s for clarity
    fig1, _ = plot_time_domain(
        t_on,
        x_on,
        title="On-bin tone at 10.0 Hz (time domain)",
        fig_id="Figure B.1a",
        t_range=t_show,
    )
    fig2, _ = plot_time_domain(
        t_off,
        x_off,
        title="Off-bin tone at 10.5 Hz (time domain)",
        fig_id="Figure B.1b",
        t_range=t_show,
    )

    # --- Dual-stack spectrum plots ---
    f_show = (0, 30)  # zoom to 0–30 Hz
    fig3, _ = plot_dual_stack_spectrum(
        freqs_on,
        P_on,
        title="On-bin tone (10.0 Hz) — no leakage",
        fig_id="Figure B.2a",
        f_range=f_show,
    )
    fig4, _ = plot_dual_stack_spectrum(
        freqs_off,
        P_off,
        title="Off-bin tone (10.5 Hz) — maximum leakage",
        fig_id="Figure B.2b",
        f_range=f_show,
    )

    return {"delta_f": delta_f, "N": N, "figs": [fig1, fig2, fig3, fig4]}


# ============================================================
# 3. Experiment B: Dual-tone chord — resolution
# ============================================================
def experiment_dual_tone():
    """
    DFT of a two-tone signal. Tests whether the bin spacing resolves
    the two components.
    """
    # --- Generate signal ---
    x, n, t = make_mixed_tones(
        [F1_CHORD, F2_CHORD], amplitudes=[A, A], duration=DURATION_LAB
    )

    N = len(x)  # number of samples
    delta_f = FS / N  # bin spacing
    separation = abs(F2_CHORD - F1_CHORD)  # tone separation (Hz)
    print(f"\nDual-tone: f1={F1_CHORD} Hz, f2={F2_CHORD} Hz")
    print(f"Separation = {separation} Hz, Δf = {delta_f:.6f} Hz")
    print(f"Separation / Δf = {separation/delta_f:.0f} bins — should be resolvable")

    # --- Compute DFT ---
    freqs, X = compute_dft(x)  # DFT of chord
    P = compute_power(X, N)  # power spectrum

    # --- Time-domain (first) ---
    fig1, _ = plot_time_domain(
        t,
        x,
        title="Dual-tone chord (time domain)",
        fig_id="Figure B.3a",
        t_range=(0, 0.5),
    )

    # --- Dual-stack spectrum ---
    fig2, _ = plot_dual_stack_spectrum(
        freqs,
        P,
        title="Dual-tone chord spectrum",
        fig_id="Figure B.3b",
        f_range=(0, 30),
    )

    return {"separation": separation, "delta_f": delta_f, "figs": [fig1, fig2]}


# ============================================================
# 4. Experiment C: Zero-padding — bins vs. resolution
# ============================================================
def experiment_zero_padding():
    """
    Zero-pad a short segment and show that bin count increases
    but frequency resolution does not.
    """
    # --- Two tones BELOW the resolution limit ---
    short_dur = 1.0  # 1-second segment → Δf_min = 1.0 Hz
    f1_zp, f2_zp = 10.0, 10.5  # separation = 0.5 Hz < Δf_min = 1.0 Hz
    x, n, t = make_mixed_tones([f1_zp, f2_zp], amplitudes=[1.0, 1.0], duration=short_dur)
    N_orig = len(x)  # original length
    N_padded = N_orig * ZERO_PAD_FACTOR  # zero-padded length

    delta_f_orig = FS / N_orig  # original bin spacing
    delta_f_padded = FS / N_padded  # padded bin spacing
    separation = f2_zp - f1_zp  # 0.5 Hz
    print("\nZero-padding experiment:")
    print(f"Original: N={N_orig}, Δf={delta_f_orig:.2f} Hz")
    print(f"Padded:   N={N_padded}, Δf={delta_f_padded:.4f} Hz")
    print(f"Resolution (determined by original N): Δf_min = {delta_f_orig:.2f} Hz")
    print(f"Tones at {f1_zp} Hz and {f2_zp} Hz - separation = {separation} Hz < Δf_min = {delta_f_orig:.2f} Hz")

    # --- Compute DFTs ---
    freqs_orig, X_orig = compute_dft(x, n_fft=N_orig)  # no zero-padding
    freqs_pad, X_pad = compute_dft(x, n_fft=N_padded)  # zero-padded
    P_orig = compute_power(X_orig, N_orig)  # power (original)
    P_pad = compute_power(X_pad, N_orig)  # power (padded, same normalization)

    # --- Dual-stack: original vs. zero-padded ---
    f_show = (5, 16)  # zoom to region of interest
    fig1, _ = plot_dual_stack_spectrum(
        freqs_orig,
        P_orig,
        title=f"10 Hz + 10.5 Hz, original (no zero-padding) - NOT resolved",
        fig_id="Figure B.4a",
        f_range=f_show,
    )
    fig2, _ = plot_dual_stack_spectrum(
        freqs_pad,
        P_pad,
        title=f"10 Hz + 10.5 Hz, zero-padded x{ZERO_PAD_FACTOR} - still NOT resolved",
        fig_id="Figure B.4b",
        f_range=f_show,
    )

    return {
        "delta_f_orig": delta_f_orig,
        "delta_f_padded": delta_f_padded,
        "figs": [fig1, fig2],
    }


# ============================================================
# 5. Run all experiments
# ============================================================
def run_lab1(save=False):
    """Run all Lab 1 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 1: The DFT of Basic Signals")
    print("=" * 60)

    results_a = experiment_on_vs_off_bin()  # Experiment A
    results_b = experiment_dual_tone()  # Experiment B
    results_c = experiment_zero_padding()  # Experiment C

    if save:  # save all figures
        all_figs = results_a["figs"] + results_b["figs"] + results_c["figs"]
        for i, fig in enumerate(all_figs, start=1):
            save_figure(fig, lab_number=1, fig_id=f"{i:02d}")

    plt.show()  # display all figures
    return results_a, results_b, results_c


if __name__ == "__main__":
    run_lab1(save=True)
