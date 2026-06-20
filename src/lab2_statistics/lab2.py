"""
Lab 2: Statistics on a Noisy Signal (↔ A.4)
---------------------------------------------
Tests: exponential distribution of bin power, noise floor estimation,
       spectral detection via threshold γ, Welch averaging.
"""
import numpy as np                                        # numerical computing
import matplotlib.pyplot as plt                           # plotting
from scipy import signal as sp_signal                     # Welch's method
import sys, os                                            # path manipulation
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (                                  # project utilities
    FS, DURATION, SEED, DPI,
    make_tone, make_noise, make_time_axis,
    plot_time_domain, plot_dual_stack_spectrum, save_figure,
)
from src.common.config import COLORMAP                    # colormap for histograms

# ============================================================
# Parameters (Table B.2 in the report)
# ============================================================
F0 = 10.0                 # tone frequency (Hz)
A_TONE = 0.5              # tone amplitude (deliberately small)
SIGMA_NOISE = 1.0         # noise standard deviation
DURATION_LAB = DURATION   # signal duration (s)
GAMMA_VALUES = [3.0, 4.6, 6.9]  # detection thresholds (Table A.3)

# Welch parameters
WELCH_SEGMENTS = [         # (segment_length_in_seconds, label)
    (DURATION, "Full (single periodogram)"),
    (20.0, "20 s segments"),
    (5.0, "5 s segments"),
    (2.0, "2 s segments"),
]


# ============================================================
# 1. Helper: compute periodogram
# ============================================================
def compute_periodogram(x, fs=FS):
    """
    Compute the one-sided periodogram: S[k] = |X[k]|² / N (Equation A.29).

    Parameters
    ----------
    x  : ndarray — input signal
    fs : int — sampling frequency (Hz)

    Returns
    -------
    freqs : ndarray — frequency axis (Hz), positive half
    Sxx   : ndarray — power spectral density estimate
    """
    N = len(x)                                            # signal length
    X = np.fft.fft(x)                                     # full DFT
    freqs = np.fft.fftfreq(N, d=1/fs)                     # frequency axis
    pos = freqs >= 0                                      # positive frequencies
    Sxx = np.abs(X[pos])**2 / N                           # periodogram (Eq. A.29)
    return freqs[pos], Sxx


# ============================================================
# 2. Experiment A: Bin distributions under noise
# ============================================================
def experiment_bin_distributions():
    """
    Generate pure noise, compute DFT, examine magnitude, phase,
    and power distributions of the bins.
    """
    # --- Generate pure noise ---
    x_noise, n, t = make_noise(sigma=SIGMA_NOISE, duration=DURATION_LAB, seed=SEED)

    N = len(x_noise)                                      # number of samples
    X = np.fft.fft(x_noise)                               # full DFT
    freqs = np.fft.fftfreq(N, d=1/FS)                     # frequency axis
    pos = (freqs > 0) & (freqs < FS/2)                    # exclude DC and Nyquist

    magnitudes = np.abs(X[pos])                           # |X[k]|
    phases = np.angle(X[pos])                             # ∠X[k]
    powers = np.abs(X[pos])**2                            # |X[k]|²

    # --- Expected values ---
    expected_power_mean = N * SIGMA_NOISE**2               # E[|X[k]|²] = Nσ² (Eq. A.26)
    print(f"N = {N}, σ = {SIGMA_NOISE}")
    print(f"Expected mean bin power: Nσ² = {expected_power_mean:.1f}")
    print(f"Measured mean bin power: {np.mean(powers):.1f}")

    # --- Histograms ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))       # three histograms side by side

    # Magnitude (Rayleigh distribution)
    axes[0].hist(magnitudes, bins=100, density=True, alpha=0.7, label="Measured")
    axes[0].set_xlabel("Magnitude |X[k]|")                # x label with units
    axes[0].set_ylabel("Probability density")             # y label
    axes[0].set_title("Magnitude distribution (Rayleigh)")
    axes[0].legend()

    # Phase (uniform on (-π, π])
    axes[1].hist(phases, bins=100, density=True, alpha=0.7, label="Measured")
    axes[1].axhline(1/(2*np.pi), color="red", ls="--", label="Uniform 1/(2π)")
    axes[1].set_xlabel("Phase ∠X[k] (radians)")           # x label with units
    axes[1].set_ylabel("Probability density")             # y label
    axes[1].set_title("Phase distribution (Uniform)")
    axes[1].legend()

    # Power (exponential distribution)
    axes[2].hist(powers, bins=100, density=True, alpha=0.7, label="Measured")
    lam = 1 / expected_power_mean                          # exponential rate (Eq. A.25)
    p_range = np.linspace(0, np.percentile(powers, 99), 200)
    axes[2].plot(p_range, lam * np.exp(-lam * p_range), "r--", label="Exponential fit")
    axes[2].set_xlabel("Power |X[k]|²")                   # x label with units
    axes[2].set_ylabel("Probability density")             # y label
    axes[2].set_title("Power distribution (Exponential)")
    axes[2].legend()

    fig.suptitle("Figure B.5 — DFT bin distributions under white Gaussian noise",
                 fontsize=11, fontstyle="italic")
    fig.set_dpi(DPI)                                      # project DPI
    fig.tight_layout()                                    # prevent clipping
    return fig


# ============================================================
# 3. Experiment B: Spectral detection of a tone in noise
# ============================================================
def experiment_spectral_detection():
    """
    Bury a tone in noise, estimate the noise floor from the spectrum,
    and apply detection thresholds from Table A.3.
    """
    # --- Generate tone + noise ---
    x_tone, _, _ = make_tone(F0, A=A_TONE, duration=DURATION_LAB)
    x_noise, _, _ = make_noise(sigma=SIGMA_NOISE, duration=DURATION_LAB, seed=SEED)
    x = x_tone + x_noise                                  # combined signal
    _, t = make_time_axis(duration=DURATION_LAB)

    N = len(x)                                            # number of samples

    # --- Time-domain (always first) ---
    fig_td, _ = plot_time_domain(t, x, title="Tone buried in noise (time domain)",
                                  fig_id="Figure B.6", t_range=(0, 2))

    # --- Periodogram ---
    freqs, Sxx = compute_periodogram(x)                   # periodogram of noisy signal

    # --- Noise floor estimation (median, Equation A.27) ---
    noise_floor = np.median(Sxx)                          # robust estimator
    print(f"\nSpectral detection:")
    print(f"Noise floor (median): {noise_floor:.2f}")

    # --- Detection thresholds (Equation A.28) ---
    fig_det, (ax_lin, ax_db) = plot_dual_stack_spectrum(
        freqs, Sxx,
        title="Tone-in-noise periodogram with detection thresholds",
        fig_id="Figure B.7", f_range=(0, 30)
    )

    colors = ["green", "orange", "red"]                   # threshold colors
    for gamma, color in zip(GAMMA_VALUES, colors):
        threshold = gamma * noise_floor                   # γ × noise floor
        false_alarm = np.exp(-gamma)                      # Pr(exceed | noise only)
        label = f"γ={gamma} (P_fa={false_alarm:.3f})"     # label with false-alarm prob
        ax_lin.axhline(threshold, color=color, ls="--", lw=1, label=label)
        ax_db.axhline(10*np.log10(max(threshold, 1e-20)),  # same threshold in dB
                      color=color, ls="--", lw=1, label=label)
        print(f"  γ={gamma}: threshold={threshold:.2f}, P_fa={false_alarm:.4f}")

    ax_lin.legend(fontsize=7, loc="upper right")          # legend on linear panel
    ax_db.legend(fontsize=7, loc="upper right")           # legend on dB panel

    # --- Check: does the tone bin exceed the threshold? ---
    bin_f0 = np.argmin(np.abs(freqs - F0))                # bin nearest to f0
    power_at_f0 = Sxx[bin_f0]                             # power at tone bin
    ratio = power_at_f0 / noise_floor                     # ratio to noise floor
    print(f"  Power at f0={F0} Hz: {power_at_f0:.2f}, ratio to floor: {ratio:.1f}")
    print(f"  Detected at γ=3? {'Yes' if ratio > 3 else 'No'}")
    print(f"  Detected at γ=6.9? {'Yes' if ratio > 6.9 else 'No'}")

    return fig_td, fig_det


# ============================================================
# 4. Experiment C: Welch averaging — variance reduction
# ============================================================
def experiment_welch_averaging():
    """
    Apply Welch's method with different segment lengths to the same
    tone-in-noise signal. Show variance reduction vs. resolution tradeoff.
    """
    # --- Generate signal (same as Experiment B) ---
    x_tone, _, _ = make_tone(F0, A=A_TONE, duration=DURATION_LAB)
    x_noise, _, _ = make_noise(sigma=SIGMA_NOISE, duration=DURATION_LAB, seed=SEED)
    x = x_tone + x_noise                                  # combined signal

    n_panels = len(WELCH_SEGMENTS)                        # number of segment lengths
    fig, axes = plt.subplots(n_panels, 2, figsize=(14, 3*n_panels),
                              sharex="col")                # columns share x-axis

    for i, (seg_dur, label) in enumerate(WELCH_SEGMENTS):
        nperseg = int(seg_dur * FS)                       # segment length in samples
        nperseg = min(nperseg, len(x))                    # cap at signal length
        noverlap = nperseg // 2                           # 50% overlap (Hann COLA)

        # --- Welch PSD ---
        f_welch, Pxx = sp_signal.welch(x, fs=FS,         # scipy Welch implementation
                                        nperseg=nperseg,
                                        noverlap=noverlap,
                                        window="hann")
        n_segs = (len(x) - nperseg) // (nperseg - noverlap) + 1  # number of segments
        delta_f = FS / nperseg                            # frequency resolution

        # --- Linear scale (left column) ---
        axes[i, 0].plot(f_welch, Pxx, linewidth=0.5)     # linear PSD
        axes[i, 0].set_ylabel("PSD (linear)")            # y label
        axes[i, 0].set_title(f"{label} — L={n_segs}, Δf={delta_f:.3f} Hz",
                              fontsize=9)
        axes[i, 0].grid(True, alpha=0.3)                  # light grid

        # --- dB scale (right column) ---
        Pxx_db = 10 * np.log10(np.maximum(Pxx, 1e-20))   # dB conversion
        axes[i, 1].plot(f_welch, Pxx_db, linewidth=0.5)  # dB PSD
        axes[i, 1].set_ylabel("PSD (dB)")                # y label
        axes[i, 1].set_title(f"{label} — L={n_segs}, Δf={delta_f:.3f} Hz",
                              fontsize=9)
        axes[i, 1].grid(True, alpha=0.3)                  # light grid

        print(f"\nWelch: {label}")
        print(f"  Segment length: {nperseg} samples ({seg_dur} s)")
        print(f"  Segments L: {n_segs}")
        print(f"  Δf: {delta_f:.4f} Hz")
        print(f"  Relative variance: {1/n_segs:.4f}")

    axes[-1, 0].set_xlabel("Frequency (Hz)")              # x label on bottom
    axes[-1, 1].set_xlabel("Frequency (Hz)")              # x label on bottom
    axes[0, 0].set_xlim(0, 30)                            # zoom to 0–30 Hz
    axes[0, 1].set_xlim(0, 30)                            # zoom to 0–30 Hz

    fig.suptitle("Figure B.8 — Welch averaging: variance vs. resolution tradeoff",
                 fontsize=11, fontstyle="italic")
    fig.set_dpi(DPI)                                      # project DPI
    fig.tight_layout()                                    # prevent clipping
    return fig


# ============================================================
# 5. Run all experiments
# ============================================================
def run_lab2(save=False):
    """Run all Lab 2 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 2: Statistics on a Noisy Signal")
    print("=" * 60)

    fig_dist = experiment_bin_distributions()              # Experiment A
    fig_td, fig_det = experiment_spectral_detection()      # Experiment B
    fig_welch = experiment_welch_averaging()                # Experiment C

    if save:                                              # save all figures
        all_figs = [fig_dist, fig_td, fig_det, fig_welch]
        for i, fig in enumerate(all_figs, start=1):
            save_figure(fig, lab_number=2, fig_id=f"{i:02d}")

    plt.show()                                            # display all figures
    return fig_dist, fig_td, fig_det, fig_welch


if __name__ == "__main__":
    run_lab2(save=True)
