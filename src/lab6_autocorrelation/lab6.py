"""
Lab 6: Autocorrelation of a Noisy Signal (↔ A.6)
--------------------------------------------------
Tests: autocorrelation as periodicity detector, Wiener-Khinchin theorem,
       phase-blindness.

Figures B.31–B.36 in the report.
"""

import numpy as np  # numerical computing
import matplotlib.pyplot as plt  # plotting
import sys
import os  # path manipulation

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (  # project utilities
    FS,
    DURATION,
    SEED,
    DPI,
    make_tone,
    make_noise,
    make_time_axis,
    plot_time_domain,
    plot_dual_stack_spectrum,
    save_figure,
)

# ============================================================
# Parameters (Table B.9 in the report)
# ============================================================
F0 = 10.0  # tone frequency (Hz)
A_TONE = 0.5  # tone amplitude (small, buried in noise)
SIGMA_NOISE = 1.0  # noise standard deviation
DURATION_LAB = 60.0  # 60 s — shorter for faster autocorrelation computation
LAB_NUMBER = 6  # lab number for figure saving
FIG_START = 31  # first figure number (B.31)


# ============================================================
# 1. Compute autocorrelation from definition
# ============================================================
def compute_autocorrelation(x):
    """
    Compute the autocorrelation r[l] = Σ x[n] x[n-l] (Equation A.44).
    Uses np.correlate for efficiency; returns positive lags only.

    Parameters
    ----------
    x : ndarray — input signal

    Returns
    -------
    lags   : ndarray — lag indices [0, 1, ..., N-1]
    r      : ndarray — autocorrelation values at each lag
    """
    N = len(x)  # signal length
    r_full = np.correlate(x, x, mode="full")  # full autocorrelation [-N+1, ..., N-1]
    r = r_full[N - 1:]  # positive lags only [0, 1, ..., N-1]
    lags = np.arange(len(r))  # lag indices
    return lags, r


# ============================================================
# 2. Experiment A: Periodicity detection in noise
# ============================================================
def experiment_periodicity():
    """
    Autocorrelation of a tone buried in noise.
    The tone is invisible in time domain but periodic peaks
    emerge in the autocorrelation at multiples of the tone period.
    """
    print("\nExperiment A: Periodicity detection")
    print("-" * 40)

    # --- Generate tone + noise ---
    x_tone, _, _ = make_tone(F0, A=A_TONE, duration=DURATION_LAB)
    x_noise, _, _ = make_noise(sigma=SIGMA_NOISE, duration=DURATION_LAB, seed=SEED)
    x = x_tone + x_noise  # combined signal
    _, t = make_time_axis(duration=DURATION_LAB)

    N = len(x)  # number of samples
    period_samples = int(FS / F0)  # expected period: fs/f0 = 25 samples

    # --- Time-domain plot (always first) ---
    fig_td, _ = plot_time_domain(
        t, x,
        title=f"Tone (f₀={F0} Hz, A={A_TONE}) buried in noise (σ={SIGMA_NOISE})",
        fig_id=f"Figure B.{FIG_START}",
        t_range=(0, 2),
    )

    # --- Compute autocorrelation ---
    lags, r = compute_autocorrelation(x)
    lag_time = lags / FS  # convert lags to seconds

    # --- Verify r[0] = total energy (Equation A.45) ---
    energy_time = np.sum(x**2)  # energy from time domain
    energy_autocorr = r[0]  # energy from autocorrelation
    print(f"  N = {N}, period = {period_samples} samples ({1/F0:.2f} s)")
    print(f"  r[0] = {energy_autocorr:.2f}, Σ|x[n]|² = {energy_time:.2f} "
          f"(match: {np.isclose(energy_autocorr, energy_time)})")

    # --- Plot autocorrelation (zoomed to first ~10 periods) ---
    max_lag = period_samples * 10  # show 10 periods
    fig_ac, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(12, 6))

    # Top: full range (shows lag-0 spike)
    ax_full.plot(lag_time[:N // 2], r[:N // 2], linewidth=0.3)
    ax_full.set_xlabel("Lag (s)")
    ax_full.set_ylabel("r[l]")
    ax_full.set_title(f"Figure B.{FIG_START + 1}", loc="left",
                       fontsize=9, fontstyle="italic")
    ax_full.set_title("Autocorrelation — full range (lag-0 spike dominates)")
    ax_full.grid(True, alpha=0.3)

    # Bottom: zoomed (periodic peaks visible)
    ax_zoom.plot(lag_time[:max_lag], r[:max_lag], linewidth=0.8)
    ax_zoom.set_xlabel("Lag (s)")
    ax_zoom.set_ylabel("r[l]")
    ax_zoom.set_title("Autocorrelation — zoomed (periodic peaks at multiples of 1/f₀)")
    ax_zoom.grid(True, alpha=0.3)

    # Mark expected period multiples
    for k in range(1, 11):
        lag_k = k * period_samples
        if lag_k < max_lag:
            ax_zoom.axvline(lag_k / FS, color="red", ls="--", lw=0.5, alpha=0.6)
    ax_zoom.annotate(
        f"Period = {period_samples} samples = {1/F0:.2f} s",
        xy=(period_samples / FS, r[period_samples]),
        xytext=(period_samples / FS + 0.05, r[period_samples] * 1.3),
        fontsize=8, color="red",
        arrowprops=dict(arrowstyle="->", color="red"),
    )

    fig_ac.set_dpi(DPI)
    fig_ac.tight_layout()

    return fig_td, fig_ac, lags, r, x


# ============================================================
# 3. Experiment B: Wiener-Khinchin verification
# ============================================================
def experiment_wiener_khinchin(x, lags, r):
    """
    Verify that DFT of autocorrelation = power spectrum |X[k]|².
    Equation (A.47).
    """
    print("\nExperiment B: Wiener-Khinchin verification")
    print("-" * 40)

    N = len(x)  # signal length

    # --- Method 1: |X[k]|² directly ---
    X = np.fft.fft(x)  # DFT of signal
    power_direct = np.abs(X)**2  # |X[k]|²

    # --- Method 2: DFT of autocorrelation ---
    # The autocorrelation r[l] for l=0..N-1 is one-sided.
    # For Wiener-Khinchin, we need the periodic autocorrelation:
    # r_periodic[l] = Σ x[n] x[(n+l) mod N]  which equals IFFT{|X[k]|²}
    # Equivalently: |X[k]|² = FFT{r_periodic[l]}
    # We can compute r_periodic via: IFFT(|FFT(x)|²) then FFT back.
    # But the direct way: use the fact that for a finite signal,
    # the circular autocorrelation via FFT is:
    r_periodic = np.fft.ifft(power_direct).real  # circular autocorrelation
    power_from_autocorr = np.abs(np.fft.fft(r_periodic))  # back to |X[k]|²

    # --- Compare ---
    freqs = np.fft.fftfreq(N, d=1/FS)
    pos = freqs >= 0

    max_error = np.max(np.abs(power_direct[pos] - power_from_autocorr[pos]))
    rel_error = max_error / np.max(power_direct[pos])
    print(f"  Max absolute error: {max_error:.6f}")
    print(f"  Relative error: {rel_error:.2e}")
    print(f"  Match: {rel_error < 1e-10}")

    # --- Plot overlay ---
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    # Linear scale
    ax_lin.plot(freqs[pos], power_direct[pos] / N, linewidth=0.8,
                label="|X[k]|²/N (direct)", alpha=0.8)
    ax_lin.plot(freqs[pos], power_from_autocorr[pos] / N, linewidth=0.8,
                label="DFT{r[l]}/N (Wiener-Khinchin)", ls="--", alpha=0.8)
    ax_lin.set_ylabel("Power")
    ax_lin.set_title(f"Figure B.{FIG_START + 2}", loc="left",
                      fontsize=9, fontstyle="italic")
    ax_lin.set_title("Wiener-Khinchin: |X[k]|² vs DFT{r[l]}")
    ax_lin.legend(fontsize=8)
    ax_lin.set_xlim(0, 30)
    ax_lin.grid(True, alpha=0.3)

    # dB scale
    P_direct_db = 10 * np.log10(np.maximum(power_direct[pos] / N, 1e-20))
    P_autocorr_db = 10 * np.log10(np.maximum(power_from_autocorr[pos] / N, 1e-20))
    ax_db.plot(freqs[pos], P_direct_db, linewidth=0.8,
               label="|X[k]|²/N (direct)", alpha=0.8)
    ax_db.plot(freqs[pos], P_autocorr_db, linewidth=0.8,
               label="DFT{r[l]}/N (Wiener-Khinchin)", ls="--", alpha=0.8)
    ax_db.set_ylabel("Power (dB)")
    ax_db.set_xlabel("Frequency (Hz)")
    ax_db.legend(fontsize=8)
    ax_db.set_xlim(0, 30)
    ax_db.grid(True, alpha=0.3)

    fig.set_dpi(DPI)
    fig.tight_layout()

    return fig


# ============================================================
# 4. Experiment C: Phase-blindness
# ============================================================
def experiment_phase_blindness():
    """
    Two tones with different phases produce identical autocorrelations
    and power spectra — proving autocorrelation discards phase.
    """
    print("\nExperiment C: Phase-blindness")
    print("-" * 40)

    dur = 10.0  # short signal for clarity

    # --- Two tones: same frequency, different phase ---
    x1, _, t = make_tone(F0, A=1.0, phi=0.0, duration=dur)  # starts at peak
    x2, _, _ = make_tone(F0, A=1.0, phi=np.pi, duration=dur)  # starts at trough

    # --- Time domain: they look different ---
    fig_td, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 4), sharex=True)
    ax1.plot(t[:500], x1[:500], linewidth=0.8)  # first 2 seconds
    ax1.set_ylabel("Amplitude")
    ax1.set_title(f"Figure B.{FIG_START + 3}", loc="left",
                   fontsize=9, fontstyle="italic")
    ax1.set_title("φ = 0 (starts at peak)")
    ax1.grid(True, alpha=0.3)
    ax2.plot(t[:500], x2[:500], linewidth=0.8, color="tab:orange")
    ax2.set_ylabel("Amplitude")
    ax2.set_xlabel("Time (s)")
    ax2.set_title("φ = π (starts at trough)")
    ax2.grid(True, alpha=0.3)
    fig_td.set_dpi(DPI)
    fig_td.tight_layout()

    # --- Autocorrelations: they are identical ---
    lags1, r1 = compute_autocorrelation(x1)
    lags2, r2 = compute_autocorrelation(x2)
    max_lag = 200  # show ~8 periods

    fig_ac, ax = plt.subplots(1, 1, figsize=(12, 4))
    ax.plot(lags1[:max_lag] / FS, r1[:max_lag], linewidth=1.0,
            label="φ = 0", alpha=0.8)
    ax.plot(lags2[:max_lag] / FS, r2[:max_lag], linewidth=1.0,
            label="φ = π", ls="--", alpha=0.8)
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("r[l]")
    ax.set_title(f"Figure B.{FIG_START + 4}", loc="left",
                  fontsize=9, fontstyle="italic")
    ax.set_title("Autocorrelation: φ = 0 vs φ = π (identical)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig_ac.set_dpi(DPI)
    fig_ac.tight_layout()

    # --- Verify numerically ---
    diff = np.max(np.abs(r1 - r2))
    print(f"  Max |r1[l] - r2[l]|: {diff:.2e}")
    print(f"  Identical: {diff < 1e-6}")

    # --- Power spectra: also identical ---
    X1 = np.fft.fft(x1)
    X2 = np.fft.fft(x2)
    P1 = np.abs(X1)**2
    P2 = np.abs(X2)**2
    power_diff = np.max(np.abs(P1 - P2))
    print(f"  Max ||X1[k]|² - |X2[k]|²|: {power_diff:.2e}")
    print(f"  Power spectra identical: {power_diff < 1e-6}")

    return fig_td, fig_ac


# ============================================================
# 5. Run all experiments
# ============================================================
def run_lab6(save=False):
    """Run all Lab 6 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 6: Autocorrelation of a Noisy Signal")
    print("=" * 60)

    fig_td_a, fig_ac_a, lags, r, x = experiment_periodicity()  # Experiment A
    fig_wk = experiment_wiener_khinchin(x, lags, r)  # Experiment B
    fig_td_c, fig_ac_c = experiment_phase_blindness()  # Experiment C

    all_figs = [
        (fig_td_a, False),  # B.31: tone-in-noise time domain
        (fig_ac_a, False),  # B.32: autocorrelation (full + zoomed)
        (fig_wk, False),  # B.33: Wiener-Khinchin overlay
        (fig_td_c, False),  # B.34: phase comparison time domain
        (fig_ac_c, False),  # B.35: phase comparison autocorrelation
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
    run_lab6(save=True)
