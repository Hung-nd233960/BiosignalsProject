"""
Lab 3: Windowing and the Dirichlet Kernel (↔ A.3)
---------------------------------------------------
Derives the DFT of the rectangular window, analyzes its properties,
and extends to Hann, Hamming, Blackman via the cosine-sum identity.

Figures B.9–B.14 in the report.
"""
import numpy as np                                        # numerical computing
import matplotlib.pyplot as plt                           # plotting
from scipy.signal import argrelextrema                    # local maxima finder
from scipy.stats import linregress                        # linear regression
import sys, os                                            # path manipulation
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import DPI, save_figure                   # project utilities
from src.common.windows import rectangular, hann, hamming, blackman

# ============================================================
# Parameters (Table B.4 in the report)
# ============================================================
M = 256                   # typical EEG STFT window: 256 samples ≈ 1.024 s at 250 Hz
N_DENSE = M * 2048        # dense sampling for smooth continuous-looking curves
LAB_NUMBER = 3            # lab number for figure saving
FIG_START = 9             # first figure number (B.9), continuing from Lab 2's B.8


# ============================================================
# Core: compute the normalized Dirichlet kernel magnitude
# ============================================================
def dirichlet_kernel(omega, M):
    """
    Normalized Dirichlet kernel magnitude:
      |D(ω)| = |sin(ωM/2) / sin(ω/2)| / M

    Parameters
    ----------
    omega : ndarray — angular frequency in radians
    M     : int — window length

    Returns
    -------
    D : ndarray — normalized magnitude (peak = 1.0 at ω = 0)
    """
    numerator = np.sin(omega * M / 2)                     # sin(ωM/2)
    denominator = np.sin(omega / 2)                       # sin(ω/2)
    with np.errstate(divide="ignore", invalid="ignore"):  # handle 0/0 at ω=0
        D = np.abs(numerator / denominator) / M           # |sin(ωM/2)/sin(ω/2)| / M
    D[np.isnan(D)] = 1.0                                  # limit at ω=0 is 1.0
    return D


def window_spectrum_dense(win_func, M, N_dense=N_DENSE):
    """
    Compute the magnitude spectrum of a window at dense frequency sampling.
    Returns continuous-looking curve in bins (not Hz).

    Parameters
    ----------
    win_func : callable — window function (e.g. hann)
    M        : int — window length
    N_dense  : int — DFT length (high zero-pad for smooth curve)

    Returns
    -------
    bins  : ndarray — frequency axis in bins (centered, -M/2 to M/2)
    W_mag : ndarray — normalized magnitude spectrum (peak = 1.0)
    """
    w = win_func(M)                                       # generate window
    W = np.fft.fft(w, n=N_dense)                          # zero-padded DFT
    W_shifted = np.fft.fftshift(W)                        # center at DC
    W_mag = np.abs(W_shifted)                             # magnitude
    W_mag = W_mag / W_mag.max()                           # normalize peak to 1.0
    bins = (np.arange(N_dense) - N_dense // 2) / (N_dense / M)  # bin axis
    return bins, W_mag


def find_sidelobe_maxima(bins, W_mag, main_lobe_bins=1.5):
    """
    Find side-lobe maxima positions and values.

    Parameters
    ----------
    bins           : ndarray — frequency axis in bins
    W_mag          : ndarray — normalized magnitude
    main_lobe_bins : float — half-width of main lobe to exclude

    Returns
    -------
    peak_bins   : ndarray — bin positions of side-lobe maxima (positive side only)
    peak_values : ndarray — magnitude values at those maxima
    """
    right_half = bins > main_lobe_bins                    # positive side, outside main lobe
    indices = np.where(right_half)[0]                     # indices in that region
    W_right = W_mag[indices]                              # magnitude in that region
    bins_right = bins[indices]                            # bins in that region

    local_max_idx = argrelextrema(W_right, np.greater, order=5)[0]  # local maxima
    peak_bins = bins_right[local_max_idx]                 # bin positions
    peak_values = W_right[local_max_idx]                  # magnitude values
    return peak_bins, peak_values


# ============================================================
# Figure B.9 — Dirichlet kernel anatomy (zoomed to side lobes)
# ============================================================
def figure_b9_dirichlet_anatomy():
    """Desmos-style graph of the Dirichlet kernel, annotated."""
    bins, D_mag = window_spectrum_dense(rectangular, M)

    null_bins = np.arange(1, 8)                           # nulls at integer bins
    peak_bins, peak_values = find_sidelobe_maxima(bins, D_mag, main_lobe_bins=1.2)
    midpoints = np.array([k + 0.5 for k in range(1, 7)])  # (2k+1)/2 rule

    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    ax.set_facecolor("white")

    plot_range = (bins >= -10) & (bins <= 10)
    ax.plot(bins[plot_range], D_mag[plot_range],
            color="#2563eb", linewidth=1.5, zorder=3)

    for k in null_bins:                                   # annotate nulls
        ax.axvline(k, color="#94a3b8", linewidth=0.5, linestyle=":", zorder=1)
        ax.axvline(-k, color="#94a3b8", linewidth=0.5, linestyle=":", zorder=1)

    for i, (b, v) in enumerate(zip(peak_bins[:6], peak_values[:6])):
        ax.plot(b, v, "o", color="#dc2626", markersize=5, zorder=4)
        if i < 3:
            ax.annotate(f"({b:.2f}, {v:.4f})",
                        xy=(b, v), xytext=(b + 0.3, v + 0.02),
                        fontsize=7, color="#dc2626",
                        arrowprops=dict(arrowstyle="-", color="#dc2626", lw=0.5))

    for mp in midpoints[:3]:                              # midpoint approximations
        ax.axvline(mp, color="#f59e0b", linewidth=0.5, linestyle="--", alpha=0.7)

    ax.annotate("Main lobe", xy=(0, 0.85), fontsize=10, ha="center",
                color="#2563eb", fontstyle="italic")
    ax.annotate("← side lobes →", xy=(4, 0.08), fontsize=9, ha="center",
                color="#dc2626", fontstyle="italic")

    ax.set_xlabel("Frequency (bins)", fontsize=11)
    ax.set_ylabel("Normalized magnitude $|W(\\omega)| / M$", fontsize=11)
    ax.set_title(f"Figure B.9 — Dirichlet kernel anatomy ($M = {M}$)",
                 fontsize=12, fontstyle="italic", loc="left")
    ax.set_xlim(-10, 10)
    ax.set_ylim(-0.01, 0.35)
    ax.grid(True, alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig


# ============================================================
# Figure B.10 — Dirichlet kernel with envelope (full range)
# ============================================================
def figure_b10_dirichlet_full():
    """Full-range Dirichlet kernel showing the main lobe and envelope."""
    bins, D_mag = window_spectrum_dense(rectangular, M)

    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    ax.set_facecolor("white")

    plot_range = (bins >= -10) & (bins <= 10)
    ax.plot(bins[plot_range], D_mag[plot_range],
            color="#2563eb", linewidth=1.5, zorder=3)

    env_bins = np.linspace(0.5, 10, 500)                  # avoid zero
    env_omega = 2 * np.pi * env_bins / M                  # angular frequency
    envelope = 1.0 / (M * np.abs(np.sin(env_omega / 2)))  # exact envelope
    ax.plot(env_bins, envelope, "--", color="#dc2626",
            linewidth=1.0, alpha=0.7, label="Envelope: $1/(M\\,|\\sin(\\omega/2)|)$")
    ax.plot(-env_bins, envelope, "--", color="#dc2626",
            linewidth=1.0, alpha=0.7)

    ax.set_xlabel("Frequency (bins)", fontsize=11)
    ax.set_ylabel("Normalized magnitude", fontsize=11)
    ax.set_title(f"Figure B.10 — Dirichlet kernel with envelope ($M = {M}$)",
                 fontsize=12, fontstyle="italic", loc="left")
    ax.set_xlim(-10, 10)
    ax.set_ylim(-0.02, 1.05)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig


# ============================================================
# Figure B.11 — First side-lobe analysis
# ============================================================
def figure_b11_first_sidelobe():
    """Zoom into the first side lobe: true max vs. k=1.5 approximation."""
    bins, D_mag = window_spectrum_dense(rectangular, M)
    peak_bins, peak_values = find_sidelobe_maxima(bins, D_mag, main_lobe_bins=1.2)

    true_bin = peak_bins[0]                               # actual position
    true_val = peak_values[0]                             # actual value
    true_dB = 20 * np.log10(true_val)                     # dB relative to main lobe

    approx_bin = 1.5                                      # midpoint approximation
    approx_omega = 2 * np.pi * approx_bin / M             # angular frequency
    approx_val = dirichlet_kernel(np.array([approx_omega]), M)[0]
    approx_dB = 20 * np.log10(max(approx_val, 1e-20))

    print(f"First side-lobe strength (M={M}):")
    print(f"  True maximum:     bin={true_bin:.3f}, magnitude={true_val:.5f}, {true_dB:.1f} dB")
    print(f"  k=1.5 approx:     bin={approx_bin:.3f}, magnitude={approx_val:.5f}, {approx_dB:.1f} dB")
    print(f"  Textbook:         -13.0 dB (asymptotic, 2/(3π) = {2/(3*np.pi):.5f})")

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.set_facecolor("white")

    zoom = (bins >= 0.8) & (bins <= 3.5)
    ax.plot(bins[zoom], D_mag[zoom], color="#2563eb", linewidth=1.5, zorder=3)

    ax.plot(true_bin, true_val, "o", color="#dc2626", markersize=8, zorder=4)
    ax.annotate(f"True max: ({true_bin:.3f}, {true_val:.5f})\n= {true_dB:.1f} dB",
                xy=(true_bin, true_val),
                xytext=(true_bin + 0.3, true_val + 0.01),
                fontsize=9, color="#dc2626",
                arrowprops=dict(arrowstyle="->", color="#dc2626"))

    ax.plot(approx_bin, approx_val, "s", color="#f59e0b", markersize=8, zorder=4)
    ax.annotate(f"k=1.5 approx: {approx_val:.5f}\n= {approx_dB:.1f} dB",
                xy=(approx_bin, approx_val),
                xytext=(approx_bin + 0.3, approx_val - 0.02),
                fontsize=9, color="#f59e0b",
                arrowprops=dict(arrowstyle="->", color="#f59e0b"))

    ax.axvline(1, color="#94a3b8", linewidth=0.5, linestyle=":", label="Nulls")
    ax.axvline(2, color="#94a3b8", linewidth=0.5, linestyle=":")

    ax.set_xlabel("Frequency (bins)", fontsize=11)
    ax.set_ylabel("Normalized magnitude", fontsize=11)
    ax.set_title(f"Figure B.11 — First side-lobe analysis ($M = {M}$)",
                 fontsize=12, fontstyle="italic", loc="left")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig, {"true_bin": true_bin, "true_val": true_val, "true_dB": true_dB,
                 "approx_val": approx_val, "approx_dB": approx_dB}


# ============================================================
# Figure B.12 — Side-lobe decay regression
# ============================================================
def figure_b12_decay_regression():
    """Collect side-lobe maxima, fit regression, measure decay rate."""
    bins, D_mag = window_spectrum_dense(rectangular, M)
    peak_bins, peak_values = find_sidelobe_maxima(bins, D_mag, main_lobe_bins=1.2)

    n_lobes = min(6, len(peak_bins))                      # use up to 6 side lobes
    pk_b = peak_bins[:n_lobes]                            # bin positions
    pk_v = peak_values[:n_lobes]                          # magnitude values
    pk_dB = 20 * np.log10(pk_v)                           # dB values

    crude_approx_bins = np.arange(1, n_lobes + 1) + 0.5  # midpoint approximation

    print(f"\nSide-lobe decay analysis (M={M}):")
    for i in range(n_lobes):
        env_approx = 2.0 / ((2*(i+1) + 1) * np.pi)       # (2k'+1)/2 rule
        env_dB = 20 * np.log10(env_approx)
        print(f"  Lobe {i+1}: actual bin={pk_b[i]:.3f}, mag={pk_v[i]:.5f} ({pk_dB[i]:.1f} dB)"
              f"  |  approx: {env_approx:.5f} ({env_dB:.1f} dB)")

    log_bins = np.log10(pk_b)                             # log of bin positions
    log_mag = np.log10(pk_v)                              # log of magnitudes
    slope, intercept, r_value, _, _ = linregress(log_bins, log_mag)
    dB_per_octave = slope * 20 * np.log10(2)

    log_crude = np.log10(crude_approx_bins)
    slope_c, _, r_c, _, _ = linregress(log_crude, log_mag)
    dB_oct_crude = slope_c * 20 * np.log10(2)

    print(f"\n  Regression (accurate k'):")
    print(f"    Slope: {slope:.3f}, R²: {r_value**2:.5f}, Decay: {dB_per_octave:.1f} dB/oct")
    print(f"  Regression (crude midpoints):")
    print(f"    Slope: {slope_c:.3f}, R²: {r_c**2:.5f}, Decay: {dB_oct_crude:.1f} dB/oct")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    zoom = (bins >= 0) & (bins <= pk_b[-1] + 2)
    ax1.plot(bins[zoom], D_mag[zoom], color="#2563eb", linewidth=1.0)
    ax1.plot(pk_b, pk_v, "o", color="#dc2626", markersize=6, zorder=4,
             label="Side-lobe maxima")
    for i, (b, v) in enumerate(zip(pk_b, pk_v)):
        ax1.annotate(f"{i+1}", xy=(b, v), xytext=(b, v + 0.005),
                     fontsize=8, ha="center", color="#dc2626")
    ax1.set_xlabel("Frequency (bins)", fontsize=11)
    ax1.set_ylabel("Normalized magnitude", fontsize=11)
    ax1.set_title("Side-lobe maxima", fontsize=10, fontstyle="italic", loc="left")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.15)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2.plot(pk_b, pk_v, "o", color="#dc2626", markersize=6, zorder=4,
             label="Measured maxima")
    fit_x = np.linspace(pk_b[0] * 0.8, pk_b[-1] * 1.2, 100)
    fit_y = 10**(intercept) * fit_x**slope
    ax2.plot(fit_x, fit_y, "--", color="#059669", linewidth=1.5,
             label=f"Fit: slope={slope:.2f}, R²={r_value**2:.4f}\n"
                   f"= {dB_per_octave:.1f} dB/oct")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Bin position (log scale)", fontsize=11)
    ax2.set_ylabel("Magnitude (log scale)", fontsize=11)
    ax2.set_title("Log-log regression → decay rate", fontsize=10,
                  fontstyle="italic", loc="left")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.15, which="both")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle(f"Figure B.12 — Side-lobe decay analysis ($M = {M}$)",
                 fontsize=12, fontstyle="italic")
    fig.set_dpi(DPI)
    fig.tight_layout()

    table_data = []
    for i in range(n_lobes):
        env_approx = 2.0 / ((2*(i+1) + 1) * np.pi)
        table_data.append({
            "lobe": i + 1,
            "bin_actual": pk_b[i], "bin_midpoint": crude_approx_bins[i],
            "mag_actual": pk_v[i], "mag_approx": env_approx,
            "dB_actual": pk_dB[i], "dB_approx": 20 * np.log10(env_approx),
        })

    return fig, {
        "slope": slope, "r_squared": r_value**2,
        "dB_per_octave": dB_per_octave,
        "slope_crude": slope_c, "r_squared_crude": r_c**2,
        "dB_per_octave_crude": dB_oct_crude,
        "table": table_data,
    }


# ============================================================
# Figure B.13 — All window spectra comparison (linear, zoomed)
# ============================================================
def figure_b13_comparison():
    """Overlay normalized |W(ω)| for all four windows, linear scale, zoomed."""
    windows = [
        ("Rectangular", rectangular, "#2563eb", 1),
        ("Hann",        hann,        "#dc2626", 2),
        ("Hamming",     hamming,     "#059669", 2),
        ("Blackman",    blackman,    "#f59e0b", 3),
    ]

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    ax.set_facecolor("white")

    for name, win_func, color, beta in windows:
        bins, W_mag = window_spectrum_dense(win_func, M)
        plot_range = (bins >= -8) & (bins <= 8)
        ax.plot(bins[plot_range], W_mag[plot_range],
                color=color, linewidth=1.5, label=f"{name} (β={beta})", zorder=3)
        ax.axvline(beta, color=color, linewidth=0.5, linestyle=":", alpha=0.5)
        ax.axvline(-beta, color=color, linewidth=0.5, linestyle=":", alpha=0.5)

    ax.set_xlabel("Frequency (bins)", fontsize=11)
    ax.set_ylabel("Normalized magnitude", fontsize=11)
    ax.set_title(f"Figure B.13 — Window spectra comparison, linear scale ($M = {M}$)",
                 fontsize=12, fontstyle="italic", loc="left")
    ax.set_xlim(-8, 8)
    ax.set_ylim(-0.01, 0.25)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig


# ============================================================
# Figure B.14 — All window spectra comparison (linear, full)
# ============================================================
def figure_b14_comparison_full():
    """Same comparison but showing the full main lobe."""
    windows = [
        ("Rectangular", rectangular, "#2563eb", 1),
        ("Hann",        hann,        "#dc2626", 2),
        ("Hamming",     hamming,     "#059669", 2),
        ("Blackman",    blackman,    "#f59e0b", 3),
    ]

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    ax.set_facecolor("white")

    for name, win_func, color, beta in windows:
        bins, W_mag = window_spectrum_dense(win_func, M)
        plot_range = (bins >= -8) & (bins <= 8)
        ax.plot(bins[plot_range], W_mag[plot_range],
                color=color, linewidth=1.5, label=f"{name} (β={beta})", zorder=3)

    ax.set_xlabel("Frequency (bins)", fontsize=11)
    ax.set_ylabel("Normalized magnitude", fontsize=11)
    ax.set_title(f"Figure B.14 — Window spectra comparison, full range ($M = {M}$)",
                 fontsize=12, fontstyle="italic", loc="left")
    ax.set_xlim(-8, 8)
    ax.set_ylim(-0.02, 1.05)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig


# ============================================================
# Run all experiments
# ============================================================
def run_lab3(save=False):
    """Run all Lab 3 experiments and optionally save figures."""
    print("=" * 60)
    print(f"Lab 3: Windowing and the Dirichlet Kernel (M = {M})")
    print("=" * 60)

    fig_anat = figure_b9_dirichlet_anatomy()               # Figure B.9
    fig_full = figure_b10_dirichlet_full()                  # Figure B.10
    fig_sl, sl_data = figure_b11_first_sidelobe()           # Figure B.11
    fig_decay, decay_data = figure_b12_decay_regression()   # Figure B.12
    fig_comp = figure_b13_comparison()                      # Figure B.13
    fig_comp_full = figure_b14_comparison_full()            # Figure B.14

    if save:
        all_figs = [fig_anat, fig_full, fig_sl, fig_decay, fig_comp, fig_comp_full]
        for i, fig in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"                # B.09 through B.14
            save_figure(fig, lab_number=LAB_NUMBER, fig_id=fig_id)

    plt.show()
    return {
        "sidelobe": sl_data,
        "decay": decay_data,
    }


if __name__ == "__main__":
    run_lab3(save=True)
