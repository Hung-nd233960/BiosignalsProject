"""
Lab 4: The STFT of a Fluctuating Signal (↔ A.5)
-------------------------------------------------
Tests: Heisenberg uncertainty tradeoff, hop size / overlap / COLA,
       multi-scale limitation of the STFT.

Figures B.15–B.22 in the report.
"""

import numpy as np  # numerical computing
import matplotlib.pyplot as plt  # plotting
from scipy import signal as sp_signal  # STFT via spectrogram
import sys
import os  # path manipulation

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.common import (  # project utilities
    FS,
    DPI,
    COLORMAP,
    make_chirp,
    make_transient,
    plot_time_domain,
    save_figure,
)

# ============================================================
# Parameters (Table B.7 in the report)
# ============================================================
F0_CHIRP = 5.0  # chirp start frequency (Hz)
F1_CHIRP = 45.0  # chirp end frequency (Hz)
DURATION_LAB = 120.0  # 120 s — shorter so chirp diagonal is visible at all zooms
MU_CHIRP = (F1_CHIRP - F0_CHIRP) / DURATION_LAB  # 0.333 Hz/s
LAB_NUMBER = 4  # lab number for figure saving
FIG_START = 15  # first figure number (B.15)

# Window lengths for Heisenberg sweep (Experiment A)
WINDOW_LENGTHS = [125, 250, 500, 1250]  # samples (0.5, 1, 2, 5 seconds)

# Overlap settings for Experiment B
OVERLAP_FRACTIONS = [0.0, 0.5, 0.75]  # 0%, 50%, 75%

# Multi-scale signal (Experiment C)
BURST_CENTER_S = 60.0  # burst at midpoint (60 s)
BURST_SIGMA_S = 0.5  # burst width: 0.5 s (short transient)
BURST_FREQ = 10.0  # burst carrier: alpha band (10 Hz)
BURST_AMP = 3.0  # burst amplitude (larger than chirp for visibility)


# ============================================================
# 1. Helper: compute and plot spectrogram (dual-stack)
# ============================================================
def compute_stft(x, fs=FS, nperseg=256, noverlap=None, window="hann"):
    """
    Compute the STFT spectrogram using scipy.signal.spectrogram.

    Parameters
    ----------
    x        : ndarray — input signal
    fs       : int — sampling frequency (Hz)
    nperseg  : int — window length M (samples)
    noverlap : int — overlap in samples (default: nperseg // 2)
    window   : str — window type

    Returns
    -------
    t_stft : ndarray — time axis (s)
    f_stft : ndarray — frequency axis (Hz)
    Sxx    : 2D ndarray — power spectrogram (freq × time)
    """
    if noverlap is None:
        noverlap = nperseg // 2  # default 50% overlap
    f_stft, t_stft, Sxx = sp_signal.spectrogram(
        x,
        fs=fs,
        nperseg=nperseg,  # window length
        noverlap=noverlap,  # overlap
        window=window,  # window type
    )
    return t_stft, f_stft, Sxx


def plot_spectrogram_dualstack(
    t_stft,
    f_stft,
    Sxx,
    title=None,
    fig_id=None,
    f_range=(0, 30),
    t_range=None,
    cmap=COLORMAP,
):
    """
    Dual-stack spectrogram: linear on top, dB on bottom.

    Parameters
    ----------
    t_stft  : ndarray — time axis (s)
    f_stft  : ndarray — frequency axis (Hz)
    Sxx     : 2D ndarray — power spectrogram
    title   : str — overall title
    fig_id  : str — figure number label
    f_range : tuple — frequency zoom (Hz)
    t_range : tuple — time zoom (s)
    cmap    : str — colormap

    Returns
    -------
    fig, (ax_lin, ax_db) : matplotlib Figure and Axes
    """
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(14, 7), sharex=True, sharey=True)

    # --- Top: linear scale ---
    im_lin = ax_lin.pcolormesh(t_stft, f_stft, Sxx, shading="gouraud", cmap=cmap)
    ax_lin.set_ylabel("Frequency (Hz)")
    fig.colorbar(im_lin, ax=ax_lin, label="Power (linear)")
    if title:
        ax_lin.set_title(title)
    if fig_id:
        ax_lin.set_title(fig_id, loc="left", fontsize=9, fontstyle="italic")

    # --- Bottom: dB scale ---
    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-20))
    im_db = ax_db.pcolormesh(t_stft, f_stft, Sxx_db, shading="gouraud", cmap=cmap)
    ax_db.set_ylabel("Frequency (Hz)")
    ax_db.set_xlabel("Time (s)")
    fig.colorbar(im_db, ax=ax_db, label="Power (dB)")

    if f_range:
        ax_lin.set_ylim(f_range)
    if t_range:
        ax_lin.set_xlim(t_range)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig, (ax_lin, ax_db)


# ============================================================
# 2. Experiment A: Heisenberg tradeoff — window length sweep
# ============================================================
def experiment_heisenberg():
    """
    Same chirp, multiple window lengths.
    Demonstrates: short window = good time / poor freq, long = opposite.
    Δt·Δf = β is constant.
    """
    print("\nExperiment A: Heisenberg tradeoff")
    print("-" * 40)

    # --- Generate chirp ---
    x, n, t = make_chirp(F0_CHIRP, MU_CHIRP, A=1.0, duration=DURATION_LAB)

    # --- Time-domain plot (always first) ---
    fig_td, _ = plot_time_domain(
        t,
        x,
        title=f"Linear chirp {F0_CHIRP}→{F1_CHIRP} Hz over {DURATION_LAB:.0f} s (time domain)",
        fig_id=f"Figure B.{FIG_START}",
        t_range=(0, 5),
    )

    # --- Spectrograms at multiple window lengths ---
    figs_stft = []
    beta = 2  # Hann window β
    for i, M in enumerate(WINDOW_LENGTHS):
        dt = M / FS  # time resolution (s)
        df = beta * FS / M  # freq resolution (Hz)
        product = dt * df  # should equal β = 2
        noverlap = M // 2  # 50% overlap (COLA for Hann)
        hop = M - noverlap  # hop size

        print(
            f"  M={M} ({dt:.2f} s): Δt={dt:.2f} s, Δf={df:.2f} Hz, "
            f"Δt·Δf={product:.1f}, hop={hop}"
        )

        t_stft, f_stft, Sxx = compute_stft(x, nperseg=M, noverlap=noverlap)

        fig_id = f"Figure B.{FIG_START + 1 + i}"
        fig, _ = plot_spectrogram_dualstack(
            t_stft,
            f_stft,
            Sxx,
            title=f"Chirp spectrogram — M={M} ({dt:.1f} s), "
            f"Δt={dt:.2f} s, Δf={df:.2f} Hz",
            fig_id=fig_id,
            f_range=(0, F1_CHIRP + 5),  # show full chirp range + margin
            t_range=(t_stft[0], t_stft[-1]),  # trim to actual STFT range
        )
        figs_stft.append(fig)

    return fig_td, figs_stft


# ============================================================
# 3. Experiment B: Overlap and the tapering problem
# ============================================================
def experiment_overlap():
    """
    Same chirp + a transient at a segment boundary.
    Shows: 0% overlap loses edge samples, 50% recovers them, 75% interpolates.
    """
    print("\nExperiment B: Overlap and tapering")
    print("-" * 40)

    # --- Generate chirp + transient ---
    x_chirp, n, t = make_chirp(F0_CHIRP, MU_CHIRP, A=1.0, duration=DURATION_LAB)
    n0 = int(BURST_CENTER_S * FS)  # center sample
    sigma_t = int(BURST_SIGMA_S * FS)  # width in samples
    x_burst, _, _ = make_transient(
        n0, sigma_t, f0=BURST_FREQ, A=BURST_AMP, duration=DURATION_LAB
    )
    x = x_chirp + x_burst  # combined signal

    # --- Burst true extent for reference lines ---
    burst_start = BURST_CENTER_S - 2 * BURST_SIGMA_S  # 59 s
    burst_end = BURST_CENTER_S + 2 * BURST_SIGMA_S  # 61 s

    M = 256  # fixed window length
    figs = []

    for overlap_frac in OVERLAP_FRACTIONS:
        noverlap = int(M * overlap_frac)  # overlap in samples
        hop = M - noverlap  # hop size
        n_cols = (len(x) - M) // hop + 1  # number of STFT columns

        print(
            f"  Overlap={overlap_frac*100:.0f}%: hop={hop}, "
            f"columns={n_cols}, segments/sample={'%.1f' % (M / hop if hop > 0 else 0)}"
        )

        t_stft, f_stft, Sxx = compute_stft(x, nperseg=M, noverlap=noverlap)

        fig_idx = len(figs)
        fig_id = f"Figure B.{FIG_START + 5 + fig_idx}"
        fig, (ax_lin, ax_db) = plot_spectrogram_dualstack(
            t_stft,
            f_stft,
            Sxx,
            title=f"Overlap={overlap_frac*100:.0f}% (hop={hop}), " f"M={M}",
            fig_id=fig_id,
            f_range=(0, 30),
            t_range=(55, 65),  # zoom near the burst
        )
        for ax in (ax_lin, ax_db):  # add burst extent lines
            ax.axvline(burst_start, color="white", ls="--", lw=1, alpha=0.8)
            ax.axvline(burst_end, color="white", ls="--", lw=1, alpha=0.8)
        ax_lin.annotate(
            "← true burst extent →",
            xy=((burst_start + burst_end) / 2, 28),
            ha="center",
            fontsize=8,
            color="white",
        )
        figs.append(fig)

    return figs


# ============================================================
# 4. Experiment C: Multi-scale limitation
# ============================================================
def experiment_multiscale():
    """
    Chirp + short alpha burst. No single window captures both:
    short M sees burst, smears chirp; long M sharpens chirp, smears burst.
    """
    print("\nExperiment C: Multi-scale limitation")
    print("-" * 40)

    # --- Generate chirp + alpha burst ---
    x_chirp, n, t = make_chirp(F0_CHIRP, MU_CHIRP, A=1.0, duration=DURATION_LAB)
    n0 = int(BURST_CENTER_S * FS)  # burst center
    sigma_t = int(BURST_SIGMA_S * FS)  # burst width
    x_burst, _, _ = make_transient(
        n0, sigma_t, f0=BURST_FREQ, A=BURST_AMP, duration=DURATION_LAB
    )
    x = x_chirp + x_burst  # combined signal

    # --- Time-domain plot zoomed to burst region ---
    fig_td, _ = plot_time_domain(
        t,
        x,
        title="Chirp + alpha burst (time domain, zoomed)",
        fig_id=f"Figure B.{FIG_START + 8}",
        t_range=(58, 62),
    )

    # --- Burst true extent: ±2σ around center ---
    burst_start = BURST_CENTER_S - 2 * BURST_SIGMA_S  # 59 s
    burst_end = BURST_CENTER_S + 2 * BURST_SIGMA_S  # 61 s

    # --- Short window: sees burst, smears chirp ---
    M_short = 125  # 0.5 s
    t_s, f_s, Sxx_s = compute_stft(x, nperseg=M_short, noverlap=M_short // 2)
    fig_short, (ax_lin_s, ax_db_s) = plot_spectrogram_dualstack(
        t_s,
        f_s,
        Sxx_s,
        title=f"Short window M={M_short} ({M_short/FS:.1f} s) — "
        f"burst visible, chirp smeared",
        fig_id=f"Figure B.{FIG_START + 9}",
        f_range=(0, 30),
        t_range=(55, 65),
    )
    for ax in (ax_lin_s, ax_db_s):  # add burst extent lines
        ax.axvline(burst_start, color="white", ls="--", lw=1, alpha=0.8)
        ax.axvline(burst_end, color="white", ls="--", lw=1, alpha=0.8)
    ax_lin_s.annotate(
        "← true burst extent →",  # label
        xy=((burst_start + burst_end) / 2, 28),
        ha="center",
        fontsize=8,
        color="white",
    )

    # --- Long window: sharpens chirp, smears burst ---
    M_long = 1250  # 5 s
    t_l, f_l, Sxx_l = compute_stft(x, nperseg=M_long, noverlap=M_long // 2)
    fig_long, (ax_lin_l, ax_db_l) = plot_spectrogram_dualstack(
        t_l,
        f_l,
        Sxx_l,
        title=f"Long window M={M_long} ({M_long/FS:.1f} s) — "
        f"chirp sharp, burst smeared",
        fig_id=f"Figure B.{FIG_START + 10}",
        f_range=(0, 30),
        t_range=(55, 65),
    )
    for ax in (ax_lin_l, ax_db_l):  # add burst extent lines
        ax.axvline(burst_start, color="white", ls="--", lw=1, alpha=0.8)
        ax.axvline(burst_end, color="white", ls="--", lw=1, alpha=0.8)
    ax_lin_l.annotate(
        "← true burst extent →",
        xy=((burst_start + burst_end) / 2, 28),
        ha="center",
        fontsize=8,
        color="white",
    )

    dt_short = M_short / FS  # time resolution
    df_short = 2 * FS / M_short  # freq resolution (β=2)
    dt_long = M_long / FS
    df_long = 2 * FS / M_long

    print(
        f"  Short M={M_short}: Δt={dt_short:.2f} s, Δf={df_short:.2f} Hz, "
        f"Δt·Δf={dt_short*df_short:.1f}"
    )
    print(
        f"  Long  M={M_long}: Δt={dt_long:.2f} s, Δf={df_long:.2f} Hz, "
        f"Δt·Δf={dt_long*df_long:.1f}"
    )
    print(
        "  Neither captures both: burst needs Δt≤1 s, "
        "chirp needs Δf≤1 Hz → impossible with Δt·Δf=2"
    )

    return fig_td, fig_short, fig_long


# ============================================================
# 5. Run all experiments
# ============================================================
def run_lab4(save=False):
    """Run all Lab 4 experiments and optionally save figures."""
    print("=" * 60)
    print("Lab 4: The STFT of a Fluctuating Signal")
    print("=" * 60)

    fig_td_a, figs_heisenberg = experiment_heisenberg()  # Experiment A
    figs_overlap = experiment_overlap()  # Experiment B
    fig_td_c, fig_short, fig_long = experiment_multiscale()  # Experiment C

    # --- Tag each figure: (fig, raster_only?) ---
    all_figs = [
        (fig_td_a, False),  # B.15: chirp time domain (line plot)
        *[(f, True) for f in figs_heisenberg],  # B.16–B.19: spectrograms (raster only)
        *[(f, True) for f in figs_overlap],  # B.20–B.22: spectrograms (raster only)
        (fig_td_c, False),  # B.23: multi-scale time domain (line plot)
        (fig_short, True),  # B.24: short window spectrogram
        (fig_long, True),  # B.25: long window spectrogram
    ]

    if save:
        for i, (fig, raster_only) in enumerate(all_figs):
            fig_id = f"{FIG_START + i:02d}"
            save_figure(
                fig, lab_number=LAB_NUMBER, fig_id=fig_id, raster_only=raster_only
            )

    plt.show()
    return [fig for fig, _ in all_figs]


if __name__ == "__main__":
    run_lab4(save=True)
