"""
Model signal generators — one function per Appendix A archetype.
All signals satisfy the Volume B lab constraints:
  - frequencies below F_MAX (100 Hz)
  - sampling at FS (250 Hz default)
  - duration at least DURATION (1200 s default)
"""
import numpy as np
from .config import FS, DURATION, SEED, F_MAX


def make_time_axis(duration=DURATION, fs=FS):
    """Return the sample indices n and the corresponding time axis t (seconds)."""
    N = int(duration * fs)           # total number of samples
    n = np.arange(N)                 # sample indices [0, 1, ..., N-1]
    t = n / fs                       # time in seconds
    return n, t


# --- AA.1: Single tone (Equation AA.1) ---
def make_tone(f0, A=1.0, phi=0.0, duration=DURATION, fs=FS):
    """
    x[n] = A * cos(2π f0 n / fs + φ)

    Parameters
    ----------
    f0  : float — frequency in Hz (must be < F_MAX)
    A   : float — amplitude
    phi : float — initial phase in radians
    duration : float — signal duration in seconds
    fs  : int — sampling frequency in Hz

    Returns
    -------
    x : ndarray — the signal samples
    n : ndarray — the sample indices
    t : ndarray — the time axis in seconds
    """
    assert f0 < F_MAX, f"f0={f0} Hz exceeds F_MAX={F_MAX} Hz"
    n, t = make_time_axis(duration, fs)    # generate time axis
    x = A * np.cos(2 * np.pi * f0 * n / fs + phi)  # Equation (AA.1)
    return x, n, t


# --- AA.2: Mixed tones (Equation AA.2) ---
def make_mixed_tones(freqs, amplitudes=None, phases=None, duration=DURATION, fs=FS):
    """
    x[n] = Σ A_i * cos(2π f_i n / fs + φ_i)

    Parameters
    ----------
    freqs      : list of float — frequencies in Hz (each must be < F_MAX)
    amplitudes : list of float — amplitudes (default: all 1.0)
    phases     : list of float — initial phases in radians (default: all 0.0)
    duration   : float — signal duration in seconds
    fs         : int — sampling frequency in Hz

    Returns
    -------
    x : ndarray — the signal samples
    n : ndarray — the sample indices
    t : ndarray — the time axis in seconds
    """
    K = len(freqs)                                        # number of components
    if amplitudes is None:
        amplitudes = [1.0] * K                            # default: unit amplitude
    if phases is None:
        phases = [0.0] * K                                # default: zero phase
    for f in freqs:
        assert f < F_MAX, f"f={f} Hz exceeds F_MAX={F_MAX} Hz"
    n, t = make_time_axis(duration, fs)                   # generate time axis
    x = np.zeros(len(n))                                  # initialize signal
    for i in range(K):                                    # sum each component
        x += amplitudes[i] * np.cos(2 * np.pi * freqs[i] * n / fs + phases[i])
    return x, n, t


# --- AA.3: Chirp (Equations AA.3, AA.4) ---
def make_chirp(f0, mu, A=1.0, phi=0.0, duration=DURATION, fs=FS):
    """
    x[n] = A * cos(2π/fs * (f0*n + μ/2 * n²/fs) + φ)
    Instantaneous frequency: f_inst[n] = f0 + μ * n / fs

    Parameters
    ----------
    f0  : float — starting frequency in Hz
    mu  : float — chirp rate in Hz/s
    A   : float — amplitude
    phi : float — initial phase in radians
    duration : float — signal duration in seconds
    fs  : int — sampling frequency in Hz

    Returns
    -------
    x : ndarray — the signal samples
    n : ndarray — the sample indices
    t : ndarray — the time axis in seconds
    """
    n, t = make_time_axis(duration, fs)                   # generate time axis
    f_end = f0 + mu * duration                            # final frequency
    assert f0 < F_MAX, f"f0={f0} Hz exceeds F_MAX={F_MAX} Hz"
    assert f_end < F_MAX, f"f_end={f_end:.1f} Hz exceeds F_MAX={F_MAX} Hz"
    phase = 2 * np.pi / fs * (f0 * n + mu / 2 * n**2 / fs) + phi  # Equation (AA.3)
    x = A * np.cos(phase)                                 # generate signal
    return x, n, t


# --- AA.4: Multi-chirp (Equation AA.5) ---
def make_multi_chirp(f0s, mus, amplitudes=None, phases=None, duration=DURATION, fs=FS):
    """
    x[n] = Σ A_i * cos(2π/fs * (f0_i*n + μ_i/2 * n²/fs) + φ_i)

    Parameters
    ----------
    f0s        : list of float — starting frequencies in Hz
    mus        : list of float — chirp rates in Hz/s
    amplitudes : list of float — amplitudes (default: all 1.0)
    phases     : list of float — initial phases in radians (default: all 0.0)
    duration   : float — signal duration in seconds
    fs         : int — sampling frequency in Hz

    Returns
    -------
    x : ndarray — the signal samples
    n : ndarray — the sample indices
    t : ndarray — the time axis in seconds
    """
    K = len(f0s)                                          # number of chirp components
    if amplitudes is None:
        amplitudes = [1.0] * K                            # default: unit amplitude
    if phases is None:
        phases = [0.0] * K                                # default: zero phase
    n, t = make_time_axis(duration, fs)                   # generate time axis
    x = np.zeros(len(n))                                  # initialize signal
    for i in range(K):                                    # sum each chirp component
        f_end = f0s[i] + mus[i] * duration                # check endpoint
        assert f0s[i] < F_MAX, f"f0={f0s[i]} Hz exceeds F_MAX={F_MAX} Hz"
        assert f_end < F_MAX, f"f_end={f_end:.1f} Hz exceeds F_MAX={F_MAX} Hz"
        phase = 2 * np.pi / fs * (f0s[i] * n + mus[i] / 2 * n**2 / fs) + phases[i]
        x += amplitudes[i] * np.cos(phase)                # accumulate
    return x, n, t


# --- AA.5: Transient / pulse (Equations AA.6, AA.7) ---
def make_transient(n0, sigma_t, f0=0.0, A=1.0, phi=0.0, duration=DURATION, fs=FS):
    """
    Gaussian-enveloped burst:
      x[n] = A * exp(-(n - n0)² / (2 σ_t²)) * cos(2π f0 n / fs + φ)
    If f0 = 0, produces a pure Gaussian pulse (no carrier).

    Parameters
    ----------
    n0      : int — center sample index of the transient
    sigma_t : float — width of the Gaussian envelope in samples
    f0      : float — carrier frequency in Hz (0 for baseband pulse)
    A       : float — peak amplitude
    phi     : float — carrier phase in radians
    duration : float — signal duration in seconds
    fs      : int — sampling frequency in Hz

    Returns
    -------
    x : ndarray — the signal samples
    n : ndarray — the sample indices
    t : ndarray — the time axis in seconds
    """
    if f0 > 0:
        assert f0 < F_MAX, f"f0={f0} Hz exceeds F_MAX={F_MAX} Hz"
    n, t = make_time_axis(duration, fs)                   # generate time axis
    envelope = A * np.exp(-((n - n0) ** 2) / (2 * sigma_t ** 2))  # Gaussian envelope
    if f0 > 0:
        carrier = np.cos(2 * np.pi * f0 * n / fs + phi)  # carrier oscillation
    else:
        carrier = np.ones(len(n))                         # no carrier (baseband)
    x = envelope * carrier                                # Equation (AA.6)
    return x, n, t


# --- AA.6: Noise (Equation AA.8) ---
def make_noise(sigma=1.0, duration=DURATION, fs=FS, seed=SEED):
    """
    x[n] = η[n],  η ~ N(0, σ²), i.i.d.

    Parameters
    ----------
    sigma    : float — standard deviation of the noise
    duration : float — signal duration in seconds
    fs       : int — sampling frequency in Hz
    seed     : int — random seed for reproducibility

    Returns
    -------
    x : ndarray — the noise samples
    n : ndarray — the sample indices
    t : ndarray — the time axis in seconds
    """
    rng = np.random.default_rng(seed)                     # seeded RNG for reproducibility
    n, t = make_time_axis(duration, fs)                   # generate time axis
    x = rng.normal(loc=0.0, scale=sigma, size=len(n))     # Equation (AA.8)
    return x, n, t
