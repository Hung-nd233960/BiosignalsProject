"""
Window functions — hand-built from the cosine-sum formula (Equation A.16)
plus the Gaussian window (Equation A.22).
"""
import numpy as np


# --- Cosine-sum general form (Equation A.16) ---
def _cosine_sum(N, coeffs):
    """
    General cosine-sum window: w[n] = Σ (-1)^p * a_p * cos(2πpn/N)

    Parameters
    ----------
    N      : int — window length in samples
    coeffs : list of float — coefficients [a_0, a_1, ..., a_P]

    Returns
    -------
    w : ndarray — window samples of length N
    """
    n = np.arange(N)                                      # sample indices [0, N-1]
    w = np.zeros(N)                                       # initialize window
    for p, a_p in enumerate(coeffs):                      # sum each cosine term
        sign = (-1) ** p                                  # alternating sign
        w += sign * a_p * np.cos(2 * np.pi * p * n / N)  # Equation (A.16)
    return w


# --- Rectangular (Equation A.17) ---
def rectangular(N):
    """w[n] = 1 for all n. β = 1, peak side-lobe = −13 dB."""
    return np.ones(N)                                     # all ones


# --- Hann (Equation A.18) ---
def hann(N):
    """w[n] = 0.5 − 0.5·cos(2πn/N). β = 2, peak side-lobe = −31.5 dB."""
    return _cosine_sum(N, [0.5, 0.5])                     # a0=0.5, a1=0.5


# --- Hamming (Equation A.19) ---
def hamming(N):
    """w[n] = 0.54 − 0.46·cos(2πn/N). β = 2, peak side-lobe = −42.7 dB."""
    return _cosine_sum(N, [0.54, 0.46])                   # a0=0.54, a1=0.46


# --- Blackman (Equation A.20) ---
def blackman(N):
    """w[n] = 0.42 − 0.5·cos(2πn/N) + 0.08·cos(4πn/N). β = 3, peak side-lobe = −58 dB."""
    return _cosine_sum(N, [0.42, 0.5, 0.08])              # a0=0.42, a1=0.5, a2=0.08


# --- Gaussian (Equation A.22) ---
def gaussian(N, sigma_ratio=0.4):
    """
    w[n] = exp(−0.5 * ((n − (N−1)/2) / σ)²)
    σ = sigma_ratio * (N−1)/2

    Parameters
    ----------
    N           : int — window length in samples
    sigma_ratio : float — σ as a fraction of the half-length (default 0.4)

    Returns
    -------
    w : ndarray — window samples of length N
    """
    n = np.arange(N)                                      # sample indices
    center = (N - 1) / 2                                  # window center
    sigma = sigma_ratio * center                          # σ in samples
    w = np.exp(-0.5 * ((n - center) / sigma) ** 2)        # Equation (A.22)
    return w


# --- Window metadata for tables and labels ---
WINDOW_REGISTRY = {
    "rectangular": {"func": rectangular, "beta": 1, "peak_sidelobe_dB": -13.0},
    "hann":        {"func": hann,        "beta": 2, "peak_sidelobe_dB": -31.5},
    "hamming":     {"func": hamming,     "beta": 2, "peak_sidelobe_dB": -42.7},
    "blackman":    {"func": blackman,    "beta": 3, "peak_sidelobe_dB": -58.0},
}
