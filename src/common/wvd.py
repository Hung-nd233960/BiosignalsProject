"""
Wigner-Ville Distribution (WVD) and Smoothed Pseudo Wigner-Ville Distribution (SPWVD)
implementations from scratch in Python, optimized using NumPy and SciPy.
Includes interpolation-by-2 to resolve the half-sample time-frequency aliasing.
"""

import numpy as np
from scipy.signal import hilbert, resample
from scipy.signal.windows import hann

def wigner_ville(x, fs, n_fft=None):
    """
    Compute the discrete Wigner-Ville Distribution (WVD) of a signal.
    Uses the analytic signal (via Hilbert transform) to suppress the DC/Nyquist self-ghosts,
    and interpolates the signal by a factor of 2 to prevent frequency aliasing.

    Parameters
    ----------
    x      : ndarray — 1D input signal
    fs     : float — sampling frequency in Hz
    n_fft  : int — number of frequency bins (default: length of x)

    Returns
    -------
    wvd    : ndarray — WVD matrix of shape (n_fft // 2, N), real-valued
    t      : ndarray — time axis in seconds
    freqs  : ndarray — frequency axis in Hz (0 to fs/2)
    """
    # 1. Convert to analytic signal to eliminate negative frequencies
    z = hilbert(x)
    N = len(z)
    if n_fft is None:
        n_fft = N
        
    # 2. Interpolate by 2 using FFT resample (zero-padding in frequency domain)
    # This doubles the sampling rate to 2*fs, which makes the Nyquist frequency 
    # of the lag autocorrelation (spacing 2m) equal to fs/2.
    z_interp = resample(z, 2 * N)
    
    wvd_interp = np.zeros((n_fft, 2 * N))
    
    # 3. Compute WVD on the interpolated signal
    for n in range(2 * N):
        L = min(n, 2 * N - 1 - n)
        # Limit lag to n_fft // 2 - 1
        L = min(L, n_fft // 2 - 1)
        
        r = np.zeros(n_fft, dtype=complex)
        r[0] = z_interp[n] * np.conj(z_interp[n])
        for m in range(1, L + 1):
            val = z_interp[n + m] * np.conj(z_interp[n - m])
            r[m] = val
            r[n_fft - m] = np.conj(val)
            
        wvd_interp[:, n] = 2.0 * np.real(np.fft.fft(r))
        
    # 4. Decimate the time axis by 2 to return to the original time resolution
    wvd_dec = wvd_interp[:, ::2]
    
    # The first half of the FFT bins covers 0 to fs/2
    freqs = np.linspace(0, fs / 2, n_fft // 2)
    t = np.arange(N) / fs
    
    return wvd_dec[:n_fft // 2, :], t, freqs


def smoothed_pseudo_wigner_ville(x, fs, h=None, g=None, n_fft=None):
    """
    Compute the Smoothed Pseudo Wigner-Ville Distribution (SPWVD) of a signal.
    Applies independent time smoothing (window g) and frequency/lag smoothing (window h).
    Uses the analytic signal and interpolates by 2 to prevent aliasing.

    Parameters
    ----------
    x      : ndarray — 1D input signal
    fs     : float — sampling frequency in Hz
    h      : ndarray or int — frequency smoothing (lag) window or its length (defaults to Hann window, length 51)
    g      : ndarray or int — time smoothing window or its length (defaults to Hann window, length 11)
    n_fft  : int — number of frequency bins (default: max(256, 2 * len(h) - 1))

    Returns
    -------
    spwvd  : ndarray — SPWVD matrix of shape (n_fft // 2, N), real-valued
    t      : ndarray — time axis in seconds
    freqs  : ndarray — frequency axis in Hz (0 to fs/2)
    """
    # 1. Convert to analytic signal
    z = hilbert(x)
    N = len(z)
    
    # Resolve window h (lag window)
    if h is None:
        h_len_orig = 51
        h_orig = hann(h_len_orig, sym=True)
    elif isinstance(h, int):
        h_len_orig = h
        h_orig = hann(h_len_orig, sym=True)
    else:
        h_orig = h
        h_len_orig = len(h)
        
    # Resolve window g (time window)
    if g is None:
        g_len_orig = 11
        g_orig = hann(g_len_orig, sym=True)
    elif isinstance(g, int):
        g_len_orig = g
        g_orig = hann(g_len_orig, sym=True)
    else:
        g_orig = g
        g_len_orig = len(g)

    # 2. Interpolate z by 2
    z_interp = resample(z, 2 * N)
    
    # 3. Create windows on the resampled grid (which has double the samples)
    if len(g_orig) > 1:
        # Scale time window g by factor of 2 in length (odd length)
        g_len_resampled = 2 * g_len_orig - 1
        g_resampled = resample(g_orig, g_len_resampled)
        # Normalize resampled time window so it sums to 1
        g_resampled = g_resampled / np.sum(g_resampled)
    else:
        # No time smoothing (PWVD case)
        g_resampled = g_orig
        g_len_resampled = len(g_resampled)
        
    h_len_resampled = 2 * h_len_orig - 1
    h_resampled = resample(h_orig, h_len_resampled)
    
    M = (h_len_resampled - 1) // 2  # half-length of resampled lag window
    P = (g_len_resampled - 1) // 2  # half-length of resampled time window
    
    if n_fft is None:
        n_fft = max(h_len_resampled, 256)
        
    spwvd_interp = np.zeros((n_fft, 2 * N))
    
    # Pad interpolated signal to handle boundary conditions easily
    pad_len = M + P
    z_pad = np.pad(z_interp, pad_len, mode='constant', constant_values=0)
    
    # 4. Compute SPWVD on resampled signal
    for n in range(2 * N):
        n_idx = n + pad_len
        r = np.zeros(n_fft, dtype=complex)
        
        # Lag m = 0
        sum_val = 0.0
        for p in range(-P, P + 1):
            sum_val += g_resampled[P + p] * z_pad[n_idx + p] * np.conj(z_pad[n_idx + p])
        r[0] = h_resampled[M] * sum_val
        
        # Lags m = 1 ... M
        # Limit lag to n_fft // 2 - 1
        L_max = min(M, n_fft // 2 - 1)
        for m in range(1, L_max + 1):
            sum_val = 0.0
            for p in range(-P, P + 1):
                sum_val += g_resampled[P + p] * z_pad[n_idx + p + m] * np.conj(z_pad[n_idx + p - m])
            val = h_resampled[M + m] * sum_val
            r[m] = val
            r[n_fft - m] = np.conj(val)
            
        spwvd_interp[:, n] = 2.0 * np.real(np.fft.fft(r))
        
    # 5. Decimate the time axis by 2 to return to the original time resolution
    spwvd_dec = spwvd_interp[:, ::2]
    
    freqs = np.linspace(0, fs / 2, n_fft // 2)
    t = np.arange(N) / fs
    
    return spwvd_dec[:n_fft // 2, :], t, freqs
