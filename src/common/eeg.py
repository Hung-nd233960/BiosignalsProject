"""
EEG utilities for Volume C — loading, preprocessing, and band-power analysis.
All EEG I/O goes through this module. Never call MNE directly in lab code.
"""

import numpy as np
from scipy import signal as sp_signal
from .config import FS, EEG_BANDS, DATA_DIR, SEED
from .windows import hann


def load_eeg(filepath, channels=None, tmin=None, tmax=None, preload=True):
    """
    Load an EEG recording via MNE and return raw data in µV.

    Parameters
    ----------
    filepath  : str — path to the EEG file (any format MNE supports)
    channels  : list of str — channel names to select (None = all)
    tmin      : float — start time in seconds (None = beginning)
    tmax      : float — end time in seconds (None = end)
    preload   : bool — load data into memory

    Returns
    -------
    data : ndarray — shape (n_channels, n_samples), in µV
    ch_names : list of str — channel names
    fs : float — sampling frequency in Hz
    times : ndarray — time axis in seconds
    """
    import mne                                            # import here to keep MNE optional for Volume B

    raw = mne.io.read_raw(filepath, preload=preload)      # load via MNE
    if channels is not None:
        raw = raw.pick(channels)                          # select channels

    fs = raw.info["sfreq"]                                # sampling frequency
    data = raw.get_data(tmin=tmin, tmax=tmax)             # shape: (n_channels, n_samples)
    data = data * 1e6                                     # V → µV

    n_samples = data.shape[1]                             # number of samples
    t_start = tmin if tmin is not None else 0.0           # start time
    times = t_start + np.arange(n_samples) / fs           # time axis in seconds

    ch_names = raw.ch_names                               # channel names

    return data, ch_names, fs, times


def compute_band_power(x, fs=FS, nperseg=None, window="hann", overlap_frac=0.5,
                       bands=None):
    """
    Compute power in each EEG band via Welch's method.

    Parameters
    ----------
    x             : ndarray — 1D signal in µV
    fs            : float — sampling frequency in Hz
    nperseg       : int — segment length in samples (default: 5 s worth)
    window        : str — window type for Welch
    overlap_frac  : float — overlap as fraction of nperseg (default: 0.5 = 50%)
    bands         : dict — band definitions {name: (f_low, f_high)} (default: EEG_BANDS)

    Returns
    -------
    band_power : dict — {band_name: power_in_µV²} for each band
    freqs      : ndarray — frequency axis from Welch (Hz)
    psd        : ndarray — full PSD in µV²/Hz
    """
    if bands is None:
        bands = EEG_BANDS                                 # default: standard EEG bands
    if nperseg is None:
        nperseg = int(5.0 * fs)                           # default: 5 s segments (Δf = 0.2 Hz)

    noverlap = int(nperseg * overlap_frac)                # overlap in samples

    freqs, psd = sp_signal.welch(
        x, fs=fs,
        nperseg=nperseg,
        noverlap=noverlap,
        window=window,
    )                                                     # PSD in µV²/Hz (if input is µV)

    band_power = {}                                       # power per band
    for name, (f_low, f_high) in bands.items():
        band_mask = (freqs >= f_low) & (freqs < f_high)   # bins in this band
        df = freqs[1] - freqs[0]                          # frequency resolution
        band_power[name] = np.sum(psd[band_mask]) * df    # integrate PSD → µV²

    return band_power, freqs, psd


def get_channel_data(data, ch_names, channel):
    """
    Extract a single channel from the multi-channel data array.

    Parameters
    ----------
    data     : ndarray — shape (n_channels, n_samples)
    ch_names : list of str — channel names
    channel  : str — name of the channel to extract

    Returns
    -------
    x : ndarray — 1D signal for the requested channel (µV)
    """
    idx = ch_names.index(channel)                         # find channel index
    return data[idx]                                      # return 1D array
