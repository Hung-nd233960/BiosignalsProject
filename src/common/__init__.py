from .config import *
from .signals import (
    make_tone,
    make_mixed_tones,
    make_chirp,
    make_multi_chirp,
    make_transient,
    make_noise,
    make_time_axis,
)
from .windows import rectangular, hann, hamming, blackman, gaussian
from .plotting import (
    plot_time_domain,
    plot_dual_stack_spectrum,
    plot_spectrogram,
    save_figure,
)
from .eeg import load_eeg, compute_band_power, get_channel_data
