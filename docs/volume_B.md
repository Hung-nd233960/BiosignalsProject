# Volume B — Application and Derivation (Labs)

> Each lab derives the mathematics, implements it in code, and verifies the theory graphically.
> Labs are paired to the theory sections they test.
> All model signals satisfy the EEG-realism constraints: frequencies below 100 Hz, duration ≥ 1200 s, sampling at 250 Hz.

---

## Basis Functions and Infrastructure

Before any lab, we establish the shared code that every experiment imports. This infrastructure enforces the project standards from CLAUDE.md — EEG-realism constraints, reproducibility, plotting conventions — in one place, so no lab needs to redefine them.

All code runs in the `biosignals` conda environment (Python ≥ 3.11, NumPy, SciPy, Matplotlib). The environment is defined in `environment.yml` at the project root.

### Constants (`src/common/config.py`)

Every lab imports its parameters from this file. No magic numbers in lab code.

```python
import os

# --- Sampling and signal parameters (Volume B lab constraints) ---
FS = 250              # sampling frequency (Hz) — default for all labs
DURATION = 1200       # minimum signal duration (s) — 20 minutes
N_SAMPLES = FS * DURATION  # total samples at default fs and duration
F_MAX = 100           # maximum allowed signal frequency (Hz)

# --- Reproducibility ---
SEED = 42             # fixed random seed for all noise generation

# --- Figure export ---
DPI = 300                         # minimum render DPI
FIGURE_FORMAT = "png"             # export format: "png", "pdf", or "svg"
COLORMAP = "viridis"              # default perceptually uniform colormap
COLORMAP_DIVERGING = "inferno"    # alternative for diverging data
FORBIDDEN_COLORMAPS = {"jet", "rainbow", "hsv"}  # never use these

# --- Paths (relative to project root) ---
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GRAPH_DIR = os.path.join(PROJECT_ROOT, "results", "graphs")
```

`FS`, `DURATION`, and `F_MAX` enforce the EEG-realism constraints from the Volume B header. `SEED` ensures all noise is reproducible. `DPI`, `COLORMAP`, and `FORBIDDEN_COLORMAPS` enforce the CLAUDE.md graph standards.

### Signal Generators (`src/common/signals.py`)

One function per Appendix A archetype. Each generator enforces the frequency constraint (`assert f0 < F_MAX`) and returns three arrays: the signal `x`, sample indices `n`, and time axis `t`.

**Time axis:**

```python
import numpy as np
from .config import FS, DURATION, SEED, F_MAX

def make_time_axis(duration=DURATION, fs=FS):
    N = int(duration * fs)           # total number of samples
    n = np.arange(N)                 # sample indices [0, 1, ..., N-1]
    t = n / fs                       # time in seconds
    return n, t
```

**Single tone** (Equation (AA.1)): $x[n] = A \cos(2\pi f_0 n / f_s + \phi)$

```python
def make_tone(f0, A=1.0, phi=0.0, duration=DURATION, fs=FS):
    assert f0 < F_MAX, f"f0={f0} Hz exceeds F_MAX={F_MAX} Hz"
    n, t = make_time_axis(duration, fs)    # generate time axis
    x = A * np.cos(2 * np.pi * f0 * n / fs + phi)  # Equation (AA.1)
    return x, n, t
```

**Mixed tones** (Equation (AA.2)): $x[n] = \sum A_i \cos(2\pi f_i n / f_s + \phi_i)$

```python
def make_mixed_tones(freqs, amplitudes=None, phases=None, duration=DURATION, fs=FS):
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
```

**Noise** (Equation (AA.8)): $x[n] = \eta[n], \quad \eta \sim \mathcal{N}(0, \sigma^2)$

```python
def make_noise(sigma=1.0, duration=DURATION, fs=FS, seed=SEED):
    rng = np.random.default_rng(seed)                     # seeded RNG for reproducibility
    n, t = make_time_axis(duration, fs)                   # generate time axis
    x = rng.normal(loc=0.0, scale=sigma, size=len(n))     # Equation (AA.8)
    return x, n, t
```

The remaining generators (`make_chirp`, `make_multi_chirp`, `make_transient`) follow the same pattern, implementing Equations (AA.3)–(AA.7). Full source: `src/common/signals.py`.

### Window Functions (`src/common/windows.py`)

Hand-built from the cosine-sum formula (Equation (A.16)), not from library calls:

$$
w[n] = \sum_{p=0}^{P} (-1)^p \, a_p \cos\!\left(\frac{2\pi p n}{N}\right) \tag{A.16}
$$

```python
import numpy as np

def _cosine_sum(N, coeffs):
    """General cosine-sum window: w[n] = Σ (-1)^p * a_p * cos(2πpn/N)"""
    n = np.arange(N)                                      # sample indices [0, N-1]
    w = np.zeros(N)                                       # initialize window
    for p, a_p in enumerate(coeffs):                      # sum each cosine term
        sign = (-1) ** p                                  # alternating sign
        w += sign * a_p * np.cos(2 * np.pi * p * n / N)  # Equation (A.16)
    return w

def rectangular(N):
    return np.ones(N)                                     # all ones

def hann(N):
    return _cosine_sum(N, [0.5, 0.5])                     # a0=0.5, a1=0.5

def hamming(N):
    return _cosine_sum(N, [0.54, 0.46])                   # a0=0.54, a1=0.46

def blackman(N):
    return _cosine_sum(N, [0.42, 0.5, 0.08])              # a0=0.42, a1=0.5, a2=0.08

def gaussian(N, sigma_ratio=0.4):
    n = np.arange(N)                                      # sample indices
    center = (N - 1) / 2                                  # window center
    sigma = sigma_ratio * center                          # σ in samples
    w = np.exp(-0.5 * ((n - center) / sigma) ** 2)        # Equation (A.22)
    return w
```

Each window is a direct implementation of its formula from Table A.1 / Equations (A.17)–(A.22). The `_cosine_sum` helper makes the general structure explicit: provide the coefficients, get the window.

### Plotting Utilities (`src/common/plotting.py`)

These functions enforce the CLAUDE.md graph standards: 300 DPI, axis labels with physical units, the dual-stack rule (linear first, dB second), and forbidden colormap rejection.

**Time-domain plot** — always the first plot for any signal:

```python
import matplotlib.pyplot as plt
from .config import DPI, COLORMAP, FIGURE_FORMAT, GRAPH_DIR, FORBIDDEN_COLORMAPS

def plot_time_domain(t, x, xlabel="Time (s)", ylabel="Amplitude",
                     title=None, fig_id=None, t_range=None):
    fig, ax = plt.subplots(1, 1, figsize=(12, 3))         # wide, short aspect
    if t_range is not None:                                # zoom to range if specified
        mask = (t >= t_range[0]) & (t <= t_range[1])       # boolean mask
        ax.plot(t[mask], x[mask], linewidth=0.5)           # plot zoomed segment
    else:
        ax.plot(t, x, linewidth=0.5)                       # plot full signal
    ax.set_xlabel(xlabel)                                  # x label with units
    ax.set_ylabel(ylabel)                                  # y label with units
    if fig_id:
        ax.set_title(fig_id, loc="left", fontsize=9, fontstyle="italic")
    ax.grid(True, alpha=0.3)                               # light grid
    fig.set_dpi(DPI)                                       # enforce minimum DPI
    fig.tight_layout()                                     # prevent label clipping
    return fig, ax
```

**Dual-stack spectrum plot** — linear scale on top (primary, physical units), dB scale on bottom (secondary, dynamic range):

```python
def plot_dual_stack_spectrum(freqs, power, xlabel="Frequency (Hz)",
                             ylabel_linear="Power", ylabel_db="Power (dB)",
                             title=None, fig_id=None, f_range=None):
    fig, (ax_lin, ax_db) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    # --- Top: linear scale (primary) ---
    ax_lin.plot(freqs, power, linewidth=0.5)               # linear power
    ax_lin.set_ylabel(ylabel_linear)                       # label with physical units
    ax_lin.grid(True, alpha=0.3)                           # light grid

    # --- Bottom: log/dB scale (secondary) ---
    power_db = 10 * np.log10(np.maximum(power, 1e-20))     # dB, floor at -200 dB
    ax_db.plot(freqs, power_db, linewidth=0.5)             # dB power
    ax_db.set_ylabel(ylabel_db)                            # label with dB units
    ax_db.set_xlabel(xlabel)                               # frequency label on bottom
    ax_db.grid(True, alpha=0.3)                            # light grid

    if f_range is not None:
        ax_lin.set_xlim(f_range)                           # zoom (applies to both)

    fig.set_dpi(DPI)
    fig.tight_layout()
    return fig, (ax_lin, ax_db)
```

The spectrogram dual-stack (`plot_spectrogram`) follows the same pattern for 2D time-frequency data, with `pcolormesh` and `viridis` colormap. It also validates colormaps against the forbidden list (`jet`, `rainbow`, `hsv`). Full source: `src/common/plotting.py`.

---
---

## B.1 — Lab 1: The DFT of Basic Signals  *(↔ A.1, A.2)*

### Introduction

The DFT (Section A.2) transforms a finite discrete signal into a finite set of frequency-domain coefficients — the bins. Section A.2.1 established that the DFT is the DTFT sampled at $N$ equally spaced frequencies; Section A.2.2 defined a bin as one sample of that continuous spectrum; Section A.2.3 drew the distinction between resolution (determined by signal duration) and bin count (determined by DFT length).

This lab tests those claims directly. We construct a single tone and a dual-tone chord, compute their DFTs, and verify:

- The bin spacing matches Equation (A.6): $\Delta f = f_s / N$.
- An on-bin tone produces a single spike with no leakage.
- Two tones separated by more than $\Delta f$ are resolved as distinct peaks.
- Zero-padding increases bin count but does not improve resolution.

### Setup

**Model signal formulas.**

Single tone (Equation (AA.1)):

$$
x[n] = A \cos\!\left(\frac{2\pi f_0}{f_s} n\right), \qquad n = 0, 1, \ldots, N-1 \tag{B.1}
$$

Dual-tone chord (Equation (AA.2)):

$$
x[n] = A_1 \cos\!\left(\frac{2\pi f_1}{f_s} n\right) + A_2 \cos\!\left(\frac{2\pi f_2}{f_s} n\right) \tag{B.2}
$$

**Lab-specific helper functions.**

The DFT computation and power spectrum extraction:

```python
import numpy as np
from src.common import FS, DURATION, make_tone, make_mixed_tones
from src.common import plot_time_domain, plot_dual_stack_spectrum

def compute_dft(x, fs=FS, n_fft=None):
    if n_fft is None:
        n_fft = len(x)                                    # DFT length = signal length
    X_full = np.fft.fft(x, n=n_fft)                       # full DFT (Equation A.3)
    freqs_full = np.fft.fftfreq(n_fft, d=1/fs)            # frequency axis in Hz
    pos = freqs_full >= 0                                  # positive frequencies only
    return freqs_full[pos], X_full[pos]                    # return positive half

def compute_power(X, N):
    return np.abs(X)**2 / N                                # P[k] = |X[k]|²/N (Eq. A.29)
```

**Experiment A — On-bin vs. off-bin.** A 1-second segment ($N = 250$, $\Delta f = 1.0$ Hz) so that leakage is visible. On-bin: $f_0 = 10.0$ Hz (bin 10, integer). Off-bin: $f_0 = 10.5$ Hz (bin 10.5, maximally between two bins — worst-case leakage).

```python
# --- Generate signals (1-second segment: Δf = 1 Hz, leakage visible) ---
x_on, n_on, t_on = make_tone(10.0, A=1.0, duration=1.0)   # on-bin tone
x_off, n_off, t_off = make_tone(10.5, A=1.0, duration=1.0) # maximally off-bin

# --- Compute DFTs and power spectra ---
freqs_on, X_on = compute_dft(x_on)                        # DFT of on-bin tone
freqs_off, X_off = compute_dft(x_off)                     # DFT of off-bin tone
P_on = compute_power(X_on, len(x_on))                     # power spectrum (on-bin)
P_off = compute_power(X_off, len(x_off))                  # power spectrum (off-bin)

# --- Plot: time domain first, then dual-stack spectrum ---
plot_time_domain(t_on, x_on, fig_id="Figure B.1a", t_range=(0, 0.5))
plot_dual_stack_spectrum(freqs_on, P_on, fig_id="Figure B.2a", f_range=(0, 30))
plot_dual_stack_spectrum(freqs_off, P_off, fig_id="Figure B.2b", f_range=(0, 30))
```

**Experiment B — Dual-tone chord.** $f_1 = 10.0$ Hz, $f_2 = 12.0$ Hz.

```python
x, n, t = make_mixed_tones([10.0, 12.0], amplitudes=[1.0, 1.0], duration=1200)
freqs, X = compute_dft(x)                                 # DFT of chord
P = compute_power(X, len(x))                              # power spectrum

plot_time_domain(t, x, fig_id="Figure B.3a", t_range=(0, 0.5))
plot_dual_stack_spectrum(freqs, P, fig_id="Figure B.3b", f_range=(0, 30))
```

**Experiment C — Zero-padding.** 1-second segment, two tones at 10 Hz and 11 Hz, 4× zero-padding.

```python
x, n, t = make_mixed_tones([10.0, 11.0], amplitudes=[1.0, 1.0], duration=1.0)
N_orig = len(x)                                           # 250 samples
N_padded = N_orig * 4                                     # 1000 samples (zero-padded)

freqs_orig, X_orig = compute_dft(x, n_fft=N_orig)        # no zero-padding
freqs_pad, X_pad = compute_dft(x, n_fft=N_padded)        # zero-padded

plot_dual_stack_spectrum(freqs_orig, compute_power(X_orig, N_orig),
                         fig_id="Figure B.4a", f_range=(5, 16))
plot_dual_stack_spectrum(freqs_pad, compute_power(X_pad, N_orig),
                         fig_id="Figure B.4b", f_range=(5, 16))
```

### Parameters

**Table B.1 — Lab 1 parameters**

| Parameter | Experiment A | Experiment B | Experiment C |
| --- | --- | --- | --- |
| $f_s$ (Hz) | 250 | 250 | 250 |
| Duration (s) | 1.0 | 1200 | 1.0 |
| $N$ (samples) | 250 | 300 000 | 250 |
| $\Delta f$ (Hz) | 1.0 | 0.000833 | 1.0 |
| $f_0$ or $f_1$ (Hz) | 10.0 / 10.5 | 10.0 | 10.0 |
| $f_2$ (Hz) | — | 12.0 | 11.0 |
| Amplitude $A$ | 1.0 | 1.0 | 1.0 |
| Zero-pad factor | — | — | 4× |

### Results

**Experiment A** — Figure B.1 shows the time-domain waveforms of both tones over the first 0.5 s. Both are clean cosines, indistinguishable by eye.

![Figure B.1a — On-bin tone, time domain](../results/graphs/lab1/figure_B_01.png)

![Figure B.1b — Off-bin tone, time domain](../results/graphs/lab1/figure_B_02.png)

Figure B.2 shows the dual-stack power spectra — and the contrast is dramatic. The on-bin tone at 10.0 Hz (Figure B.2a) produces a single spike with zero leakage: the dB panel shows −200 dB (numerical floor) at all other bins. The off-bin tone at 10.5 Hz (Figure B.2b) shows **maximum leakage**: the tone's energy is split between bins 10 and 11 (neither captures it fully), and the side lobes spread power across the entire spectrum. The dB panel never drops below −25 dB — energy is everywhere.

![Figure B.2a — On-bin tone spectrum, no leakage](../results/graphs/lab1/figure_B_03.png)

![Figure B.2b — Off-bin tone spectrum, maximum leakage](../results/graphs/lab1/figure_B_04.png)

This is the Dirichlet kernel (Lab 3, Equation (B.11)) in action. At 10.5 Hz, the tone falls exactly at the midpoint between two bins ($f_0 / \Delta f = 10.5$, half-integer). The DFT evaluates the DTFT at the bin frequencies, and none of them align with the tone — every bin sees the tone through a side lobe. The 1-second duration ($N = 250$, $\Delta f = 1.0$ Hz) makes this effect maximally visible; at the full 1200-second duration, the bin grid is so fine ($\Delta f = 0.000833$ Hz) that nearly every frequency lands on a bin and leakage vanishes. This is why windowing (Section A.3, Lab 3) exists: to suppress the side lobes that cause this leakage.

**Experiment B** — Figure B.3 shows the dual-tone chord. The time-domain plot (Figure B.3a) shows the expected beat pattern. The spectrum (Figure B.3b) shows two clean spikes at 10 Hz and 12 Hz, fully resolved. The separation (2.0 Hz) is $2400 \times \Delta f$.

![Figure B.3a — Dual-tone chord, time domain](../results/graphs/lab1/figure_B_05.png)

![Figure B.3b — Dual-tone chord spectrum](../results/graphs/lab1/figure_B_06.png)

**Experiment C** — Figure B.4 compares the original and zero-padded spectra. The original ($N = 250$, $\Delta f = 1.0$ Hz, Figure B.4a) shows two tones at the resolution limit — barely distinguishable. The zero-padded spectrum ($N = 1000$, $\Delta f = 0.25$ Hz, Figure B.4b) shows the same two lobes sampled more densely — smoother, but no sharper. Zero-padding interpolated the same DTFT curve; it did not resolve the two tones.

![Figure B.4a — Original, no zero-padding](../results/graphs/lab1/figure_B_07.png)

![Figure B.4b — Zero-padded 4×](../results/graphs/lab1/figure_B_08.png)

### Verification

| Prediction (Volume A) | Measured | Confirmed? |
| --- | --- | --- |
| $\Delta f = f_s/N = 250/250 = 1.0$ Hz (Eq. A.6) | 1.0 Hz | Yes |
| On-bin tone (10.0 Hz) → single spike, no leakage | Zero leakage (−200 dB floor) | Yes |
| Off-bin tone (10.5 Hz) → energy spread across all bins | Leakage across entire band (never below −25 dB) | Yes |
| Two tones at 2.0 Hz separation → resolved | Two distinct spikes | Yes |
| Zero-padding (4×) changes $\Delta f$ but not resolution | Bin spacing reduced, resolution unchanged | Yes |

### Conclusion

The DFT behaves as Section A.2 predicts. Bin spacing is $f_s / N$; on-bin tones produce zero leakage. The off-bin experiment makes the cost of the rectangular window unmistakable: a single tone at 10.5 Hz — maximally between two bins — leaks energy across every bin in the spectrum. This is the Dirichlet kernel's side-lobe structure (Appendix B) made visible. The zero-padding experiment confirms that more bins do not mean more resolution. Windowing (Appendix B, Lab 3) is the remedy for leakage.

---

## B.2 — Lab 2: Statistics on a Noisy Signal  *(↔ A.4)*

### Introduction

Section A.4 established that DFT bins under noise are random variables: the power at each bin follows an exponential distribution (Equation (A.25)), the noise floor can be estimated from the spectrum itself (Equation (A.27)), and detection uses a threshold derived from the exponential model (Equation (A.28)). Section A.4.3 showed that the periodogram's variance does not decrease with $N$, and Section A.4.4 introduced Welch's method as the fix.

This lab tests:

- DFT bins of pure noise follow the predicted distributions (Rayleigh magnitude, uniform phase, exponential power).
- A tone buried in noise is detected using the spectral threshold $\gamma$, with false-alarm probability $e^{-\gamma}$.
- Welch averaging reduces variance at the cost of frequency resolution.

### Setup

**Model signals.**

Pure noise (Equation (AA.8)):

$$
x[n] = \eta[n], \qquad \eta[n] \sim \mathcal{N}(0, \sigma^2) \tag{B.3}
$$

Tone buried in noise:

$$
x[n] = A \cos\!\left(\frac{2\pi f_0}{f_s} n\right) + \eta[n] \tag{B.4}
$$

**Experiment A — Bin distributions.** Compute the DFT of pure noise and histogram the magnitude, phase, and power of all bins.

```python
from src.common import FS, DURATION, SEED, make_noise

x_noise, n, t = make_noise(sigma=1.0, duration=1200, seed=42)

N = len(x_noise)                                          # 300,000 samples
X = np.fft.fft(x_noise)                                   # full DFT
freqs = np.fft.fftfreq(N, d=1/FS)                         # frequency axis
pos = (freqs > 0) & (freqs < FS/2)                        # exclude DC and Nyquist

magnitudes = np.abs(X[pos])                               # |X[k]| → Rayleigh
phases = np.angle(X[pos])                                 # ∠X[k] → Uniform(-π, π]
powers = np.abs(X[pos])**2                                # |X[k]|² → Exponential
```

**Experiment B — Spectral detection.** Bury a tone ($A = 0.5$) in noise ($\sigma = 1.0$), estimate the noise floor, and apply thresholds.

```python
from src.common import make_tone

x_tone, _, _ = make_tone(10.0, A=0.5, duration=1200)
x_noise, _, _ = make_noise(sigma=1.0, duration=1200, seed=42)
x = x_tone + x_noise                                      # combined signal

# --- Periodogram ---
freqs, Sxx = compute_periodogram(x)                       # |X[k]|²/N

# --- Noise floor (median estimator, Equation A.27) ---
noise_floor = np.median(Sxx)

# --- Detection thresholds (Equation A.28) ---
for gamma in [3.0, 4.6, 6.9]:
    threshold = gamma * noise_floor                       # γ × noise floor
    false_alarm = np.exp(-gamma)                          # Pr(exceed | noise only)
```

Where `compute_periodogram` is:

```python
def compute_periodogram(x, fs=FS):
    N = len(x)                                            # signal length
    X = np.fft.fft(x)                                     # full DFT
    freqs = np.fft.fftfreq(N, d=1/fs)                     # frequency axis
    pos = freqs >= 0                                      # positive frequencies
    Sxx = np.abs(X[pos])**2 / N                           # periodogram (Eq. A.29)
    return freqs[pos], Sxx
```

**Experiment C — Welch averaging.** Apply `scipy.signal.welch` at four segment lengths.

```python
from scipy import signal as sp_signal

for seg_dur in [1200, 20.0, 5.0, 2.0]:                   # segment durations (s)
    nperseg = int(seg_dur * FS)                           # segment length in samples
    nperseg = min(nperseg, len(x))                        # cap at signal length
    noverlap = nperseg // 2                               # 50% overlap (Hann COLA)

    f_welch, Pxx = sp_signal.welch(x, fs=FS,             # scipy Welch implementation
                                    nperseg=nperseg,
                                    noverlap=noverlap,
                                    window="hann")
```

### Parameters

**Table B.2 — Lab 2 parameters**

| Parameter | Value |
| --- | --- |
| $f_s$ (Hz) | 250 |
| Duration (s) | 1200 |
| $N$ (samples) | 300 000 |
| Tone frequency $f_0$ (Hz) | 10.0 |
| Tone amplitude $A$ | 0.5 |
| Noise $\sigma$ | 1.0 |
| Seed | 42 |
| Detection thresholds $\gamma$ | 3.0, 4.6, 6.9 |
| Welch segments (s) | 1200, 20, 5, 2 |
| Welch window | Hann |
| Welch overlap | 50% |

### Results

**Experiment A** — Figure B.5 shows three histograms of DFT bin statistics for pure noise:

![Figure B.5 — DFT bin distributions under white Gaussian noise](../results/graphs/lab2/figure_B_01.png)

- *Magnitude* $|X[k]|$: Rayleigh distribution — zero at origin, single peak, long tail.
- *Phase* $\angle X[k]$: uniform on $(-\pi, \pi]$. The theoretical line at $1/(2\pi) \approx 0.159$ matches the histogram.
- *Power* $|X[k]|^2$: exponential distribution. The theoretical curve ($\lambda = 1/(N\sigma^2)$) overlays the histogram closely. Measured mean: 300,413. Predicted: $N\sigma^2 = 300\,000$. Deviation: 0.14%.

**Experiment B** — Figure B.6 shows the time domain of the tone-in-noise signal over the first 2 seconds. The tone ($A = 0.5$) is invisible — buried in noise ($\sigma = 1.0$). Time-domain inspection cannot detect it.

![Figure B.6 — Tone buried in noise, time domain](../results/graphs/lab2/figure_B_02.png)

Figure B.7 shows the periodogram with detection thresholds. The tone at 10 Hz produces a power of 18,569 — a ratio of 26,788× the noise floor. Detected at all three thresholds:

![Figure B.7 — Periodogram with detection thresholds](../results/graphs/lab2/figure_B_03.png)

| Threshold $\gamma$ | $P_{fa} = e^{-\gamma}$ | Threshold value | Detected? |
| --- | --- | --- | --- |
| 3.0 | 0.050 | 2.08 | Yes (ratio = 26,788) |
| 4.6 | 0.010 | 3.19 | Yes |
| 6.9 | 0.001 | 4.78 | Yes |

**Experiment C** — Figure B.8 shows the Welch progression (dual-stack: linear left, dB right):

![Figure B.8 — Welch averaging progression](../results/graphs/lab2/figure_B_04.png)

| Segment | $L$ | $\Delta f$ (Hz) | Relative variance | Spectrum appearance |
| --- | --- | --- | --- | --- |
| Full (1200 s) | 1 | 0.0008 | 1.000 | Ragged, ±10 dB fluctuations |
| 20 s | 119 | 0.050 | 0.008 | Smooth, clean tone peak |
| 5 s | 479 | 0.200 | 0.002 | Very smooth, flat noise floor |
| 2 s | 1199 | 0.500 | 0.001 | Smoothest, but tone peak is wide |

### Verification

| Prediction (Volume A) | Measured | Confirmed? |
| --- | --- | --- |
| $\mathbb{E}[\|X[k]\|^2] = N\sigma^2 = 300\,000$ (Eq. A.26) | 300,413 | Yes (0.14% deviation) |
| Phase uniform on $(-\pi, \pi]$ (A.4.1) | Uniform histogram | Yes |
| Power exponential with mean $N\sigma^2$ (Eq. A.25) | Exponential fit matches | Yes |
| $\gamma = 3.0 \Rightarrow P_{fa} = 0.050$ (Eq. A.28) | Tone detected (ratio 26,788) | Yes |
| Welch variance $\propto 1/L$ (Eq. A.32) | Progressive smoothing | Yes |
| Welch resolution degrades to $\beta \cdot f_s/M$ (Eq. A.21) | $\Delta f$ increases with shorter $M$ | Yes |

### Conclusion

The spectral statistics framework from Section A.4 holds. Bin power under noise follows the exponential distribution. The noise floor estimated from the spectrum itself (median) enables detection thresholds grounded in probability, not arbitrary σ rules. The tone at 10 Hz — invisible in the time domain — is detected with a ratio exceeding 26,000 in the frequency domain.

Welch's method demonstrates the resolution-variance tradeoff: 5-second segments ($\Delta f = 0.2$ Hz, $L = 479$) produce a smooth spectrum that resolves all EEG bands while keeping variance below 1%.

---

## Appendix B — Window Derivations  *(↔ A.3)*

> All derivations use $M = 256$ samples (≈ 1.024 s at $f_s = 250$ Hz — the typical EEG STFT epoch length). All spectra are normalized by $M$ so that the main-lobe peak is 1. Graphs are rendered at high zero-pad ($2048 \times M$) for visual smoothness. Code and figures: `src/appendix_b/appendix_b.py`.

---

### The Dirichlet Kernel

The DFT of the rectangular window $w[n] = 1$ for $n = 0, 1, \ldots, M-1$ is computed from the finite geometric series. Starting from the definition:

$$
W(\omega) = \sum_{n=0}^{M-1} e^{-j\omega n} \tag{B.5}
$$

This is a geometric series with ratio $r = e^{-j\omega}$:

$$
W(\omega) = \frac{1 - e^{-j\omega M}}{1 - e^{-j\omega}} \tag{B.6}
$$

To extract the magnitude, we factor out the half-angle from both numerator and denominator. From the numerator:

$$
1 - e^{-j\omega M} = e^{-j\omega M/2}\left(e^{j\omega M/2} - e^{-j\omega M/2}\right) = e^{-j\omega M/2} \cdot 2j\sin\!\left(\frac{\omega M}{2}\right) \tag{B.7}
$$

From the denominator:

$$
1 - e^{-j\omega} = e^{-j\omega/2}\left(e^{j\omega/2} - e^{-j\omega/2}\right) = e^{-j\omega/2} \cdot 2j\sin\!\left(\frac{\omega}{2}\right) \tag{B.8}
$$

Dividing Equation (B.7) by Equation (B.8):

$$
W(\omega) = \frac{e^{-j\omega M/2} \cdot 2j\sin(\omega M/2)}{e^{-j\omega/2} \cdot 2j\sin(\omega/2)} = e^{-j\omega(M-1)/2} \cdot \frac{\sin(\omega M/2)}{\sin(\omega/2)} \tag{B.9}
$$

The factor $e^{-j\omega(M-1)/2}$ is a linear phase (time shift to the window center). It affects phase but not magnitude. The magnitude is the **Dirichlet kernel**:

$$
|W(\omega)| = \left|\frac{\sin(\omega M/2)}{\sin(\omega/2)}\right| \tag{B.10}
$$

**Normalization.** At $\omega = 0$, both sine functions vanish, but L'Hôpital's rule (or direct substitution via $\lim_{x \to 0} \sin(Mx)/\sin(x) = M$) gives $|W(0)| = M$. We normalize by $M$ so that the peak equals 1:

$$
D(\omega) = \frac{1}{M}\left|\frac{\sin(\omega M/2)}{\sin(\omega/2)}\right| \tag{B.11}
$$

Equation (B.11) is the normalized Dirichlet kernel. It is the spectral footprint of the rectangular window: every tone in a rectangularly-windowed DFT is convolved with this shape (Section A.3.2).

**Code:**

```python
import numpy as np

def dirichlet_kernel(omega, M):
    """Normalized Dirichlet kernel: |sin(ωM/2) / sin(ω/2)| / M"""
    numerator = np.sin(omega * M / 2)                     # sin(ωM/2)
    denominator = np.sin(omega / 2)                       # sin(ω/2)
    with np.errstate(divide="ignore", invalid="ignore"):  # handle 0/0 at ω=0
        D = np.abs(numerator / denominator) / M           # Equation (B.11)
    D[np.isnan(D)] = 1.0                                  # limit at ω=0 is 1.0
    return D
```

---

### Anatomy of the Dirichlet Kernel

Figure B.9 shows the normalized Dirichlet kernel $D(\omega)$ at $M = 256$, plotted as a continuous function of frequency in bins.

![Figure B.9 — Dirichlet kernel anatomy](../results/graphs/lab3/figure_B_09.png)

**Nulls.** The numerator $\sin(\omega M/2)$ vanishes when $\omega M/2 = k\pi$ for integer $k \neq 0$, i.e. at:

$$
\omega_k = \frac{2\pi k}{M}, \qquad k = \pm 1, \pm 2, \ldots \tag{B.12}
$$

In bin units ($\text{bin} = \omega M / (2\pi)$), the nulls fall at **integer bins**: $k = \pm 1, \pm 2, \ldots$. These are visible as the zero-crossings in Figure B.9.

**Main lobe.** The central peak between the first nulls at $k = -1$ and $k = +1$ is the **main lobe**. Its width is 2 bins (from $-1$ to $+1$). This is the narrowest possible main lobe — rectangular pays for it with the highest side lobes.

**Side lobes.** Between each pair of adjacent nulls lies a **side lobe** — a local maximum of $D(\omega)$. Figure B.9 annotates the first three side-lobe maxima with their positions and magnitudes.

**The skew observation.** The side-lobe maxima are **not** centered between the nulls. The first maximum is at bin 1.430, not 1.500. This is because $D(\omega) = |\sin(\omega M/2) / \sin(\omega/2)|$ is not a pure sinusoid — it is a ratio of two sines with different frequencies. The denominator $1/\sin(\omega/2)$ is a monotonically decreasing envelope that pulls each maximum slightly toward the origin (toward the main lobe).

The **midpoint approximation** $\omega \approx (2k+1)\pi/M$ (i.e. bin $\approx k + 0.5$) is commonly used and close, but not exact. This matters for the decay-rate analysis below, where we compare actual maxima positions against the approximation.

**Envelope.** Figure B.10 shows the kernel with the envelope $1/(M \cdot |\sin(\omega/2)|)$ overlaid.

![Figure B.10 — Dirichlet kernel with envelope](../results/graphs/lab3/figure_B_10.png)

For large $\omega$ (far from the main lobe), $\sin(\omega/2) \approx \omega/2$ breaks down, but the envelope still tracks the side-lobe peaks accurately. The side lobes touch the envelope because $|\sin(\omega M/2)|$ reaches 1 near (but not exactly at) each maximum.

**Code for dense spectrum computation:**

```python
M = 256                                                   # EEG-typical window length
N_DENSE = M * 2048                                        # high zero-pad for smoothness

def window_spectrum_dense(win_func, M, N_dense=N_DENSE):
    w = win_func(M)                                       # generate window
    W = np.fft.fft(w, n=N_dense)                          # zero-padded DFT
    W_shifted = np.fft.fftshift(W)                        # center at DC
    W_mag = np.abs(W_shifted)                             # magnitude
    W_mag = W_mag / W_mag.max()                           # normalize peak to 1.0
    bins = (np.arange(N_dense) - N_dense // 2) / (N_dense / M)
    return bins, W_mag
```

---

### Properties of the Rectangular Window

#### First side-lobe strength

The first side lobe is the tallest, and its level relative to the main lobe determines the worst-case leakage from an off-bin tone.

**Method (a): visual measurement from the graph.**

Figure B.11 zooms into the first side lobe.

![Figure B.11 — First side-lobe analysis](../results/graphs/lab3/figure_B_11.png)

The true maximum is located by finding the local peak of the computed spectrum:

```python
from scipy.signal import argrelextrema

def find_sidelobe_maxima(bins, W_mag, main_lobe_bins=1.5):
    right_half = bins > main_lobe_bins                    # positive side, outside main lobe
    indices = np.where(right_half)[0]                     # indices in that region
    W_right = W_mag[indices]                              # magnitude in that region
    bins_right = bins[indices]                            # bins in that region
    local_max_idx = argrelextrema(W_right, np.greater, order=5)[0]
    return bins_right[local_max_idx], W_right[local_max_idx]
```

Result: the first side-lobe maximum is at **bin 1.430** with normalized magnitude **0.21724**, corresponding to $20\log_{10}(0.21724) = -13.3$ dB.

**Method (b): analytical approximation via $k = 1.5$.**

Between the nulls at bin 1 and bin 2, the midpoint approximation places the maximum at bin 1.5, i.e. $\omega = 3\pi/M$. The numerator $|\sin(\omega M/2)| = |\sin(3\pi/2)| = 1$, and the denominator is $M \cdot |\sin(3\pi/(2M))|$. For large $M$, $\sin(3\pi/(2M)) \approx 3\pi/(2M)$, giving:

$$
D\!\left(\frac{3\pi}{M}\right) \approx \frac{1}{M \cdot \frac{3\pi}{2M}} = \frac{2}{3\pi} \approx 0.21221 \tag{B.13}
$$

In dB: $20\log_{10}(2/(3\pi)) = -13.5$ dB.

**Comparison:**

| Method | Bin position | Magnitude | dB |
| --- | --- | --- | --- |
| (a) True maximum | 1.430 | 0.21724 | −13.3 |
| (b) $k = 1.5$ approximation | 1.500 | 0.21222 | −13.5 |
| Textbook (asymptotic) | — | $2/(3\pi)$ | −13.0 |

The three values agree within 0.5 dB. The textbook value of −13.0 dB is the asymptotic limit for $M \to \infty$; at $M = 256$ the actual maximum is slightly higher (−13.3 dB) because the discrete kernel has not fully converged. The $k = 1.5$ approximation underestimates slightly because the true maximum is skewed toward the main lobe (bin 1.430, not 1.500).

#### Decay rate per octave

**Method (a): analytical approximation.**

The $k$-th side-lobe peak is near $\omega \approx (2k+1)\pi/M$, where the numerator $|\sin(\omega M/2)| \approx 1$. The envelope is then:

$$
D_k \approx \frac{1}{M \cdot \sin\!\left(\frac{(2k+1)\pi}{2M}\right)} \tag{B.14}
$$

For large $M$ and moderate $k$, $\sin(x) \approx x$, so:

$$
D_k \approx \frac{2}{(2k+1)\pi} \tag{B.15}
$$

Equation (B.15) decays as $1/k$, which is $1/\omega$ since $\omega \propto k$. In dB, doubling $\omega$ (one octave) reduces the level by $20\log_{10}(2) = 6.02$ dB. This is the **6 dB/octave** rolloff ($1/\omega^1$ decay) characteristic of the rectangular window, caused by the discontinuity in the window value at the edges.

**Method (b): statistical measurement via regression.**

We collect the actual positions and magnitudes of the first 6 side-lobe maxima and fit a line on the log-log plot.

**Table B.4 — Side-lobe maxima of the Dirichlet kernel ($M = 256$)**

| Lobe | Actual bin | Midpoint approx | Actual magnitude | $(2k+1)\pi$ approx | Actual dB | Approx dB |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.430 | 1.5 | 0.21724 | 0.21221 | −13.3 | −13.5 |
| 2 | 2.459 | 2.5 | 0.12839 | 0.12732 | −17.8 | −17.9 |
| 3 | 3.471 | 3.5 | 0.09135 | 0.09095 | −20.8 | −20.8 |
| 4 | 4.478 | 4.5 | 0.07095 | 0.07074 | −23.0 | −23.0 |
| 5 | 5.481 | 5.5 | 0.05802 | 0.05787 | −24.7 | −24.8 |
| 6 | 6.484 | 6.5 | 0.04908 | 0.04897 | −26.2 | −26.2 |

Table B.4 confirms that the $(2k+1)\pi$ approximation (Equation (B.15)) is excellent for $k \geq 2$ — within 0.2 dB. The deviation is largest at $k = 1$ (the first side lobe), where the skew effect is strongest.

**Regression code:**

```python
from scipy.stats import linregress

log_bins = np.log10(peak_bins)                            # log of bin positions
log_mag = np.log10(peak_values)                           # log of magnitudes
slope, intercept, r_value, _, _ = linregress(log_bins, log_mag)
dB_per_octave = slope * 20 * np.log10(2)                  # convert to dB/octave
```

**Results** (Figure B.12):

![Figure B.12 — Side-lobe decay analysis](../results/graphs/lab3/figure_B_12.png)

| Regression method | Slope | R² | dB/octave |
| --- | --- | --- | --- |
| Accurate $k'$ (actual maxima positions) | −0.985 | 0.99997 | −5.9 |
| Crude integer $k$ (midpoint approximation) | −1.014 | 0.99996 | −6.1 |
| Theoretical | −1.000 | — | −6.0 |

Both regressions confirm the $1/\omega$ decay to within 0.1 dB/octave. The R² values (> 0.9999) indicate near-perfect fit to the power-law model. The crude midpoint approximation is marginally closer to the theoretical slope, but both are effectively exact at this precision.

#### From Rectangular to Hann, Hamming, Blackman

The cosine-sum windows (Equation (A.16)) can be decomposed into shifted rectangular windows using the identity:

$$
\cos\!\left(\frac{2\pi p n}{M}\right) = \frac{1}{2}\left(e^{j 2\pi p n / M} + e^{-j 2\pi p n / M}\right) \tag{B.16}
$$

Multiplying the rectangular window $w_{\text{rect}}[n] = 1$ by $\cos(2\pi p n / M)$ is equivalent to shifting the Dirichlet kernel by $\pm p$ bins in frequency:

$$
\text{DFT}\!\left\{w_{\text{rect}}[n] \cdot \cos\!\left(\frac{2\pi p n}{M}\right)\right\} = \frac{1}{2}\left[W(\omega - 2\pi p/M) + W(\omega + 2\pi p/M)\right] \tag{B.17}
$$

where $W(\omega)$ is the rectangular window DFT from Equation (B.9).

**Hann** ($a_0 = 0.5, a_1 = 0.5$):

$$
w_{\text{Hann}}[n] = 0.5 - 0.5\cos\!\left(\frac{2\pi n}{M}\right) = 0.5 \cdot w_{\text{rect}}[n] - 0.5\cos\!\left(\frac{2\pi n}{M}\right) \cdot w_{\text{rect}}[n]
$$

Taking the DFT and applying Equation (B.17) with $p = 1$:

$$
W_{\text{Hann}}(\omega) = 0.5\,W(\omega) - 0.25\,W(\omega - 2\pi/M) - 0.25\,W(\omega + 2\pi/M) \tag{B.18}
$$

Equation (B.18) is the Hann window spectrum: one Dirichlet kernel centered at the origin (weight 0.5) minus two copies shifted by $\pm 1$ bin (weight 0.25 each). The negative signs arise from the $(-1)^p$ alternation in Equation (A.16).

**Expanding to the pure sine form.** Equation (B.18) is conceptually clear — three shifted copies — but to analyze cancellation quantitatively, we substitute the sin/sin expression for each shifted kernel and simplify.

Each shifted kernel has the form (from Equation (B.9)):

$$
W(\omega - 2\pi p/M) = e^{-j(\omega - 2\pi p/M)(M-1)/2} \cdot \frac{\sin\!\left(\frac{(\omega - 2\pi p/M) M}{2}\right)}{\sin\!\left(\frac{\omega - 2\pi p/M}{2}\right)}
$$

The key simplification is in the numerator. For any integer $p$:

$$
\sin\!\left(\frac{(\omega - 2\pi p/M) M}{2}\right) = \sin\!\left(\frac{\omega M}{2} - p\pi\right) = (-1)^p \sin\!\left(\frac{\omega M}{2}\right) \tag{B.19}
$$

Equation (B.19) is the critical identity: shifting by $p$ bins multiplies the numerator by $(-1)^p$ but otherwise leaves it unchanged. The $\sin(\omega M/2)$ factor is **shared by all shifted kernels**.

Substituting into Equation (B.18) for **Hann** ($p = 0$ and $p = \pm 1$, so $(-1)^0 = 1$ and $(-1)^1 = -1$):

$$
W_{\text{Hann}}(\omega) = e^{-j\phi(\omega)} \cdot \sin\!\left(\frac{\omega M}{2}\right) \cdot \left[\frac{0.5}{\sin(\omega/2)} + \frac{0.25}{\sin(\omega/2 - \pi/M)} + \frac{0.25}{\sin(\omega/2 + \pi/M)}\right] \tag{B.20}
$$

where $e^{-j\phi(\omega)}$ collects the phase terms (which do not affect the magnitude). The negative signs from the $(-1)^p$ factor cancel the negative signs from the $-0.25$ coefficients in Equation (B.18), so **all three terms in the bracket are positive**.

The normalized magnitude is:

$$
\frac{|W_{\text{Hann}}(\omega)|}{M} = \frac{|\sin(\omega M/2)|}{M} \cdot \left[\frac{0.5}{\sin(\omega/2)} + \frac{0.25}{\sin(\omega/2 - \pi/M)} + \frac{0.25}{\sin(\omega/2 + \pi/M)}\right] \tag{B.21}
$$

Equation (B.21) is the **pure sine form** of the Hann window spectrum. It is a single expression: one shared numerator $\sin(\omega M/2)$ multiplied by a sum of three $1/\sin$ terms. The structure is transparent:

- **Zeros:** the shared numerator $\sin(\omega M/2) = 0$ at $\omega = 2\pi k/M$ (integer bins), same as rectangular — but the denominator terms also vanish at $\omega = 2\pi/M$ (bin 1) and $\omega = -2\pi/M$ (bin $-1$), creating 0/0 forms. By L'Hôpital, these evaluate to finite values — they become part of the main lobe rather than nulls. The first true null is pushed out to bin 2, doubling the main-lobe width from 2 bins (rectangular) to 4 bins (Hann).

- **Side-lobe cancellation:** in the side-lobe region (bins > 2), all three $1/\sin$ terms are nonzero and positive. Their sum is much smaller than the single $1/\sin(\omega/2)$ term of rectangular because the shifted terms $1/\sin(\omega/2 \pm \pi/M)$ partially cancel the central term. The cancellation is not exact — a residual remains — but it reduces the side-lobe level from −13.3 dB to −31.5 dB.

The same expansion for **Hamming** ($a_0 = 0.54, a_1 = 0.46$):

$$
\frac{|W_{\text{Hamming}}(\omega)|}{M} = \frac{|\sin(\omega M/2)|}{M} \cdot \left[\frac{0.54}{\sin(\omega/2)} + \frac{0.23}{\sin(\omega/2 - \pi/M)} + \frac{0.23}{\sin(\omega/2 + \pi/M)}\right] \tag{B.22}
$$

The structure is identical to Hann — same shared numerator, same three $1/\sin$ terms — but the coefficients are different: $0.54$ on the central term vs. $0.23$ on the shifted terms (compare Hann's $0.5$ and $0.25$). These coefficients were chosen to minimize the peak side-lobe level. The result is −42.7 dB, deeper than Hann's −31.5 dB.

However, Hamming's coefficients do not sum to zero at the edges: $w[0] = 0.54 - 0.46 = 0.08 \neq 0$. This means the window has a value discontinuity at its boundaries. In the sine form, this manifests as the three $1/\sin$ terms not cancelling to higher order at large $\omega$ — the residual decays as $1/\omega^1$ (6 dB/oct), same as rectangular, even though the nearest side lobes are much lower.

For **Blackman** ($a_0 = 0.42, a_1 = 0.5, a_2 = 0.08$), the same expansion gives five $1/\sin$ terms. The $p = 2$ shifts use Equation (B.19) with $(-1)^2 = +1$:

$$
\frac{|W_{\text{Blackman}}(\omega)|}{M} = \frac{|\sin(\omega M/2)|}{M} \cdot \left[\frac{0.42}{\sin(\omega/2)} + \frac{0.25}{\sin(\omega/2 - \pi/M)} + \frac{0.25}{\sin(\omega/2 + \pi/M)} + \frac{0.04}{\sin(\omega/2 - 2\pi/M)} + \frac{0.04}{\sin(\omega/2 + 2\pi/M)}\right] \tag{B.23}
$$

Five terms, one shared numerator. The denominator has 0/0 forms at bins $\pm 1$ and $\pm 2$, pushing the first true null to bin 3 and widening the main lobe to 6 bins. The two extra terms provide a second level of cancellation in the side-lobe region, driving the peak side-lobe to −58 dB. Blackman goes to zero at the edges ($0.42 - 0.50 + 0.08 = 0$), so the five-term residual decays as $1/\omega^3$ (18 dB/oct).

**Summary: the pure sine forms.**

**Table B.5 — Pure sine form expressions**

All windows share the numerator $|\sin(\omega M/2)|$. The bracket $[\cdots]$ is the weighted sum of $1/\sin$ terms. Let $\alpha = \omega/2$ for compactness.

**Rectangular** (Equation (B.11), 1 term):

$$
\frac{|W_{\text{rect}}(\omega)|}{M} = \frac{|\sin(\alpha M)|}{M \cdot |\sin(\alpha)|} \tag{B.24}
$$

**Hann** (Equation (B.21), 3 terms):

$$
\frac{|W_{\text{Hann}}(\omega)|}{M} = \frac{|\sin(\alpha M)|}{M} \left[\frac{0.5}{\sin(\alpha)} + \frac{0.25}{\sin(\alpha - \pi/M)} + \frac{0.25}{\sin(\alpha + \pi/M)}\right] \tag{B.25}
$$

**Hamming** (Equation (B.22), 3 terms):

$$
\frac{|W_{\text{Ham}}(\omega)|}{M} = \frac{|\sin(\alpha M)|}{M} \left[\frac{0.54}{\sin(\alpha)} + \frac{0.23}{\sin(\alpha - \pi/M)} + \frac{0.23}{\sin(\alpha + \pi/M)}\right] \tag{B.26}
$$

**Blackman** (Equation (B.23), 5 terms):

$$
\frac{|W_{\text{Blk}}(\omega)|}{M} = \frac{|\sin(\alpha M)|}{M} \left[\frac{0.42}{\sin(\alpha)} + \frac{0.25}{\sin(\alpha - \pi/M)} + \frac{0.25}{\sin(\alpha + \pi/M)} + \frac{0.04}{\sin(\alpha - 2\pi/M)} + \frac{0.04}{\sin(\alpha + 2\pi/M)}\right] \tag{B.27}
$$

Every window in the cosine-sum family is: one shared $\sin(\omega M/2)$ numerator (which sets the nulls), multiplied by a weighted sum of $1/\sin$ terms with shifted arguments (which sets the side-lobe structure). The window design problem reduces to choosing the weights that produce the desired cancellation pattern.

**Why Hamming's rolloff is slow despite its low side lobes.** From Table B.5, Hann and Hamming have the same three-term structure. At large $\omega$, each $1/\sin$ term behaves as $1/(\omega/2 + \text{shift}) \approx 2/\omega$. The three terms combine as:

For Hann: $0.5 + 0.25 + 0.25 = 1.0$ — but the signs of the shifted terms (after accounting for the full complex expression) produce cancellation at order $1/\omega$ and $1/\omega^2$, leaving a residual at $1/\omega^3$. This happens because $w[0] = w[M-1] = 0$, so the window and its first derivative vanish at the edges.

For Hamming: the same three terms, but with coefficients $0.54 + 0.23 + 0.23 = 1.0$. The coefficients are tuned to cancel a specific side-lobe peak, not to cancel the $1/\omega$ term of the asymptotic expansion. The $1/\omega$ term survives because $w[0] = 0.08 \neq 0$. Result: the nearest side lobes are lower (−42.7 dB vs. −31.5 dB), but the far side lobes decay at only 6 dB/oct instead of 18 dB/oct.

This is the precise mechanism behind the rolloff rule from Section A.3.4: the rolloff rate is determined by edge smoothness, which in the pure sine form manifests as the degree of asymptotic cancellation among the $1/\sin$ terms.

#### Comparative graph (linear scale)

Figure B.13 overlays the normalized spectra $D(\omega)$ of all four windows on one plot in linear scale, zoomed to the side-lobe region ($|D| \leq 0.25$). Figure B.14 shows the same comparison at full scale.

![Figure B.13 — Window spectra comparison, linear scale (zoomed)](../results/graphs/lab3/figure_B_13.png)

![Figure B.14 — Window spectra comparison, full range](../results/graphs/lab3/figure_B_14.png)

The tradeoff is immediately visible:

- **Rectangular** has the narrowest main lobe (2 bins) but the tallest side lobes (≈ 0.22).
- **Hann** and **Hamming** have equal main-lobe widths (4 bins). Hamming's side lobes are lower in the nearest lobes but do not decay as fast.
- **Blackman** has the widest main lobe (6 bins) but side lobes so small they are barely visible on the linear scale.

This is the main-lobe/side-lobe tradeoff from Section A.3.4, now derived from first principles: each window is a weighted sum of shifted Dirichlet kernels, and the weights determine how much cancellation occurs in the side-lobe region.

**Summary of derived properties:**

**Table B.6 — Window properties derived in Lab 3**

| Window | Formula (shifted Dirichlet kernels) | Main-lobe width (bins) | Peak side-lobe (dB) | Rolloff | Edge value |
| --- | --- | --- | --- | --- | --- |
| Rectangular | $D(\omega)$ | 2 | −13.3 | 6 dB/oct | 1 |
| Hann | $0.5D - 0.25D_{\pm 1}$ | 4 | −31.5 | 18 dB/oct | 0 |
| Hamming | $0.54D - 0.23D_{\pm 1}$ | 4 | −42.7 | 6 dB/oct | 0.08 |
| Blackman | $0.42D - 0.25D_{\pm 1} + 0.04D_{\pm 2}$ | 6 | −58 | 18 dB/oct | 0 |

These are the numbers that Lab 3 will confirm experimentally with the two-tone resolution test.

---

## B.3 — Lab 3: Two-Tone Resolution Test  *(↔ A.3, Appendix B)*

### Introduction

Appendix B derived the window spectra from first principles and established their main-lobe widths (Table B.6). Section A.3.5 stated the resolution limit $\Delta f_{\min} \approx \beta \cdot f_s/N$, where $\beta$ is the main-lobe half-width in bins. This lab confirms the practical consequence: can two tones be distinguished at the predicted separation?

### Setup

Two tones at $f_1 = 10$ Hz and $f_2 = 10 + \Delta$ Hz in a 5-second segment ($N = 1250$, $\Delta f = 0.20$ Hz). The separation $\Delta$ is swept from 0.05 Hz to 3.0 Hz. For each separation and each window, the DFT is computed and a dip-detection criterion determines resolvability: two tones are resolved if the valley between their peaks falls below 80% of the average peak height.

```python
import numpy as np
from src.common import FS, make_mixed_tones
from src.common.windows import rectangular, hann, blackman

test_dur = 5.0                                            # 5-second segment
N_test = int(test_dur * FS)                               # 1250 samples
delta_f = FS / N_test                                     # 0.20 Hz bin spacing

for sep in np.arange(0.05, 3.0, 0.05):                   # sweep separations
    x, _, _ = make_mixed_tones([10.0, 10.0 + sep], duration=test_dur)

    for name, win_func, beta in [
        ("Rectangular", rectangular, 1),
        ("Hann",        hann,        2),
        ("Blackman",    blackman,    3),
    ]:
        w = win_func(N_test)                              # generate window
        x_w = x * w                                       # apply window
        X_w = np.fft.fft(x_w)                             # DFT
        P = np.abs(X_w)**2                                # power

        # --- Resolvability: is there a dip between the two peaks? ---
        freqs = np.fft.fftfreq(N_test, d=1/FS)            # frequency axis
        pos = freqs >= 0                                  # positive half
        f_pos, P_pos = freqs[pos], P[pos]                 # positive frequencies
        band = (f_pos >= 10 - 1) & (f_pos <= 10 + sep + 1)  # region around tones
        P_band, f_band = P_pos[band], f_pos[band]        # power in region
        peak1 = np.argmin(np.abs(f_band - 10.0))         # bin nearest f1
        peak2 = np.argmin(np.abs(f_band - (10.0 + sep))) # bin nearest f2
        valley = np.min(P_band[peak1:peak2+1])            # minimum between peaks
        peak_avg = (P_band[peak1] + P_band[peak2]) / 2   # average peak height
        resolved = valley < 0.8 * peak_avg                # dip below 80%?
```

Full source: `src/lab3_windows/lab3.py`.

### Parameters

**Table B.3 — Lab 3 parameters**

| Parameter | Value |
| --- | --- |
| $f_s$ (Hz) | 250 |
| Duration (s) | 5.0 |
| $N$ (samples) | 1250 |
| $\Delta f$ (Hz) | 0.20 |
| $f_1$ (Hz) | 10.0 |
| $f_2$ (Hz) | $10.0 + \Delta$ (sweep) |
| Separation sweep $\Delta$ (Hz) | 0.05 to 3.0, step 0.05 |
| Windows tested | Rectangular, Hann, Blackman |
| Resolvability criterion | Valley < 80% of average peak |

### Results

Figure B.11 shows resolvability vs. tone separation for each window. The step-function transition from "not resolved" to "resolved" matches the predicted $\Delta f_{\min} = \beta \cdot f_s/N$:

| Window | $\beta$ (from Table B.6) | Predicted $\Delta f_{\min}$ (Hz) | Measured transition (Hz) |
| --- | --- | --- | --- |
| Rectangular | 1 | $1 \times 0.20 = 0.20$ | ≈ 0.20 |
| Hann | 2 | $2 \times 0.20 = 0.40$ | ≈ 0.40 |
| Blackman | 3 | $3 \times 0.20 = 0.60$ | ≈ 0.60 |

### Verification

| Prediction (Appendix B / Eq. A.21) | Measured | Confirmed? |
| --- | --- | --- |
| Rectangular: $\Delta f_{\min} = 0.20$ Hz ($\beta = 1$) | ≈ 0.20 Hz | Yes |
| Hann: $\Delta f_{\min} = 0.40$ Hz ($\beta = 2$) | ≈ 0.40 Hz | Yes |
| Blackman: $\Delta f_{\min} = 0.60$ Hz ($\beta = 3$) | ≈ 0.60 Hz | Yes |

### Conclusion

The resolution limit $\Delta f_{\min} \approx \beta \cdot f_s / N$ from Appendix B is confirmed experimentally. The main-lobe widths derived from the Dirichlet kernel decomposition (Table B.6) correctly predict the minimum separation at which two tones are distinguishable.

For EEG at $f_s = 250$ Hz: a 1-second window ($M = 256$) gives $\Delta f = 0.977$ Hz. With a Hann window ($\beta = 2$), $\Delta f_{\min} = 1.95$ Hz — sufficient to separate the standard EEG bands (narrowest gap: δ–θ at 4 Hz). With 1200-second lab signals, the resolution margin is enormous (Table A.2). The window choice for EEG is driven by side-lobe suppression, not resolution — and Appendix B has derived exactly why.

---

---

## B.4 — Lab 4: The STFT of a Fluctuating Signal  *(↔ A.5)*

### Introduction

The DFT (Lab 1) gives frequency content but discards all timing. Welch's method (Lab 2) averages over time explicitly. For any signal whose frequency content changes — an EEG rhythm that comes and goes, a chirp that sweeps — we need both axes simultaneously. The STFT (Section A.5) provides this: it slides a windowed DFT across time and keeps each segment's spectrum indexed by position. The result is the **spectrogram** — the first usable time-frequency representation in this report.

This lab tests:

- The STFT definition and the spectrogram as a 2D time-frequency plot (A.5.1).
- The **Heisenberg uncertainty tradeoff**: $\Delta t \cdot \Delta f = \beta$ is constant; window length $M$ slides between time and frequency resolution (A.5.2).
- **Overlap and the tapering problem**: how the window's taper suppresses edge samples, and how overlap (COLA) recovers them (A.5.3).
- The **multi-scale limitation**: no single STFT window can capture both a narrow-band chirp and a short transient simultaneously.

### Setup

**Model signal.** A linear chirp sweeping from $f_0 = 5$ Hz to $f_1 = 45$ Hz over 120 seconds (Equation (AA.3)):

$$
x[n] = A \cos\!\left(\frac{2\pi}{f_s}\left(f_0 n + \frac{\mu}{2}\frac{n^2}{f_s}\right)\right) \tag{B.28}
$$

with chirp rate $\mu = (f_1 - f_0)/T = 40/120 = 0.333$ Hz/s. The instantaneous frequency at time $t$ is $f_{\text{inst}}(t) = 5 + 0.333\,t$ Hz. The signal duration is 120 s (shorter than the 1200 s lab default so the chirp diagonal is visible at all zoom levels).

For Experiments B and C, an alpha burst (10 Hz Gaussian-enveloped transient, Equation (AA.6)) is added at $t = 60$ s:

$$
x_{\text{burst}}[n] = A_b \exp\!\left(-\frac{(n - n_0)^2}{2\sigma_t^2}\right) \cos\!\left(\frac{2\pi f_b}{f_s} n\right) \tag{B.29}
$$

with $A_b = 3$, $f_b = 10$ Hz, $\sigma_t = 125$ samples (0.5 s).

**STFT computation:**

```python
import numpy as np
from scipy import signal as sp_signal
from src.common import FS, make_chirp, make_transient, plot_time_domain, save_figure

def compute_stft(x, fs=FS, nperseg=256, noverlap=None, window="hann"):
    if noverlap is None:
        noverlap = nperseg // 2                           # default 50% overlap
    f_stft, t_stft, Sxx = sp_signal.spectrogram(
        x, fs=fs, nperseg=nperseg,                        # window length
        noverlap=noverlap,                                # overlap
        window=window,                                    # window type
    )
    return t_stft, f_stft, Sxx
```

**Experiment A — Heisenberg tradeoff.** The chirp is analyzed with Hann windows at four lengths ($M = 125, 250, 500, 1250$ samples), all with 50% overlap.

```python
F0_CHIRP = 5.0                                            # start frequency (Hz)
F1_CHIRP = 45.0                                           # end frequency (Hz)
DURATION_LAB = 120.0                                      # signal duration (s)
MU_CHIRP = (F1_CHIRP - F0_CHIRP) / DURATION_LAB           # 0.333 Hz/s

x, n, t = make_chirp(F0_CHIRP, MU_CHIRP, A=1.0, duration=DURATION_LAB)

for M in [125, 250, 500, 1250]:                           # window length sweep
    dt = M / FS                                           # time resolution (s)
    df = 2 * FS / M                                       # freq resolution (Hz), β=2
    noverlap = M // 2                                     # 50% overlap (Hann COLA)
    t_stft, f_stft, Sxx = compute_stft(x, nperseg=M, noverlap=noverlap)
    # ... plot dual-stack spectrogram
```

**Experiment B — Overlap and tapering.** The chirp + burst signal is analyzed with a fixed $M = 256$ at three overlap settings (0%, 50%, 75%), zoomed to the burst region (55–65 s). White dashed lines mark the true burst extent (±2σ).

```python
x_chirp, _, _ = make_chirp(F0_CHIRP, MU_CHIRP, A=1.0, duration=DURATION_LAB)
n0 = int(60.0 * FS)                                      # burst at t=60 s
sigma_t = int(0.5 * FS)                                   # σ = 0.5 s
x_burst, _, _ = make_transient(n0, sigma_t, f0=10.0, A=3.0, duration=DURATION_LAB)
x = x_chirp + x_burst                                    # combined signal

M = 256                                                   # fixed window
for overlap_frac in [0.0, 0.5, 0.75]:                    # overlap sweep
    noverlap = int(M * overlap_frac)
    t_stft, f_stft, Sxx = compute_stft(x, nperseg=M, noverlap=noverlap)
    # ... plot with burst reference lines at t = 59 s and t = 61 s
```

**Experiment C — Multi-scale limitation.** The same chirp + burst signal, analyzed with a short window ($M = 125$, 0.5 s) and a long window ($M = 1250$, 5 s), zoomed to the burst region. White dashed lines mark the true burst extent.

Full source: `src/lab4_stft/lab4.py`.

### Parameters

**Table B.7 — Lab 4 parameters**

| Parameter | Experiment A | Experiment B | Experiment C |
| --- | --- | --- | --- |
| $f_s$ (Hz) | 250 | 250 | 250 |
| Duration (s) | 120 | 120 | 120 |
| Chirp $f_0$ (Hz) | 5.0 | 5.0 | 5.0 |
| Chirp $f_1$ (Hz) | 45.0 | 45.0 | 45.0 |
| Chirp rate $\mu$ (Hz/s) | 0.333 | 0.333 | 0.333 |
| Burst center (s) | — | 60.0 | 60.0 |
| Burst $\sigma_t$ (s) | — | 0.5 | 0.5 |
| Burst freq (Hz) | — | 10.0 | 10.0 |
| Burst amplitude | — | 3.0 | 3.0 |
| Window | Hann | Hann | Hann |
| Window lengths $M$ | 125, 250, 500, 1250 | 256 | 125, 1250 |
| Overlap | 50% | 0%, 50%, 75% | 50% |

### Results

**Experiment A — Heisenberg tradeoff.**

Figure B.15 shows the chirp in the time domain (first 5 seconds).

![Figure B.15 — Linear chirp 5→45 Hz, time domain](../results/graphs/lab4/figure_B_15.png)

Figures B.16–B.19 show the spectrogram of the same chirp at four window lengths. The diagonal sweeps from 5 Hz to 45 Hz over 120 seconds. The key observation: the diagonal's **thickness** changes with $M$, but the $\Delta t \cdot \Delta f$ product is constant at $\beta = 2$.

![Figure B.16 — M=125 (0.5 s): thick diagonal, Δf=4.0 Hz](../results/graphs/lab4/figure_B_16.png)

![Figure B.17 — M=250 (1.0 s): moderate thickness, Δf=2.0 Hz](../results/graphs/lab4/figure_B_17.png)

![Figure B.18 — M=500 (2.0 s): thin diagonal, Δf=1.0 Hz](../results/graphs/lab4/figure_B_18.png)

![Figure B.19 — M=1250 (5.0 s): very thin diagonal, Δf=0.4 Hz, but staircase steps in time](../results/graphs/lab4/figure_B_19.png)

| Window $M$ (samples) | $\Delta t$ (s) | $\Delta f$ (Hz) | $\Delta t \cdot \Delta f$ | Diagonal appearance |
| --- | --- | --- | --- | --- |
| 125 (0.5 s) | 0.50 | 4.00 | 2.0 | Thick, fuzzy — good time steps |
| 250 (1.0 s) | 1.00 | 2.00 | 2.0 | Moderate thickness |
| 500 (2.0 s) | 2.00 | 1.00 | 2.0 | Thin — frequency well resolved |
| 1250 (5.0 s) | 5.00 | 0.40 | 2.0 | Very thin, but 5 s time steps |

Every row has $\Delta t \cdot \Delta f = 2.0$ (Hann's $\beta$). The slider moves, the area doesn't shrink. This is the uncertainty principle (Equation (A.40)) made visible.

**Experiment B — Overlap and tapering.**

Figures B.20–B.22 show the same chirp + burst signal analyzed with $M = 256$ at three overlap levels, zoomed to the burst region. White dashed lines mark the true burst extent (±2σ = 59–61 s).

![Figure B.20 — 0% overlap: burst visible but gaps from tapering](../results/graphs/lab4/figure_B_20.png)

![Figure B.21 — 50% overlap: COLA satisfied, uniform coverage](../results/graphs/lab4/figure_B_21.png)

![Figure B.22 — 75% overlap: smoother but no sharper than 50%](../results/graphs/lab4/figure_B_22.png)

| Overlap | Hop $H$ | Columns | Segments per sample | Observation |
| --- | --- | --- | --- | --- |
| 0% | 256 | 117 | 1.0 | Burst shows gaps; edge samples suppressed by taper |
| 50% | 128 | 233 | 2.0 | COLA satisfied; burst cleanly captured within reference lines |
| 75% | 64 | 465 | 4.0 | Smoother time axis, but no finer resolution than 50% |

At 0% overlap, the Hann window multiplies edge samples by zero — features at segment boundaries can be lost (Section A.5.3). At 50%, the COLA condition (Equation (A.42)) is satisfied: every sample receives equal total weight, and the burst is captured completely. At 75%, the spectrogram has more columns (finer time grid) but no additional resolution — the same distinction as zero-padding (Section A.2.3), now on the time axis.

**Experiment C — Multi-scale limitation.**

Figure B.23 shows the chirp + alpha burst in the time domain, zoomed to the burst region.

![Figure B.23 — Chirp + alpha burst, time domain (zoomed)](../results/graphs/lab4/figure_B_23.png)

Figures B.24–B.25 show the same signal analyzed with a short window and a long window. White dashed lines mark the true burst extent.

![Figure B.24 — Short window M=125 (0.5 s): burst localized, chirp smeared](../results/graphs/lab4/figure_B_24.png)

![Figure B.25 — Long window M=1250 (5.0 s): chirp sharp, burst smeared far beyond true extent](../results/graphs/lab4/figure_B_25.png)

| Window | $\Delta t$ (s) | $\Delta f$ (Hz) | Burst | Chirp |
| --- | --- | --- | --- | --- |
| Short ($M = 125$) | 0.50 | 4.00 | Localized within reference lines | Smeared into a broad band |
| Long ($M = 1250$) | 5.00 | 0.40 | Smeared far beyond reference lines | Sharp, thin diagonal |

The short window captures the burst correctly (energy stays within the white dashed lines) but smears the chirp into a thick band. The long window sharpens the chirp into a thin diagonal but smears the burst across the entire 10-second view — far beyond its true 2-second extent. No single window captures both: the burst needs $\Delta t \leq 1$ s, the chirp needs $\Delta f \leq 1$ Hz, but $\Delta t \cdot \Delta f = 2$ means you cannot have both simultaneously.

### Verification

| Prediction (Volume A) | Measured | Confirmed? |
| --- | --- | --- |
| $\Delta t \cdot \Delta f = \beta = 2$ for Hann (Eq. A.39) | 2.0 at all four $M$ values | Yes |
| COLA at 50% overlap for Hann (Eq. A.42) | Uniform coverage, no gaps | Yes |
| Overlap beyond COLA adds columns, not resolution (A.5.3) | 75% smoother but not sharper than 50% | Yes |
| Multi-scale signals cannot be captured by a single $M$ | Short $M$ smears chirp, long $M$ smears burst | Yes |

### Conclusion

The STFT is the first tool in this report that answers "what frequency is present at what time." The Heisenberg tradeoff is real and inescapable: $\Delta t \cdot \Delta f = \beta$ is constant across all window lengths (Figures B.16–B.19). Overlap solves the tapering problem (Figure B.20 vs B.21), and the COLA condition guarantees uniform sample coverage.

The multi-scale experiment (Figures B.24–B.25) reveals the STFT's fundamental limitation: it forces a single choice of $M$ for the entire signal. When the signal contains features at different time-frequency scales — a narrow-band rhythm and a short transient, as EEG often does — no single window captures both. This limitation motivates the Wigner-Ville Distribution (Lab 7), which is not bound by the uncertainty principle in the same way.

---

*Next: B.5 — Two-Tone Resolution on the Spectrogram. The resolution limit $\Delta f_{\min} \approx \beta \cdot f_s / M$ from Lab 3 is confirmed visually: two stationary tones on a spectrogram — too close and they merge into one line, far enough and they split into two.*
