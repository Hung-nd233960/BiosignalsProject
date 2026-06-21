# From the DFT to the SPWVD: Time-Frequency Analysis Applied to Neonatal EEG

## Volume A - Theoretical Background

**Author:** Nguyen Duc Hung - 20233960

### What this report is about

This report studies the fundamentals of digital signal processing in the frequency and time-frequency domains, starting from the Discrete Fourier Transform (DFT) and building through the Short-Time Fourier Transform (STFT), autocorrelation, and the Wigner-Ville Distribution (WVD) to the Smoothed Pseudo Wigner-Ville Distribution (SPWVD). The theory is developed in Volume A, validated on model signals in Volume B (8 labs + appendices), and applied to a real neonatal EEG recording in Volume C.

### What this volume covers

Volume A develops the theoretical foundation used throughout the report:

- **A.1** Signal theory: sampling, discrete frequency, Nyquist, signal classification, energy/power, the signal as a vector
- **A.2** The DFT: bins, resolution, zero-padding, Parseval's theorem
- **A.3** Windowing: leakage, the cosine-sum family, the resolution limit, the M vs M-1 convention
- **A.4** Statistics: bin distributions under noise, the noise floor, spectral detection, Welch's method
- **A.5** The STFT: the spectrogram, the Heisenberg uncertainty principle, overlap and COLA
- **A.6** Autocorrelation: periodicity detection, Wiener-Khinchin, phase-blindness, cross-correlation
- **A.7** The WVD: instantaneous autocorrelation, chirp sharpness, cross-terms, the analytic signal
- **A.8** The PWVD and SPWVD: lag window, time window, the two-knob tradeoff, Cohen's class
- **Appendix A** Signal taxonomy: the six archetypes and their behaviour under each transform

All formulas are discrete-time and discrete-frequency. The continuous-time form is mentioned for motivation where necessary, but the working formula is always the discrete version.

## A.1 Signal Theory

### A.1.1 What a Discrete Signal Means

A continuous physical signal - the voltage on a scalp electrode, the pressure in a microphone - exists for all time $t$. We cannot store a continuous function. Sampling produces a discrete signal by recording the value at regularly spaced instants:

$$
x[n] = x_c(n / f_s), \qquad n = 0, 1, \ldots, N-1 \tag{A.1}
$$

where $x_c(t)$ is the continuous original, $f_s$ is the sampling rate in Hz (samples per second), and $N$ is the total number of samples collected.

What is **kept**: the signal's value at each sample point. If the signal is band-limited (no frequency content above $f_s/2$), these values are sufficient to reconstruct the continuous original exactly (Shannon's sampling theorem).

What is **lost**: everything between samples. If the signal contains energy above $f_s/2$ - which real signals always do, to some extent - that energy folds back into the representable band as aliasing (Section A.1.3).

The notation: $x[n]$ with square brackets denotes a discrete signal. The index $n$ is a dimensionless integer. Physical time is recovered as $t_n = n / f_s$ seconds.

### A.1.2 Discrete Frequency and the Unit Circle

Sampling discretizes time. What happens to frequency?

A continuous sinusoid $\cos(2\pi f t)$ sampled at rate $f_s$ becomes:

$$
x[n] = \cos(2\pi f \cdot n / f_s) = \cos(\omega n) \tag{A.2}
$$

where $\omega = 2\pi f / f_s$ is the **normalized angular frequency** in radians per sample. This is the natural frequency unit for discrete signals: it measures how much phase the sinusoid accumulates per sample.

The complex exponential $e^{j\omega n}$ traces the unit circle in the complex plane. One full revolution ($\omega = 2\pi$) corresponds to one cycle per sample, which is frequency $f = f_s$. Half a revolution ($\omega = \pi$) corresponds to $f = f_s/2$.

A critical consequence: frequencies that differ by $f_s$ produce identical samples:

$$
\cos(2\pi f \cdot n / f_s) = \cos(2\pi (f + f_s) \cdot n / f_s) \tag{A.3}
$$

because $\cos(\omega n + 2\pi n) = \cos(\omega n)$ for integer $n$. The frequency axis is **periodic** with period $f_s$. We only need to examine $f \in [0, f_s)$ - or equivalently $\omega \in [0, 2\pi)$ - to see everything the discrete signal can represent.

### A.1.3 Nyquist and Aliasing

From Equation (A.3), frequencies $f$ and $f_s - f$ produce the same samples when viewed as cosines: $\cos(2\pi f n/f_s) = \cos(2\pi(f_s - f)n/f_s)$ (since cosine is even and $\cos(2\pi n - \theta) = \cos(\theta)$). They are **aliases** of each other - indistinguishable after sampling.

The **Nyquist frequency** is the boundary:

$$
f_N = \frac{f_s}{2} \tag{A.4}
$$

Frequencies below $f_N$ have a unique representation. Frequencies above $f_N$ fold back: a tone at $f > f_N$ appears at its alias $f_s - f < f_N$. For example, an 80 Hz tone sampled at 200 Hz appears as a 120 Hz tone - wait, $200 - 80 = 120 > f_N = 100$. More precisely: the alias wraps modulo $f_s$, and the apparent frequency in $[0, f_N]$ is $|f \mod f_s - f_s \cdot \text{round}(f/f_s)|$. For 80 Hz at $f_s = 200$: the tone is below $f_N = 100$, so no aliasing. For 130 Hz at $f_s = 200$: the alias is $200 - 130 = 70$ Hz.

**Anti-aliasing in practice.** Before sampling, an analog low-pass filter removes energy above $f_N$. For EEG: the dataset used in Volume C has $f_s = 200$ Hz ($f_N = 100$ Hz) and is filtered at 50 Hz - a conservative cutoff well below Nyquist, providing a safety margin.

### A.1.4 Sampling Count vs. Sampling Duration

The sampling rate $f_s$ and the number of samples $N$ are two independent choices. They control different things:

| Knob | What it sets | What it does NOT set |
| --- | --- | --- |
| $f_s$ (sampling rate) | Nyquist frequency $f_N = f_s/2$ (bandwidth) | Frequency resolution |
| $N$ (sample count) | Duration $T = N/f_s$; frequency resolution $\Delta f = 1/T$ | Nyquist frequency |

Increasing $f_s$ extends the bandwidth (higher frequencies become representable) but does not improve the ability to distinguish two nearby frequencies. Increasing $N$ (at the same $f_s$) extends the observation time, which improves frequency resolution (derived in Section A.2.2) but does not change the Nyquist limit.

For the EEG dataset: $f_s = 200$ Hz gives bandwidth to 100 Hz. $N = 228\,000$ gives $T = 1140$ s (19 minutes) and $\Delta f = 1/1140 = 0.00088$ Hz. These are independent facts about the recording.

### A.1.5 Signal Classification

Not all signals behave the same way under analysis. Two independent classifications determine which tools are appropriate:

**Deterministic vs. random.**

- A **deterministic** signal is fully specified by a formula: $x[n] = A\cos(2\pi f_0 n / f_s + \phi)$. Given the parameters, every sample is known in advance. Its DFT is exact and reproducible.
- A **random** (stochastic) signal is drawn from a probability distribution. Each realization is different. White noise $\eta[n] \sim \mathcal{N}(0, \sigma^2)$ is the canonical example. Its DFT is a random variable - each bin has a magnitude and phase drawn from a distribution, not a fixed value. This is why spectral analysis requires statistics (Section A.4).

**Stationary vs. non-stationary.**

- A **stationary** signal has statistical properties that do not change over time. A constant-frequency tone, or white noise with fixed variance, are stationary. Any segment of the signal looks statistically the same as any other.
- A **non-stationary** signal changes over time. A chirp (frequency sweeps), an EEG rhythm that comes and goes, a transient artifact - these have time-varying frequency content. The global DFT averages over the entire duration and loses the timing information. This is why the STFT (Section A.5) and WVD (Section A.7) exist.

These two axes are independent. The six signal archetypes from Appendix A map onto them:

**Table A.1 - Signal classification and appropriate tools**

| Signal type | Example archetypes | Primary tool |
| --- | --- | --- |
| Deterministic + stationary | Single tone, mixed tones | DFT (A.2) |
| Deterministic + non-stationary | Chirp, transient | STFT (A.5), WVD (A.7) |
| Random + stationary | White noise | DFT + statistics (A.4) |
| Random + non-stationary | EEG | Full toolbox: STFT + Welch + WVD/SPWVD |

EEG sits in the hardest quadrant: random and non-stationary. This is why the report builds the full chain of tools - each one handles a piece of the problem.

### A.1.6 Energy and Power

Before we can discuss "how much energy is at 10 Hz," we need to define energy in the discrete domain.

**Signal energy** is the sum of squared magnitudes:

$$
E = \sum_{n=0}^{N-1} |x[n]|^2 \tag{A.5}
$$

For a real signal in uV, $E$ has units of uV$^2$ $\cdot$ samples. It measures the total "weight" of the signal.

**Signal power** is energy per sample:

$$
P = \frac{E}{N} = \frac{1}{N}\sum_{n=0}^{N-1} |x[n]|^2 \tag{A.6}
$$

Power is normalized so it does not grow with signal length. A 10-second recording of a tone has more energy than a 1-second recording, but the same power.

These definitions connect forward:

- **Parseval's theorem** (Section A.2.4) states that the total energy $E$ can be computed equivalently from the DFT coefficients: $E = \frac{1}{N}\sum_k |X[k]|^2$. This requires $E$ to be defined first.
- **Power spectral density** (Section A.4) distributes the power across frequency bins. Each bin's PSD value has units of uV$^2$/Hz.
- **Band-power analysis** (Volume C) integrates PSD over an EEG frequency band to get power in uV$^2$. The chain starts here.

### A.1.7 The Signal as a Vector

A length-$N$ discrete signal $x[n]$ can be written as a vector:

$$
\mathbf{x} = (x[0], \; x[1], \; \ldots, \; x[N-1]) \in \mathbb{R}^N \qquad \text{(A.7)}
$$

This is not a metaphor. The signal literally is a point in $N$-dimensional space. Each sample $x[n]$ is one coordinate.

**Inner product.** The similarity between two signals is measured by:

$$
\langle \mathbf{x}, \mathbf{y} \rangle = \sum_{n=0}^{N-1} x[n] \, y^*[n] \tag{A.8}
$$

where $y^*$ denotes the complex conjugate (for real signals, $y^* = y$). When $\langle \mathbf{x}, \mathbf{y} \rangle = 0$, the signals are **orthogonal** - they share no structure, they are independent directions in $\mathbb{R}^N$.

**Norm.** The length of the signal vector is:

$$
\|\mathbf{x}\| = \sqrt{\langle \mathbf{x}, \mathbf{x} \rangle} = \sqrt{\sum_{n=0}^{N-1} |x[n]|^2} = \sqrt{E} \tag{A.9}
$$

The norm squared is the energy from Equation (A.5). Energy is the squared length of the signal vector.

**Why this matters.** The DFT basis vectors (Section A.2.4) are orthonormal: they have unit norm and are mutually orthogonal. The DFT is a **rotation** in $\mathbb{R}^N$ - a change from the sample-value coordinates to the frequency coordinates. Rotations preserve length, so $\|\mathbf{x}\|$ is the same in both bases. That is Parseval's theorem: energy conservation is not a special property of the DFT, it is a geometric consequence of orthonormality.

## A.2 The DTFT and the DFT

### A.2.1 From the DTFT to the DFT

We begin with what the reader already knows. Given a finite-length discrete-time signal $x[n]$ of $N$ samples, the **Discrete-Time Fourier Transform** is defined as:

$$
X(\omega) = \sum_{n=0}^{N-1} x[n] \, e^{-j\omega n} \tag{A.10}
$$

where $\omega$ is the **normalized angular frequency** in radians per sample, and $X(\omega)$ is a continuous, periodic function of $\omega$ with period $2\pi$.

Equation (A.10) is the formula from a Signals & Systems course. It is correct, complete, and - in practice - unusable.

The problem is not on the time side. We already discretized time when we sampled: the signal is a finite list of $N$ numbers. The problem is on the **frequency side**: $X(\omega)$ is defined for every real value of $\omega$ in $[0, 2\pi)$. That is an uncountably infinite set of values. No computer stores a continuous function. No plot has infinite resolution. To actually *use* the frequency-domain representation, we must **sample it**.

The question is: at which frequencies do we evaluate $X(\omega)$?

For a signal whose frequency content is known in advance - say, a single sinusoid at exactly 10 Hz - we could evaluate $X(\omega)$ at just that one frequency and be done. But for a **general signal** - one whose spectral content we do not know beforehand - we have no such shortcut. We need a **systematic grid** of frequency points that covers the full band, so that no feature is missed regardless of where it falls.

The simplest such grid is $N$ equally spaced points around the unit circle:

$$
\omega_k = \frac{2\pi k}{N}, \qquad k = 0, 1, \ldots, N-1 \tag{A.11}
$$

Substituting Equation (A.11) into the DTFT (A.10) gives the **Discrete Fourier Transform**:

$$
X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j 2\pi k n / N}, \qquad k = 0, 1, \ldots, N-1 \tag{A.12}
$$

Equation (A.12) is the DTFT evaluated at $N$ specific frequencies. Nothing more, nothing less. It does not introduce new mathematics - it inherits everything from the DTFT. What it introduces is a **finite, computable representation** of the frequency content.

The inverse transform in Equation (A.13) confirms the relationship is exact - $N$ samples in time, $N$ samples in frequency, no information lost:

$$
x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \, e^{j 2\pi k n / N}, \qquad n = 0, 1, \ldots, N-1 \tag{A.13}
$$

The DTFT and the DFT are not two different transforms. The DFT is the DTFT made computable by choosing where to look.

### A.2.2 What a Bin Is

Each value $X[k]$ in the DFT output is called a **bin**. The word is worth understanding precisely, because much of the confusion in spectral analysis comes from treating bins as something they are not.

**A bin is one sample of the continuous DTFT, taken at a specific frequency.**

The DTFT $X(\omega)$ is a continuous curve. The DFT samples this curve at $N$ points. Each sample is a bin. Between the bins, the DTFT still exists - it has values there - but the DFT does not report them.

**Bin spacing in normalized frequency:**

The distance between adjacent bins on the $\omega$-axis is:

$$
\Delta\omega = \frac{2\pi}{N} \tag{A.14}
$$

**Bin spacing in physical frequency:**

To convert to Hz, recall that $\omega = 2\pi f / f_s$, where $f$ is the physical frequency in Hz and $f_s$ is the sampling rate. Applying this to Equation (A.14), the spacing between adjacent bins in Hz is:

$$
\Delta f = \frac{f_s}{N} \tag{A.15}
$$

Equation (A.15) is one of the most important numbers in spectral analysis. It depends on **two quantities only**: the sampling rate $f_s$ and the number of samples $N$. Since the signal duration is $T = N / f_s$, we can equivalently write Equation (A.15) as:

$$
\Delta f = \frac{1}{T} \tag{A.16}
$$

The bin spacing is the reciprocal of the observation time. A 1-second signal gives 1 Hz bin spacing. A 20-second signal gives 0.05 Hz bin spacing. This is not a design choice - it is a consequence of sampling the DTFT at $N$ equally spaced points.

**Which physical frequency does bin $k$ represent?**

$$
f_k = k \cdot \Delta f = \frac{k \cdot f_s}{N} \tag{A.17}
$$

For $k = 0$: the DC component (zero frequency).
For $k = N/2$ (if $N$ is even): the Nyquist frequency $f_s / 2$.
For $k > N/2$: these bins correspond to negative frequencies, aliased into the range $(f_s/2, \, f_s)$. For real-valued signals, they are conjugate mirrors of the bins below $N/2$ and carry no new information.

### A.2.3 Resolution vs. Bin Count

These two concepts are routinely conflated. They are not the same thing.

**Bin count** is the number of DFT output values: $N$ bins from a length-$N$ signal. It determines how densely the DTFT is sampled.

**Frequency resolution** is the ability to distinguish two nearby frequency components in the signal. It is determined by the **signal duration** $T = N / f_s$, and its fundamental limit is given by Equation (A.18):

$$
\Delta f_{\min} = \frac{1}{T} = \frac{f_s}{N} \tag{A.18}
$$

Two tones separated by less than $\Delta f_{\min}$ cannot be resolved - not because we lack bins, but because the DTFT itself merges them into a single lobe. The underlying continuous spectrum does not contain two distinct peaks.

**Zero-padding** is the operation that exposes this distinction. If we take our $N$-sample signal and append zeros to create a length-$M$ sequence ($M > N$), then compute a length-$M$ DFT:

- The **bin count** increases from $N$ to $M$. The DTFT is now sampled more finely.
- The **bin spacing** decreases from $f_s/N$ to $f_s/M$. We see more points on the same curve.
- The **frequency resolution does not change**. The DTFT curve itself - its shape, its lobes, its ability to separate two tones - is determined entirely by the original $N$ samples. Zero-padding samples the same curve more densely; it does not sharpen it.

Zero-padding is interpolation in the frequency domain: more points on the same curve. It makes plots smoother and peak locations easier to read, but it does not reveal detail that the signal duration did not capture.

To genuinely improve resolution, there is only one option: **record more signal** - increase $N$ at the same $f_s$ (equivalently, increase $T$).

### A.2.4 The DFT as an Orthonormal Basis; Parseval's Theorem

The DFT has a clean linear-algebraic interpretation. Define the $N$ basis vectors:

$$
\phi_k[n] = \frac{1}{\sqrt{N}} \, e^{j 2\pi k n / N}, \qquad n = 0, 1, \ldots, N-1 \tag{A.19}
$$

These vectors are **orthonormal**: for any two basis indices $k$ and $m$,

$$
\sum_{n=0}^{N-1} \phi_k[n] \, \phi_m^*[n] = \begin{cases} 1 & \text{if } k = m \\ 0 & \text{if } k \neq m \end{cases} \tag{A.20}
$$

The DFT decomposes $x[n]$ into coefficients on this basis. Each bin $X[k]$ measures how much of the signal aligns with the complex exponential at frequency $\omega_k$. The inverse DFT reconstructs the signal as a weighted sum of these basis vectors. This is a change of basis in $\mathbb{C}^N$ - reversible, lossless, and exact.

**Parseval's theorem** (Equation (A.21)) states that total energy is preserved across the transform:

$$
\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |X[k]|^2 \tag{A.21}
$$

The left side is the total energy computed in the time domain. The right side is the total energy computed from the frequency-domain coefficients, scaled by $1/N$. They are equal.

This is not approximate. It is a direct consequence of the orthonormality of the DFT basis - a unitary transform preserves inner products and, therefore, norms.

The practical consequence: by Equation (A.21), if we want to know how much energy a signal carries in a particular frequency band, we can sum $|X[k]|^2 / N$ over the bins in that band. The answer is exact, not an estimate. This is the foundation of **band-power analysis**, which Volume C will apply to EEG frequency bands (δ, θ, α, β, γ).

## A.3 Leakage and Windowing

### A.3.1 The Hidden Rectangular Window

The DFT in Equation (A.12) sums from $n = 0$ to $n = N - 1$. This means any signal that exists outside that range is silently set to zero. We can make this explicit by defining a **rectangular window**:

$$
w_{\text{rect}}[n] = \begin{cases} 1 & \text{if } 0 \leq n \leq N-1 \\ 0 & \text{otherwise} \end{cases} \tag{A.22}
$$

The DFT of a signal $x[n]$ is actually the DFT of the product $x[n] \cdot w_{\text{rect}}[n]$. We never asked for this multiplication - it is a consequence of having a finite number of samples.

The rectangular window has a sharp edge: at $n = 0$ the signal abruptly begins, and at $n = N - 1$ it abruptly ends. If the signal is not perfectly periodic within the $N$-sample frame - meaning if the last sample does not connect smoothly back to the first - then this abrupt cut creates a discontinuity at the boundary. The DFT, which implicitly treats the signal as periodic, sees this discontinuity as part of the signal.

This is the origin of **spectral leakage**: the discontinuity introduces frequency content that was not present in the original signal.

### A.3.2 Leakage: The Window's Spectrum Convolved onto the Signal

Multiplication in the time domain is convolution in the frequency domain. The DTFT of the windowed signal $x[n] \cdot w[n]$ is:

$$
X_w(\omega) = \frac{1}{2\pi} \int_{-\pi}^{\pi} X(\theta) \, W(\omega - \theta) \, d\theta \tag{A.23}
$$

where $W(\omega)$ is the DTFT of the window itself. Each frequency component in $X(\theta)$ gets smeared by the shape of $W(\omega)$.

In discrete terms, this means: the DFT does not show the signal's true spectrum. It shows the signal's spectrum convolved with the window's spectrum. Every tone in the signal is replaced by a copy of $W(\omega)$ centered at that tone's frequency.

The DTFT of the rectangular window is the **Dirichlet kernel** (derived in Appendix B.1):

$$
W_{\text{rect}}(\omega) = e^{-j\omega(N-1)/2} \cdot \frac{\sin(\omega N / 2)}{\sin(\omega / 2)} \tag{A.24}
$$

Equation (A.24) has a tall main lobe centered at $\omega = 0$ and a series of side lobes that decay slowly on either side. These side lobes are the source of leakage: energy from a tone at one frequency spills into bins at other frequencies through the side lobes of $W_{\text{rect}}(\omega)$.

**On-bin vs. off-bin placement.** If a tone's frequency falls exactly on a DFT bin (i.e. $f = k \cdot f_s / N$ for integer $k$), the main lobe of $W(\omega)$ lands exactly on that bin and all other bins fall on the nulls of the Dirichlet kernel. Leakage vanishes - by coincidence of alignment, not by design. If the tone is even slightly off-bin, the nulls no longer align with the bins, and leakage appears across the entire spectrum. This sensitivity to placement is a fundamental fragility of the rectangular window.

### A.3.3 The Cosine-Sum Window Family

The goal of a non-rectangular window is to reduce side-lobe levels - accepting a wider main lobe in exchange. The most widely used family is the **cosine-sum** windows, defined by the general formula:

$$
w[n] = \sum_{p=0}^{P} (-1)^p \, a_p \cos\!\left(\frac{2\pi p n}{N}\right), \qquad n = 0, 1, \ldots, N-1 \tag{A.25}
$$

Each member of the family is determined by its coefficients $a_0, a_1, \ldots, a_P$.

**The $N$ vs $N-1$ convention.** Equation (A.25) divides by $N$ in the cosine argument. Many textbooks write $N-1$ instead:

$$
w_{\text{sym}}[n] = \sum_{p=0}^{P} (-1)^p \, a_p \cos\!\left(\frac{2\pi p n}{N-1}\right) \tag{A.16b}
$$

The two conventions exist because they serve different purposes.

**Symmetric** ($N-1$, Equation (A.16b)): the window is symmetric around its center, with $w[0] = w[N-1]$ exactly (for Hann: both are exactly zero). This is the form found in most textbooks because it is mathematically clean to define - a single arch from zero to zero - and is the natural choice for **filter design**, where the window is applied once and symmetry matters for linear-phase responses.

Used by:

- `numpy.hanning(N)` - docs state the formula as $w(n) = 0.5 - 0.5\cos(2\pi n / (M-1))$
- `scipy.signal.windows.hann(N, sym=True)` - default is `sym=True`
- **MATLAB** `hann(N)` and `hamming(N)` - default is `'symmetric'`. MATLAB's docs define $w(n) = 0.5(1 - \cos(2\pi n / N))$ with window length $L = N + 1$, which is equivalent to dividing by $M-1$.

**Periodic** ($N$, Equation (A.25)): the window is one period of a length-$N$ periodic sequence. At the edges, $w[0] = 0$ but $w[N-1] \neq 0$ (for Hann at $N = 256$: $w[N-1] \approx 0.00015$). This is equivalent to generating a length-$(N+1)$ symmetric window and dropping the last sample. This is the natural choice for **DFT-based spectral analysis**, and the reason is tied to how the DFT treats its input.

Used by:

- `scipy.signal.windows.hann(N, sym=False)`
- `scipy.signal.get_window('hann', N, fftbins=True)` - default is `fftbins=True`
- `scipy.signal.spectrogram` - passes string windows to `get_window` with `fftbins=True` ("DFT-even by default", per the scipy docs)
- **MATLAB** `hann(N, 'periodic')` and `hamming(N, 'periodic')` - generates $N+1$ symmetric points, drops the last. MATLAB's docs: "the periodic option is useful for spectral analysis because the DFT assumes periodicity."

**Why spectral analysis uses periodic.** The DFT implicitly treats its input as periodic - sample $N$ is the same as sample $0$. When a symmetric window is applied and the DFT repeats the windowed signal, the boundary looks like:

$$
\ldots, \; w[N{-}2] \cdot x[N{-}2], \; \underbrace{w[N{-}1] \cdot x[N{-}1]}_{= 0}, \; \underbrace{w[0] \cdot x[0]}_{= 0}, \; w[1] \cdot x[1], \; \ldots
$$

Two consecutive zeros at every period boundary - the last sample of one period and the first sample of the next are both multiplied by zero. One of those zeros is redundant: it wastes a sample without adding information.

The periodic convention avoids this. With $w[N{-}1] \approx 0.00015$ (not exactly zero), the periodic repetition gives:

$$
\ldots, \; w[N{-}2] \cdot x[N{-}2], \; \underbrace{w[N{-}1] \cdot x[N{-}1]}_{\approx 0}, \; \underbrace{w[0] \cdot x[0]}_{= 0}, \; w[1] \cdot x[1], \; \ldots
$$

Only one zero at the boundary. The window tiles cleanly without doubling.

This is why scipy's spectral functions (`spectrogram`, `welch`, `stft`) default to periodic windows - they call `get_window` with `fftbins=True`. The scipy docs call this "DFT-even." Meanwhile, `numpy.hanning` and scipy's default `sym=True` produce symmetric windows designed for filter applications, not DFT analysis.

**The practical difference is negligible.** At $N = 256$:

- Edge value difference: $0.00015$ (0.015% of peak).
- Peak side-lobe difference: $0.11$ dB.
- Maximum sample-by-sample difference: $0.0094$ (less than 1%).

No measurement or plot in this report would change if we switched conventions. The difference vanishes as $N$ grows.

**The mathematical difference matters.** The periodic convention ($N$) has a critical algebraic property: when the cosine-sum window's DFT is decomposed into shifted Dirichlet kernels (Lab 3), the shift by $p$ bins introduces a numerator:

$$
\sin\!\left(\frac{\omega M}{2} - p\pi\right) = (-1)^p \sin\!\left(\frac{\omega M}{2}\right)
$$

This factorization requires the cosine argument to divide by $N$ (so the shift $2\pi p / N$ produces exactly $p\pi$ in the numerator after multiplication by $N/2$). With $N-1$, the product $p\pi \cdot N/(N-1)$ is not an integer multiple of $\pi$, and the clean $(-1)^p$ extraction fails - the closed-form decomposition into shifted Dirichlet kernels cannot be carried out.

**Table A.2 - Library defaults for window conventions**

| Library / Function | Default | Symmetric ($N-1$) | Periodic ($N$) |
| --- | --- | --- | --- |
| `numpy.hanning(N)` | Symmetric | only option | - |
| `scipy.signal.windows.hann(N)` | Symmetric (`sym=True`) | `sym=True` | `sym=False` |
| `scipy.signal.get_window('hann', N)` | Periodic (`fftbins=True`) | `fftbins=False` | `fftbins=True` |
| `scipy.signal.spectrogram` | Periodic | - | via `get_window` |
| MATLAB `hann(N)` | Symmetric (`'symmetric'`) | `'symmetric'` | `'periodic'` |
| MATLAB `hamming(N)` | Symmetric (`'symmetric'`) | `'symmetric'` | `'periodic'` |

General-purpose functions default to symmetric; spectral analysis functions default to periodic. Both conventions converge as $N \to \infty$ (Appendix B).

**This report uses the periodic convention ($N$) throughout** - in all formulas, all code (`src/common/windows.py`), and all library calls. This matches scipy's spectral analysis defaults and enables the closed-form derivations in Lab 3. Lab 3 confirms numerically that both conventions produce the same window properties at $M = 256$ to within 0.11 dB.

Table A.3 lists the four standard members, their coefficients, and their spectral properties. The edge values in Table A.3 refer to the periodic convention.

**Table A.3 - Cosine-sum window family: coefficients and spectral properties**

| Property | Rectangular | Hann | Hamming | Blackman |
| --- | --- | --- | --- | --- |
| Order $P$ | 0 | 1 | 1 | 2 |
| $a_0$ | 1 | 0.5 | 0.54 | 0.42 |
| $a_1$ | - | 0.5 | 0.46 | 0.5 |
| $a_2$ | - | - | - | 0.08 |
| Edge value $w[0] = w[N{-}1]$ | 1 | 0 | 0.08 | 0 |
| Main-lobe width (bins) | 2 | 4 | 4 | 6 |
| $\beta$ (resolution factor) | 1 | 2 | 2 | 3 |
| Peak side-lobe (dB) | −13 | −31.5 | −42.7 | −58 |
| Side-lobe rolloff | 6 dB/oct | 18 dB/oct | 6 dB/oct | 18 dB/oct |
| Rolloff order | $1/\omega^1$ | $1/\omega^3$ | $1/\omega^1$ | $1/\omega^3$ |

The individual formulas, expanded from Equation (A.25), are:

**Rectangular** ($P = 0$):

$$
w_{\text{rect}}[n] = 1 \tag{A.26}
$$

**Hann** ($P = 1$):

$$
w_{\text{Hann}}[n] = 0.5 - 0.5\cos\!\left(\frac{2\pi n}{N}\right) \tag{A.27}
$$

**Hamming** ($P = 1$):

$$
w_{\text{Hamming}}[n] = 0.54 - 0.46\cos\!\left(\frac{2\pi n}{N}\right) \tag{A.28}
$$

**Blackman** ($P = 2$):

$$
w_{\text{Blackman}}[n] = 0.42 - 0.5\cos\!\left(\frac{2\pi n}{N}\right) + 0.08\cos\!\left(\frac{4\pi n}{N}\right) \tag{A.29}
$$

The design logic behind each window is visible in Table A.3:

- **Rectangular** is the baseline. No shaping, narrowest main lobe, but the worst side lobes (−13 dB) and the slowest rolloff ($1/\omega$). Every off-bin tone leaks broadly.
- **Hann** goes to zero at both edges ($w[0] = w[N{-}1] = 0$). This eliminates the boundary discontinuity in both the window value and its first derivative, producing an $18$ dB/oct rolloff. The cost is a main lobe twice as wide as rectangular.
- **Hamming** is optimized to minimize the peak side-lobe level (−42.7 dB) rather than the rolloff rate. It does not go to zero at the edges ($w[0] = 0.08$), so the boundary discontinuity remains in the value - hence the rolloff rate stays at $6$ dB/oct, same as rectangular. But the specific choice of $a_0 = 0.54$, $a_1 = 0.46$ places a cancellation that pushes the nearest side lobes far below those of Hann.
- **Blackman** uses a second cosine term ($P = 2$) and goes to zero at both edges. The result is the lowest peak side-lobe (−58 dB) with $18$ dB/oct rolloff, but the widest main lobe (6 bins). It trades the most resolution for the cleanest spectrum.

### A.3.4 Main-Lobe vs. Side-Lobe Tradeoff; the Rolloff Rule

Table A.3 shows the central tradeoff of window design: **narrower main lobe ↔ higher side lobes**. No window escapes this. A narrow main lobe preserves frequency resolution but lets leakage through; a wide main lobe suppresses leakage but blurs nearby tones together.

The **side-lobe rolloff rate** is governed by the smoothness of the window at its edges:

- If $w[n]$ has a **discontinuity in value** at the boundary (rectangular, Hamming), the rolloff is $1/\omega^1$ - i.e. 6 dB per octave.
- If $w[n]$ goes to zero at both edges and has a **continuous first derivative** (Hann, Blackman), the rolloff is $1/\omega^3$ - i.e. 18 dB per octave.
- In general, each additional order of continuity at the boundary adds $1/\omega^2$ to the rolloff (12 dB/oct per derivative that vanishes).

This is the $1/\omega^p$ rule referenced in the table of contents. It connects a time-domain property (edge smoothness) to a frequency-domain consequence (how fast distant leakage falls off). The rule explains why Hamming, despite its low peak side-lobe, has the same rolloff as rectangular: its edge value is $0.08$, not $0$, so the value discontinuity persists.

### A.3.5 The Resolution Limit

In Section A.2.3, we stated that the frequency resolution is $\Delta f_{\min} = f_s / N$. That was for the rectangular window. With a general window, the resolution limit becomes:

$$
\Delta f_{\min} \approx \beta \cdot \frac{f_s}{N} \tag{A.30}
$$

where $\beta$ is the **main-lobe width factor** - the half-width of the main lobe measured in bins. The values of $\beta$ for the standard windows are listed in Table A.3.

Equation (A.30) is the practical resolution limit. Two tones separated by less than $\beta \cdot f_s / N$ in frequency will merge into a single lobe regardless of any other processing. The rectangular window ($\beta = 1$) gives the best resolution; Blackman ($\beta = 3$) gives the worst. This is not a deficiency of Blackman - it is the price of its −58 dB side-lobe suppression.

For the EEG context: at $f_s = 250$ Hz with a $T = 20$ s window ($N = 5000$), the bin spacing is $\Delta f = 0.05$ Hz. The resolution limits are:

**Table A.4 - Resolution limits for standard windows at $f_s = 250$ Hz, $N = 5000$**

| Window | $\beta$ | $\Delta f_{\min}$ (Hz) | Can resolve adjacent EEG bands? |
| --- | --- | --- | --- |
| Rectangular | 1 | 0.05 | Yes - all standard bands |
| Hann | 2 | 0.10 | Yes - all standard bands |
| Hamming | 2 | 0.10 | Yes - all standard bands |
| Blackman | 3 | 0.15 | Yes - all standard bands |

Table A.4 shows that with a 20-second observation window, even the widest window (Blackman) resolves 0.15 Hz - far finer than the narrowest EEG band gap (δ-θ boundary at 4 Hz). The lab constraint of $T \geq 1200$ s makes this margin even more comfortable. The window choice for EEG is therefore driven entirely by side-lobe suppression, not by resolution.

### A.3.6 The Gaussian Window (deferred)

The cosine-sum family is not the only option. The **Gaussian window**:

$$
w_{\text{Gauss}}[n] = \exp\!\left(-\frac{1}{2}\left(\frac{n - (N-1)/2}{\sigma}\right)^2\right), \qquad n = 0, 1, \ldots, N-1 \tag{A.31}
$$

where $\sigma$ controls the width, has a unique property: its Fourier transform is also a Gaussian. This makes it the **minimum-uncertainty window** - the one that achieves equality in the time-frequency uncertainty bound $\Delta t \cdot \Delta f \geq 1/(4\pi)$.

We flag it here because it will reappear in A.8, where the SPWVD's smoothing windows are chosen. The Gaussian's minimum-uncertainty property makes it the natural choice for time-frequency smoothing, where the tradeoff between time and frequency resolution must be managed per axis. The full discussion is deferred to that section.

## A.4 Statistics and the DFT

### A.4.1 Each Bin as a Random Variable

Consider a signal that is pure white Gaussian noise: $x[n] \sim \mathcal{N}(0, \sigma_x^2)$, with each sample independent. What does its DFT look like?

Each DFT bin is a weighted sum of $N$ independent Gaussian samples (Equation (A.12)). By linearity, $X[k]$ is itself a complex Gaussian random variable. Its real and imaginary parts are independent, each with variance $N\sigma_x^2 / 2$:

$$
\text{Re}\{X[k]\}, \; \text{Im}\{X[k]\} \sim \mathcal{N}\!\left(0, \; \frac{N\sigma_x^2}{2}\right) \tag{A.32}
$$

The **magnitude** $|X[k]|$ follows a Rayleigh distribution, and the **phase** $\angle X[k]$ is uniformly distributed on $(-\pi, \pi]$. The phase carries no information about the noise - it is purely random, uniformly spread.

The quantity we care about for spectral analysis is the **power** at each bin:

$$
P[k] = |X[k]|^2 = \text{Re}\{X[k]\}^2 + \text{Im}\{X[k]\}^2 \tag{A.33}
$$

Since $P[k]$ is the sum of squares of two independent Gaussians, it follows an **exponential distribution**:

$$
P[k] \sim \text{Exponential}\!\left(\lambda = \frac{1}{N\sigma_x^2}\right) \tag{A.34}
$$

with mean $\mathbb{E}[P[k]] = N\sigma_x^2$ and - crucially - standard deviation also equal to $N\sigma_x^2$. The coefficient of variation is 1: the standard deviation equals the mean. This means a single bin's power estimate is inherently noisy, with 100% relative uncertainty, regardless of how many samples $N$ we use.

For comparison, a deterministic signal (e.g. a pure tone) concentrates its energy in a few bins, leaving most bins near zero. This produces CV $\gg$ 1 - the standard deviation far exceeds the mean because a few bins are very large while the majority are very small. The CV of the bin power distribution therefore distinguishes noise-like spectra (CV $\approx$ 1, energy spread across all bins) from spectra with concentrated features (CV $\gg$ 1, energy in a few bins). This property is used in Volume C, Section C.4.

Now add a deterministic tone to the noise: $x[n] = A\cos(2\pi f_0 n / f_s) + \eta[n]$. At the bin $k_0$ nearest to $f_0$, the DFT output is a complex Gaussian shifted by the tone's contribution. The power $P[k_0]$ now follows a **non-central chi-squared distribution** - its mean is elevated above the noise-only level by the tone's energy. At all other bins (away from $f_0$), the distribution remains exponential as in Equation (A.34).

This is the statistical model. A tone reveals itself as a bin whose power is improbably large under the null hypothesis of noise alone.

### A.4.2 The Noise Floor and Spectral Detection

The **noise floor** is the expected power level of bins that contain only noise. From Equation (A.34), it is:

$$
P_{\text{floor}} = \mathbb{E}[P[k]] = N\sigma_x^2 \tag{A.35}
$$

In practice, we do not know $\sigma_x^2$ in advance. We estimate $P_{\text{floor}}$ from the spectrum itself. A robust estimator is the **median** of the bin powers:

$$
\hat{P}_{\text{floor}} = \text{median}\{P[0], P[1], \ldots, P[N/2]\} \tag{A.36}
$$

The median is preferred over the mean because a few bins containing genuine signal components inflate the mean but barely affect the median.

**Detection criterion.** Under the exponential model in Equation (A.34), the probability that a noise-only bin exceeds a threshold $\gamma \cdot \hat{P}_{\text{floor}}$ is:

$$
\Pr\!\left(P[k] > \gamma \cdot P_{\text{floor}}\right) = e^{-\gamma} \tag{A.37}
$$

Equation (A.37) is the detection rule, derived from the spectral model - not from a generic "$2\sigma$" rule. The threshold $\gamma$ is chosen to control the false-alarm rate:

**Table A.5 - Detection thresholds from the exponential noise model**

| Threshold $\gamma$ | False-alarm probability $e^{-\gamma}$ | Interpretation |
| --- | --- | --- |
| 3 | 0.050 | 5% of noise-only bins exceed this |
| 4.6 | 0.010 | 1% false-alarm rate |
| 6.9 | 0.001 | 0.1% false-alarm rate |

The key difference from generic statistics: the threshold is not "$2\sigma$ because that is common." It is $\gamma \cdot \hat{P}_{\text{floor}}$ because the exponential distribution of bin power under noise gives us Equation (A.37), and we choose $\gamma$ to match the false-alarm probability we can tolerate. The noise floor $\hat{P}_{\text{floor}}$ is estimated from the spectrum, so the threshold adapts to the actual noise level in the recording - not to an assumed distribution of time-domain samples.

This also explains why crude time-domain thresholding (e.g. "remove samples beyond $2\sigma$") is problematic: it operates on samples, not on frequency content. Removing a time-domain sample affects every bin in the DFT (by Equation (A.12), every bin is a sum over all samples). The energy removed is spread unpredictably across the spectrum, breaking the Parseval relationship in Equation (A.21) in ways that are difficult to characterize.

### A.4.3 The Periodogram Variance Problem

The **periodogram** is the name for the power spectrum estimated from a single DFT:

$$
\hat{S}[k] = \frac{1}{N} |X[k]|^2 \tag{A.38}
$$

Equation (A.38) is the natural estimator: compute the DFT, take the magnitude squared, scale by $1/N$. It is an **unbiased** estimator of the true power spectrum - its expected value equals the true spectral density at each bin.

But it is not a **consistent** estimator. This is the surprising and counterintuitive result:

$$
\text{Var}\!\left(\hat{S}[k]\right) \approx S[k]^2 \tag{A.39}
$$

The variance of the periodogram at bin $k$ is proportional to the square of the true spectral value - and this does not decrease as $N$ increases. Doubling the signal length $N$ halves the bin spacing (Equation (A.15)) and doubles the number of bins, but each bin's power estimate remains just as noisy.

This is a direct consequence of the exponential distribution in Equation (A.34): the coefficient of variation (standard deviation / mean) is always 1, independent of $N$. More data gives us finer frequency spacing, not more reliable power estimates.

In plain terms: a single DFT of a long signal gives a ragged, noisy spectrum, and making the signal longer does not smooth it. The raggedness is statistical, not a failure of implementation. Any single realization of a noisy signal produces a periodogram that fluctuates wildly around the true spectral shape.

### A.4.4 Averaging to Reduce Variance: Welch's Method

The periodogram variance problem has a direct solution: **average multiple independent periodograms**. If we have $L$ independent estimates $\hat{S}_1[k], \hat{S}_2[k], \ldots, \hat{S}_L[k]$, the averaged estimate:

$$
\bar{S}[k] = \frac{1}{L} \sum_{\ell=1}^{L} \hat{S}_\ell[k] \tag{A.40}
$$

has variance reduced by a factor of $L$:

$$
\text{Var}\!\left(\bar{S}[k]\right) = \frac{1}{L} \, S[k]^2 \tag{A.41}
$$

**Welch's method** obtains the $L$ independent estimates from a single recording by dividing the signal into overlapping segments, windowing each segment, computing its periodogram, and averaging:

1. Divide the $N$-sample signal into segments of length $M$, with overlap of $D$ samples.
2. Apply a window $w[n]$ (from Section A.3) to each segment.
3. Compute the periodogram of each windowed segment.
4. Average all $L$ periodograms.

The number of segments is:

$$
L = \left\lfloor \frac{N - M}{M - D} \right\rfloor + 1 \tag{A.42}
$$

The tradeoff is explicit and quantitative:

- Each segment has $M < N$ samples, so the bin spacing per segment is $f_s / M$ - coarser than the $f_s / N$ achievable from the full signal.
- The frequency resolution degrades from $\Delta f = f_s / N$ to $\Delta f = \beta \cdot f_s / M$ (Equation (A.30)).
- The variance decreases by the factor $1/L$.

**Table A.6 - Welch's method: the resolution-variance tradeoff at $f_s = 250$ Hz, $N = 300\,000$ (1200 s), Hann window ($\beta = 2$), 50% overlap**

| Segment length $M$ | Segments $L$ | $\Delta f$ (Hz) | Relative variance |
| --- | --- | --- | --- |
| 300 000 (full) | 1 | 0.0017 | 1.00 |
| 30 000 (120 s) | 19 | 0.017 | 0.053 |
| 5 000 (20 s) | 119 | 0.10 | 0.0084 |
| 1 250 (5 s) | 479 | 0.40 | 0.0021 |
| 500 (2 s) | 1 199 | 1.0 | 0.00083 |

Table A.6 shows the tradeoff for our lab signal parameters. A single periodogram of the full 1200-second signal gives the finest frequency resolution (0.0017 Hz) but with 100% relative variance - the spectrum is ragged and unreliable. At the other extreme, 2-second segments give 1199 averages and a smooth spectrum, but the resolution is only 1 Hz - too coarse to separate neighboring EEG sub-bands.

The practical choice for EEG depends on what needs to be resolved. For separating the standard bands (δ, θ, α, β, γ), the boundaries are at 4, 8, 13, and 30 Hz - separations of at least 4 Hz. A segment length of $M = 1250$ (5 s, $\Delta f = 0.40$ Hz) resolves these boundaries comfortably with nearly 500 averages. For finer structure within a band - say, distinguishing 9 Hz from 11 Hz alpha - longer segments are needed ($M = 5000$, giving $\Delta f = 0.10$ Hz).

Welch's method does not bypass the uncertainty in each bin; it trades frequency detail for statistical reliability. The total energy is still conserved (Equation (A.21) applies to each segment), and the noise floor estimated from the averaged spectrum is far more stable than from a single periodogram. Detection via Equation (A.37) applied to a Welch-averaged spectrum is correspondingly more reliable.

## A.5 The STFT and Spectrograms

### A.5.1 The Windowed DFT Slid in Time

The STFT is built from exactly the components already developed: a window (Section A.3), a DFT (Section A.2), and the sliding-segment idea from Welch (Section A.4.4). The difference is that Welch averages the segment spectra into one; the STFT **keeps each segment spectrum indexed by its position in time**.

Given a signal $x[n]$, a window $w[n]$ of length $M$, and a hop size $H$ (the number of samples the window advances per step), the STFT is:

$$
X[m, k] = \sum_{n=0}^{M-1} x[n + mH] \, w[n] \, e^{-j 2\pi k n / M}, \qquad k = 0, 1, \ldots, M-1 \tag{A.43}
$$

where $m$ is the **segment index** (the time position) and $k$ is the **frequency bin** within that segment. Each value $X[m, k]$ is one DFT bin of one windowed segment.

The output of Equation (A.43) is a **two-dimensional grid** - the **time-frequency plane**:

- The **time axis** is discrete, with spacing $\Delta t_{\text{grid}} = H / f_s$ seconds between columns.
- The **frequency axis** is discrete, with spacing $\Delta f = f_s / M$ Hz between rows (Equation (A.15), applied to each segment of length $M$).
- Each grid cell contains a complex number $X[m, k]$.

The **spectrogram** is the squared magnitude of the STFT:

$$
S[m, k] = |X[m, k]|^2 \tag{A.44}
$$

Equation (A.44) is the power at time step $m$ and frequency bin $k$. It is the object we plot: a 2D image with time on the horizontal axis, frequency on the vertical axis, and color representing power.

The physical frequency and time of each cell are:

$$
f_k = \frac{k \cdot f_s}{M} \quad \text{(Hz)}, \qquad t_m = \frac{m \cdot H}{f_s} \tag{A.45}
$$

The connection to Welch is now explicit. Welch computes $S[m, k]$ for all $m$ and $k$, then averages over $m$: $\bar{S}[k] = \frac{1}{L}\sum_m S[m,k]$. The STFT keeps the $m$ axis. Welch answers "what frequencies are present on average?"; the STFT answers "what frequencies are present, and when?"

### A.5.2 The Uncertainty Principle

This is the central tradeoff of the entire report.

The window length $M$ controls both time resolution and frequency resolution, but in opposite directions:

**Time resolution.** The STFT localizes a signal event to the window that contains it. The finest time localization is the window duration:

$$
\Delta t = \frac{M}{f_s} \tag{A.46}
$$

A short window (small $M$) gives fine time resolution: a 1-second burst of alpha rhythm lands in one or two columns. A long window (large $M$) smears the burst across a single column that also contains whatever came before and after.

**Frequency resolution.** Within each window, the DFT resolves frequencies to (from Equation (A.30)):

$$
\Delta f = \beta \cdot \frac{f_s}{M} \tag{A.47}
$$

A long window (large $M$) gives fine frequency resolution: a 10 Hz tone and a 10.5 Hz tone appear as separate peaks. A short window (small $M$) merges them into one lobe.

Equations (A.46) and (A.47) compete. Their product:

$$
\Delta t \cdot \Delta f = \beta \tag{A.48}
$$

For the rectangular window ($\beta = 1$), this is already $\geq 1$. For Hann ($\beta = 2$), it is $\geq 2$. The fundamental lower bound, achievable only by the Gaussian window (Equation (A.31)), is:

$$
\Delta t \cdot \Delta f \geq \frac{1}{4\pi} \tag{A.49}
$$

Equation (A.49) is the **Heisenberg-Gabor uncertainty principle** for time-frequency analysis. It is not a limitation of the STFT algorithm. It is a mathematical property of Fourier pairs: a function and its Fourier transform cannot both be arbitrarily concentrated. No linear time-frequency representation - no matter how cleverly constructed - escapes this bound.

**What the uncertainty principle means in practice.** The window length $M$ is a single slider between two extremes. You do not choose "good" or "bad" - you choose *which axis to sacrifice*:

**Table A.7 - The uncertainty tradeoff at $f_s = 250$ Hz, Hann window ($\beta = 2$)**

| Window $M$ (samples) | Duration $\Delta t$ (s) | Freq. resolution $\Delta f$ (Hz) | $\Delta t \cdot \Delta f$ | What it can resolve | What it cannot resolve |
| --- | --- | --- | --- | --- | --- |
| 125 (0.5 s) | 0.5 | 4.0 | 2.0 | Sub-second transients, blinks | Cannot separate α (8-13) from θ (4-8) cleanly |
| 250 (1.0 s) | 1.0 | 2.0 | 2.0 | Second-scale events; α vs. θ bands | Cannot resolve structure within a band |
| 500 (2.0 s) | 2.0 | 1.0 | 2.0 | Individual Hz within bands | Cannot localize events shorter than 2 s |
| 1 250 (5.0 s) | 5.0 | 0.4 | 2.0 | Fine spectral detail (9 vs. 11 Hz) | Cannot localize events shorter than 5 s |
| 5 000 (20.0 s) | 20.0 | 0.1 | 2.0 | Sub-Hz spectral structure | Cannot track any temporal variation faster than 20 s |

Table A.7 makes the tradeoff concrete for EEG. Every row has the same $\Delta t \cdot \Delta f = 2.0$ (the Hann window's $\beta$). Moving down the table buys frequency detail and costs time detail, in exact proportion. No row is "best" - each is suited to a different question about the signal.

The last row illustrates the extreme: a 20-second window resolves 0.1 Hz, but the spectrogram has one column per 20 seconds - it cannot track any rhythm that changes faster than that. At the other extreme, a 0.5-second window catches sub-second transients but cannot distinguish alpha from theta. The STFT forces you to decide what matters before you compute.

This is the limitation that motivates the second half of the report. The WVD (Section A.7) is not a linear representation - it is quadratic - and therefore is not bound by Equation (A.49) in the same way. It achieves sharper localization on both axes simultaneously, at a different cost (cross-terms). The SPWVD (Section A.8) manages that cost with two independent smoothing knobs. But the uncertainty principle remains the reference: any improvement the WVD family offers is measured against the STFT's bound.

### A.5.3 Hop Size, Overlap, and the COLA Condition

The STFT in Equation (A.43) has two parameters beyond the window itself: the **hop size** $H$ and, equivalently, the **overlap** $D = M - H$. These are often treated as minor implementation details. They are not. Without proper overlap, the window's tapering creates a systematic problem: **samples near segment boundaries are suppressed or lost**.

#### The tapering problem

Every non-rectangular window tapers toward zero at its edges. From Table A.3, a Hann window satisfies $w[0] = w[N-1] = 0$: the first and last samples of each segment are multiplied by zero. A Blackman window does the same. Even a Hamming window, which does not reach zero, reduces edge samples to 8% of their original amplitude.

Consider what happens when adjacent segments do not overlap ($H = M$, zero overlap). Each sample of the original signal appears in exactly one segment. But within that segment, the window assigns it a weight:

- Samples near the center receive weight $\approx 1$.
- Samples near the edges receive weight $\approx 0$ (for Hann/Blackman) or weight $\approx 0.08$ (for Hamming).

The result: samples that happen to fall near a segment boundary are effectively erased from the analysis. Their contribution to every frequency bin in that segment is multiplied by a number close to zero. This is not a minor effect - for a Hann window, the samples at positions $n = 0$ and $n = M - 1$ are multiplied by exactly zero. They contribute nothing.

If the signal has an important feature (a spike, a burst onset, a transient) at a segment boundary, that feature is suppressed. Not by filtering, not by noise, but by the window's taper. The feature exists in the data; the STFT simply does not see it.

#### How overlap solves the tapering problem

Overlap means that adjacent segments share some samples. With hop size $H < M$, each sample appears in multiple segments. In each segment, it sits at a different position relative to the window center, and therefore receives a different weight.

The critical question is: **what is the total weight a sample receives, summed across all segments it appears in?**

For a sample at position $n$ in the original signal, this total weight is:

$$
W_{\text{total}}[n] = \sum_{m} w[n - mH] \tag{A.50}
$$

where the sum runs over all segment indices $m$ for which the sample falls within the window. If $W_{\text{total}}[n]$ varies with $n$, then different samples receive different total representation in the STFT - some parts of the signal are emphasized, others are suppressed. This is an amplitude modulation imposed by the analysis, not by the signal.

**Example: Hann window with no overlap ($H = M$).** Each sample appears in exactly one segment. $W_{\text{total}}[n] = w[n \bmod M]$, which varies from 0 at the edges to 1 at the center. Samples at segment boundaries are lost.

**Example: Hann window with 50% overlap ($H = M/2$).** Each sample appears in exactly two segments. In one segment, it sits at position $p$; in the adjacent segment, it sits at position $p + M/2$. The Hann window satisfies:

$$
w_{\text{Hann}}[p] + w_{\text{Hann}}[p + M/2] = 1 \tag{A.51}
$$

Equation (A.51) states that the two overlapping Hann windows sum to a constant. Every sample receives the same total weight of 1, regardless of where it falls relative to segment boundaries. The tapering problem is completely eliminated: no sample is suppressed, no sample is over-represented.

#### The COLA condition

Equation (A.51) is a specific instance of the **Constant Overlap-Add (COLA) condition**:

$$
\sum_{m} w[n - mH] = C \qquad \text{(A.52)}
$$

where $C$ is a constant (typically normalized to 1). When COLA is satisfied:

1. **Every sample receives equal total weight.** The analysis does not favor any time position over another.
2. **The STFT is invertible.** The original signal can be perfectly reconstructed from the STFT coefficients by overlap-adding the inverse-DFT of each segment. This means no information is lost - the time-frequency plane is a complete, lossless representation.
3. **Parseval's theorem extends to the STFT.** The total energy computed from the spectrogram equals the total signal energy (with appropriate normalization).

**Table A.8 - COLA-satisfying overlap for standard windows**

| Window | Minimum COLA overlap | Hop size $H$ | Segments per sample |
| --- | --- | --- | --- |
| Rectangular | 0% | $M$ | 1 |
| Hann | 50% | $M/2$ | 2 |
| Hamming | 50% | $M/2$ | 2 (approximate COLA) |
| Blackman | 67% | $M/3$ | 3 |

Table A.8 shows the minimum overlap each window needs to satisfy COLA. The rectangular window satisfies COLA trivially at any overlap (including zero) because it does not taper - but it has the worst spectral leakage (Section A.3). Hann at 50% is the standard practical choice: it satisfies COLA exactly, requires only 2× the computation of zero-overlap, and has −31.5 dB peak side-lobes.

Hamming at 50% satisfies COLA only approximately - the total weight oscillates slightly around a constant because $w[0] = 0.08 \neq 0$. In practice the deviation is small (< 0.2%), and 50% overlap is used regardless. Blackman requires 67% overlap (each sample appears in 3 segments) to achieve COLA, which increases computation by 3× relative to zero-overlap.

#### Oversampling in time: more columns, not more resolution

Increasing the overlap beyond the COLA minimum (e.g. using 75% overlap with a Hann window instead of 50%) produces more time-axis columns in the spectrogram. This is the **time-axis analogue of zero-padding** (Section A.2.3):

- More columns = finer time-axis grid = smoother-looking spectrogram.
- **Time resolution does not improve.** The resolution is set by the window length $M$ (Equation (A.46)), not by the hop size. A smaller hop interpolates the same information onto a denser grid, but does not reveal temporal detail that the window duration cannot capture.

The parallel is exact: zero-padding adds frequency bins without adding frequency resolution; overlap beyond COLA adds time columns without adding time resolution. Both are interpolation - useful for visualization, but not for resolving features.

### A.5.4 Reading a Spectrogram; the Cost of Each Parameter Choice

The spectrogram in Equation (A.44) is a 2D image. Reading it requires knowing what signal features look like on the time-frequency plane. The mapping follows directly from the signal taxonomy in Appendix A:

**Table A.9 - Signal archetypes on the spectrogram**

| Archetype | Time-domain appearance | Spectrogram signature |
| --- | --- | --- |
| Stationary tone | Constant-frequency sinusoid | Horizontal line at $f_0$ |
| Mixed tones | Sum of sinusoids | Parallel horizontal lines |
| Chirp | Frequency-swept sinusoid | Diagonal line (slope = sweep rate) |
| Transient / pulse | Short burst or impulse | Vertical stripe (broadband, localized in time) |
| Noise | Random, no structure | Uniform speckle across the plane |

Table A.9 is the key to interpretation: horizontal = frequency-stable, vertical = time-localized, diagonal = frequency-changing. Any real signal - including EEG - is a mixture of these archetypes, and the spectrogram decomposes it into their superposition on the time-frequency plane.

**The cost of each parameter choice.** Choosing the STFT parameters means deciding in advance which features matter most. There is no universal setting. Table A.10 summarizes the practical consequences:

**Table A.10 - STFT parameter choices and their consequences**

| Parameter | Increasing it... | Improves... | Degrades... |
| --- | --- | --- | --- |
| Window length $M$ | Longer segments | Frequency resolution (Eq. A.38) | Time resolution (Eq. A.37) |
| Hop size $H$ | Fewer, sparser columns | Computation speed | Time-axis smoothness (not resolution) |
| Overlap $M - H$ | More columns, better COLA | Sample coverage; visual smoothness | Computation cost; redundancy |
| Window taper (e.g. Rect → Hann → Blackman) | Stronger taper | Side-lobe suppression (Table A.3) | Main-lobe width (frequency resolution) |

**The honest conclusion.** The STFT is the first tool in this report that can answer "what frequency is present at what time." It is the workhorse of practical time-frequency analysis, and for many signals it is sufficient. But the uncertainty principle (Equation (A.49)) means it always blurs something: either time or frequency, controlled by $M$, with no way to have both.

For a signal whose features all live at roughly the same time-frequency scale - say, an EEG recording dominated by a steady alpha rhythm - one window length captures everything. But EEG often contains features at multiple scales simultaneously: a narrow-band alpha rhythm (needs long $M$ to resolve) and a short transient artifact (needs short $M$ to localize). No single STFT window captures both. Two separate STFTs with different $M$ can show each feature, but there is no principled way to merge them into a single representation.

This is the limitation that the Wigner-Ville Distribution addresses. The WVD (Section A.7) is not bound by the uncertainty principle in the same way because it is quadratic, not linear. It achieves sharper localization on both axes - but at a cost (cross-terms) that requires its own management (Sections A.7-A.8). Before reaching the WVD, we need one more tool: the autocorrelation function (Section A.6), which provides the bridge.

## A.6 Autocorrelation

### A.6.1 Definition; the Periodicity Detector

The **autocorrelation** of a discrete signal $x[n]$ of length $N$ is defined as:

$$
r[l] = \sum_{n=0}^{N-1-|l|} x[n] \, x^*[n - l], \qquad l = -(N-1), \ldots, 0, \ldots, N-1 \tag{A.53}
$$

where $l$ is the **lag** (the shift in samples) and $x^*$ denotes the complex conjugate. For real-valued signals, $x^* = x$ and the conjugate can be dropped.

Equation (A.53) compares the signal with a shifted copy of itself. At each lag $l$, it multiplies the signal sample-by-sample with the version shifted by $l$ positions, and sums. If the signal and its shifted copy are similar (aligned peaks, aligned troughs), the sum is large. If they are dissimilar (peaks meeting troughs), the sum is small or negative.

Two properties follow directly from the definition:

**Property 1: $r[0]$ is the total signal energy.**

$$
r[0] = \sum_{n=0}^{N-1} |x[n]|^2 \tag{A.54}
$$

This is the same quantity that appears on the left side of Parseval's theorem (Equation (A.21)). The autocorrelation at lag zero is the energy - the maximum value of $r[l]$, since no shift can produce greater alignment than zero shift.

**Property 2: Symmetry.**

$$
r[-l] = r^*[l] \tag{A.55}
$$

For real signals, this simplifies to $r[-l] = r[l]$: the autocorrelation is an even function. Negative lags carry no new information beyond positive lags.

**Why "periodicity detector"?** If $x[n]$ contains a periodic component with period $P$ samples, then $x[n]$ and $x[n - P]$ are nearly identical - they are the same oscillation, shifted by one full cycle. The product $x[n] \cdot x[n - P]$ is large and positive at every $n$, so $r[P]$ is large. Likewise $r[2P]$, $r[3P]$, etc. The autocorrelation of a periodic signal is itself periodic at the same period, producing peaks at lags that are multiples of $P$.

This is what makes autocorrelation useful in the presence of noise. Noise has no periodic structure, so its autocorrelation is large only at $l = 0$ (its energy) and fluctuates near zero for all other lags. A tone buried in noise produces an autocorrelation that is noisy near $l = 0$ but shows clear periodic peaks at lags $l = P, 2P, 3P, \ldots$ - the periodicity survives even when the tone is invisible in the time-domain waveform.

### A.6.2 The Wiener-Khinchin Theorem

The **Wiener-Khinchin theorem** states that the autocorrelation and the power spectrum are a Fourier pair:

$$
|X[k]|^2 = \sum_{l=-(N-1)}^{N-1} r[l] \, e^{-j 2\pi k l / N} \tag{A.56}
$$

The DFT of the autocorrelation $r[l]$ is the power spectrum $|X[k]|^2$. Equivalently, the autocorrelation is the inverse DFT of the power spectrum.

This is not an approximation. It follows from substituting the definition of $r[l]$ (Equation (A.53)) and rearranging. The proof in the discrete case is short:

Starting from $r[l] = \sum_n x[n] \, x^*[n-l]$ and taking the DFT over $l$, the double sum separates into two independent sums - one giving $X[k]$ and the other giving $X^*[k]$ - whose product is $|X[k]|^2$. The details are algebraic bookkeeping; the result is exact.

The theorem connects two views of the same information:

- The **power spectrum** $|X[k]|^2$ tells you how much energy is at each frequency.
- The **autocorrelation** $r[l]$ tells you how self-similar the signal is at each lag.

They are not two different measurements. They are the same measurement expressed in two domains: one in frequency, one in lag. Knowing either one fully determines the other.

This connection is the reason Section A.6 exists. The WVD (Section A.7) will be defined as the Fourier transform of an **instantaneous** autocorrelation - an autocorrelation computed at each time position. The Wiener-Khinchin theorem is the static version of this idea: a global autocorrelation (no time index) whose Fourier transform is the global power spectrum (no time index). The WVD generalizes both sides to include time.

### A.6.3 Why Autocorrelation Discards Phase

The power spectrum $|X[k]|^2$ contains no phase information. Writing $X[k] = |X[k]| \, e^{j\phi[k]}$, the power spectrum is $|X[k]|^2$ - the phase $\phi[k]$ is gone. This means:

- A cosine starting at $t = 0$ and a cosine starting at $t = 0.5$ s have **identical** power spectra (same frequency, same amplitude, different phase).
- A signal and any time-shifted version of itself have **identical** autocorrelations and power spectra.

This is not a defect - it is inherent in what autocorrelation measures. Autocorrelation asks: "does the signal repeat itself after $l$ samples?" The answer depends on the signal's periodicity and amplitude, not on when it started. Phase encodes the "when"; autocorrelation is blind to it.

**The consequence for time-frequency analysis.** The power spectrum (and therefore the autocorrelation) tells you *which* frequencies are present and *how strong* they are, but not *when* they occur. This is why the global DFT power spectrum cannot distinguish a 10 Hz tone that plays for the entire recording from a 10 Hz burst that lasts one second - both deposit the same energy into the same bin, just spread differently in time.

The STFT (Section A.5) solved this by windowing: compute the power spectrum locally in time, so "when" is recovered at the cost of frequency resolution. But we now see the problem from the autocorrelation side: the global autocorrelation $r[l]$ averages over all time, so time information is lost.

The WVD's approach (Section A.7) will be different. Instead of windowing the signal and then computing the autocorrelation, it computes an **autocorrelation at each time instant** - the instantaneous autocorrelation $r_n[l] = x[n + l/2] \, x^*[n - l/2]$ - and then takes the Fourier transform over $l$. This preserves the time index $n$ through the entire computation. The Wiener-Khinchin theorem (Equation (A.56)) becomes a time-indexed family of Fourier pairs: one instantaneous autocorrelation and one instantaneous power spectrum at each $n$. That is the WVD.

### A.6.4 Autocorrelation Signatures of the Signal Archetypes

Before moving to the WVD, it is useful to see what autocorrelation does to each signal archetype. These signatures carry directly into the WVD, where they appear along the lag axis at each time position.

**Table A.11 - Autocorrelation signatures of the signal archetypes**

| Archetype | Signal form | Autocorrelation $r[l]$ | Key feature |
| --- | --- | --- | --- |
| Single tone ($f_0$) | $A\cos(2\pi f_0 n / f_s)$ | Cosine at the same frequency: $\propto \cos(2\pi f_0 l / f_s)$ | Periodic, undamped - frequency preserved, phase lost |
| Mixed tones ($f_1, f_2$) | $A_1\cos(\ldots) + A_2\cos(\ldots)$ | Sum of cosines at $f_1$ and $f_2$ | Each tone contributes independently; no cross-terms |
| Impulse | $\delta[n - n_0]$ | $\delta[l]$ (spike at lag zero only) | No self-similarity at any nonzero lag |
| White noise | $\eta[n] \sim \mathcal{N}(0, \sigma^2)$ | $\approx N\sigma^2 \cdot \delta[l]$ | Energy spike at $l = 0$; near-zero for $l \neq 0$ |
| Tone in noise | $A\cos(\ldots) + \eta[n]$ | Cosine at $f_0$ riding on a noise spike at $l = 0$ | Periodicity emerges at lags $l \gg 0$ where noise vanishes |

Table A.11 shows why autocorrelation is a periodicity detector (Section A.6.1). The tone and mixed-tone archetypes produce periodic autocorrelations - the frequency information survives, but the phase (start time) does not. Noise contributes only at $l = 0$. A tone buried in noise is invisible in the time domain but visible in the autocorrelation at lags beyond the noise spike.

Two observations from Table A.11 set up the WVD:

1. **No cross-terms for mixed tones.** The autocorrelation of $x_1[n] + x_2[n]$ contains $r_{x_1}[l] + r_{x_2}[l]$ plus cross-terms $r_{x_1 x_2}[l]$. For uncorrelated tones, the cross-terms vanish. This will **not** be the case for the WVD, where the quadratic structure produces cross-terms that do not vanish - the central complication of Sections A.7-A.8.

2. **No time index.** Every entry in Table A.11 is a single function of lag $l$, averaged over all time. A chirp (frequency sweeping over time) would show a blurred autocorrelation - an average of cosines at many frequencies - losing the sweep information. The WVD fixes this by computing the autocorrelation at each time instant separately.

### A.6.5 Cross-Correlation

Autocorrelation compares a signal with shifted copies of **itself**. Cross-correlation compares two **different** signals:

$$
r_{xy}[l] = \sum_{n=0}^{N-1-|l|} x[n] \cdot y^*[n - l] \tag{A.57}
$$

The formula is identical to Equation (A.53) except that $x$ is compared to $y$, not to itself. If $y = x$, Equation (A.57) reduces to autocorrelation.

**What it measures.** At each lag $l$, cross-correlation asks: "how similar are $x$ and $y$ when $y$ is shifted by $l$ samples?" A large $r_{xy}[l]$ means the two signals share structure at that lag. A near-zero $r_{xy}[l]$ for all $l$ means the signals are unrelated.

**Cross-spectral density.** The Wiener-Khinchin theorem (Equation (A.56)) generalizes: the DFT of the cross-correlation is the **cross-spectrum**:

$$
S_{xy}[k] = X[k] \cdot Y^*[k] = \text{DFT}\{r_{xy}[l]\} \tag{A.58}
$$

where $X[k]$ and $Y[k]$ are the DFTs of $x$ and $y$. The cross-spectrum is complex - its magnitude measures how much shared power exists at each frequency, and its phase measures the timing relationship.

**Normalized cross-correlation.** To get a dimensionless measure of similarity, normalize by the energies:

$$
\rho_{xy} = \frac{r_{xy}[0]}{\sqrt{r_{xx}[0] \cdot r_{yy}[0]}} = \frac{\sum_n x[n] \cdot y[n]}{\sqrt{\sum_n |x[n]|^2 \cdot \sum_n |y[n]|^2}} \tag{A.59}
$$

This is the **Pearson correlation coefficient** at zero lag. It ranges from $-1$ (perfectly anti-correlated) through $0$ (uncorrelated) to $+1$ (identical up to scale). In practice, `np.corrcoef(x, y)[0, 1]` computes this directly.

**Why this matters for EEG.** An EEG recording has multiple channels. Cross-correlation between an EEG channel (CZ) and an auxiliary channel (ECG) answers: "does the heartbeat appear in the brain signal?" If $|\rho_{xy}| \approx 0$, there is no contamination. If $|\rho_{xy}|$ is significant, the ECG is leaking into the EEG - and the cross-spectrum $S_{xy}[k]$ shows at which frequencies the leakage occurs. This is used in Volume C, Section C.4.

## A.7 The Wigner-Ville Distribution

### A.7.1 Instantaneous Autocorrelation and Definition

The global autocorrelation $r[l]$ (Section A.6) compares a signal with its shifted copy over the entire time record, collapsing the time index. To construct a time-varying spectrum, we need a localized measure of self-similarity: the **instantaneous autocorrelation function (IAF)**.

A naive discretization would evaluate $x[n + m/2] \cdot x^*[n - m/2]$, which requires half-integer sample values when $m$ is odd. To stay on the integer sample grid, the **discrete-time IAF** uses a double-lag step:

$$
R_x[n, m] = x[n + m] \cdot x^*[n - m] \tag{A.60}
$$

where $n$ is the discrete time sample and $m$ is the lag. The **Discrete-Time Wigner-Ville Distribution (DT-WVD)** is the DFT of this sequence with respect to $m$:

$$
W_x[n, k] = 2 \sum_{m=-M}^{M} x[n + m] \, x^*[n - m] \, e^{-j \frac{4\pi k m}{N_f}} \tag{A.61}
$$

where $M$ is the maximum lag and $N_f$ is the number of frequency bins. The factor of 2 is a convention to align the energy scale.

Equation (A.61) is the working formula used in all code. It is a discrete sum over lags, evaluated at each time sample $n$ - no integrals, no continuous variables.

**The frequency aliasing problem.** The double-lag step $x[n+m] \cdot x^*[n-m]$ means the effective sampling interval of the autocorrelation sequence is $2\Delta t$, halving the Nyquist frequency from $f_s/2$ to $f_s/4$. To prevent aliasing, the input signal must be **interpolated by a factor of 2** before computing the WVD. In the interpolated signal (sampled at $2f_s$), the lag step restores the Nyquist limit to $f_s/2$.

### A.7.2 Sharpness of the Single Chirp

The principal advantage of the WVD over the STFT is its ability to bypass the Heisenberg uncertainty principle (Equation (A.49)) for single-component signals. We demonstrate this on a discrete linear chirp (analytic signal):

$$
x[n] = A \, e^{j 2\pi \left(\frac{f_0 n}{f_s} + \frac{\mu n^2}{2 f_s^2}\right)} \tag{A.62}
$$

where $f_0$ is the starting frequency and $\mu$ is the chirp rate in Hz/s. The instantaneous frequency at sample $n$ is $f_{\text{inst}}[n] = f_0 + \mu n / f_s$.

**Computing the IAF.** Substituting Equation (A.62) into Equation (A.60):

$$
R_x[n, m] = x[n+m] \cdot x^*[n-m] = A^2 \, e^{j(\phi(n+m) - \phi(n-m))}
$$

Expanding the phase difference:

$$
\phi(n+m) - \phi(n-m) = \frac{2\pi f_0}{f_s}\left[(n+m) - (n-m)\right] + \frac{\pi\mu}{f_s^2}\left[(n+m)^2 - (n-m)^2\right]
$$

$$
= \frac{4\pi f_0 m}{f_s} + \frac{4\pi\mu n m}{f_s^2} = \frac{4\pi m}{f_s}\left(f_0 + \frac{\mu n}{f_s}\right) = \frac{4\pi m}{f_s} f_{\text{inst}}[n] \tag{A.63}
$$

The IAF simplifies to a complex exponential in $m$ whose frequency depends on the instantaneous frequency at time $n$:

$$
R_x[n, m] = A^2 \, e^{j \frac{4\pi m}{f_s} f_{\text{inst}}[n]} \tag{A.64}
$$

**Computing the DT-WVD.** Substituting Equation (A.64) into Equation (A.61):

$$
W_x[n, k] = 2A^2 \sum_{m=-M}^{M} e^{j 4\pi m \left(\frac{f_{\text{inst}}[n]}{f_s} - \frac{k}{N_f}\right)} \tag{A.65}
$$

Equation (A.65) is a **Dirichlet kernel** - the same structure as the rectangular window DFT derived in Lab 3. It peaks sharply when the exponent vanishes:

$$
\frac{f_{\text{inst}}[n]}{f_s} = \frac{k}{N_f} \quad \Rightarrow \quad k = \frac{N_f \cdot f_{\text{inst}}[n]}{f_s} \tag{A.66}
$$

At this bin, the sum evaluates to $2A^2(2M+1)$. Away from it, the terms oscillate and largely cancel. The peak width is approximately $1/(2M+1)$ bins - the longer the lag range $M$, the sharper the peak. As $M \to \infty$, the peak approaches a delta function.

This is the discrete equivalent of the continuous result $W_x(t,f) = A^2 \delta(f - f_{\text{inst}}(t))$: the DT-WVD of a linear chirp is a **sharp peak** that tracks the instantaneous frequency exactly, with resolution that improves with $M$ and is not bound by the STFT's uncertainty tradeoff. No integrals, no delta functions assumed - just a finite sum that produces a sharp peak.

### A.7.3 The Quadratic Nature and Cross-Terms

The WVD's sharpness comes from its **quadratic (bilinear) structure**: the signal is multiplied by itself. This introduces a severe complication for multi-component signals. The WVD is non-linear; superposition does not apply.

Let $x[n] = x_1[n] + x_2[n]$. The IAF expands as:

$$
R_x[n, m] = R_{x_1}[n, m] + R_{x_2}[n, m] + x_1[n+m] \cdot x_2^*[n-m] + x_2[n+m] \cdot x_1^*[n-m] \tag{A.67}
$$

Taking the DFT over $m$:

$$
W_x[n, k] = W_{x_1}[n, k] + W_{x_2}[n, k] + 2\,\text{Re}\{W_{x_1 x_2}[n, k]\} \tag{A.68}
$$

The third term is the **cross-term**. For two discrete tones $x_1[n] = A_1 e^{j(2\pi f_1 n/f_s + \phi_1)}$ and $x_2[n] = A_2 e^{j(2\pi f_2 n/f_s + \phi_2)}$, the cross-term evaluates to:

$$
2\,\text{Re}\{W_{x_1 x_2}[n, k]\} = 2A_1 A_2 \cos\!\left(\frac{4\pi(f_1 - f_2)n}{f_s} + (\phi_1 - \phi_2)\right) \cdot D_k\!\left(\frac{f_1 + f_2}{2}\right) \tag{A.69}
$$

where $D_k$ denotes the Dirichlet kernel centered at the midpoint frequency $(f_1 + f_2)/2$. This defines the three fundamental properties of WVD cross-terms:

1. **Midpoint location:** the cross-term peaks at $f_c = (f_1 + f_2)/2$, exactly halfway between the two components.
2. **Temporal oscillation:** it oscillates in time with frequency $|f_1 - f_2|$ (the difference frequency).
3. **Amplitude:** the peak amplitude is $2A_1 A_2$ - twice the product of the component amplitudes. If $A_1 = A_2 = 1$, the cross-term is twice as large as the auto-terms.

For $K$ components, there are $K$ auto-terms and $K(K-1)/2$ cross-terms. The cross-terms grow quadratically and can dominate the representation.

### A.7.4 Duality of Cross-Terms (Tones vs. Impulses)

Cross-terms exhibit a strict duality between the time and frequency domains:

**Table A.12 - Duality of WVD cross-terms**

| Separation | Midpoint location | Oscillation direction | Oscillation rate |
| --- | --- | --- | --- |
| Two tones ($\Delta f$) | Midpoint frequency $(f_1+f_2)/2$ | Time axis | $\Delta f$ Hz |
| Two impulses ($\Delta t$) | Midpoint time $(t_1+t_2)/2$ | Frequency axis | $\Delta t$ s |

For tones separated in frequency, the cross-term oscillates in time. For impulses separated in time, the cross-term oscillates in frequency. For a general time-frequency separation, the cross-term is centered at the joint midpoint and oscillates perpendicular to the line connecting the two components.

### A.7.5 Analytic Signal (Hilbert Transform) and DC Self-Ghost

A real-valued tone $x[n] = A\cos(2\pi f_0 n / f_s)$ consists of two complex exponentials at $+f_0$ and $-f_0$. The WVD treats these as two components and generates a cross-term between them:

- Midpoint: $(f_0 + (-f_0))/2 = 0$ Hz (DC)
- Oscillation: $|f_0 - (-f_0)| = 2f_0$ Hz

This is the **DC self-ghost** - an oscillating artifact at 0 Hz that corrupts the low-frequency plane. In EEG, where delta-band activity (0.5-4 Hz) is clinically significant, this would be catastrophic.

The solution: convert the real signal to its **analytic signal** before computing the WVD:

$$
z[n] = x[n] + j\,\mathcal{H}\{x[n]\} \tag{A.70}
$$

where $\mathcal{H}\{\cdot\}$ is the Hilbert transform. The analytic signal retains only positive frequencies, eliminating the $-f_0$ component and preventing the DC self-ghost. This step is mandatory for any practical WVD implementation.

## A.8 The PWVD and the SPWVD

### A.8.1 The Lag Window $h$ (Frequency Smoothing)

To suppress cross-terms, the **Pseudo Wigner-Ville Distribution (PWVD)** windows the IAF in the lag domain with a symmetric window $h[m]$:

$$
PW_x[n, k] = 2 \sum_{m=-M}^{M} h[m] \, x[n+m] \, x^*[n-m] \, e^{-j \frac{4\pi k m}{N_f}} \tag{A.71}
$$

Multiplication by $h[m]$ in the lag domain is convolution in frequency. The PWVD smooths the WVD along the frequency axis, suppressing cross-terms that oscillate in the frequency direction (those from components separated in time, per Table A.12).

### A.8.2 Limitation: Time-Oscillating Ghosts

The PWVD performs no smoothing along the time axis. Cross-terms that oscillate in time (from components separated in frequency, per Table A.12) survive the PWVD completely intact. For EEG with simultaneous rhythms in different bands, this is a severe limitation.

### A.8.3 The Time Window $g$ (Time Smoothing)

The **Smoothed Pseudo Wigner-Ville Distribution (SPWVD)** adds a second window $g[p]$ that averages across adjacent time steps:

$$
SPW_x[n, k] = 2 \sum_{m=-M}^{M} h[m] \left(\sum_{p=-P}^{P} g[p] \, x[n-p+m] \, x^*[n-p-m]\right) e^{-j \frac{4\pi k m}{N_f}} \tag{A.72}
$$

The double summation provides **two independent smoothing knobs**:

- Inner sum: convolves the IAF in time with $g[p]$, suppressing time-oscillating cross-terms.
- Outer sum: windows in lag with $h[m]$ and computes the DFT, suppressing frequency-oscillating cross-terms.

### A.8.4 The Smoothing Tradeoff

Smoothing suppresses cross-terms but also smears genuine features. The SPWVD's advantage over the STFT: the two smoothing axes are **decoupled**:

- Time resolution is set by $g[p]$ only.
- Frequency resolution is set by $h[m]$ only.

The STFT ties both to a single window length $M$. The SPWVD allows independent tuning.

**The minimum smoothing rule.** To suppress a cross-term, the window must average at least one full cycle of its oscillation:

To suppress a time-oscillating ghost from components separated by $\Delta f$ Hz, the time window duration must satisfy:

$$
T_g \geq \frac{1}{\Delta f} \tag{A.73}
$$

To suppress a frequency-oscillating ghost from components separated by $\Delta t$ s, the lag window bandwidth must satisfy:

$$
B_H \geq \frac{1}{\Delta t} \tag{A.74}
$$

Hann or Gaussian windows are preferred for both $g$ and $h$ - their transforms decay rapidly without strong side lobes.

### A.8.5 Cohen's Class: A Unifying Framework

The WVD, PWVD, SPWVD, and spectrogram are all members of **Cohen's class** of quadratic time-frequency distributions. Any shift-invariant distribution can be written as:

$$
C_x[n, k] = \sum_p \sum_m \Pi[p, m] \, W_x[n-p, k-m] \tag{A.75}
$$

where $\Pi[p, m]$ is a 2D smoothing kernel applied to the raw WVD:

| Distribution | Kernel $\Pi$ | Resolution | Cross-terms |
| --- | --- | --- | --- |
| Raw WVD | $\delta[p]\,\delta[m]$ (no smoothing) | Maximum | Maximum |
| Spectrogram | WVD of the STFT window (coupled) | Heisenberg-limited | Fully suppressed |
| SPWVD | $g[p] \cdot H[m]$ (separable, independent) | Tunable per axis | Controllable |

The raw WVD and the spectrogram are the two extremes of a single family. The SPWVD sits between them, offering a configurable bridge: trade resolution for ghost suppression in a controlled, axis-independent manner. Wavelets and other quadratic distributions also belong to Cohen's class but use different kernel shapes - a direction beyond the scope of this report.

# Appendix A - Signal Taxonomy

*This appendix defines the six signal archetypes that serve as the reference grid for the entire report. Every lab in Volume B constructs model signals from these archetypes; every analysis in Volume C maps real EEG features back to them. Each archetype is presented with its discrete-time formula, key properties, common real-world sources, the DSP techniques suited to it, and its role in the signal-processing narrative of this report. The summary table at the end collects their behaviour under each transform.*

## AA.1 Single Tone

**Mathematical formula.**

$$
x[n] = A \cos\!\left(\frac{2\pi f_0}{f_s} n + \phi\right), \qquad n = 0, 1, \ldots, N-1 \tag{AA.1}
$$

where $A$ is the amplitude, $f_0$ is the frequency in Hz, $f_s$ is the sampling rate, and $\phi$ is the initial phase.

**Properties.**

- **Stationary.** The frequency content does not change over time. Every segment of the signal has the same spectrum.
- **Narrowband.** All energy is concentrated at a single frequency $f_0$ (and its mirror at $-f_0$ for a real signal).
- **Periodic.** The signal repeats with period $P = f_s / f_0$ samples. If $P$ is not an integer, the signal is not exactly periodic within any finite window - the source of off-bin leakage (Section A.3.2).
- **Deterministic.** Fully specified by three parameters ($A$, $f_0$, $\phi$); no randomness.

**DFT behaviour.** The DFT of Equation (AA.1) places energy in the bin(s) nearest to $f_0$. If $f_0$ falls exactly on a bin ($f_0 = k \cdot f_s / N$ for integer $k$), the result is a single spike at bin $k$ with magnitude $AN/2$. If $f_0$ is off-bin, the energy leaks across neighbouring bins through the window's side lobes (Section A.3.2).

**Common sources.**

- EEG: a dominant alpha rhythm (≈10 Hz) in a resting, eyes-closed recording can approximate a single tone over short segments.
- Power-line interference: 50 Hz (or 60 Hz) contamination appears as a persistent single tone in many biomedical recordings.
- Calibration signals: a known-frequency test tone used to verify equipment.

**Techniques suited to this archetype.**

- The **DFT** (Section A.2) is sufficient. A single tone is the archetype the DFT is designed for - it maps directly to one bin (or a small cluster under leakage).
- **Windowing** (Section A.3) matters only for the off-bin case, where it controls side-lobe leakage.
- The **STFT** adds no useful information for a stationary tone - every column looks the same.

**Role in signal processing.** The single tone is the base case. Every transform in this report is easiest to understand when applied to a single tone first. It is the signal for which the DFT is exact, the STFT is redundant, the WVD is clean (except for the DC ghost from real-valued signals), and the SPWVD is trivial. It is also the building block: every other archetype is either a sum of tones, a tone whose parameters change, or the absence of tonal structure.

## AA.2 Mixed Tones

**Mathematical formula.**

$$
x[n] = \sum_{i=1}^{K} A_i \cos\!\left(\frac{2\pi f_i}{f_s} n + \phi_i\right) \tag{AA.2}
$$

where $K$ is the number of components, each with its own amplitude $A_i$, frequency $f_i$, and phase $\phi_i$.

**Properties.**

- **Stationary.** Like the single tone, the frequency content is constant over time.
- **Multi-narrowband.** Energy is concentrated at $K$ discrete frequencies. The spectrum has $K$ spikes (assuming they are resolvable).
- **Resolvability condition.** Two tones at $f_1$ and $f_2$ are distinguishable in the DFT only if $|f_1 - f_2| > \Delta f_{\min} = \beta \cdot f_s / N$ (Equation (A.30)). If they are closer, the window's main lobes overlap and the two tones merge into one.
- **Superposition.** The DFT is linear: $\text{DFT}(x_1 + x_2) = \text{DFT}(x_1) + \text{DFT}(x_2)$. Each tone contributes independently to the spectrum.

**DFT behaviour.** Each tone produces a spike at its corresponding bin. The spikes do not interact in the DFT (linearity). The challenge is purely resolution: are the tones far enough apart for the chosen window to separate them?

**Common sources.**

- EEG: simultaneous rhythms in different bands - alpha (10 Hz) and beta (20 Hz) coexisting in a recording. At the scalp, multiple neural oscillators contribute simultaneously.
- Music: a chord is a sum of tones.
- DTMF (telephone signalling): each key is encoded as a pair of tones.

**Techniques suited to this archetype.**

- The **DFT** resolves them if the separation exceeds $\Delta f_{\min}$. Window choice (Table A.3) determines the resolution limit.
- **Welch's method** (Section A.4.4) gives a statistically reliable estimate of the power at each tone.
- The **STFT** is again redundant for stationary mixed tones - every column is the same.

**Role in signal processing.** Mixed tones test the resolution limit. They are the archetype where the window choice in Section A.3 becomes consequential: rectangular resolves the closest pairs, but leaks the most; Blackman leaks the least, but merges pairs that rectangular could separate. In the WVD, mixed tones introduce **cross-terms** - spurious energy at the midpoint frequency $(f_1 + f_2)/2$ - which do not appear in the DFT. This is the first signal that breaks the WVD and motivates the SPWVD.

## AA.3 Chirp (Fluctuating Tone)

**Mathematical formula (linear chirp).**

$$
x[n] = A \cos\!\left(\frac{2\pi}{f_s}\left(f_0 n + \frac{\mu}{2} \frac{n^2}{f_s}\right) + \phi\right) \tag{AA.3}
$$

where $f_0$ is the starting frequency, $\mu$ is the **chirp rate** in Hz/s (the rate at which the instantaneous frequency changes), and the instantaneous frequency at sample $n$ is:

$$
f_{\text{inst}}[n] = f_0 + \mu \cdot \frac{n}{f_s} \tag{AA.4}
$$

**Properties.**

- **Non-stationary.** The frequency changes over time - the defining property. The signal is not the same in every segment.
- **Narrowband at each instant.** At any given time, the chirp has a well-defined instantaneous frequency. But over the full duration, the frequency sweeps a range, so the global spectrum is broadband.
- **Continuous phase.** The phase function $2\pi(f_0 n + \mu n^2 / (2 f_s)) / f_s$ is smooth - there are no discontinuities. This makes the chirp well-behaved (per the CLAUDE.md signal constraints).

**DFT behaviour.** The global DFT spreads energy across all bins between $f_0$ and $f_0 + \mu T$, producing a smeared band rather than a spike. The DFT cannot tell you *when* each frequency was present - it averages over the entire signal.

**STFT behaviour.** The chirp appears as a **diagonal line** on the spectrogram, with slope equal to the chirp rate $\mu$. However, the line is blurred: the uncertainty principle (Equation (A.49)) smears it along both axes. A short window gives a sharp time position but a thick frequency spread; a long window gives a thin frequency trace but poor time localization.

**WVD behaviour.** The WVD of a single chirp produces a **razor-sharp diagonal** - a delta function along the instantaneous frequency trajectory $f_{\text{inst}}[n]$. No blur, no smearing. This is the WVD's signature result: for a single-component signal, it achieves perfect time-frequency localization.

**Common sources.**

- EEG: frequency modulation within rhythmic activity - an alpha rhythm whose frequency drifts from 9 Hz to 11 Hz over several seconds.
- Neural oscillatory bursts that spin up (increasing frequency) or spin down (decreasing frequency).
- Radar and sonar: chirp pulses are standard for range-Doppler estimation.

**Techniques suited to this archetype.**

- The **STFT** shows the sweep direction and approximate trajectory, limited by the uncertainty principle.
- The **WVD** gives the exact trajectory for a single chirp - the clearest demonstration of the WVD's advantage over the STFT.
- The **DFT** is not useful - it collapses the time axis and smears the chirp into a featureless band.

**Role in signal processing.** The chirp is the archetype that justifies the WVD. It is the signal for which the STFT's blur is most visibly a limitation, and the WVD's sharpness is most visibly a gain. Every comparison between the STFT and WVD in this report starts with the chirp. It is also the archetype closest to real EEG non-stationarity: neural rhythms are rarely fixed-frequency; they modulate.

## AA.4 Multi-Chirp (Mixed Fluctuating Tones)

**Mathematical formula.**

$$
x[n] = \sum_{i=1}^{K} A_i \cos\!\left(\frac{2\pi}{f_s}\left(f_{0,i} \, n + \frac{\mu_i}{2} \frac{n^2}{f_s}\right) + \phi_i\right) \tag{AA.5}
$$

Each component has its own starting frequency $f_{0,i}$, chirp rate $\mu_i$, amplitude $A_i$, and phase $\phi_i$.

**Properties.**

- **Non-stationary and multi-component.** Multiple frequency trajectories evolving simultaneously.
- **The hardest case for time-frequency analysis.** No single linear method (DFT, STFT) can resolve all components sharply, and the WVD's quadratic nature generates cross-terms between every pair of components.

**DFT behaviour.** Overlapping smeared bands - worse than a single chirp, since the swept ranges may overlap in frequency.

**STFT behaviour.** Multiple blurred diagonals. If the chirps cross in time-frequency, the crossing region is particularly difficult to interpret.

**WVD behaviour.** Each component appears as a sharp diagonal (as in Section AA.3), but between every pair of components $(i, j)$, a **cross-term** appears:

- Located at the **midpoint** in both time and frequency between the two components.
- Oscillating at the **difference frequency** $|f_i(t) - f_j(t)|$.
- The cross-term carries as much energy as the auto-terms (the genuine components). It is not a small perturbation.

For $K$ components, there are $K$ auto-terms and $K(K-1)/2$ cross-terms. The cross-terms grow quadratically with the number of components and can dominate the representation.

**Common sources.**

- EEG: simultaneous rhythms in different bands, each with its own frequency modulation - e.g. a drifting alpha and a drifting beta coexisting.
- Multi-speaker audio, multi-target radar, overlapping neural oscillatory bursts.

**Techniques suited to this archetype.**

- The **STFT** shows blurred but cross-term-free diagonals (linearity prevents cross-terms).
- The **WVD** shows sharp diagonals buried in cross-term interference.
- The **SPWVD** (Section A.8) suppresses cross-terms with two independent smoothing knobs while preserving most of the WVD's sharpness - this is the archetype the SPWVD is designed for.

**Role in signal processing.** The multi-chirp is the archetype that breaks the raw WVD. It is the reason the WVD alone is not a practical tool, and the reason the PWVD and SPWVD exist. Every real multi-component signal - including EEG - behaves like a multi-chirp at some level, so understanding cross-terms is not optional.

## AA.5 Transient / Pulse

**Mathematical formula (Gaussian-enveloped burst).**

$$
x[n] = A \exp\!\left(-\frac{(n - n_0)^2}{2\sigma_t^2}\right) \cos\!\left(\frac{2\pi f_0}{f_s} n + \phi\right) \tag{AA.6}
$$

where $n_0$ is the center position, $\sigma_t$ controls the duration (in samples), and $f_0$ is the carrier frequency. The limiting case $\sigma_t \to 0$ approaches an impulse $\delta[n - n_0]$.

A pure impulse (no carrier) is:

$$
x[n] = A \, \delta[n - n_0] \tag{AA.7}
$$

**Properties.**

- **Time-localized.** The energy is concentrated around $n_0$. The shorter the transient, the more localized.
- **Broadband.** By the uncertainty principle (Equation (A.49)), a short signal must have a wide bandwidth. A pure impulse has a flat spectrum - equal energy at all frequencies.
- **Non-stationary.** The signal is not the same in every segment - it exists only briefly.
- **The dual of a tone.** A tone is perfectly localized in frequency and completely delocalized in time. A transient is the opposite: localized in time, spread in frequency. They are the two endpoints of the uncertainty tradeoff.

**DFT behaviour.** A broadband smear - the DFT spreads the transient's energy across many bins. For a pure impulse, $|X[k]| = A$ for all $k$ (flat magnitude spectrum). The DFT cannot tell you *when* the transient occurred.

**STFT behaviour.** A **vertical stripe** in the spectrogram - energy at all frequencies, concentrated in the time column(s) that contain the transient. A short STFT window localizes the stripe tightly; a long window smears it across multiple columns.

**WVD behaviour.** A sharp vertical stripe (auto-term), plus **broadband ghosts** if any other component is present. For a single impulse in isolation, the WVD is clean. For an impulse plus a tone, a cross-term oscillates at the tone's frequency, centered between the impulse time and the tone's presence.

**Common sources.**

- EEG: blink artifacts (sharp transients from eye movement), epileptiform spikes, K-complexes in sleep EEG, electrode pops.
- Biomedical: QRS complex in ECG (a sharp, brief waveform).
- General: clicks, impacts, switching transients.

**Techniques suited to this archetype.**

- The **STFT** with a short window gives the best time localization - but at the cost of frequency resolution (which is irrelevant for a broadband event).
- **Statistical detection** (Section A.4.2): a transient creates a sudden energy spike that exceeds the noise floor across many bins simultaneously. Detecting it is a time-domain energy test, not a frequency-domain peak test.
- The **autocorrelation** of a transient is a spike at lag zero (Equation (A.54)) - no periodicity.

**Role in signal processing.** The transient is the opposite extreme from the tone. It tests time localization rather than frequency resolution. In EEG, transient detection (Part C.4) is primarily about identifying artifacts (blinks, electrode pops) and clinically relevant events (spikes). The transient archetype is also the one where the STFT's short-window mode works best, and where the WVD's cross-terms are most problematic (broadband ghosts from interaction with any coexisting tone or rhythm).

## AA.6 Noise

**Mathematical formula (white Gaussian noise).**

$$
x[n] = \eta[n], \qquad \eta[n] \sim \mathcal{N}(0, \sigma^2) \tag{AA.8}
$$

Each sample is drawn independently from a Gaussian distribution with mean zero and variance $\sigma^2$.

**Properties.**

- **Stationary** (in the statistical sense). The distribution of any sample is the same at all times. But any individual realization looks irregular and non-repeating.
- **Broadband.** The expected power spectrum is flat: $\mathbb{E}[|X[k]|^2] = N\sigma^2$ for all $k$ (Equation (A.35)). Energy is spread equally across all frequencies.
- **No deterministic structure.** No periodicity, no frequency trajectory, no time-localized event. Noise is defined by its statistical properties, not by a formula for $x[n]$.
- **Phase is uniformly random.** Each bin's phase $\angle X[k]$ is independent and uniform on $(-\pi, \pi]$ (Section A.4.1).

**DFT behaviour.** A ragged, flat spectrum: each bin's magnitude fluctuates around the same mean $N\sigma^2$, with the exponential distribution described in Equation (A.34). No bin stands out systematically - but any single realization has random peaks and valleys that can be mistaken for tonal features if statistics are not applied.

**STFT behaviour.** Uniform speckle across the entire time-frequency plane. No horizontal lines, no vertical stripes, no diagonals - just random variation at each $(m, k)$ cell. This is the visual baseline: any structure that appears on top of the speckle is a candidate for a genuine signal feature.

**WVD behaviour.** Speckle with self-interference. The WVD of noise is noisier than the STFT of noise because the quadratic structure creates cross-terms between different noise components. The WVD of noise has higher variance than the spectrogram of noise - the sharpness that helps for tones and chirps hurts for noise.

**SPWVD behaviour.** The smoothing windows reduce the variance of the noise floor but do not remove it. Noise is the one archetype that no amount of time-frequency smoothing can clean up - it has no structure to preserve.

**Common sources.**

- Thermal noise in amplifiers and electrodes (Johnson-Nyquist noise).
- Biological background activity in EEG: the aggregate electrical activity of millions of neurons that do not belong to the rhythm under study.
- Quantization noise from analog-to-digital conversion.
- EMG contamination: muscle activity produces broadband noise-like signals at higher frequencies.

**Techniques suited to this archetype.**

- **Statistical characterization** (Section A.4): the noise floor, the exponential distribution of bin power, and the detection threshold from Equation (A.37).
- **Averaging** (Welch's method, Section A.4.4): reduces variance of the spectral estimate, revealing the true noise level.
- **Autocorrelation** (Section A.6): white noise has $r[l] \approx 0$ for $l \neq 0$, which distinguishes it from any periodic or quasi-periodic signal.

**Role in signal processing.** Noise is not a signal to analyze - it is the **background against which signals are detected**. Every detection, estimation, and resolution result in this report (Equation (A.37), Equation (A.39), Table A.5) is stated relative to the noise. Without a noise model, "this bin has power $P$" is meaningless; with the exponential model from Section A.4.1, it becomes "this bin's power exceeds the noise floor by a factor of $\gamma$ with false-alarm probability $e^{-\gamma}$." Noise defines the floor; signals are what rises above it.

## Summary Table

The following table collects the behaviour of all six archetypes under each transform. Every cell has been justified in the subsections above and in the corresponding theory sections of Volume A. This table is referenced throughout Volume B (labs) and Volume C (EEG application).

**Table AA.1 - Signal archetypes × transform behaviour**

| Archetype | Time-domain form | DFT | STFT | WVD | SPWVD |
| --- | --- | --- | --- | --- | --- |
| Single tone (AA.1) | $A\cos(2\pi f_0 n / f_s + \phi)$ | Single spike at $f_0$ | Horizontal line | Sharp horizontal line (+ DC ghost if real) | Clean line |
| Mixed tones (AA.2) | $\sum A_i \cos(\ldots)$ | $K$ spikes | Parallel horizontal lines | Lines + midpoint cross-term ghosts | Lines; ghosts suppressed |
| Chirp (AA.3) | Swept-frequency sinusoid | Smeared band | Blurred diagonal | Razor-sharp diagonal | Sharp diagonal |
| Multi-chirp (AA.4) | Sum of chirps | Overlapping smears | Crossing blurred diagonals | Sharp diagonals + moving cross-term ghosts | Diagonals; ghosts suppressed |
| Transient / pulse (AA.5) | Short burst or impulse | Flat (broadband) | Vertical stripe | Sharp vertical stripe (+ broadband ghosts) | Vertical stripe; ghosts suppressed |
| Noise (AA.6) | Random process | Flat, random phase | Uniform speckle | Speckle (+ self-interference) | Speckle; not removed by smoothing |

**Reading the table.** Each row is an archetype; each column is a transform. Moving left to right across the columns traces the progression of this report: from the global DFT (no time information) to the STFT (time-frequency, limited by uncertainty) to the WVD (sharp but cross-term-corrupted) to the SPWVD (sharp with cross-terms managed). The diagonal of the table - the cells where each archetype meets its best-suited transform - is where the analysis produces the clearest result. Off-diagonal cells show where a transform struggles with a particular archetype.

**The two failure modes.** The table reveals two systematic failures:

1. **The STFT always blurs.** Every STFT cell uses words like "blurred" or gives less detail than the WVD cell in the same row. This is the uncertainty principle (Equation (A.49)) manifesting across all archetypes.

2. **The WVD always adds ghosts for multi-component signals.** Every WVD cell for a multi-component archetype (mixed tones, multi-chirp, transient coexisting with other components) mentions cross-terms. This is the quadratic nature (Section A.7.3) manifesting across all archetypes.

The SPWVD column shows the compromise: it inherits most of the WVD's sharpness while suppressing most of the cross-terms, at the cost of some smoothing. The one archetype it cannot help is noise - there is no structure to preserve.
