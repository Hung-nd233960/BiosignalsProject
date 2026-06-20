# Volume A — Theoretical Background

## A.2 The DTFT and the DFT

*Gap: Signals & Systems gives you the DTFT as a continuous function of frequency. What it never makes explicit is that you cannot store, compute, or display a continuous function — so the frequency axis must be discretized too. The DFT is that discretization, and a "bin" is the name for each sample on the frequency grid.*

---

### A.2.1 From the DTFT to the DFT

We begin with what the reader already knows. Given a finite-length discrete-time signal $x[n]$ of $N$ samples, the **Discrete-Time Fourier Transform** is defined as:

$$
X(\omega) = \sum_{n=0}^{N-1} x[n] \, e^{-j\omega n} \tag{A.1}
$$

where $\omega$ is the **normalized angular frequency** in radians per sample, and $X(\omega)$ is a continuous, periodic function of $\omega$ with period $2\pi$.

Equation (A.1) is the formula from a Signals & Systems course. It is correct, complete, and — in practice — unusable.

The problem is not on the time side. We already discretized time when we sampled: the signal is a finite list of $N$ numbers. The problem is on the **frequency side**: $X(\omega)$ is defined for every real value of $\omega$ in $[0, 2\pi)$. That is an uncountably infinite set of values. No computer stores a continuous function. No plot has infinite resolution. To actually *use* the frequency-domain representation, we must **sample it**.

The question is: at which frequencies do we evaluate $X(\omega)$?

For a signal whose frequency content is known in advance — say, a single sinusoid at exactly 10 Hz — we could evaluate $X(\omega)$ at just that one frequency and be done. But for a **general signal** — one whose spectral content we do not know beforehand — we have no such shortcut. We need a **systematic grid** of frequency points that covers the full band, so that no feature is missed regardless of where it falls.

The simplest such grid is $N$ equally spaced points around the unit circle:

$$
\omega_k = \frac{2\pi k}{N}, \qquad k = 0, 1, \ldots, N-1 \tag{A.2}
$$

Substituting Equation (A.2) into the DTFT (A.1) gives the **Discrete Fourier Transform**:

$$
X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j 2\pi k n / N}, \qquad k = 0, 1, \ldots, N-1 \tag{A.3}
$$

Equation (A.3) is the DTFT evaluated at $N$ specific frequencies. Nothing more, nothing less. It does not introduce new mathematics — it inherits everything from the DTFT. What it introduces is a **finite, computable representation** of the frequency content.

The inverse transform in Equation (A.4) confirms the relationship is exact — $N$ samples in time, $N$ samples in frequency, no information lost:

$$
x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \, e^{j 2\pi k n / N}, \qquad n = 0, 1, \ldots, N-1 \tag{A.4}
$$

The DTFT and the DFT are not two different transforms. The DFT is the DTFT made computable by choosing where to look.

---

### A.2.2 What a Bin Is

Each value $X[k]$ in the DFT output is called a **bin**. The word is worth understanding precisely, because much of the confusion in spectral analysis comes from treating bins as something they are not.

**A bin is one sample of the continuous DTFT, taken at a specific frequency.**

The DTFT $X(\omega)$ is a continuous curve. The DFT samples this curve at $N$ points. Each sample is a bin. Between the bins, the DTFT still exists — it has values there — but the DFT does not report them.

**Bin spacing in normalized frequency:**

The distance between adjacent bins on the $\omega$-axis is:

$$
\Delta\omega = \frac{2\pi}{N} \quad \text{(radians per sample)} \tag{A.5}
$$

**Bin spacing in physical frequency:**

To convert to Hz, recall that $\omega = 2\pi f / f_s$, where $f$ is the physical frequency in Hz and $f_s$ is the sampling rate. Applying this to Equation (A.5), the spacing between adjacent bins in Hz is:

$$
\Delta f = \frac{f_s}{N} \quad \text{(Hz)} \tag{A.6}
$$

Equation (A.6) is one of the most important numbers in spectral analysis. It depends on **two quantities only**: the sampling rate $f_s$ and the number of samples $N$. Since the signal duration is $T = N / f_s$, we can equivalently write Equation (A.6) as:

$$
\Delta f = \frac{1}{T} \quad \text{(Hz)} \tag{A.7}
$$

The bin spacing is the reciprocal of the observation time. A 1-second signal gives 1 Hz bin spacing. A 20-second signal gives 0.05 Hz bin spacing. This is not a design choice — it is a consequence of sampling the DTFT at $N$ equally spaced points.

**Which physical frequency does bin $k$ represent?**

$$
f_k = k \cdot \Delta f = \frac{k \cdot f_s}{N} \quad \text{(Hz)} \tag{A.8}
$$

For $k = 0$: the DC component (zero frequency).
For $k = N/2$ (if $N$ is even): the Nyquist frequency $f_s / 2$.
For $k > N/2$: these bins correspond to negative frequencies, aliased into the range $(f_s/2, \, f_s)$. For real-valued signals, they are conjugate mirrors of the bins below $N/2$ and carry no new information.

---

### A.2.3 Resolution vs. Bin Count

These two concepts are routinely conflated. They are not the same thing.

**Bin count** is the number of DFT output values: $N$ bins from a length-$N$ signal. It determines how densely the DTFT is sampled.

**Frequency resolution** is the ability to distinguish two nearby frequency components in the signal. It is determined by the **signal duration** $T = N / f_s$, and its fundamental limit is given by Equation (A.9):

$$
\Delta f_{\min} = \frac{1}{T} = \frac{f_s}{N} \quad \text{(Hz)} \tag{A.9}
$$

Two tones separated by less than $\Delta f_{\min}$ cannot be resolved — not because we lack bins, but because the DTFT itself merges them into a single lobe. The underlying continuous spectrum does not contain two distinct peaks.

**Zero-padding** is the operation that exposes this distinction. If we take our $N$-sample signal and append zeros to create a length-$M$ sequence ($M > N$), then compute a length-$M$ DFT:

- The **bin count** increases from $N$ to $M$. The DTFT is now sampled more finely.
- The **bin spacing** decreases from $f_s/N$ to $f_s/M$. We see more points on the same curve.
- The **frequency resolution does not change**. The DTFT curve itself — its shape, its lobes, its ability to separate two tones — is determined entirely by the original $N$ samples. Zero-padding samples the same curve more densely; it does not sharpen it.

Zero-padding is interpolation in the frequency domain: more points on the same curve. It makes plots smoother and peak locations easier to read, but it does not reveal detail that the signal duration did not capture.

To genuinely improve resolution, there is only one option: **record more signal** — increase $N$ at the same $f_s$ (equivalently, increase $T$).

---

### A.2.4 The DFT as an Orthonormal Basis; Parseval's Theorem

The DFT has a clean linear-algebraic interpretation. Define the $N$ basis vectors:

$$
\phi_k[n] = \frac{1}{\sqrt{N}} \, e^{j 2\pi k n / N}, \qquad n = 0, 1, \ldots, N-1 \tag{A.10}
$$

These vectors are **orthonormal**: for any two basis indices $k$ and $m$,

$$
\sum_{n=0}^{N-1} \phi_k[n] \, \phi_m^*[n] = \begin{cases} 1 & \text{if } k = m \\ 0 & \text{if } k \neq m \end{cases} \tag{A.11}
$$

The DFT decomposes $x[n]$ into coefficients on this basis. Each bin $X[k]$ measures how much of the signal aligns with the complex exponential at frequency $\omega_k$. The inverse DFT reconstructs the signal as a weighted sum of these basis vectors. This is a change of basis in $\mathbb{C}^N$ — reversible, lossless, and exact.

**Parseval's theorem** (Equation (A.12)) states that total energy is preserved across the transform:

$$
\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |X[k]|^2 \tag{A.12}
$$

The left side is the total energy computed in the time domain. The right side is the total energy computed from the frequency-domain coefficients, scaled by $1/N$. They are equal.

This is not approximate. It is a direct consequence of the orthonormality of the DFT basis — a unitary transform preserves inner products and, therefore, norms.

The practical consequence: by Equation (A.12), if we want to know how much energy a signal carries in a particular frequency band, we can sum $|X[k]|^2 / N$ over the bins in that band. The answer is exact, not an estimate. This is the foundation of **band-power analysis**, which Volume C will apply to EEG frequency bands (δ, θ, α, β, γ).

---

---

## A.3 Leakage and Windowing

*Gap: Signals & Systems never mentions windowing. The DFT's summation limits silently impose a rectangular window on every signal, and the artifacts this creates — spectral leakage — are invisible unless you know to look for them. This section makes the window explicit and shows how to choose a better one.*

---

### A.3.1 The Hidden Rectangular Window

The DFT in Equation (A.3) sums from $n = 0$ to $n = N - 1$. This means any signal that exists outside that range is silently set to zero. We can make this explicit by defining a **rectangular window**:

$$
w_{\text{rect}}[n] = \begin{cases} 1 & \text{if } 0 \leq n \leq N-1 \\ 0 & \text{otherwise} \end{cases} \tag{A.13}
$$

The DFT of a signal $x[n]$ is actually the DFT of the product $x[n] \cdot w_{\text{rect}}[n]$. We never asked for this multiplication — it is a consequence of having a finite number of samples.

The rectangular window has a sharp edge: at $n = 0$ the signal abruptly begins, and at $n = N - 1$ it abruptly ends. If the signal is not perfectly periodic within the $N$-sample frame — meaning if the last sample does not connect smoothly back to the first — then this abrupt cut creates a discontinuity at the boundary. The DFT, which implicitly treats the signal as periodic, sees this discontinuity as part of the signal.

This is the origin of **spectral leakage**: the discontinuity introduces frequency content that was not present in the original signal.

---

### A.3.2 Leakage: The Window's Spectrum Convolved onto the Signal

Multiplication in the time domain is convolution in the frequency domain. The DTFT of the windowed signal $x[n] \cdot w[n]$ is:

$$
X_w(\omega) = \frac{1}{2\pi} \int_{-\pi}^{\pi} X(\theta) \, W(\omega - \theta) \, d\theta \tag{A.14}
$$

where $W(\omega)$ is the DTFT of the window itself. Each frequency component in $X(\theta)$ gets smeared by the shape of $W(\omega)$.

In discrete terms, this means: the DFT does not show the signal's true spectrum. It shows the signal's spectrum convolved with the window's spectrum. Every tone in the signal is replaced by a copy of $W(\omega)$ centered at that tone's frequency.

The DTFT of the rectangular window is the **Dirichlet kernel** (derived in Appendix B.1):

$$
W_{\text{rect}}(\omega) = e^{-j\omega(N-1)/2} \cdot \frac{\sin(\omega N / 2)}{\sin(\omega / 2)} \tag{A.15}
$$

Equation (A.15) has a tall main lobe centered at $\omega = 0$ and a series of side lobes that decay slowly on either side. These side lobes are the source of leakage: energy from a tone at one frequency spills into bins at other frequencies through the side lobes of $W_{\text{rect}}(\omega)$.

**On-bin vs. off-bin placement.** If a tone's frequency falls exactly on a DFT bin (i.e. $f = k \cdot f_s / N$ for integer $k$), the main lobe of $W(\omega)$ lands exactly on that bin and all other bins fall on the nulls of the Dirichlet kernel. Leakage vanishes — by coincidence of alignment, not by design. If the tone is even slightly off-bin, the nulls no longer align with the bins, and leakage appears across the entire spectrum. This sensitivity to placement is a fundamental fragility of the rectangular window.

---

### A.3.3 The Cosine-Sum Window Family

The goal of a non-rectangular window is to reduce side-lobe levels — accepting a wider main lobe in exchange. The most widely used family is the **cosine-sum** windows, defined by the general formula:

$$
w[n] = \sum_{p=0}^{P} (-1)^p \, a_p \cos\!\left(\frac{2\pi p n}{N}\right), \qquad n = 0, 1, \ldots, N-1 \tag{A.16}
$$

Each member of the family is determined by its coefficients $a_0, a_1, \ldots, a_P$. Table A.1 lists the four standard members, their coefficients, and their spectral properties.

**Table A.1 — Cosine-sum window family: coefficients and spectral properties**

| Property | Rectangular | Hann | Hamming | Blackman |
| --- | --- | --- | --- | --- |
| Order $P$ | 0 | 1 | 1 | 2 |
| $a_0$ | 1 | 0.5 | 0.54 | 0.42 |
| $a_1$ | — | 0.5 | 0.46 | 0.5 |
| $a_2$ | — | — | — | 0.08 |
| Edge value $w[0] = w[N{-}1]$ | 1 | 0 | 0.08 | 0 |
| Main-lobe width (bins) | 2 | 4 | 4 | 6 |
| $\beta$ (resolution factor) | 1 | 2 | 2 | 3 |
| Peak side-lobe (dB) | −13 | −31.5 | −42.7 | −58 |
| Side-lobe rolloff | 6 dB/oct | 18 dB/oct | 6 dB/oct | 18 dB/oct |
| Rolloff order | $1/\omega^1$ | $1/\omega^3$ | $1/\omega^1$ | $1/\omega^3$ |

The individual formulas, expanded from Equation (A.16), are:

**Rectangular** ($P = 0$):

$$
w_{\text{rect}}[n] = 1 \tag{A.17}
$$

**Hann** ($P = 1$):

$$
w_{\text{Hann}}[n] = 0.5 - 0.5\cos\!\left(\frac{2\pi n}{N}\right) \tag{A.18}
$$

**Hamming** ($P = 1$):

$$
w_{\text{Hamming}}[n] = 0.54 - 0.46\cos\!\left(\frac{2\pi n}{N}\right) \tag{A.19}
$$

**Blackman** ($P = 2$):

$$
w_{\text{Blackman}}[n] = 0.42 - 0.5\cos\!\left(\frac{2\pi n}{N}\right) + 0.08\cos\!\left(\frac{4\pi n}{N}\right) \tag{A.20}
$$

The design logic behind each window is visible in Table A.1:

- **Rectangular** is the baseline. No shaping, narrowest main lobe, but the worst side lobes (−13 dB) and the slowest rolloff ($1/\omega$). Every off-bin tone leaks broadly.
- **Hann** goes to zero at both edges ($w[0] = w[N{-}1] = 0$). This eliminates the boundary discontinuity in both the window value and its first derivative, producing an $18$ dB/oct rolloff. The cost is a main lobe twice as wide as rectangular.
- **Hamming** is optimized to minimize the peak side-lobe level (−42.7 dB) rather than the rolloff rate. It does not go to zero at the edges ($w[0] = 0.08$), so the boundary discontinuity remains in the value — hence the rolloff rate stays at $6$ dB/oct, same as rectangular. But the specific choice of $a_0 = 0.54$, $a_1 = 0.46$ places a cancellation that pushes the nearest side lobes far below those of Hann.
- **Blackman** uses a second cosine term ($P = 2$) and goes to zero at both edges. The result is the lowest peak side-lobe (−58 dB) with $18$ dB/oct rolloff, but the widest main lobe (6 bins). It trades the most resolution for the cleanest spectrum.

---

### A.3.4 Main-Lobe vs. Side-Lobe Tradeoff; the Rolloff Rule

Table A.1 shows the central tradeoff of window design: **narrower main lobe ↔ higher side lobes**. No window escapes this. A narrow main lobe preserves frequency resolution but lets leakage through; a wide main lobe suppresses leakage but blurs nearby tones together.

The **side-lobe rolloff rate** is governed by the smoothness of the window at its edges:

- If $w[n]$ has a **discontinuity in value** at the boundary (rectangular, Hamming), the rolloff is $1/\omega^1$ — i.e. 6 dB per octave.
- If $w[n]$ goes to zero at both edges and has a **continuous first derivative** (Hann, Blackman), the rolloff is $1/\omega^3$ — i.e. 18 dB per octave.
- In general, each additional order of continuity at the boundary adds $1/\omega^2$ to the rolloff (12 dB/oct per derivative that vanishes).

This is the $1/\omega^p$ rule referenced in the table of contents. It connects a time-domain property (edge smoothness) to a frequency-domain consequence (how fast distant leakage falls off). The rule explains why Hamming, despite its low peak side-lobe, has the same rolloff as rectangular: its edge value is $0.08$, not $0$, so the value discontinuity persists.

---

### A.3.5 The Resolution Limit

In Section A.2.3, we stated that the frequency resolution is $\Delta f_{\min} = f_s / N$. That was for the rectangular window. With a general window, the resolution limit becomes:

$$
\Delta f_{\min} \approx \beta \cdot \frac{f_s}{N} \tag{A.21}
$$

where $\beta$ is the **main-lobe width factor** — the half-width of the main lobe measured in bins. The values of $\beta$ for the standard windows are listed in Table A.1.

Equation (A.21) is the practical resolution limit. Two tones separated by less than $\beta \cdot f_s / N$ in frequency will merge into a single lobe regardless of any other processing. The rectangular window ($\beta = 1$) gives the best resolution; Blackman ($\beta = 3$) gives the worst. This is not a deficiency of Blackman — it is the price of its −58 dB side-lobe suppression.

For the EEG context: at $f_s = 250$ Hz with a $T = 20$ s window ($N = 5000$), the bin spacing is $\Delta f = 0.05$ Hz. The resolution limits are:

**Table A.2 — Resolution limits for standard windows at $f_s = 250$ Hz, $N = 5000$**

| Window | $\beta$ | $\Delta f_{\min}$ (Hz) | Can resolve adjacent EEG bands? |
| --- | --- | --- | --- |
| Rectangular | 1 | 0.05 | Yes — all standard bands |
| Hann | 2 | 0.10 | Yes — all standard bands |
| Hamming | 2 | 0.10 | Yes — all standard bands |
| Blackman | 3 | 0.15 | Yes — all standard bands |

Table A.2 shows that with a 20-second observation window, even the widest window (Blackman) resolves 0.15 Hz — far finer than the narrowest EEG band gap (δ–θ boundary at 4 Hz). The lab constraint of $T \geq 1200$ s makes this margin even more comfortable. The window choice for EEG is therefore driven entirely by side-lobe suppression, not by resolution.

---

### A.3.6 The Gaussian Window (deferred)

The cosine-sum family is not the only option. The **Gaussian window**:

$$
w_{\text{Gauss}}[n] = \exp\!\left(-\frac{1}{2}\left(\frac{n - (N-1)/2}{\sigma}\right)^2\right), \qquad n = 0, 1, \ldots, N-1 \tag{A.22}
$$

where $\sigma$ controls the width, has a unique property: its Fourier transform is also a Gaussian. This makes it the **minimum-uncertainty window** — the one that achieves equality in the time-frequency uncertainty bound $\Delta t \cdot \Delta f \geq 1/(4\pi)$.

We flag it here because it will reappear in A.8, where the SPWVD's smoothing windows are chosen. The Gaussian's minimum-uncertainty property makes it the natural choice for time-frequency smoothing, where the tradeoff between time and frequency resolution must be managed per axis. The full discussion is deferred to that section.

---

---

## A.4 Statistics and the DFT

*Gap: Signals & Systems treats the DFT output as a deterministic, exact description of the signal. In practice, every recorded signal contains noise, and each DFT bin becomes a random variable. This section develops the statistical framework needed to answer the only question that matters: "is this spectral feature real, or is it noise?" — answered from inside the spectral framework, not by applying generic statistics to the time-domain signal.*

---

### A.4.1 Each Bin as a Random Variable

Consider a signal that is pure white Gaussian noise: $x[n] \sim \mathcal{N}(0, \sigma_x^2)$, with each sample independent. What does its DFT look like?

Each DFT bin is a weighted sum of $N$ independent Gaussian samples (Equation (A.3)). By linearity, $X[k]$ is itself a complex Gaussian random variable. Its real and imaginary parts are independent, each with variance $N\sigma_x^2 / 2$:

$$
\text{Re}\{X[k]\}, \; \text{Im}\{X[k]\} \sim \mathcal{N}\!\left(0, \; \frac{N\sigma_x^2}{2}\right) \tag{A.23}
$$

The **magnitude** $|X[k]|$ follows a Rayleigh distribution, and the **phase** $\angle X[k]$ is uniformly distributed on $(-\pi, \pi]$. The phase carries no information about the noise — it is purely random, uniformly spread.

The quantity we care about for spectral analysis is the **power** at each bin:

$$
P[k] = |X[k]|^2 = \text{Re}\{X[k]\}^2 + \text{Im}\{X[k]\}^2 \tag{A.24}
$$

Since $P[k]$ is the sum of squares of two independent Gaussians, it follows an **exponential distribution**:

$$
P[k] \sim \text{Exponential}\!\left(\lambda = \frac{1}{N\sigma_x^2}\right) \tag{A.25}
$$

with mean $\mathbb{E}[P[k]] = N\sigma_x^2$ and — crucially — standard deviation also equal to $N\sigma_x^2$. The coefficient of variation is 1: the standard deviation equals the mean. This means a single bin's power estimate is inherently noisy, with 100% relative uncertainty, regardless of how many samples $N$ we use.

Now add a deterministic tone to the noise: $x[n] = A\cos(2\pi f_0 n / f_s) + \eta[n]$. At the bin $k_0$ nearest to $f_0$, the DFT output is a complex Gaussian shifted by the tone's contribution. The power $P[k_0]$ now follows a **non-central chi-squared distribution** — its mean is elevated above the noise-only level by the tone's energy. At all other bins (away from $f_0$), the distribution remains exponential as in Equation (A.25).

This is the statistical model. A tone reveals itself as a bin whose power is improbably large under the null hypothesis of noise alone.

---

### A.4.2 The Noise Floor and Spectral Detection

The **noise floor** is the expected power level of bins that contain only noise. From Equation (A.25), it is:

$$
P_{\text{floor}} = \mathbb{E}[P[k]] = N\sigma_x^2 \tag{A.26}
$$

In practice, we do not know $\sigma_x^2$ in advance. We estimate $P_{\text{floor}}$ from the spectrum itself. A robust estimator is the **median** of the bin powers:

$$
\hat{P}_{\text{floor}} = \text{median}\{P[0], P[1], \ldots, P[N/2]\} \tag{A.27}
$$

The median is preferred over the mean because a few bins containing genuine signal components inflate the mean but barely affect the median.

**Detection criterion.** Under the exponential model in Equation (A.25), the probability that a noise-only bin exceeds a threshold $\gamma \cdot \hat{P}_{\text{floor}}$ is:

$$
\Pr\!\left(P[k] > \gamma \cdot P_{\text{floor}}\right) = e^{-\gamma} \tag{A.28}
$$

Equation (A.28) is the detection rule, derived from the spectral model — not from a generic "$2\sigma$" rule. The threshold $\gamma$ is chosen to control the false-alarm rate:

**Table A.3 — Detection thresholds from the exponential noise model**

| Threshold $\gamma$ | False-alarm probability $e^{-\gamma}$ | Interpretation |
| --- | --- | --- |
| 3 | 0.050 | 5% of noise-only bins exceed this |
| 4.6 | 0.010 | 1% false-alarm rate |
| 6.9 | 0.001 | 0.1% false-alarm rate |

The key difference from generic statistics: the threshold is not "$2\sigma$ because that is common." It is $\gamma \cdot \hat{P}_{\text{floor}}$ because the exponential distribution of bin power under noise gives us Equation (A.28), and we choose $\gamma$ to match the false-alarm probability we can tolerate. The noise floor $\hat{P}_{\text{floor}}$ is estimated from the spectrum, so the threshold adapts to the actual noise level in the recording — not to an assumed distribution of time-domain samples.

This also explains why crude time-domain thresholding (e.g. "remove samples beyond $2\sigma$") is problematic: it operates on samples, not on frequency content. Removing a time-domain sample affects every bin in the DFT (by Equation (A.3), every bin is a sum over all samples). The energy removed is spread unpredictably across the spectrum, breaking the Parseval relationship in Equation (A.12) in ways that are difficult to characterize.

---

### A.4.3 The Periodogram Variance Problem

The **periodogram** is the name for the power spectrum estimated from a single DFT:

$$
\hat{S}[k] = \frac{1}{N} |X[k]|^2 \tag{A.29}
$$

Equation (A.29) is the natural estimator: compute the DFT, take the magnitude squared, scale by $1/N$. It is an **unbiased** estimator of the true power spectrum — its expected value equals the true spectral density at each bin.

But it is not a **consistent** estimator. This is the surprising and counterintuitive result:

$$
\text{Var}\!\left(\hat{S}[k]\right) \approx S[k]^2 \tag{A.30}
$$

The variance of the periodogram at bin $k$ is proportional to the square of the true spectral value — and this does not decrease as $N$ increases. Doubling the signal length $N$ halves the bin spacing (Equation (A.6)) and doubles the number of bins, but each bin's power estimate remains just as noisy.

This is a direct consequence of the exponential distribution in Equation (A.25): the coefficient of variation (standard deviation / mean) is always 1, independent of $N$. More data gives us finer frequency spacing, not more reliable power estimates.

In plain terms: a single DFT of a long signal gives a ragged, noisy spectrum, and making the signal longer does not smooth it. The raggedness is statistical, not a failure of implementation. Any single realization of a noisy signal produces a periodogram that fluctuates wildly around the true spectral shape.

---

### A.4.4 Averaging to Reduce Variance: Welch's Method

The periodogram variance problem has a direct solution: **average multiple independent periodograms**. If we have $L$ independent estimates $\hat{S}_1[k], \hat{S}_2[k], \ldots, \hat{S}_L[k]$, the averaged estimate:

$$
\bar{S}[k] = \frac{1}{L} \sum_{\ell=1}^{L} \hat{S}_\ell[k] \tag{A.31}
$$

has variance reduced by a factor of $L$:

$$
\text{Var}\!\left(\bar{S}[k]\right) = \frac{1}{L} \, S[k]^2 \tag{A.32}
$$

**Welch's method** obtains the $L$ independent estimates from a single recording by dividing the signal into overlapping segments, windowing each segment, computing its periodogram, and averaging:

1. Divide the $N$-sample signal into segments of length $M$, with overlap of $D$ samples.
2. Apply a window $w[n]$ (from Section A.3) to each segment.
3. Compute the periodogram of each windowed segment.
4. Average all $L$ periodograms.

The number of segments is:

$$
L = \left\lfloor \frac{N - M}{M - D} \right\rfloor + 1 \tag{A.33}
$$

The tradeoff is explicit and quantitative:

- Each segment has $M < N$ samples, so the bin spacing per segment is $f_s / M$ — coarser than the $f_s / N$ achievable from the full signal.
- The frequency resolution degrades from $\Delta f = f_s / N$ to $\Delta f = \beta \cdot f_s / M$ (Equation (A.21)).
- The variance decreases by the factor $1/L$.

**Table A.4 — Welch's method: the resolution-variance tradeoff at $f_s = 250$ Hz, $N = 300\,000$ (1200 s), Hann window ($\beta = 2$), 50% overlap**

| Segment length $M$ | Segments $L$ | $\Delta f$ (Hz) | Relative variance |
| --- | --- | --- | --- |
| 300 000 (full) | 1 | 0.0017 | 1.00 |
| 30 000 (120 s) | 19 | 0.017 | 0.053 |
| 5 000 (20 s) | 119 | 0.10 | 0.0084 |
| 1 250 (5 s) | 479 | 0.40 | 0.0021 |
| 500 (2 s) | 1 199 | 1.0 | 0.00083 |

Table A.4 shows the tradeoff for our lab signal parameters. A single periodogram of the full 1200-second signal gives the finest frequency resolution (0.0017 Hz) but with 100% relative variance — the spectrum is ragged and unreliable. At the other extreme, 2-second segments give 1199 averages and a smooth spectrum, but the resolution is only 1 Hz — too coarse to separate neighboring EEG sub-bands.

The practical choice for EEG depends on what needs to be resolved. For separating the standard bands (δ, θ, α, β, γ), the boundaries are at 4, 8, 13, and 30 Hz — separations of at least 4 Hz. A segment length of $M = 1250$ (5 s, $\Delta f = 0.40$ Hz) resolves these boundaries comfortably with nearly 500 averages. For finer structure within a band — say, distinguishing 9 Hz from 11 Hz alpha — longer segments are needed ($M = 5000$, giving $\Delta f = 0.10$ Hz).

Welch's method does not bypass the uncertainty in each bin; it trades frequency detail for statistical reliability. The total energy is still conserved (Equation (A.12) applies to each segment), and the noise floor estimated from the averaged spectrum is far more stable than from a single periodogram. Detection via Equation (A.28) applied to a Welch-averaged spectrum is correspondingly more reliable.

---

---

## A.5 The STFT and Spectrograms

*Gap: Signals & Systems never introduces the time-frequency plane. The DFT gives frequency content but discards all timing; Welch's method (A.4.4) averages over time explicitly. For any signal whose frequency content changes — an EEG rhythm that comes and goes, a chirp that sweeps — we need both axes simultaneously. The STFT is the first tool that provides this, and the spectrogram is the first usable time-frequency representation in this report.*

---

### A.5.1 The Windowed DFT Slid in Time

The STFT is built from exactly the components already developed: a window (Section A.3), a DFT (Section A.2), and the sliding-segment idea from Welch (Section A.4.4). The difference is that Welch averages the segment spectra into one; the STFT **keeps each segment spectrum indexed by its position in time**.

Given a signal $x[n]$, a window $w[n]$ of length $M$, and a hop size $H$ (the number of samples the window advances per step), the STFT is:

$$
X[m, k] = \sum_{n=0}^{M-1} x[n + mH] \, w[n] \, e^{-j 2\pi k n / M}, \qquad k = 0, 1, \ldots, M-1 \tag{A.34}
$$

where $m$ is the **segment index** (the time position) and $k$ is the **frequency bin** within that segment. Each value $X[m, k]$ is one DFT bin of one windowed segment.

The output of Equation (A.34) is a **two-dimensional grid** — the **time-frequency plane**:

- The **time axis** is discrete, with spacing $\Delta t_{\text{grid}} = H / f_s$ seconds between columns.
- The **frequency axis** is discrete, with spacing $\Delta f = f_s / M$ Hz between rows (Equation (A.6), applied to each segment of length $M$).
- Each grid cell contains a complex number $X[m, k]$.

The **spectrogram** is the squared magnitude of the STFT:

$$
S[m, k] = |X[m, k]|^2 \tag{A.35}
$$

Equation (A.35) is the power at time step $m$ and frequency bin $k$. It is the object we plot: a 2D image with time on the horizontal axis, frequency on the vertical axis, and color representing power.

The physical frequency and time of each cell are:

$$
f_k = \frac{k \cdot f_s}{M} \quad \text{(Hz)}, \qquad t_m = \frac{m \cdot H}{f_s} \quad \text{(s)} \tag{A.36}
$$

The connection to Welch is now explicit. Welch computes $S[m, k]$ for all $m$ and $k$, then averages over $m$: $\bar{S}[k] = \frac{1}{L}\sum_m S[m,k]$. The STFT keeps the $m$ axis. Welch answers "what frequencies are present on average?"; the STFT answers "what frequencies are present, and when?"

---

### A.5.2 The Uncertainty Principle

This is the central tradeoff of the entire report.

The window length $M$ controls both time resolution and frequency resolution, but in opposite directions:

**Time resolution.** The STFT localizes a signal event to the window that contains it. The finest time localization is the window duration:

$$
\Delta t = \frac{M}{f_s} \quad \text{(s)} \tag{A.37}
$$

A short window (small $M$) gives fine time resolution: a 1-second burst of alpha rhythm lands in one or two columns. A long window (large $M$) smears the burst across a single column that also contains whatever came before and after.

**Frequency resolution.** Within each window, the DFT resolves frequencies to (from Equation (A.21)):

$$
\Delta f = \beta \cdot \frac{f_s}{M} \quad \text{(Hz)} \tag{A.38}
$$

A long window (large $M$) gives fine frequency resolution: a 10 Hz tone and a 10.5 Hz tone appear as separate peaks. A short window (small $M$) merges them into one lobe.

Equations (A.37) and (A.38) compete. Their product:

$$
\Delta t \cdot \Delta f = \beta \tag{A.39}
$$

For the rectangular window ($\beta = 1$), this is already $\geq 1$. For Hann ($\beta = 2$), it is $\geq 2$. The fundamental lower bound, achievable only by the Gaussian window (Equation (A.22)), is:

$$
\Delta t \cdot \Delta f \geq \frac{1}{4\pi} \tag{A.40}
$$

Equation (A.40) is the **Heisenberg–Gabor uncertainty principle** for time-frequency analysis. It is not a limitation of the STFT algorithm. It is a mathematical property of Fourier pairs: a function and its Fourier transform cannot both be arbitrarily concentrated. No linear time-frequency representation — no matter how cleverly constructed — escapes this bound.

**What the uncertainty principle means in practice.** The window length $M$ is a single slider between two extremes. You do not choose "good" or "bad" — you choose *which axis to sacrifice*:

**Table A.5 — The uncertainty tradeoff at $f_s = 250$ Hz, Hann window ($\beta = 2$)**

| Window $M$ (samples) | Duration $\Delta t$ (s) | Freq. resolution $\Delta f$ (Hz) | $\Delta t \cdot \Delta f$ | What it can resolve | What it cannot resolve |
| --- | --- | --- | --- | --- | --- |
| 125 (0.5 s) | 0.5 | 4.0 | 2.0 | Sub-second transients, blinks | Cannot separate α (8–13) from θ (4–8) cleanly |
| 250 (1.0 s) | 1.0 | 2.0 | 2.0 | Second-scale events; α vs. θ bands | Cannot resolve structure within a band |
| 500 (2.0 s) | 2.0 | 1.0 | 2.0 | Individual Hz within bands | Cannot localize events shorter than 2 s |
| 1 250 (5.0 s) | 5.0 | 0.4 | 2.0 | Fine spectral detail (9 vs. 11 Hz) | Cannot localize events shorter than 5 s |
| 5 000 (20.0 s) | 20.0 | 0.1 | 2.0 | Sub-Hz spectral structure | Cannot track any temporal variation faster than 20 s |

Table A.5 makes the tradeoff concrete for EEG. Every row has the same $\Delta t \cdot \Delta f = 2.0$ (the Hann window's $\beta$). Moving down the table buys frequency detail and costs time detail, in exact proportion. No row is "best" — each is suited to a different question about the signal.

The last row illustrates the extreme: a 20-second window resolves 0.1 Hz, but the spectrogram has one column per 20 seconds — it cannot track any rhythm that changes faster than that. At the other extreme, a 0.5-second window catches sub-second transients but cannot distinguish alpha from theta. The STFT forces you to decide what matters before you compute.

This is the limitation that motivates the second half of the report. The WVD (Section A.7) is not a linear representation — it is quadratic — and therefore is not bound by Equation (A.40) in the same way. It achieves sharper localization on both axes simultaneously, at a different cost (cross-terms). The SPWVD (Section A.8) manages that cost with two independent smoothing knobs. But the uncertainty principle remains the reference: any improvement the WVD family offers is measured against the STFT's bound.

---

### A.5.3 Hop Size, Overlap, and the COLA Condition

The STFT in Equation (A.34) has two parameters beyond the window itself: the **hop size** $H$ and, equivalently, the **overlap** $D = M - H$. These are often treated as minor implementation details. They are not. Without proper overlap, the window's tapering creates a systematic problem: **samples near segment boundaries are suppressed or lost**.

#### The tapering problem

Every non-rectangular window tapers toward zero at its edges. From Table A.1, a Hann window satisfies $w[0] = w[N-1] = 0$: the first and last samples of each segment are multiplied by zero. A Blackman window does the same. Even a Hamming window, which does not reach zero, reduces edge samples to 8% of their original amplitude.

Consider what happens when adjacent segments do not overlap ($H = M$, zero overlap). Each sample of the original signal appears in exactly one segment. But within that segment, the window assigns it a weight:

- Samples near the center receive weight $\approx 1$.
- Samples near the edges receive weight $\approx 0$ (for Hann/Blackman) or weight $\approx 0.08$ (for Hamming).

The result: samples that happen to fall near a segment boundary are effectively erased from the analysis. Their contribution to every frequency bin in that segment is multiplied by a number close to zero. This is not a minor effect — for a Hann window, the samples at positions $n = 0$ and $n = M - 1$ are multiplied by exactly zero. They contribute nothing.

If the signal has an important feature (a spike, a burst onset, a transient) at a segment boundary, that feature is suppressed. Not by filtering, not by noise, but by the window's taper. The feature exists in the data; the STFT simply does not see it.

#### How overlap solves the tapering problem

Overlap means that adjacent segments share some samples. With hop size $H < M$, each sample appears in multiple segments. In each segment, it sits at a different position relative to the window center, and therefore receives a different weight.

The critical question is: **what is the total weight a sample receives, summed across all segments it appears in?**

For a sample at position $n$ in the original signal, this total weight is:

$$
W_{\text{total}}[n] = \sum_{m} w[n - mH] \tag{A.41}
$$

where the sum runs over all segment indices $m$ for which the sample falls within the window. If $W_{\text{total}}[n]$ varies with $n$, then different samples receive different total representation in the STFT — some parts of the signal are emphasized, others are suppressed. This is an amplitude modulation imposed by the analysis, not by the signal.

**Example: Hann window with no overlap ($H = M$).** Each sample appears in exactly one segment. $W_{\text{total}}[n] = w[n \bmod M]$, which varies from 0 at the edges to 1 at the center. Samples at segment boundaries are lost.

**Example: Hann window with 50% overlap ($H = M/2$).** Each sample appears in exactly two segments. In one segment, it sits at position $p$; in the adjacent segment, it sits at position $p + M/2$. The Hann window satisfies:

$$
w_{\text{Hann}}[p] + w_{\text{Hann}}[p + M/2] = 1 \quad \text{for all } p \tag{A.42}
$$

Equation (A.42) states that the two overlapping Hann windows sum to a constant. Every sample receives the same total weight of 1, regardless of where it falls relative to segment boundaries. The tapering problem is completely eliminated: no sample is suppressed, no sample is over-represented.

#### The COLA condition

Equation (A.42) is a specific instance of the **Constant Overlap-Add (COLA) condition**:

$$
\sum_{m} w[n - mH] = C \quad \text{for all } n \tag{A.43}
$$

where $C$ is a constant (typically normalized to 1). When COLA is satisfied:

1. **Every sample receives equal total weight.** The analysis does not favor any time position over another.
2. **The STFT is invertible.** The original signal can be perfectly reconstructed from the STFT coefficients by overlap-adding the inverse-DFT of each segment. This means no information is lost — the time-frequency plane is a complete, lossless representation.
3. **Parseval's theorem extends to the STFT.** The total energy computed from the spectrogram equals the total signal energy (with appropriate normalization).

**Table A.6 — COLA-satisfying overlap for standard windows**

| Window | Minimum COLA overlap | Hop size $H$ | Segments per sample |
| --- | --- | --- | --- |
| Rectangular | 0% | $M$ | 1 |
| Hann | 50% | $M/2$ | 2 |
| Hamming | 50% | $M/2$ | 2 (approximate COLA) |
| Blackman | 67% | $M/3$ | 3 |

Table A.6 shows the minimum overlap each window needs to satisfy COLA. The rectangular window satisfies COLA trivially at any overlap (including zero) because it does not taper — but it has the worst spectral leakage (Section A.3). Hann at 50% is the standard practical choice: it satisfies COLA exactly, requires only 2× the computation of zero-overlap, and has −31.5 dB peak side-lobes.

Hamming at 50% satisfies COLA only approximately — the total weight oscillates slightly around a constant because $w[0] = 0.08 \neq 0$. In practice the deviation is small (< 0.2%), and 50% overlap is used regardless. Blackman requires 67% overlap (each sample appears in 3 segments) to achieve COLA, which increases computation by 3× relative to zero-overlap.

#### Oversampling in time: more columns, not more resolution

Increasing the overlap beyond the COLA minimum (e.g. using 75% overlap with a Hann window instead of 50%) produces more time-axis columns in the spectrogram. This is the **time-axis analogue of zero-padding** (Section A.2.3):

- More columns = finer time-axis grid = smoother-looking spectrogram.
- **Time resolution does not improve.** The resolution is set by the window length $M$ (Equation (A.37)), not by the hop size. A smaller hop interpolates the same information onto a denser grid, but does not reveal temporal detail that the window duration cannot capture.

The parallel is exact: zero-padding adds frequency bins without adding frequency resolution; overlap beyond COLA adds time columns without adding time resolution. Both are interpolation — useful for visualization, but not for resolving features.

---

### A.5.4 Reading a Spectrogram; the Cost of Each Parameter Choice

The spectrogram in Equation (A.35) is a 2D image. Reading it requires knowing what signal features look like on the time-frequency plane. The mapping follows directly from the signal taxonomy in Appendix A:

**Table A.7 — Signal archetypes on the spectrogram**

| Archetype | Time-domain appearance | Spectrogram signature |
| --- | --- | --- |
| Stationary tone | Constant-frequency sinusoid | Horizontal line at $f_0$ |
| Mixed tones | Sum of sinusoids | Parallel horizontal lines |
| Chirp | Frequency-swept sinusoid | Diagonal line (slope = sweep rate) |
| Transient / pulse | Short burst or impulse | Vertical stripe (broadband, localized in time) |
| Noise | Random, no structure | Uniform speckle across the plane |

Table A.7 is the key to interpretation: horizontal = frequency-stable, vertical = time-localized, diagonal = frequency-changing. Any real signal — including EEG — is a mixture of these archetypes, and the spectrogram decomposes it into their superposition on the time-frequency plane.

**The cost of each parameter choice.** Choosing the STFT parameters means deciding in advance which features matter most. There is no universal setting. Table A.8 summarizes the practical consequences:

**Table A.8 — STFT parameter choices and their consequences**

| Parameter | Increasing it... | Improves... | Degrades... |
| --- | --- | --- | --- |
| Window length $M$ | Longer segments | Frequency resolution (Eq. A.38) | Time resolution (Eq. A.37) |
| Hop size $H$ | Fewer, sparser columns | Computation speed | Time-axis smoothness (not resolution) |
| Overlap $M - H$ | More columns, better COLA | Sample coverage; visual smoothness | Computation cost; redundancy |
| Window taper (e.g. Rect → Hann → Blackman) | Stronger taper | Side-lobe suppression (Table A.1) | Main-lobe width (frequency resolution) |

**The honest conclusion.** The STFT is the first tool in this report that can answer "what frequency is present at what time." It is the workhorse of practical time-frequency analysis, and for many signals it is sufficient. But the uncertainty principle (Equation (A.40)) means it always blurs something: either time or frequency, controlled by $M$, with no way to have both.

For a signal whose features all live at roughly the same time-frequency scale — say, an EEG recording dominated by a steady alpha rhythm — one window length captures everything. But EEG often contains features at multiple scales simultaneously: a narrow-band alpha rhythm (needs long $M$ to resolve) and a short transient artifact (needs short $M$ to localize). No single STFT window captures both. Two separate STFTs with different $M$ can show each feature, but there is no principled way to merge them into a single representation.

This is the limitation that the Wigner-Ville Distribution addresses. The WVD (Section A.7) is not bound by the uncertainty principle in the same way because it is quadratic, not linear. It achieves sharper localization on both axes — but at a cost (cross-terms) that requires its own management (Sections A.7–A.8). Before reaching the WVD, we need one more tool: the autocorrelation function (Section A.6), which provides the bridge.

---

---

## A.6 Autocorrelation

*Gap: autocorrelation is sometimes mentioned in Signals & Systems but rarely connected to spectral analysis. This section introduces it as the bridge between the DFT/STFT world (Sections A.2–A.5) and the Wigner-Ville family (Sections A.7–A.8). The two ideas planted here — the Wiener-Khinchin theorem and the phase-blindness of the power spectrum — are what make the WVD definition natural rather than arbitrary.*

---

### A.6.1 Definition; the Periodicity Detector

The **autocorrelation** of a discrete signal $x[n]$ of length $N$ is defined as:

$$
r[l] = \sum_{n=0}^{N-1-|l|} x[n] \, x^*[n - l], \qquad l = -(N-1), \ldots, 0, \ldots, N-1 \tag{A.44}
$$

where $l$ is the **lag** (the shift in samples) and $x^*$ denotes the complex conjugate. For real-valued signals, $x^* = x$ and the conjugate can be dropped.

Equation (A.44) compares the signal with a shifted copy of itself. At each lag $l$, it multiplies the signal sample-by-sample with the version shifted by $l$ positions, and sums. If the signal and its shifted copy are similar (aligned peaks, aligned troughs), the sum is large. If they are dissimilar (peaks meeting troughs), the sum is small or negative.

Two properties follow directly from the definition:

**Property 1: $r[0]$ is the total signal energy.**

$$
r[0] = \sum_{n=0}^{N-1} |x[n]|^2 \tag{A.45}
$$

This is the same quantity that appears on the left side of Parseval's theorem (Equation (A.12)). The autocorrelation at lag zero is the energy — the maximum value of $r[l]$, since no shift can produce greater alignment than zero shift.

**Property 2: Symmetry.**

$$
r[-l] = r^*[l] \tag{A.46}
$$

For real signals, this simplifies to $r[-l] = r[l]$: the autocorrelation is an even function. Negative lags carry no new information beyond positive lags.

**Why "periodicity detector"?** If $x[n]$ contains a periodic component with period $P$ samples, then $x[n]$ and $x[n - P]$ are nearly identical — they are the same oscillation, shifted by one full cycle. The product $x[n] \cdot x[n - P]$ is large and positive at every $n$, so $r[P]$ is large. Likewise $r[2P]$, $r[3P]$, etc. The autocorrelation of a periodic signal is itself periodic at the same period, producing peaks at lags that are multiples of $P$.

This is what makes autocorrelation useful in the presence of noise. Noise has no periodic structure, so its autocorrelation is large only at $l = 0$ (its energy) and fluctuates near zero for all other lags. A tone buried in noise produces an autocorrelation that is noisy near $l = 0$ but shows clear periodic peaks at lags $l = P, 2P, 3P, \ldots$ — the periodicity survives even when the tone is invisible in the time-domain waveform.

---

### A.6.2 The Wiener-Khinchin Theorem

The **Wiener-Khinchin theorem** states that the autocorrelation and the power spectrum are a Fourier pair:

$$
|X[k]|^2 = \sum_{l=-(N-1)}^{N-1} r[l] \, e^{-j 2\pi k l / N} \tag{A.47}
$$

The DFT of the autocorrelation $r[l]$ is the power spectrum $|X[k]|^2$. Equivalently, the autocorrelation is the inverse DFT of the power spectrum.

This is not an approximation. It follows from substituting the definition of $r[l]$ (Equation (A.44)) and rearranging. The proof in the discrete case is short:

Starting from $r[l] = \sum_n x[n] \, x^*[n-l]$ and taking the DFT over $l$, the double sum separates into two independent sums — one giving $X[k]$ and the other giving $X^*[k]$ — whose product is $|X[k]|^2$. The details are algebraic bookkeeping; the result is exact.

The theorem connects two views of the same information:

- The **power spectrum** $|X[k]|^2$ tells you how much energy is at each frequency.
- The **autocorrelation** $r[l]$ tells you how self-similar the signal is at each lag.

They are not two different measurements. They are the same measurement expressed in two domains: one in frequency, one in lag. Knowing either one fully determines the other.

This connection is the reason Section A.6 exists. The WVD (Section A.7) will be defined as the Fourier transform of an **instantaneous** autocorrelation — an autocorrelation computed at each time position. The Wiener-Khinchin theorem is the static version of this idea: a global autocorrelation (no time index) whose Fourier transform is the global power spectrum (no time index). The WVD generalizes both sides to include time.

---

### A.6.3 Why Autocorrelation Discards Phase

The power spectrum $|X[k]|^2$ contains no phase information. Writing $X[k] = |X[k]| \, e^{j\phi[k]}$, the power spectrum is $|X[k]|^2$ — the phase $\phi[k]$ is gone. This means:

- A cosine starting at $t = 0$ and a cosine starting at $t = 0.5$ s have **identical** power spectra (same frequency, same amplitude, different phase).
- A signal and any time-shifted version of itself have **identical** autocorrelations and power spectra.

This is not a defect — it is inherent in what autocorrelation measures. Autocorrelation asks: "does the signal repeat itself after $l$ samples?" The answer depends on the signal's periodicity and amplitude, not on when it started. Phase encodes the "when"; autocorrelation is blind to it.

**The consequence for time-frequency analysis.** The power spectrum (and therefore the autocorrelation) tells you *which* frequencies are present and *how strong* they are, but not *when* they occur. This is why the global DFT power spectrum cannot distinguish a 10 Hz tone that plays for the entire recording from a 10 Hz burst that lasts one second — both deposit the same energy into the same bin, just spread differently in time.

The STFT (Section A.5) solved this by windowing: compute the power spectrum locally in time, so "when" is recovered at the cost of frequency resolution. But we now see the problem from the autocorrelation side: the global autocorrelation $r[l]$ averages over all time, so time information is lost.

The WVD's approach (Section A.7) will be different. Instead of windowing the signal and then computing the autocorrelation, it computes an **autocorrelation at each time instant** — the instantaneous autocorrelation $r_n[l] = x[n + l/2] \, x^*[n - l/2]$ — and then takes the Fourier transform over $l$. This preserves the time index $n$ through the entire computation. The Wiener-Khinchin theorem (Equation (A.47)) becomes a time-indexed family of Fourier pairs: one instantaneous autocorrelation and one instantaneous power spectrum at each $n$. That is the WVD.

---

### A.6.4 Autocorrelation Signatures of the Signal Archetypes

Before moving to the WVD, it is useful to see what autocorrelation does to each signal archetype. These signatures carry directly into the WVD, where they appear along the lag axis at each time position.

**Table A.9 — Autocorrelation signatures of the signal archetypes**

| Archetype | Signal form | Autocorrelation $r[l]$ | Key feature |
| --- | --- | --- | --- |
| Single tone ($f_0$) | $A\cos(2\pi f_0 n / f_s)$ | Cosine at the same frequency: $\propto \cos(2\pi f_0 l / f_s)$ | Periodic, undamped — frequency preserved, phase lost |
| Mixed tones ($f_1, f_2$) | $A_1\cos(\ldots) + A_2\cos(\ldots)$ | Sum of cosines at $f_1$ and $f_2$ | Each tone contributes independently; no cross-terms |
| Impulse | $\delta[n - n_0]$ | $\delta[l]$ (spike at lag zero only) | No self-similarity at any nonzero lag |
| White noise | $\eta[n] \sim \mathcal{N}(0, \sigma^2)$ | $\approx N\sigma^2 \cdot \delta[l]$ | Energy spike at $l = 0$; near-zero for $l \neq 0$ |
| Tone in noise | $A\cos(\ldots) + \eta[n]$ | Cosine at $f_0$ riding on a noise spike at $l = 0$ | Periodicity emerges at lags $l \gg 0$ where noise vanishes |

Table A.9 shows why autocorrelation is a periodicity detector (Section A.6.1). The tone and mixed-tone archetypes produce periodic autocorrelations — the frequency information survives, but the phase (start time) does not. Noise contributes only at $l = 0$. A tone buried in noise is invisible in the time domain but visible in the autocorrelation at lags beyond the noise spike.

Two observations from Table A.9 set up the WVD:

1. **No cross-terms for mixed tones.** The autocorrelation of $x_1[n] + x_2[n]$ contains $r_{x_1}[l] + r_{x_2}[l]$ plus cross-terms $r_{x_1 x_2}[l]$. For uncorrelated tones, the cross-terms vanish. This will **not** be the case for the WVD, where the quadratic structure produces cross-terms that do not vanish — the central complication of Sections A.7–A.8.

2. **No time index.** Every entry in Table A.9 is a single function of lag $l$, averaged over all time. A chirp (frequency sweeping over time) would show a blurred autocorrelation — an average of cosines at many frequencies — losing the sweep information. The WVD fixes this by computing the autocorrelation at each time instant separately.

---

*Next: A.7 — The Wigner-Ville Distribution. The global autocorrelation becomes instantaneous; the Wiener-Khinchin Fourier transform becomes time-indexed. The result is the sharpest possible time-frequency representation of a single-component signal — and, for multi-component signals, the cross-term problem that drives the rest of the report.*

---
---

# Appendix A — Signal Taxonomy

*This appendix defines the six signal archetypes that serve as the reference grid for the entire report. Every lab in Volume B constructs model signals from these archetypes; every analysis in Volume C maps real EEG features back to them. Each archetype is presented with its discrete-time formula, key properties, common real-world sources, the DSP techniques suited to it, and its role in the signal-processing narrative of this report. The summary table at the end collects their behaviour under each transform.*

---

## AA.1 Single Tone

**Mathematical formula.**

$$
x[n] = A \cos\!\left(\frac{2\pi f_0}{f_s} n + \phi\right), \qquad n = 0, 1, \ldots, N-1 \tag{AA.1}
$$

where $A$ is the amplitude, $f_0$ is the frequency in Hz, $f_s$ is the sampling rate, and $\phi$ is the initial phase.

**Properties.**

- **Stationary.** The frequency content does not change over time. Every segment of the signal has the same spectrum.
- **Narrowband.** All energy is concentrated at a single frequency $f_0$ (and its mirror at $-f_0$ for a real signal).
- **Periodic.** The signal repeats with period $P = f_s / f_0$ samples. If $P$ is not an integer, the signal is not exactly periodic within any finite window — the source of off-bin leakage (Section A.3.2).
- **Deterministic.** Fully specified by three parameters ($A$, $f_0$, $\phi$); no randomness.

**DFT behaviour.** The DFT of Equation (AA.1) places energy in the bin(s) nearest to $f_0$. If $f_0$ falls exactly on a bin ($f_0 = k \cdot f_s / N$ for integer $k$), the result is a single spike at bin $k$ with magnitude $AN/2$. If $f_0$ is off-bin, the energy leaks across neighbouring bins through the window's side lobes (Section A.3.2).

**Common sources.**

- EEG: a dominant alpha rhythm (≈10 Hz) in a resting, eyes-closed recording can approximate a single tone over short segments.
- Power-line interference: 50 Hz (or 60 Hz) contamination appears as a persistent single tone in many biomedical recordings.
- Calibration signals: a known-frequency test tone used to verify equipment.

**Techniques suited to this archetype.**

- The **DFT** (Section A.2) is sufficient. A single tone is the archetype the DFT is designed for — it maps directly to one bin (or a small cluster under leakage).
- **Windowing** (Section A.3) matters only for the off-bin case, where it controls side-lobe leakage.
- The **STFT** adds no useful information for a stationary tone — every column looks the same.

**Role in signal processing.** The single tone is the base case. Every transform in this report is easiest to understand when applied to a single tone first. It is the signal for which the DFT is exact, the STFT is redundant, the WVD is clean (except for the DC ghost from real-valued signals), and the SPWVD is trivial. It is also the building block: every other archetype is either a sum of tones, a tone whose parameters change, or the absence of tonal structure.

---

## AA.2 Mixed Tones

**Mathematical formula.**

$$
x[n] = \sum_{i=1}^{K} A_i \cos\!\left(\frac{2\pi f_i}{f_s} n + \phi_i\right) \tag{AA.2}
$$

where $K$ is the number of components, each with its own amplitude $A_i$, frequency $f_i$, and phase $\phi_i$.

**Properties.**

- **Stationary.** Like the single tone, the frequency content is constant over time.
- **Multi-narrowband.** Energy is concentrated at $K$ discrete frequencies. The spectrum has $K$ spikes (assuming they are resolvable).
- **Resolvability condition.** Two tones at $f_1$ and $f_2$ are distinguishable in the DFT only if $|f_1 - f_2| > \Delta f_{\min} = \beta \cdot f_s / N$ (Equation (A.21)). If they are closer, the window's main lobes overlap and the two tones merge into one.
- **Superposition.** The DFT is linear: $\text{DFT}(x_1 + x_2) = \text{DFT}(x_1) + \text{DFT}(x_2)$. Each tone contributes independently to the spectrum.

**DFT behaviour.** Each tone produces a spike at its corresponding bin. The spikes do not interact in the DFT (linearity). The challenge is purely resolution: are the tones far enough apart for the chosen window to separate them?

**Common sources.**

- EEG: simultaneous rhythms in different bands — alpha (10 Hz) and beta (20 Hz) coexisting in a recording. At the scalp, multiple neural oscillators contribute simultaneously.
- Music: a chord is a sum of tones.
- DTMF (telephone signalling): each key is encoded as a pair of tones.

**Techniques suited to this archetype.**

- The **DFT** resolves them if the separation exceeds $\Delta f_{\min}$. Window choice (Table A.1) determines the resolution limit.
- **Welch's method** (Section A.4.4) gives a statistically reliable estimate of the power at each tone.
- The **STFT** is again redundant for stationary mixed tones — every column is the same.

**Role in signal processing.** Mixed tones test the resolution limit. They are the archetype where the window choice in Section A.3 becomes consequential: rectangular resolves the closest pairs, but leaks the most; Blackman leaks the least, but merges pairs that rectangular could separate. In the WVD, mixed tones introduce **cross-terms** — spurious energy at the midpoint frequency $(f_1 + f_2)/2$ — which do not appear in the DFT. This is the first signal that breaks the WVD and motivates the SPWVD.

---

## AA.3 Chirp (Fluctuating Tone)

**Mathematical formula (linear chirp).**

$$
x[n] = A \cos\!\left(\frac{2\pi}{f_s}\left(f_0 n + \frac{\mu}{2} \frac{n^2}{f_s}\right) + \phi\right) \tag{AA.3}
$$

where $f_0$ is the starting frequency, $\mu$ is the **chirp rate** in Hz/s (the rate at which the instantaneous frequency changes), and the instantaneous frequency at sample $n$ is:

$$
f_{\text{inst}}[n] = f_0 + \mu \cdot \frac{n}{f_s} \quad \text{(Hz)} \tag{AA.4}
$$

**Properties.**

- **Non-stationary.** The frequency changes over time — the defining property. The signal is not the same in every segment.
- **Narrowband at each instant.** At any given time, the chirp has a well-defined instantaneous frequency. But over the full duration, the frequency sweeps a range, so the global spectrum is broadband.
- **Continuous phase.** The phase function $2\pi(f_0 n + \mu n^2 / (2 f_s)) / f_s$ is smooth — there are no discontinuities. This makes the chirp well-behaved (per the CLAUDE.md signal constraints).

**DFT behaviour.** The global DFT spreads energy across all bins between $f_0$ and $f_0 + \mu T$, producing a smeared band rather than a spike. The DFT cannot tell you *when* each frequency was present — it averages over the entire signal.

**STFT behaviour.** The chirp appears as a **diagonal line** on the spectrogram, with slope equal to the chirp rate $\mu$. However, the line is blurred: the uncertainty principle (Equation (A.40)) smears it along both axes. A short window gives a sharp time position but a thick frequency spread; a long window gives a thin frequency trace but poor time localization.

**WVD behaviour.** The WVD of a single chirp produces a **razor-sharp diagonal** — a delta function along the instantaneous frequency trajectory $f_{\text{inst}}[n]$. No blur, no smearing. This is the WVD's signature result: for a single-component signal, it achieves perfect time-frequency localization.

**Common sources.**

- EEG: frequency modulation within rhythmic activity — an alpha rhythm whose frequency drifts from 9 Hz to 11 Hz over several seconds.
- Neural oscillatory bursts that spin up (increasing frequency) or spin down (decreasing frequency).
- Radar and sonar: chirp pulses are standard for range-Doppler estimation.

**Techniques suited to this archetype.**

- The **STFT** shows the sweep direction and approximate trajectory, limited by the uncertainty principle.
- The **WVD** gives the exact trajectory for a single chirp — the clearest demonstration of the WVD's advantage over the STFT.
- The **DFT** is not useful — it collapses the time axis and smears the chirp into a featureless band.

**Role in signal processing.** The chirp is the archetype that justifies the WVD. It is the signal for which the STFT's blur is most visibly a limitation, and the WVD's sharpness is most visibly a gain. Every comparison between the STFT and WVD in this report starts with the chirp. It is also the archetype closest to real EEG non-stationarity: neural rhythms are rarely fixed-frequency; they modulate.

---

## AA.4 Multi-Chirp (Mixed Fluctuating Tones)

**Mathematical formula.**

$$
x[n] = \sum_{i=1}^{K} A_i \cos\!\left(\frac{2\pi}{f_s}\left(f_{0,i} \, n + \frac{\mu_i}{2} \frac{n^2}{f_s}\right) + \phi_i\right) \tag{AA.5}
$$

Each component has its own starting frequency $f_{0,i}$, chirp rate $\mu_i$, amplitude $A_i$, and phase $\phi_i$.

**Properties.**

- **Non-stationary and multi-component.** Multiple frequency trajectories evolving simultaneously.
- **The hardest case for time-frequency analysis.** No single linear method (DFT, STFT) can resolve all components sharply, and the WVD's quadratic nature generates cross-terms between every pair of components.

**DFT behaviour.** Overlapping smeared bands — worse than a single chirp, since the swept ranges may overlap in frequency.

**STFT behaviour.** Multiple blurred diagonals. If the chirps cross in time-frequency, the crossing region is particularly difficult to interpret.

**WVD behaviour.** Each component appears as a sharp diagonal (as in Section AA.3), but between every pair of components $(i, j)$, a **cross-term** appears:

- Located at the **midpoint** in both time and frequency between the two components.
- Oscillating at the **difference frequency** $|f_i(t) - f_j(t)|$.
- The cross-term carries as much energy as the auto-terms (the genuine components). It is not a small perturbation.

For $K$ components, there are $K$ auto-terms and $K(K-1)/2$ cross-terms. The cross-terms grow quadratically with the number of components and can dominate the representation.

**Common sources.**

- EEG: simultaneous rhythms in different bands, each with its own frequency modulation — e.g. a drifting alpha and a drifting beta coexisting.
- Multi-speaker audio, multi-target radar, overlapping neural oscillatory bursts.

**Techniques suited to this archetype.**

- The **STFT** shows blurred but cross-term-free diagonals (linearity prevents cross-terms).
- The **WVD** shows sharp diagonals buried in cross-term interference.
- The **SPWVD** (Section A.8) suppresses cross-terms with two independent smoothing knobs while preserving most of the WVD's sharpness — this is the archetype the SPWVD is designed for.

**Role in signal processing.** The multi-chirp is the archetype that breaks the raw WVD. It is the reason the WVD alone is not a practical tool, and the reason the PWVD and SPWVD exist. Every real multi-component signal — including EEG — behaves like a multi-chirp at some level, so understanding cross-terms is not optional.

---

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
- **Broadband.** By the uncertainty principle (Equation (A.40)), a short signal must have a wide bandwidth. A pure impulse has a flat spectrum — equal energy at all frequencies.
- **Non-stationary.** The signal is not the same in every segment — it exists only briefly.
- **The dual of a tone.** A tone is perfectly localized in frequency and completely delocalized in time. A transient is the opposite: localized in time, spread in frequency. They are the two endpoints of the uncertainty tradeoff.

**DFT behaviour.** A broadband smear — the DFT spreads the transient's energy across many bins. For a pure impulse, $|X[k]| = A$ for all $k$ (flat magnitude spectrum). The DFT cannot tell you *when* the transient occurred.

**STFT behaviour.** A **vertical stripe** in the spectrogram — energy at all frequencies, concentrated in the time column(s) that contain the transient. A short STFT window localizes the stripe tightly; a long window smears it across multiple columns.

**WVD behaviour.** A sharp vertical stripe (auto-term), plus **broadband ghosts** if any other component is present. For a single impulse in isolation, the WVD is clean. For an impulse plus a tone, a cross-term oscillates at the tone's frequency, centered between the impulse time and the tone's presence.

**Common sources.**

- EEG: blink artifacts (sharp transients from eye movement), epileptiform spikes, K-complexes in sleep EEG, electrode pops.
- Biomedical: QRS complex in ECG (a sharp, brief waveform).
- General: clicks, impacts, switching transients.

**Techniques suited to this archetype.**

- The **STFT** with a short window gives the best time localization — but at the cost of frequency resolution (which is irrelevant for a broadband event).
- **Statistical detection** (Section A.4.2): a transient creates a sudden energy spike that exceeds the noise floor across many bins simultaneously. Detecting it is a time-domain energy test, not a frequency-domain peak test.
- The **autocorrelation** of a transient is a spike at lag zero (Equation (A.45)) — no periodicity.

**Role in signal processing.** The transient is the opposite extreme from the tone. It tests time localization rather than frequency resolution. In EEG, transient detection (Part C.4) is primarily about identifying artifacts (blinks, electrode pops) and clinically relevant events (spikes). The transient archetype is also the one where the STFT's short-window mode works best, and where the WVD's cross-terms are most problematic (broadband ghosts from interaction with any coexisting tone or rhythm).

---

## AA.6 Noise

**Mathematical formula (white Gaussian noise).**

$$
x[n] = \eta[n], \qquad \eta[n] \sim \mathcal{N}(0, \sigma^2), \quad \text{i.i.d.} \tag{AA.8}
$$

Each sample is drawn independently from a Gaussian distribution with mean zero and variance $\sigma^2$.

**Properties.**

- **Stationary** (in the statistical sense). The distribution of any sample is the same at all times. But any individual realization looks irregular and non-repeating.
- **Broadband.** The expected power spectrum is flat: $\mathbb{E}[|X[k]|^2] = N\sigma^2$ for all $k$ (Equation (A.26)). Energy is spread equally across all frequencies.
- **No deterministic structure.** No periodicity, no frequency trajectory, no time-localized event. Noise is defined by its statistical properties, not by a formula for $x[n]$.
- **Phase is uniformly random.** Each bin's phase $\angle X[k]$ is independent and uniform on $(-\pi, \pi]$ (Section A.4.1).

**DFT behaviour.** A ragged, flat spectrum: each bin's magnitude fluctuates around the same mean $N\sigma^2$, with the exponential distribution described in Equation (A.25). No bin stands out systematically — but any single realization has random peaks and valleys that can be mistaken for tonal features if statistics are not applied.

**STFT behaviour.** Uniform speckle across the entire time-frequency plane. No horizontal lines, no vertical stripes, no diagonals — just random variation at each $(m, k)$ cell. This is the visual baseline: any structure that appears on top of the speckle is a candidate for a genuine signal feature.

**WVD behaviour.** Speckle with self-interference. The WVD of noise is noisier than the STFT of noise because the quadratic structure creates cross-terms between different noise components. The WVD of noise has higher variance than the spectrogram of noise — the sharpness that helps for tones and chirps hurts for noise.

**SPWVD behaviour.** The smoothing windows reduce the variance of the noise floor but do not remove it. Noise is the one archetype that no amount of time-frequency smoothing can clean up — it has no structure to preserve.

**Common sources.**

- Thermal noise in amplifiers and electrodes (Johnson-Nyquist noise).
- Biological background activity in EEG: the aggregate electrical activity of millions of neurons that do not belong to the rhythm under study.
- Quantization noise from analog-to-digital conversion.
- EMG contamination: muscle activity produces broadband noise-like signals at higher frequencies.

**Techniques suited to this archetype.**

- **Statistical characterization** (Section A.4): the noise floor, the exponential distribution of bin power, and the detection threshold from Equation (A.28).
- **Averaging** (Welch's method, Section A.4.4): reduces variance of the spectral estimate, revealing the true noise level.
- **Autocorrelation** (Section A.6): white noise has $r[l] \approx 0$ for $l \neq 0$, which distinguishes it from any periodic or quasi-periodic signal.

**Role in signal processing.** Noise is not a signal to analyze — it is the **background against which signals are detected**. Every detection, estimation, and resolution result in this report (Equation (A.28), Equation (A.30), Table A.3) is stated relative to the noise. Without a noise model, "this bin has power $P$" is meaningless; with the exponential model from Section A.4.1, it becomes "this bin's power exceeds the noise floor by a factor of $\gamma$ with false-alarm probability $e^{-\gamma}$." Noise defines the floor; signals are what rises above it.

---

## Summary Table

The following table collects the behaviour of all six archetypes under each transform. Every cell has been justified in the subsections above and in the corresponding theory sections of Volume A. This table is referenced throughout Volume B (labs) and Volume C (EEG application).

**Table AA.1 — Signal archetypes × transform behaviour**

| Archetype | Time-domain form | DFT | STFT | WVD | SPWVD |
| --- | --- | --- | --- | --- | --- |
| Single tone (AA.1) | $A\cos(2\pi f_0 n / f_s + \phi)$ | Single spike at $f_0$ | Horizontal line | Sharp horizontal line (+ DC ghost if real) | Clean line |
| Mixed tones (AA.2) | $\sum A_i \cos(\ldots)$ | $K$ spikes | Parallel horizontal lines | Lines + midpoint cross-term ghosts | Lines; ghosts suppressed |
| Chirp (AA.3) | Swept-frequency sinusoid | Smeared band | Blurred diagonal | Razor-sharp diagonal | Sharp diagonal |
| Multi-chirp (AA.4) | Sum of chirps | Overlapping smears | Crossing blurred diagonals | Sharp diagonals + moving cross-term ghosts | Diagonals; ghosts suppressed |
| Transient / pulse (AA.5) | Short burst or impulse | Flat (broadband) | Vertical stripe | Sharp vertical stripe (+ broadband ghosts) | Vertical stripe; ghosts suppressed |
| Noise (AA.6) | Random process | Flat, random phase | Uniform speckle | Speckle (+ self-interference) | Speckle; not removed by smoothing |

**Reading the table.** Each row is an archetype; each column is a transform. Moving left to right across the columns traces the progression of this report: from the global DFT (no time information) to the STFT (time-frequency, limited by uncertainty) to the WVD (sharp but cross-term-corrupted) to the SPWVD (sharp with cross-terms managed). The diagonal of the table — the cells where each archetype meets its best-suited transform — is where the analysis produces the clearest result. Off-diagonal cells show where a transform struggles with a particular archetype.

**The two failure modes.** The table reveals two systematic failures:

1. **The STFT always blurs.** Every STFT cell uses words like "blurred" or gives less detail than the WVD cell in the same row. This is the uncertainty principle (Equation (A.40)) manifesting across all archetypes.

2. **The WVD always adds ghosts for multi-component signals.** Every WVD cell for a multi-component archetype (mixed tones, multi-chirp, transient coexisting with other components) mentions cross-terms. This is the quadratic nature (Section A.7.3) manifesting across all archetypes.

The SPWVD column shows the compromise: it inherits most of the WVD's sharpness while suppressing most of the cross-terms, at the cost of some smoothing. The one archetype it cannot help is noise — there is no structure to preserve.
