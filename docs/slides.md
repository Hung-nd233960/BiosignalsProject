---
marp: true
theme: default
paginate: true
math: mathjax
style: |
  section {
    font-family: 'Roboto', sans-serif;
    font-size: 24px;
  }
  code {
    font-family: 'Consolas', monospace;
  }
  h1 { font-size: 36px; }
  h2 { font-size: 30px; }
  table { font-size: 20px; }
---

# Digital Signal Processing: Time-Frequency Analysis from the DFT to the SPWVD

**Application to EEG**

Nguyen Duc Hung - 20233960

Digital Biosignal Processing - Final Report

---

## Roadmap

**Volume A** - Theory: sampling, DFT, windowing, STFT, autocorrelation, WVD, SPWVD
**Volume B** - Labs: 8 experiments confirming the theory + Appendix B (M vs M-1 convergence)
**Volume C** - Application to real EEG data

The spine: **sampling → bins → leakage → statistics → STFT/uncertainty → autocorrelation → WVD → SPWVD**

Each concept builds on the previous. The endpoint: a time-frequency tool (SPWVD) that works on real EEG.

---

<!-- _class: lead -->

# A.1 - Signal Theory

*Sampling, discrete frequency, Nyquist, and the two independent knobs.*

---

## Sampling: Two Independent Knobs

**Sampling rate** $f_s$ controls the highest frequency we can represent (Nyquist: $f_{\max} = f_s/2$).

**Sampling count** $N$ controls how long we observe (duration: $T = N/f_s$).

These are independent choices. Increasing $f_s$ does not give us more frequency resolution. Increasing $N$ does not give us a higher Nyquist frequency.

| Knob | Controls | Does NOT control |
|---|---|---|
| $f_s$ (rate) | Nyquist frequency $f_s/2$ | Frequency resolution |
| $N$ (count) | Duration $T$, resolution $\Delta f$ | Nyquist frequency |

---

<!-- _class: lead -->

# A.2 - The DTFT and the DFT

*From the continuous DTFT to the computable DFT. What a bin is.*

---

## From the DTFT to the DFT

The DTFT of a finite signal $x[n]$ is a **continuous** function of frequency:

$$X(\omega) = \sum_{n=0}^{N-1} x[n] \, e^{-j\omega n}$$

Problem: we cannot store a continuous function. Solution: **sample it** at $N$ equally spaced points:

$$X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j 2\pi k n / N}, \qquad k = 0, 1, \ldots, N-1$$

The DFT is the DTFT evaluated at $N$ specific frequencies. Each value $X[k]$ is called a **bin**.

---

## What a Bin Is

A bin is one sample of the continuous DTFT, taken at frequency $f_k = k \cdot f_s / N$.

**Bin spacing:**

$$\Delta f = \frac{f_s}{N} = \frac{1}{T} \quad \text{(Hz)}$$

The bin spacing is the reciprocal of the observation time. This is not a design choice - it is a consequence of sampling the DTFT.

- 1-second signal → $\Delta f = 1$ Hz
- 20-second signal → $\Delta f = 0.05$ Hz
- 1200-second signal → $\Delta f = 0.000833$ Hz

---

## Resolution ≠ Bin Count

**Bin count** = number of DFT outputs ($N$). Can be increased by zero-padding.

**Frequency resolution** = ability to distinguish two nearby tones. Determined by **signal duration** $T$, not DFT length.

$$\Delta f_{\min} = \frac{1}{T} = \frac{f_s}{N}$$

**Zero-padding** adds bins (finer sampling of the same curve) but does NOT add resolution (the curve's shape is unchanged). More points on the same curve ≠ sharper features.

---

<!-- _class: lead -->

# A.3 - Leakage and Windowing

*The hidden rectangular window, spectral leakage, and the cosine-sum family.*

---

## The Hidden Rectangular Window

The DFT sums from $n = 0$ to $N - 1$. Everything outside is silently set to zero.

This is a **rectangular window**: $w[n] = 1$ for $0 \leq n \leq N-1$, zero otherwise.

If the signal is not perfectly periodic within the $N$-sample frame, the abrupt truncation creates a discontinuity → **spectral leakage**.

The DFT of the rectangular window is the **Dirichlet kernel**:

$$D(\omega) = \frac{1}{M}\left|\frac{\sin(\omega M/2)}{\sin(\omega/2)}\right|$$

Its side lobes are the source of leakage.

---

## The Cosine-Sum Window Family

Trade main-lobe width for side-lobe suppression:

| Window | Main-lobe (bins) | Peak side-lobe | Rolloff | $\beta$ |
|---|---|---|---|---|
| Rectangular | 2 | -13 dB | 6 dB/oct | 1 |
| Hann | 4 | -31.5 dB | 18 dB/oct | 2 |
| Hamming | 4 | -42.7 dB | 6 dB/oct | 2 |
| Blackman | 6 | -58 dB | 18 dB/oct | 3 |

**Resolution limit** with window:

$$\Delta f_{\min} \approx \beta \cdot \frac{f_s}{N}$$

Wider main lobe = worse resolution, but much cleaner spectrum.

---

<!-- _class: lead -->

# A.4 - Statistics and the DFT

*DFT bins under noise are random variables. Detection and averaging.*

---

## Each Bin is a Random Variable

Under white Gaussian noise, the **power** at each bin follows an **exponential distribution**:

$$P[k] \sim \text{Exponential}\left(\lambda = \frac{1}{N\sigma^2}\right)$$

Key consequence: the coefficient of variation = 1. A single bin's power estimate has **100% relative uncertainty**, regardless of $N$.

**Detection:** a tone is detected if its bin power exceeds $\gamma \times$ noise floor, where:

$$\Pr(P[k] > \gamma \cdot P_{\text{floor}}) = e^{-\gamma}$$

This is derived from the spectral model, not from a generic "$2\sigma$" rule.

---

## The Periodogram Variance Problem + Welch

A single DFT of a long signal gives **fine frequency spacing** but **ragged, unreliable** power estimates. More data does not help - the variance does not decrease with $N$.

**Welch's method:** divide into $L$ overlapping segments, average their periodograms.

- Variance reduced by $1/L$
- Resolution degrades to $\beta \cdot f_s / M$ (per segment)
- The tradeoff: resolution vs. statistical reliability

| Segments | $\Delta f$ | Relative variance |
|---|---|---|
| 1 (full) | 0.0008 Hz | 1.000 |
| 119 (20 s) | 0.050 Hz | 0.008 |
| 479 (5 s) | 0.200 Hz | 0.002 |

---

<!-- _class: lead -->

# A.5 - The STFT and Spectrograms

*The windowed DFT slid in time. The uncertainty principle.*

---

## The Spectrogram

Slide a windowed DFT across time, keep each spectrum indexed by position:

$$X[m, k] = \sum_{n=0}^{M-1} x[n + mH] \, w[n] \, e^{-j 2\pi k n / M}$$

The spectrogram $S[m, k] = |X[m, k]|^2$ is a 2D image:
- Horizontal axis: **time** (spacing $H/f_s$)
- Vertical axis: **frequency** (spacing $f_s/M$)
- Color: **power**

This is the first tool that answers "what frequency is present at what time."

---

## The Uncertainty Principle

**The central tradeoff of the entire report.**

$$\Delta t \cdot \Delta f = \beta \qquad \text{(constant for a given window)}$$

| Window $M$ | $\Delta t$ (s) | $\Delta f$ (Hz) | Can resolve |
|---|---|---|---|
| 125 (0.5 s) | 0.5 | 4.0 | Sub-second transients |
| 250 (1.0 s) | 1.0 | 2.0 | Alpha vs. theta bands |
| 1250 (5.0 s) | 5.0 | 0.4 | Fine spectral detail |

Short window → good time, poor frequency.
Long window → good frequency, poor time.
**No window captures both.** This motivates the WVD.

---

## Overlap and COLA

Windows taper toward zero at the edges → edge samples are suppressed.

**Overlap** fixes this: each sample appears in multiple segments at different positions.

**COLA condition** (Hann at 50% overlap):

$$w[p] + w[p + M/2] = 1 \quad \text{for all } p$$

Every sample receives equal total weight. No information lost.

Overlap beyond COLA (e.g. 75%) adds time columns but NOT time resolution - same distinction as zero-padding on the frequency axis.

---

<!-- _class: lead -->

# A.6 - Autocorrelation

*Periodicity detection, Wiener-Khinchin, and why phase is lost.*

---

## Autocorrelation and Wiener-Khinchin

**Autocorrelation:** compare a signal with shifted copies of itself.

$$r[l] = \sum_{n} x[n] \, x^*[n - l]$$

- $r[0]$ = total signal energy (Parseval)
- Periodic signal → periodic $r[l]$ → periodicity detector

**Wiener-Khinchin theorem:** the DFT of the autocorrelation is the power spectrum.

$$|X[k]|^2 = \text{DFT}\{r[l]\}$$

Autocorrelation and power spectrum are the same information in two domains.

---

## Phase is Lost

Autocorrelation tells you **which** frequencies are present and **how strong**, but NOT **when** they occur.

Two tones with $\phi = 0$ and $\phi = \pi$: different time-domain signals, **identical** autocorrelations.

This is the limitation the WVD addresses:
- Global autocorrelation → **instantaneous** autocorrelation at each time $n$
- Wiener-Khinchin → a **time-indexed** family of Fourier pairs
- The result: the WVD = the sharpest time-frequency representation

---

<!-- _class: lead -->

# Signal Taxonomy

*The reference grid for the entire report.*

---

## Six Archetypes x Four Transforms

| Archetype | DFT | STFT | WVD | SPWVD |
|---|---|---|---|---|
| Single tone | spike | horiz. line | sharp line | clean |
| Mixed tones | spikes | parallel lines | lines + ghosts | ghosts suppressed |
| Chirp | smear | blurred diagonal | **sharp diagonal** | sharp |
| Multi-chirp | smears | blurred crossings | diagonals + ghosts | ghosts suppressed |
| Transient | flat | vertical stripe | stripe + ghosts | ghosts suppressed |
| Noise | flat | speckle | speckle | speckle (can't fix) |

Every lab in Volume B and every analysis in Volume C maps back to this table.

---

<!-- _class: lead -->

# What's Next

---

## Volume B: Confirm Everything

8 labs, each testing a specific claim from Volume A:

1. **DFT** - bins, leakage, zero-padding
2. **Statistics** - exponential distribution, Welch
3. **Windowing** - Dirichlet kernel derivation
4. **STFT** - Heisenberg tradeoff
5. **Resolution** - two-tone test on spectrograms
6. **Autocorrelation** - periodicity, Wiener-Khinchin
7. **WVD** - sharp chirp, cross-terms
8. **SPWVD** - two-knob ghost suppression

Then **Volume C**: apply it all to real EEG.
