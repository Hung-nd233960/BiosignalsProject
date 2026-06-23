
# From the DFT to the SPWVD
## Time-Frequency Analysis Applied to Neonatal EEG

**Authors:** Nguyen Duc Hung - 20233960, Bui Phuong Duy - 23233957, Tran Viet Bach - 23233954

Digital Biosignal Processing - Final Report

---

## Why Time-Frequency Analysis?

- Electroencephalography (EEG) signals are non-stationary: their frequency content changes over time.
- The classical DFT provides a global spectrum but destroys all temporal information.
- Multi-component signals (multiple rhythms, bursts, artifacts) require tools that resolve both time and frequency.
- Neonatal EEG exhibits discontinuous burst patterns that a single spectrum cannot characterize.
- The goal is a progression of tools, each fixing a specific limitation of the previous one.

---

## Three-Volume Structure

- **Volume A (Theory):** Derives the DFT, windowing, STFT, autocorrelation, WVD, and SPWVD from first principles in the discrete domain.
- **Volume B (Laboratories):** Eight labs validate each theory section on model signals with known ground truth.
- **Volume C (Application):** The complete toolchain is applied to a real neonatal EEG recording (19 minutes, 24 channels).
- All mathematical formulas are discrete-time and discrete-frequency throughout.
- The report is self-contained: every tool applied in Volume C is derived in Volume A and tested in Volume B.

---

## The Tool Progression

- **DFT:** Computes the global spectrum but provides no time information whatsoever.
- **STFT:** Adds a time axis via a sliding window but is bound by the Heisenberg uncertainty limit.
- **WVD:** Achieves sub-Heisenberg resolution for single components but generates catastrophic cross-terms for multiple components.
- **SPWVD:** Introduces two independent smoothing knobs to suppress cross-terms at controllable resolution cost.

Each tool in the progression fixes the specific limitation of its predecessor.

$$\text{DFT} \xrightarrow{\text{+time}} \text{STFT} \xrightarrow{\text{+sharpness}} \text{WVD} \xrightarrow{\text{+smoothing}} \text{SPWVD}$$

---

## Signal Taxonomy

- **Deterministic stationary:** Periodic signals (sum of sines) with unchanging frequency content over time.
- **Deterministic non-stationary:** Signals whose frequency content changes (chirps, bursts, transients).
- **Random stationary:** Stochastic processes with constant statistical properties (white noise, colored noise).
- **Random non-stationary:** Signals with time-varying statistics (the most general class; real EEG belongs here).
- The DFT is sufficient only for the first class; each additional class demands a more powerful tool.

---

## Part 1: Volume A - Theory

The theoretical foundation: six sections building from the DFT to the SPWVD.

- A.1-A.2: The DFT and its frequency grid.
- A.3: Windowing and the resolution-leakage tradeoff.
- A.4: Spectral statistics and detection.
- A.5: The STFT and Heisenberg uncertainty.
- A.6: Autocorrelation and phase-blindness.
- A.7-A.8: The WVD and SPWVD.

---

## A.1-A.2: What the DFT Computes

$$X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j2\pi kn/N} \qquad \text{(A.5)}$$

- The DFT decomposes a length-$N$ signal into $N$ complex coefficients, one per frequency bin.
- Each coefficient $X[k]$ measures the correlation between the signal and a complex exponential at frequency $k \cdot f_s / N$.
- The transform is invertible: $N$ time samples map to $N$ frequency samples with no information loss.
- The magnitude $|X[k]|$ gives the amplitude; the phase $\angle X[k]$ gives the temporal alignment.

---

## A.1-A.2: Bin Spacing and the Frequency Grid

$$\Delta f = \frac{f_s}{N} = \frac{1}{T} \qquad \text{(A.6)}$$

- The DFT evaluates the spectrum at $N$ equally spaced frequencies separated by $\Delta f = f_s / N$.
- The bin spacing depends only on the sampling rate $f_s$ and the number of samples $N$.
- Equivalently, $\Delta f = 1/T$ where $T$ is the signal duration in seconds.
- Finer bin spacing requires a longer signal, not a higher sampling rate.
- The frequency grid is fixed once the signal length is chosen.

---

## A.1-A.2: On-Bin Versus Off-Bin Frequencies

- An **on-bin** frequency is an exact multiple of $\Delta f$: it aligns perfectly with one DFT bin.
- All energy concentrates in a single bin; all other bins read exactly zero (for that component).
- An **off-bin** frequency falls between two bins: its energy leaks across all $N$ bins.
- This leakage is called **spectral leakage** and is not a numerical error but a structural property of the finite-length DFT.
- The worst case occurs at the midpoint between two bins (maximally off-bin).

---

## A.1-A.2: Zero-Padding Is Interpolation, Not Resolution

- Appending zeros to the signal before the DFT increases $N$ and decreases the bin spacing $\Delta f$.
- The spectrum appears smoother because additional frequency samples are computed between the original bins.
- However, zero-padding does not add new information about the signal.
- Two tones closer than $1/T$ Hz apart remain unresolved regardless of how many zeros are appended.
- Zero-padding is spectral interpolation; genuine resolution requires a longer observation of the signal.

---

## A.3: Why Windowing Is Needed

- The DFT implicitly assumes the signal repeats with period $N$ (periodic extension).
- If the signal does not complete an integer number of cycles within $N$ samples, the periodic extension introduces discontinuities.
- Truncation to $N$ samples is equivalent to multiplication by a rectangular window.
- The rectangular window has large side lobes ($-13$ dB first side lobe) that cause severe spectral leakage.
- Applying a tapered window reduces the discontinuity at the edges and suppresses leakage.

---

## A.3: The Dirichlet Kernel

- The DFT of the rectangular window produces the Dirichlet kernel: $W_{\text{rect}}(\omega) = e^{-j\omega(M-1)/2} \cdot \frac{\sin(\omega M/2)}{\sin(\omega/2)}$.
- The kernel has a central **main lobe** of width $4\pi/M$ (in normalized frequency).
- Flanking the main lobe are **side lobes** that decay at only $-20$ dB per octave.
- All energy that leaks into the side lobes is energy misattributed to wrong frequencies.
- The side-lobe level and decay rate are the two metrics that characterize leakage severity.

---

## A.3: Cosine-Sum Windows

- All standard windows (Hann, Hamming, Blackman) are sums of shifted cosines applied to the rectangular window.
- Each cosine term adds a pair of shifted Dirichlet kernels that destructively interfere with the side lobes.
- The Hann window uses one cosine term and achieves $-32$ dB first side lobe with $-60$ dB/octave decay.
- The Hamming window optimizes the first side lobe to $-43$ dB but retains the slow $-20$ dB/octave decay.
- The Blackman window uses two cosine terms: $-58$ dB first side lobe, $-60$ dB/octave decay, but triple main-lobe width.

---

## A.3: The Pure Sine Form

- The numerator of every cosine-sum window spectrum factors into a product of pure sine functions.
- This factorization reveals that the nulls of the window spectrum are located at integer multiples of $2\pi/M$.
- The side-lobe decay rate is determined by the order of the zero at $\omega = 0$: higher-order zeros yield faster decay.
- Hann and Blackman have third-order and fifth-order zeros respectively, explaining their $-60$ dB/octave decay.
- This derivation connects the window coefficients directly to the measurable side-lobe behavior.

---

## A.3: Window Comparison Table

| Window | First Side Lobe | Decay Rate | Main Lobe Width |
|--------|----------------|------------|-----------------|
| Rectangular | $-13$ dB | $-20$ dB/oct | $2 f_s/N$ |
| Hann | $-32$ dB | $-60$ dB/oct | $4 f_s/N$ |
| Hamming | $-43$ dB | $-20$ dB/oct | $4 f_s/N$ |
| Blackman | $-58$ dB | $-60$ dB/oct | $6 f_s/N$ |

- No window is universally best: the choice depends on whether leakage suppression or frequency resolution is more critical.
- The rectangular window provides the best resolution but the worst leakage.

---

## A.4: Bin Power Distribution Under Noise

- For white Gaussian noise, each DFT bin power $P[k] = |X[k]|^2/N$ follows an exponential distribution.
- The coefficient of variation (CV = standard deviation / mean) equals exactly 1.0 for exponentially distributed data.
- This means bin power fluctuates by 100% of its mean value, regardless of signal length.
- Increasing $N$ adds more bins but does not reduce the variance of any single bin.
- The periodogram is an unbiased but inconsistent estimator of the power spectrum.

---

## A.4: Detection Threshold and False-Alarm Probability

$$P_{fa} = e^{-\gamma} \qquad \text{(A.37)}$$

- A spectral peak is deemed genuine only if it exceeds $\gamma \cdot \hat{P}_{\text{floor}}$, where $\hat{P}_{\text{floor}}$ is the estimated noise floor.
- The probability that a noise-only bin exceeds this threshold is $P_{fa} = e^{-\gamma}$.
- For $P_{fa} = 0.01$ (1% false alarm): $\gamma = \ln(100) \approx 4.6$.
- This is a principled detection rule derived from the exponential model, not an arbitrary "$k\sigma$" heuristic.
- The threshold adapts to the actual noise level in the spectrum.

---

## A.4: Welch's Method and the Resolution-Variance Tradeoff

- Welch's method divides the signal into $K$ overlapping segments, computes the periodogram of each, and averages.
- Averaging $K$ segments reduces variance by a factor of approximately $K$ (standard deviation by $\sqrt{K}$).
- The cost: each segment has length $M < N$, so frequency resolution degrades from $f_s/N$ to $\beta \cdot f_s/M$.
- This is the fundamental **resolution-variance tradeoff**: better statistics require shorter segments, which blur frequency.
- Welch's method does not bypass uncertainty; it trades frequency detail for statistical reliability.

---

## A.5: The STFT Definition

- The Short-Time Fourier Transform applies a window $w[n]$ centered at time $n_0$, then computes the DFT of the windowed segment.
- The window slides along the signal with hop size $H$, producing a spectrum at each time step.
- The result is a two-dimensional function: magnitude versus time and frequency (the spectrogram).
- The STFT adds the time axis that the DFT lacks, enabling detection of non-stationary behavior.
- The window length $M$ is the single design parameter that controls the time-frequency tradeoff.

---

## A.5: The Heisenberg Uncertainty Principle

$$\Delta t \cdot \Delta f \geq \beta \qquad \text{(A.49)}$$

- A short window ($M$ small) gives good time resolution $\Delta t$ but poor frequency resolution $\Delta f$.
- A long window ($M$ large) gives good frequency resolution but poor time resolution.
- The product $\Delta t \cdot \Delta f$ cannot be made smaller than a constant $\beta$ that depends on the window shape.
- For the Hann window, $\beta = 2$, meaning $\Delta t \cdot \Delta f \geq 2$.
- This is not a limitation of implementation but a fundamental mathematical bound.

---

## A.5: Hop Size and the COLA Condition

- The hop size $H$ determines how frequently the window advances along the signal.
- For perfect reconstruction, the windows must satisfy the Constant Overlap-Add (COLA) condition: $\sum_m w[n - mH] = \text{const}$ for all $n$.
- Typical choices: $H = M/2$ (50% overlap) for Hann, $H = M/4$ (75% overlap) for smoother time tracking.
- Smaller hop sizes produce more time frames but increase computation proportionally.
- COLA ensures that no part of the signal is over- or under-represented in the analysis.

---

## A.5: Multi-Scale Limitation

- The STFT uses a single window length for the entire signal.
- A chirp that spans both low and high frequencies cannot be optimally resolved by one window.
- Low frequencies need a long window (many cycles required); high frequencies need a short window (rapid changes).
- This one-window limitation motivates the search for representations that decouple time and frequency resolution.
- The WVD addresses this by eliminating the analysis window entirely.

---

## A.6: Autocorrelation Definition

$$R_{xx}[l] = \sum_{n} x[n] \cdot x^*[n - l] \qquad \text{(A.54)}$$

- The autocorrelation measures the self-similarity of a signal at lag $l$.
- A periodic signal produces periodic peaks in the autocorrelation at multiples of the period.
- The autocorrelation is always symmetric: $R_{xx}[l] = R_{xx}^*[-l]$.
- It reaches its maximum at lag zero: $R_{xx}[0]$ equals the total signal energy.
- Autocorrelation detects periodicity without requiring knowledge of the signal's frequency.

---

## A.6: The Wiener-Khinchin Theorem

$$S_{xx}[k] = \text{DFT}\{R_{xx}[l]\} \qquad \text{(A.56)}$$

- The DFT of the autocorrelation equals the power spectral density (PSD).
- This theorem connects the time-domain concept of self-similarity to the frequency-domain concept of power distribution.
- The PSD is always real and non-negative, because it is the magnitude-squared spectrum.
- The autocorrelation discards phase information: signals with identical spectra but different phases have identical autocorrelations.
- Phase-blindness is the price of the elegant autocorrelation-to-spectrum relationship.

---

## A.6: Cross-Correlation

- The cross-correlation $R_{xy}[l] = \sum_n x[n] \cdot y^*[n-l]$ measures shared structure between two signals.
- A peak in $R_{xy}[l]$ at lag $l_0$ indicates that $x$ and $y$ share a component with time offset $l_0$.
- Normalizing by the geometric mean of energies gives the correlation coefficient $\rho \in [-1, 1]$.
- Cross-correlation is used in Volume C to test whether auxiliary channels (ECG, EMG) contaminate the EEG.
- A low $|\rho|$ confirms independence; a high $|\rho|$ flags potential artifact contamination.

---

## A.7: The Instantaneous Autocorrelation

$$R_x[n, m] = z[n+m] \cdot z^*[n-m] \qquad \text{(A.60)}$$

- The instantaneous autocorrelation (IAF) generalizes autocorrelation to be time-dependent.
- At each time instant $n$, it measures the signal's self-similarity over lag $m$ centered at $n$.
- The signal $z[n]$ is the analytic signal (complex-valued, obtained via the Hilbert transform).
- Unlike the global autocorrelation, the IAF preserves time localization.
- The IAF is the fundamental building block of the Wigner-Ville Distribution.

---

## A.7: The WVD Definition

$$W_x[n, k] = \sum_{m} R_x[n, m] \, e^{-j2\pi km/N_f} \qquad \text{(A.61)}$$

- The WVD computes the DFT of the instantaneous autocorrelation over the lag variable $m$.
- It maps the signal onto the time-frequency plane without any analysis window.
- No window means no Heisenberg constraint: the WVD achieves perfect time-frequency concentration for a single component.
- A linear chirp produces a razor-sharp line along its instantaneous frequency trajectory.
- The WVD satisfies the marginal properties: summing over frequency gives instantaneous power, summing over time gives the power spectrum.

---

## A.7: Single-Component Sharpness

- For a single linear chirp $z[n] = e^{j(\alpha n + \beta n^2)}$, the WVD is a Dirac delta along the instantaneous frequency $f_i[n] = \alpha + 2\beta n$.
- This is the sharpest possible time-frequency representation: zero spread in both dimensions simultaneously.
- The STFT of the same chirp produces a blurred stripe whose width is determined by the Heisenberg bound.
- The WVD bypasses Heisenberg because it is quadratic (bilinear) in the signal, not linear.
- For single-component signals, the WVD is the ideal time-frequency representation.

---

## A.7: Cross-Terms

- For a signal with two components $z = z_1 + z_2$, the WVD produces: $W_z = W_{z_1} + W_{z_2} + 2\,\text{Re}(W_{z_1 z_2})$.
- The cross-term $W_{z_1 z_2}$ appears at the **midpoint** between the two auto-terms in time-frequency.
- Cross-terms oscillate at a frequency proportional to the separation between the components.
- Cross-terms are as energetic as the auto-terms: they do not decay with distance.
- For $N$ components, there are $N(N-1)/2$ cross-terms versus only $N$ auto-terms, making the WVD unusable for complex signals.

---

## A.7: The Analytic Signal

$$z[n] = x[n] + j\,\mathcal{H}\{x[n]\} \qquad \text{(A.70)}$$

- The analytic signal removes negative-frequency content via the Hilbert transform.
- A real signal $x[n]$ has a symmetric spectrum: $X[-k] = X^*[k]$. This symmetry generates a cross-term at DC (the "DC self-ghost").
- The analytic signal eliminates the negative frequencies, removing the DC ghost from the WVD.
- The Hilbert transform is computed in practice via zeroing the negative-frequency DFT bins and inverse-transforming.
- All WVD computations in this report use the analytic signal.

---

## A.8: The PWVD (Lag Window Only)

- The Pseudo Wigner-Ville Distribution (PWVD) applies a lag window $h[m]$ to the instantaneous autocorrelation before the DFT.
- The lag window limits the range of $m$ values that contribute, smoothing the frequency axis.
- Frequency-oscillating cross-terms (from components separated in time) are suppressed.
- However, time-oscillating cross-terms (from components separated in frequency) survive because no time smoothing is applied.
- The PWVD is a half-solution: it addresses one axis but leaves the other uncontrolled.

---

## A.8: PWVD Limitation

- Two tones at different frequencies but present simultaneously generate time-oscillating cross-terms.
- These ghosts oscillate along the time axis at a rate proportional to the frequency separation.
- The lag window $h[m]$ cannot suppress them because it only smooths in frequency.
- A second smoothing operation along the time axis is required.
- This motivates the full Smoothed Pseudo Wigner-Ville Distribution (SPWVD).

---

## A.8: The SPWVD Definition

$$\text{SPWVD}_x[n, k] = \sum_m h[m] \left(\sum_p g[p] \, z[n\!+\!p\!+\!m] \, z^*[n\!+\!p\!-\!m]\right) e^{-j2\pi km/N_f} \quad \text{(A.72)}$$

- The SPWVD adds a time window $g[p]$ that averages the instantaneous autocorrelation along the time axis.
- Window $h[m]$ smooths frequency; window $g[p]$ smooths time.
- Both ghost types (frequency-oscillating and time-oscillating) are now suppressed.
- The degree of smoothing on each axis is controlled independently by the length of $h$ and $g$.
- The SPWVD is the most general member of Cohen's class that we derive and apply.

---

## A.8: Two-Knob Independence

- Increasing $|h|$ (lag window length) sharpens frequency resolution but allows time-oscillating ghosts.
- Increasing $|g|$ (time window length) sharpens time resolution but allows frequency-oscillating ghosts.
- The two knobs are independent: adjusting one does not constrain the other.
- This is fundamentally different from the STFT, where a single window length controls both axes simultaneously.
- The SPWVD decouples the time-frequency tradeoff into two separate, manageable decisions.

---

## Volume A Summary: Tool Progression

| Tool | Adds | Costs | Limitation |
|------|------|-------|------------|
| DFT | Global spectrum | No time information | Cannot detect non-stationarity |
| STFT | Time-frequency plane | Heisenberg limit | $\Delta t \cdot \Delta f \geq \beta$ |
| WVD | Sub-Heisenberg sharpness | Cross-terms | Unusable for $\geq 2$ components |
| SPWVD | Controllable smoothing | Some resolution loss | Residual ghosts on complex signals |

- Each tool fixes a specific limitation of the previous one.
- The cost column explains why the simpler tool is not simply obsolete.

---

## Cohen's Class: A Unifying Framework

- All four tools (DFT, STFT, WVD, SPWVD) are members of Cohen's class of time-frequency representations.
- Cohen's class is defined as the WVD convolved with a two-dimensional kernel function.
- The STFT spectrogram uses the WVD of the window as its kernel.
- The SPWVD uses the product $h[m] \cdot g[p]$ as a separable kernel.
- The raw WVD uses a delta kernel (no smoothing), which is why it has perfect resolution but uncontrolled interference.

---

## Part 2: Volume B - Laboratories

Eight labs validate the theory on model signals with known ground truth.

- All signals: $f_s = 250$ Hz, duration $\geq 1200$ s, all components below 100 Hz.
- Each lab asks a specific question and answers it with quantitative evidence.
- The result figures are the deliverables; verification closes the loop against Volume A predictions.

---

## Lab 1: Purpose — What Happens at Off-Bin Frequencies?

- **Question:** What happens when a signal's frequency falls between DFT bins, and can zero-padding fix it?
- The DFT evaluates the spectrum at fixed bins spaced $\Delta f = f_s / N$ apart.
- A frequency that lands exactly on a bin is captured perfectly; one that falls between bins leaks energy across all bins.
- Zero-padding increases the number of frequency samples but does not add new information about the signal.
- This lab tests the distinction between spectral interpolation and genuine frequency resolution.

---

## Lab 1: On-Bin Capture Is Perfect

![](../../results/graphs/lab1/figure_B_02.png){ width=90% }

- The 10 Hz tone at $f_s = 250$ Hz, $N = 5000$ lands exactly on bin 200 ($\Delta f = 0.05$ Hz).
- All energy concentrates in a single DFT bin; every other bin reads exactly zero.
- This is the ideal case: the DFT captures the frequency with no distortion whatsoever.

---

## Lab 1: Off-Bin Leakage Reveals the Dirichlet Kernel

![](../../results/graphs/lab1/figure_B_02.png){ width=90% }

- The 10.5 Hz tone falls between bins 209 and 210, distributing energy across all $N$ bins.
- The leakage pattern traces the Dirichlet kernel shape, confirming the theoretical prediction from Section A.2.
- This is not a numerical error but a structural property of the finite-length rectangular window.
- The magnitude of leakage (first side lobe at $-13$ dB) matches the sinc approximation exactly.

---

## Lab 1: Zero-Padding Interpolates but Does Not Resolve

![](../../results/graphs/lab1/figure_B_08.png){ width=90% }

- Zero-padding the 10 + 10.5 Hz signal produces two apparent peaks separated by 0.5 Hz.
- The peaks appear resolved because zero-padding adds frequency samples between the original bins.
- However, this is spectral interpolation: no new information about the signal has been created.
- Genuine resolution of two tones requires observing the signal for a longer duration ($T > 1/\Delta f$).

---

## Lab 1: Conclusion — Leakage Motivates Windowing

- On-bin frequencies are captured perfectly; off-bin frequencies leak with the Dirichlet kernel envelope.
- Zero-padding is interpolation, not resolution. Two tones closer than $1/T$ Hz remain fundamentally unresolved.
- The rectangular window's $-13$ dB side lobes cause severe leakage for off-bin components.
- This motivates Lab 3: applying tapered windows to suppress the side lobes and reduce leakage.
- Resolution requires longer observation of the signal, not computational tricks after the fact.

---

## Lab 2: Purpose — Distinguishing Peaks from Noise

- **Question:** How do we determine whether a spectral peak is a real signal or a random noise fluctuation?
- Under white Gaussian noise, each DFT bin power follows an exponential distribution (CV = 1.0).
- The periodogram is unbiased but inconsistent: variance does not decrease with signal length.
- A principled detection threshold is needed, derived from the noise statistics rather than arbitrary heuristics.
- Welch averaging offers a path to reduce variance at the cost of frequency resolution.

---

## Lab 2: Bin Power Matches the Exponential Distribution

![](../../results/graphs/lab2/figure_B_01.png){ width=90% }

- The histogram of bin powers from white noise follows the exponential distribution with measured CV = 1.00.
- This confirms Section A.4: bin power fluctuates by 100% of its mean regardless of signal length.
- The detection threshold $\gamma = \ln(1/P_{fa})$ derives directly from this distribution.
- At $P_{fa} = 1\%$: $\gamma = 4.6$, meaning a genuine peak must exceed 4.6 times the noise floor.

---

## Lab 2: Welch Averaging Reduces Variance

![](../../results/graphs/lab2/figure_B_04.png){ width=90% }

- Welch averaging with $K$ segments reduces the standard deviation of the noise floor by $1/\sqrt{K}$.
- The spectrum becomes visually flatter, making genuine peaks stand out above the reduced noise fluctuations.
- The cost is frequency resolution: each segment has length $M < N$, degrading resolution from $f_s/N$ to $\beta \cdot f_s/M$.
- This is the resolution-variance tradeoff: better statistics require shorter segments that blur frequency.

---

## Lab 2: Conclusion — Detection Threshold and Welch Tradeoff

- Measured CV = 1.00, confirming the exponential model for noise bin powers (Section A.4).
- The threshold $P_{fa} = e^{-\gamma}$ provides a principled detection rule: $\gamma = 4.6$ for 1% false alarm.
- Welch trades resolution for reliability: variance reduction scales as $1/\sqrt{K}$, confirmed across multiple segment counts.
- This tradeoff is applied in Volume C to justify every Welch parameter choice on real EEG.

---

## Lab 3: Purpose — Deriving Window Spectra from First Principles

- **Question:** Why do different windows have different side-lobe levels, and can we derive them analytically?
- All standard windows (Hann, Hamming, Blackman) are sums of shifted cosines applied to the rectangular window.
- The rectangular window's DFT produces the Dirichlet kernel, which serves as the fundamental building block.
- The pure sine form factorization should predict all measurable window properties from the coefficients alone.
- This lab derives and verifies the complete window comparison table from Section A.3.

---

## Lab 3: Dirichlet Kernel Anatomy

![](../../results/graphs/lab3/figure_B_11.png){ width=90% }

- The Dirichlet kernel displays the main lobe (width $4\pi/M$), side lobes (first at $-13.3$ dB), and nulls at integer multiples of $2\pi/M$.
- Side-lobe decay follows $-20$ dB/octave, determined by the first-order zero at the origin.
- The main-lobe width scales inversely with window length $M$: longer windows produce narrower main lobes.
- Every property of cosine-sum windows is inherited from this kernel through shifted superposition.

---

## Lab 3: Pure Sine Form Derives All Window Spectra

- The spectrum of any cosine-sum window factors into shifted Dirichlet kernels via Equation (B.20).
- The shifted kernels destructively interfere with the side lobes of the central kernel.
- Hann achieves $-60$ dB/octave decay through a third-order zero; Blackman achieves the same through a fifth-order zero.
- The null positions, side-lobe levels, and decay rates are all predicted analytically from the cosine coefficients alone.
- This factorization is the key insight: window design reduces to choosing coefficients that cancel side lobes.

---

## Lab 3: Window Comparison — Four Spectra Overlaid

![](../../results/graphs/lab3/figure_B_14.png){ width=90% }

- Rectangular ($-13$ dB, $-20$ dB/oct), Hann ($-32$ dB, $-60$ dB/oct), Hamming ($-43$ dB, $-20$ dB/oct), Blackman ($-58$ dB, $-60$ dB/oct).
- Hann and Blackman achieve fast decay through higher-order zeros; Hamming optimizes the first side lobe only.
- Main-lobe widths: $2f_s/M$, $4f_s/M$, $4f_s/M$, $6f_s/M$ respectively.
- No window is universally best: the choice depends on whether leakage suppression or resolution is more critical.

---

## Lab 3: Conclusion — All Properties Derived from First Principles

- All measured first side-lobe levels match predictions within 0.1 dB across all four windows.
- Decay rate regressions confirm $-20$ or $-60$ dB/octave as predicted by the zero order at DC.
- The side-lobe cancellation mechanism is demonstrated analytically and confirmed numerically.
- Lab 3 establishes the quantitative foundation for all subsequent window choices in Labs 4, 5, and Volume C.

---

## Lab 4: Purpose — Can We See Frequency Change Over Time?

- **Question:** Can we track frequency evolution on a spectrogram, and what limits the STFT's resolution?
- The STFT slides a window along the signal and computes the DFT of each segment.
- A linear chirp (5 to 45 Hz over 120 s) provides a known instantaneous frequency trajectory for validation.
- The Heisenberg uncertainty principle predicts that $\Delta t \cdot \Delta f \geq \beta$, regardless of window length.
- This lab tests whether both short and long windows yield the same Heisenberg product.

---

## Lab 4: Short Window — Good Time, Poor Frequency

![](../../results/graphs/lab4/figure_B_17.png){ width=90% }

- A short window ($M = 64$, 0.256 s) tracks the chirp trajectory closely in time.
- The frequency axis is heavily blurred: $\Delta f = 7.81$ Hz, spanning nearly the entire chirp bandwidth.
- Time localization is excellent ($\Delta t = 0.256$ s), but individual frequency components cannot be distinguished.
- The spectrogram captures the time evolution at the cost of frequency precision.

---

## Lab 4: Long Window — Good Frequency, Poor Time

![](../../results/graphs/lab4/figure_B_18.png){ width=90% }

- A long window ($M = 1024$, 4.096 s) resolves frequency precisely: $\Delta f = 0.49$ Hz.
- The chirp trajectory is smeared along the time axis by $\Delta t = 4.096$ s.
- Frequency structure is clear, but the temporal progression of the chirp is lost over multi-second intervals.
- The complementary limitation confirms that time and frequency resolution cannot be improved simultaneously.

---

## Lab 4: Heisenberg Product Is Constant

- Short window: $\Delta t = 0.256$ s, $\Delta f = 7.81$ Hz, product $\Delta t \cdot \Delta f = 2.0$.
- Long window: $\Delta t = 4.096$ s, $\Delta f = 0.49$ Hz, product $\Delta t \cdot \Delta f = 2.0$.
- Both products equal the Hann-window Heisenberg constant $\beta = 2$, confirming Equation (A.49).
- The product is invariant: improving one axis degrades the other by exactly the same factor.
- No single window can simultaneously achieve fine time and fine frequency resolution.

---

## Lab 4: Conclusion — The Spectrogram Is the Primary Tool

- The STFT successfully adds a time axis to the spectrum, enabling detection of non-stationary behavior.
- The Heisenberg limit $\Delta t \cdot \Delta f \geq 2$ is confirmed quantitatively for the Hann window.
- The one-window limitation motivates the search for representations that decouple time and frequency resolution.
- The spectrogram is now the primary analysis tool; its resolution limit is tested in Lab 5.

---

## Lab 5: Purpose — When Do Two Tones Merge on the Spectrogram?

- **Question:** At what frequency separation do two tones become unresolvable on the spectrogram?
- The theoretical resolution limit for the Hann window is $\Delta f_{\min} = \beta \cdot f_s / M$ (from Lab 3's main-lobe width).
- For $M = 500$ (2.0 s) at $f_s = 250$ Hz: $\Delta f_{\min} = 1.0$ Hz.
- Two tones at 20 Hz and 20 + $\Delta f$ Hz are generated with varying separation to test this prediction.
- This lab connects Lab 3's analytical window theory to practical spectrogram performance.

---

## Lab 5: Tones Above the Resolution Limit Are Resolved

![](../../results/graphs/lab5/figure_B_28.png){ width=90% }

- Two tones separated by more than $\Delta f_{\min} = 1.0$ Hz appear as two distinct horizontal lines on the spectrogram.
- The separation is clearly visible in both the linear and dB representations.
- Each tone occupies its own frequency band with a visible gap between them.
- The spectrogram correctly identifies two separate spectral components.

---

## Lab 5: Tones Below the Limit Merge

- Two tones separated by less than $\Delta f_{\min}$ merge into a single spectral blob that cannot be split.
- No amount of post-processing or display adjustment can recover the two separate components.
- The merged blob is wider than a single tone would produce, but the two-component structure is invisible.
- This is a fundamental resolution limit, not a display or computational limitation.

---

## Lab 5: Conclusion — Lab 3 Predictions Confirmed on Spectrograms

- Predicted resolution limit: $\Delta f_{\min} = 2 f_s / M = 1.0$ Hz (Hann window, $M = 500$).
- Measured: tones at $\Delta f = 1.5$ Hz are resolved; tones at $\Delta f = 0.5$ Hz are merged.
- The transition occurs at approximately $\Delta f = 1.0$ Hz, matching the prediction exactly.
- Lab 3's window theory (main-lobe width $= 4f_s/M$) directly predicts the spectrogram resolution limit.
- The window choice for EEG is now justified: the Hann window at $M = 500$ resolves tones 1 Hz apart.

---

## Lab 6: Purpose — Detecting Periodicity Without Phase

- **Question:** Can we detect periodic structure and shared components without requiring phase information?
- The autocorrelation measures self-similarity at each lag, revealing periodicity as repeating peaks.
- The Wiener-Khinchin theorem connects autocorrelation to the power spectrum (which discards phase).
- Phase-blindness is both a limitation and a feature: it simplifies detection at the cost of temporal alignment.
- Cross-correlation extends this to pairs of signals, quantifying shared spectral content.

---

## Lab 6: Phase-Blindness Demonstrated

![](../../results/graphs/lab6/figure_B_33.png){ width=90% }

- Two signals with identical frequency content but different phase relationships produce identical autocorrelations.
- A cosine and a sine of the same frequency are indistinguishable by their autocorrelation.
- This confirms the Wiener-Khinchin theorem: the autocorrelation maps to the power spectrum, which discards phase.
- Phase-blindness means autocorrelation detects what frequencies are present, not when they start.

---

## Lab 6: Cross-Correlation Detects Shared Components

![](../../results/graphs/lab6/figure_B_37.png){ width=90% }

- Two signals sharing a common 10 Hz component produce a cross-correlation peak with $\rho = 0.107$.
- Components present in only one signal contribute no cross-correlation peak.
- The normalized correlation coefficient quantifies the degree of shared spectral structure between any two signals.
- This technique is applied in Volume C to test whether auxiliary channels (ECG, EMG) contaminate the EEG.

---

## Lab 6: Conclusion — Autocorrelation as the WVD Building Block

- Phase-blindness is confirmed: identical spectra with different phases produce identical autocorrelations.
- Cross-correlation provides a principled measure of shared structure between signals ($\rho = 0.107$ for the shared tone).
- The instantaneous autocorrelation $R_x[n, m] = z[n+m] \cdot z^*[n-m]$ is the time-local generalization.
- This building block, when Fourier-transformed over lag $m$, produces the Wigner-Ville Distribution (Lab 7).

---

## Lab 7: Purpose — Can the WVD Bypass the Heisenberg Limit?

- **Question:** Does the WVD achieve sub-Heisenberg resolution, and what is the cost for multi-component signals?
- The WVD computes the DFT of the instantaneous autocorrelation, using no analysis window.
- Without a window, the Heisenberg constraint does not apply: theoretically perfect concentration is possible.
- However, the bilinear (quadratic) nature of the WVD generates cross-terms between every pair of components.
- This lab tests both the resolution advantage and the cross-term catastrophe on controlled signals.

---

## Lab 7: Single Chirp — STFT Blurred, WVD Razor-Sharp

![](../../results/graphs/lab7/figure_B_40.png){ width=90% }

- The STFT produces a thick, blurred diagonal stripe limited by the Heisenberg bound ($\Delta t \cdot \Delta f = 2$).
- The WVD produces a razor-sharp line tracking the instantaneous frequency with zero spread.
- For a single component, the WVD completely bypasses the Heisenberg limit, achieving perfect concentration.
- This confirms Section A.7: the WVD is the ideal time-frequency representation for single-component signals.

---

## Lab 7: Two Components — Cross-Terms at the Midpoint

![](../../results/graphs/lab7/figure_B_42.png){ width=90% }

- A chirp plus a constant-frequency tone produces a clean two-component STFT spectrogram.
- The WVD shows both components sharply but adds an oscillating interference pattern at the midpoint.
- The cross-term is located at $(t_1 + t_2)/2$, $(f_1 + f_2)/2$ and is as energetic as the auto-terms.
- The display becomes uninterpretable: one cannot distinguish real components from ghost artifacts.

---

## Lab 7: Real Versus Analytic Signal — DC Ghost Removal

![](../../results/graphs/lab7/figure_B_44.png){ width=90% }

- The real-valued signal's symmetric spectrum generates a cross-term at DC (the "self-ghost").
- Computing the analytic signal via the Hilbert transform removes negative frequencies and eliminates the DC ghost.
- The analytic signal is a prerequisite for meaningful WVD computation on any real-valued data.
- All subsequent WVD and SPWVD computations use the analytic signal exclusively.

---

## Lab 7: Conclusion — Super-Resolution Is Real but Cross-Terms Are Catastrophic

- Single chirp: the WVD peak tracks the true instantaneous frequency with zero deviation (sub-Heisenberg confirmed).
- Cross-term location: measured at the midpoint in both time and frequency, confirming Section A.7.3.
- Cross-term oscillation frequency: proportional to the component separation, as predicted analytically.
- The DC ghost is completely removed by the analytic signal transformation.
- Super-resolution is real but cross-terms make the raw WVD unusable for multi-component signals.

---

## Lab 8: Purpose — Suppressing Cross-Terms Without Losing Sharpness

- **Question:** Can the SPWVD suppress cross-terms while preserving the WVD's time-frequency concentration?
- The PWVD applies a lag window $h[m]$ to smooth the frequency axis (suppresses frequency-oscillating ghosts).
- The SPWVD adds a time window $g[p]$ to smooth the time axis (suppresses time-oscillating ghosts).
- The two windows are independent: adjusting one does not constrain the other.
- This lab tests the step-by-step progression from WVD to PWVD to SPWVD on controlled multi-component signals.

---

## Lab 8: WVD to PWVD to SPWVD Progression

![](../../results/graphs/lab8/figure_B_46.png){ width=90% }

- **WVD** (top): sharp auto-terms with heavy cross-term contamination between all component pairs.
- **PWVD** (middle): lag window $h$ smooths frequency, suppressing frequency-oscillating ghosts; time-oscillating ghosts survive.
- **SPWVD** (bottom): time window $g$ added, suppressing both ghost types with controllable resolution cost.
- The progression demonstrates that each smoothing operation targets a specific ghost type independently.

---

## Lab 8: Two Impulses — Frequency-Oscillating Ghosts Suppressed

![](../../results/graphs/lab8/figure_B_49.png){ width=90% }

- Two impulses separated in time generate frequency-oscillating cross-terms in the WVD.
- The PWVD (lag window only) suppresses these ghosts because it smooths along the frequency axis.
- This confirms the duality: time-separated components generate ghosts that oscillate along the frequency axis.
- The suppression is 15-25 dB depending on the lag window length $|h|$.

---

## Lab 8: Two-Knob Sweep — Independent Control

![](../../results/graphs/lab8/figure_B_47.png){ width=90% }

- Case 1 ($|h|=101$, $|g|=5$): strong frequency smoothing, sharp time resolution, but time-oscillating ghosts survive.
- Case 2 ($|h|=25$, $|g|=31$): both axes smoothed moderately, both ghost types suppressed.
- The two knobs are independent: changing $|h|$ does not affect time resolution, changing $|g|$ does not affect frequency resolution.
- This is fundamentally different from the STFT, where a single window length controls both axes simultaneously.

---

## Lab 8: Conclusion — Duality Confirmed, Ready for Real EEG

- The PWVD suppresses frequency-oscillating ghosts by 15-25 dB; the SPWVD achieves 20-30 dB reduction overall.
- The duality is exact: time-separated components need frequency smoothing, frequency-separated components need time smoothing.
- The two-knob tradeoff is confirmed: doubling $|h|$ halves frequency resolution, doubling $|g|$ halves time resolution.
- The SPWVD handles both ghost types with independent control. It is calibrated and ready for real EEG in Volume C.

---

## Appendix B: The M vs M-1 Convention

- The periodic window convention (denominator $M$) and the symmetric convention (denominator $M-1$) produce slightly different window shapes.
- For DFT-based spectral analysis, the periodic convention is correct because it aligns with the DFT's implicit periodicity.
- Both conventions converge as $M \rightarrow \infty$: the difference is $O(1/M)$.
- At $M = 256$ (typical for EEG): the maximum pointwise difference between conventions is less than 0.1%.
- This report uses the periodic convention throughout, consistent with the DFT framework.

---

## Appendix B: CV of Signal Archetypes

- **White noise:** CV $\approx 1.0$ (exponential distribution of bin powers).
- **Pure tone (on-bin):** CV $\approx 87$ (almost all energy in one bin, the rest near zero).
- **Linear chirp:** CV $\approx 2.3$ (energy spread across many bins but not uniformly).
- The CV is a single-number diagnostic: it distinguishes signal types from their spectral statistics alone.
- Volume C uses this diagnostic to test whether the EEG noise floor matches the exponential model.

---

## Volume B Summary

- Eight labs form a progressive chain: each lab's conclusion motivates the next lab's question.
- The key finding at each stage: leakage (Lab 1), detection (Lab 2), windowing (Lab 3), time-frequency (Lab 4), resolution (Lab 5), phase-blindness (Lab 6), super-resolution (Lab 7), smoothing (Lab 8).
- All Volume A predictions are confirmed quantitatively with measured values matching theory.
- The SPWVD emerges as the most capable tool, with independently controllable time and frequency smoothing.
- The complete toolchain is now validated and ready for application to real neonatal EEG in Volume C.

---

## Part 3: Volume C - Real EEG Application

Each section asks a question, explains why a specific tool answers it, then presents the result.

- The analysis is question-driven: C.1 triages the signal, and subsequent sections follow where the data leads.
- Every tool applied here was derived in Volume A and validated in Volume B.
- The pattern: Question, Why This Tool, Result, Answer.

---

## What Signal Are We Working With?

- **Subject:** NORB00055 (neonatal), European Data Format (EDF).
- **Sampling rate:** 200 Hz. **Duration:** 1140 s (19 minutes). **Channels:** 24 (19 EEG + 5 auxiliary).
- Data loaded via MNE-Python. All amplitudes reported in $\mu$V.
- This is a single continuous recording of a neonatal brain during quiet sleep.
- The signal characteristics are unknown at the outset; triage (C.1) determines the analysis direction.

---

## Why CZ as the Primary Channel?

- CZ sits at the vertex (top center of the scalp) in the international 10-20 electrode system.
- It receives contributions from both hemispheres with minimal regional bias.
- This makes CZ the least biased single channel for whole-brain triage.
- Other channels are used only when the analysis specifically requires it (e.g., artifact assessment in C.4).
- Every channel switch is stated explicitly with justification.

---

## What Tools Do We Bring?

- All tools are imported unchanged from `src/common/`: the same infrastructure validated in Volume B.
- The toolchain: DFT, Welch PSD, STFT, autocorrelation, cross-correlation, WVD, SPWVD.
- **No filtering is applied.** Filter theory (FIR/IIR design) was not derived in Volumes A-B.
- The only preprocessing available is **segment selection**: choosing clean, short segments for WVD analysis.
- Window functions use the same periodic convention validated in Lab 3; plotting follows the dual-stack rule.

---

## C.1: What Kind of Signal Is This EEG?

- **The question:** Before choosing analysis tools, we must classify the signal. What archetype does this neonatal EEG resemble?
- Is it stationary or non-stationary? Periodic or bursty? Narrowband or broadband?
- The triage uses two tools: the time-domain plot (visual inspection) and the Welch PSD (spectral summary).
- The classification determines the direction for all subsequent sections (C.2 through C.5).
- This is the data-first principle: the signal dictates the analysis, not a predetermined sequence.

---

## C.1 Result: Time-Domain Plot

![](../../results/graphs/volume_c/c1/figure_C_02.png){ width=90% }

- The full 19-minute CZ recording shows large-amplitude bursts separated by quiet intervals.
- Peak-to-peak amplitude reaches approximately $\pm 200$ $\mu$V during bursts.
- The signal is visibly non-stationary: activity is concentrated in discrete episodes, not continuous oscillation.
- The burst pattern is immediately apparent without any spectral analysis.

---

## C.1 Result: Welch PSD Dual-Stack

![](../../results/graphs/volume_c/c1/figure_C_05.png){ width=90% }

- The Welch PSD (linear top, dB bottom) shows overwhelming delta-band dominance.
- Power drops steeply above 4 Hz, with no visible alpha or beta peaks.
- The dB representation reveals the full dynamic range: over 40 dB between delta peak and high-frequency floor.
- The spectrum is consistent with a delta-dominated, low-frequency signal.

---

## C.1 Result: Band Power Distribution

![](../../results/graphs/volume_c/c1/figure_C_06.png){ width=90% }

- Delta (0.5-4 Hz): **91.8%** of total power.
- Theta (4-8 Hz): 5.7%. Alpha (8-13 Hz): 1.0%. Beta (13-30 Hz): 1.0%.
- The absence of alpha and beta rhythms is consistent with neonatal EEG physiology.
- A single frequency band accounts for nearly all signal energy.

---

## C.1 Triage Decision

- **Classification:** Delta-dominated, bursty, non-stationary, whole-brain synchronous.
- The multi-channel heatmap confirms that the burst pattern appears across all 19 EEG channels simultaneously.
- This signal belongs to the "random non-stationary" archetype from the signal taxonomy (Appendix A).
- **Direction set:** C.2 investigates whether delta is real oscillation or 1/f background. C.3 uses the STFT to characterize the burst pattern in time. C.4 verifies channel cleanliness before C.5 applies the WVD.
- Each subsequent section must explain why it was chosen based on what triage found.

---

## C.2: Is Delta Rhythmic Oscillation or 1/f Background?

- **The question:** The triage showed delta dominance, but is the 91.8% delta power a genuine neural oscillation or merely the steep spectral slope that all biological signals exhibit?
- **Why the DFT answers this:** The log-log PSD separates 1/f background (a straight line) from discrete peaks (deviations above the line). If delta peaks exist above the 1/f extrapolation, they represent structured oscillatory activity.
- The stationary DFT provides the highest frequency resolution on the full recording ($\Delta f = 0.4$ Hz with $M = 1000$, 5.0 s segments).
- Welch parameters justified by Lab 2: $K \approx 455$ segments yields low variance; Lab 3 validates the Hann window choice.

---

## C.2 Result: Log-Log PSD with 1/f Fit

![](../../results/graphs/volume_c/c2/figure_C_09.png){ width=90% }

- On a log-log scale, the PSD follows $1/f^{3.18}$ from 5-40 Hz.
- The slope ($-3.18$) is steeper than pink noise ($-1$) or brown noise ($-2$).
- Delta peaks at 0.4-0.6 Hz sit above the 1/f extrapolation, indicating genuine rhythmic activity beyond the background slope.

---

## C.2 Result: Delta-Band Zoom

![](../../results/graphs/volume_c/c2/figure_C_10.png){ width=90% }

- Zooming into the 0-5 Hz range reveals distinct peaks at 0.4, 0.5, and 0.6 Hz.
- These peaks rise above the 1/f fit, confirming they represent quasi-periodic delta oscillations.
- The peak structure suggests a fundamental frequency near 0.5 Hz with slight spread.
- However, the stationary DFT cannot determine whether these oscillations persist continuously or arrive in bursts.

---

## C.2 Finding

- Delta peaks exist and represent genuine oscillatory activity above the 1/f background.
- The spectral slope of $-3.18$ is steeper than typical colored noise, indicating strong low-frequency dominance.
- **The limitation:** The stationary DFT collapses all 19 minutes into a single spectrum. It cannot distinguish continuous oscillation from intermittent bursts.
- The time axis is needed to resolve this ambiguity.
- This directly motivates C.3: the STFT adds the time dimension that the DFT lacks.

---

## C.3: Does Delta Persist Continuously or Come in Bursts?

- **The question:** C.2 confirmed delta peaks exist, but are they present throughout the recording or concentrated in discrete episodes? The time-domain plot suggests bursts, but quantification requires a time-frequency tool.
- **Why the STFT answers this:** The STFT provides a time-resolved spectrum (Lab 4). By extracting delta-band power as a function of time, it can detect when delta activity is present and when it is absent.
- The STFT is the correct tool here because it adds the time axis that the stationary DFT lacks.
- The Heisenberg limit ($\Delta t \cdot \Delta f \geq 2$) constrains precision but does not prevent burst detection.

---

## C.3 Result: Full-Recording Spectrogram

![](../../results/graphs/volume_c/c3/figure_C_12.png){ width=90% }

- The full 19-minute spectrogram reveals delta bursts as bright vertical stripes concentrated below 4 Hz.
- The burst pattern is clearly non-stationary: bursts occur irregularly with quiet gaps between them.
- The STFT transforms the single-spectrum view of C.2 into a time-resolved map that confirms the bursty hypothesis.

---

## C.3 Result: Zoomed 60-Second Segment

![](../../results/graphs/volume_c/c3/figure_C_13.png){ width=90% }

- Zooming to 60 seconds reveals individual burst events with duration of 5-15 seconds each.
- Between bursts, the spectrogram shows only low-level broadband activity (the inter-burst interval).
- The burst structure is the dominant non-stationary feature of this neonatal EEG.
- Individual burst onsets and offsets are visible but blurred by the Heisenberg limit.

---

## C.3 Result: Burst Detection and Overlay

![](../../results/graphs/volume_c/c3/figure_C_16.png){ width=90% }

- A burst threshold of $2 \times$ median delta power ($= 3479$ $\mu$V$^2$) identifies burst episodes.
- **19%** of the recording is classified as burst; **81%** is quiet (inter-burst interval).
- The maximum-to-median power ratio is 17x, confirming that bursts are highly energetic events.
- Detected bursts overlay precisely on the visible high-amplitude episodes in the raw time-domain signal.

---

## C.3 Answer

- Delta activity is **discontinuous**. It arrives in bursts occupying 19% of the recording, not as continuous oscillation.
- This is the dominant non-stationary feature of this neonatal EEG signal.
- The STFT resolved the ambiguity that the stationary DFT could not: delta peaks in the spectrum come from intermittent high-power episodes, not sustained rhythmic activity.
- **The limitation:** The STFT was Heisenberg-limited in its ability to precisely localize burst onsets and offsets. The WVD/SPWVD (C.5) may sharpen this.
- Before applying the computationally expensive WVD, C.4 first verifies that the channel is clean.

---

## C.4: Is CZ Clean Enough for WVD Analysis?

- **The question:** Before applying the WVD (which is sensitive to any signal component), we must verify that CZ is not contaminated by non-brain signals (ECG, EMG, EOG).
- **Why cross-correlation answers this:** Lab 6 established that the normalized correlation coefficient $\rho$ quantifies shared spectral structure between two signals. If auxiliary channels show high $\rho$ with CZ, artifact contamination is present.
- Additionally, the CV diagnostic from Lab 2 tests whether the noise floor follows the expected exponential model.
- A clean channel means the WVD will analyze brain activity, not cardiac or muscular artifacts.

---

## C.4 Result: Cross-Correlation Versus CZ

![](../../results/graphs/volume_c/c4/figure_C_20.png){ width=90% }

- All auxiliary-vs-CZ correlation coefficients are below $\rho = 0.03$.
- This confirms that no detectable ECG, EMG, or EOG artifact contaminates the CZ channel.
- Inter-EEG correlations range from $\rho = 0.47$ to $\rho = 0.78$ (shared brain activity, expected for whole-brain synchronous bursts).
- The cross-correlation technique from Lab 6 provides a principled, quantitative artifact assessment.

---

## C.4 Result and Conclusion

- The noise floor CV = 1.11 (alpha band, 8-13 Hz), close to the ideal exponential value of 1.00.
- The 11% deviation is small; possible causes include residual physiological activity or non-Gaussianity of biological noise.
- **Verdict:** CZ is clean for WVD analysis. No auxiliary signal contamination detected.
- Artifacts were not removed from the recording. Instead, C.5 selects a clean burst segment for WVD analysis.
- This is consistent with the project methodology: segment selection is the only preprocessing available without filter theory.

---

## C.5: Can the SPWVD Sharpen the Burst Timing Beyond the STFT?

- **The question:** The STFT resolved the burst pattern but the Heisenberg limit blurred burst onsets and offsets. Can the SPWVD achieve sub-Heisenberg resolution on a selected burst, revealing finer temporal structure within the burst?
- **Why the SPWVD answers this:** Lab 7 demonstrated that the WVD bypasses Heisenberg for single components. Lab 8 showed that the SPWVD suppresses cross-terms while preserving much of the WVD's sharpness.
- The strategy: select a short, clean burst segment (few components, no artifacts) where the SPWVD can operate effectively.
- This is the payoff question of the entire tool progression.

---

## C.5: How Do We Define and Select a Burst Segment?

- Burst detection reuses the C.3 threshold ($2 \times$ median delta power).
- **First attempt:** The strongest burst (t = 842.5 s, power = 15597 $\mu$V$^2$, 8.7x median) is rejected because it contains 44 flat samples indicating amplifier saturation.
- **Fallback:** The 75th percentile burst (t = 65.0 s, power = 5435 $\mu$V$^2$, 3.0x median) is clean with no flat samples or discontinuities.

![](../../results/graphs/volume_c/c5/figure_C_23.png){ width=90% }

- The accepted segment contains approximately 3 clear delta cycles over 2 seconds. It is suitable for WVD/SPWVD analysis.

---

## C.5: What Does the STFT Show on This Segment?

![](../../results/graphs/volume_c/c5/figure_C_24.png){ width=90% }

- The STFT of the selected burst segment shows a broad delta blob, limited by Heisenberg uncertainty.
- Frequency resolution at $M = 200$ (1.0 s): $\Delta f = 2.0$ Hz, spanning the entire delta band width.
- The STFT confirms "there is delta activity" but cannot resolve fine temporal structure within the burst.
- This establishes the baseline that the WVD and SPWVD aim to improve upon.

---

## C.5: What Does the Raw WVD Show?

![](../../results/graphs/volume_c/c5/figure_C_25.png){ width=90% }

- The raw WVD of the clean burst segment shows severe cross-term contamination.
- **49% of all values are negative.** The WVD is not a valid power distribution for this signal.
- Even a 2-second segment with approximately 3 delta cycles generates enough cross-terms to corrupt the representation.
- This confirms Lab 7's prediction (Section A.7.3): the raw WVD is unusable on multi-component real signals.

---

## C.5: What Does the SPWVD Recover?

![](../../results/graphs/volume_c/c5/figure_C_26.png){ width=90% }

- The SPWVD with calibrated windows ($|h| = 65$, $|g| = 31$) produces the sharpest readable time-frequency view.
- The **linear-scale** representation shows individual delta cycles resolved in time, with frequency content concentrated below 2 Hz.
- The **dB-scale** representation still shows residual circular artifacts (cross-term remnants at reduced amplitude).
- The linear panel is the highest-quality time-frequency representation achieved in this report.

---

## C.5: Three-Way Comparison — STFT Versus WVD Versus SPWVD

![](../../results/graphs/volume_c/c5/figure_C_27.png){ width=90% }

- **STFT:** Blurred but readable, no artifacts. The safe baseline.
- **WVD:** Sharp auto-terms visible but corrupted by oscillating cross-terms. Uninterpretable.
- **SPWVD:** The sharpest readable representation, combining sub-Heisenberg resolution with interpretability.
- The progression validates the entire theoretical framework: each tool fixes the previous tool's limitation.

---

## C.5: Why Not Filter Before the WVD?

- A bandpass filter (0.5-4 Hz) would isolate the delta band, reducing the number of components entering the WVD.
- Fewer components means fewer cross-terms and a dramatically cleaner SPWVD.
- However, filter theory (FIR/IIR design) was not derived in Volumes A or B.
- Applying a tool without deriving it violates the report's foundational principle of self-containment.
- **Segment selection** is the only preprocessing available, and it produces a usable result. This is a scope boundary, not a failure.

---

## C.6: What Worked and What Did Not?

| Tool | Applied To | Finding |
| ---- | ---------- | ------- |
| DFT + Welch | Full recording | Delta dominance (91.8%), 1/f slope $= -3.18$ |
| STFT | Full recording | Burst pattern: 19% burst, 81% quiet, max/median = 17x |
| Cross-correlation | Auxiliary vs CZ | CZ is clean ($\rho < 0.03$) |
| SPWVD | 2 s burst segment | Sub-Heisenberg resolution of individual delta cycles |

- Each tool contributed a finding that the previous tool could not provide.
- The progression DFT, STFT, SPWVD provided increasingly detailed views of the same signal.

---

## C.6: The Honest Limits

- **Raw WVD on real EEG:** 49% negative values, cross-term soup, completely unreadable.
- **SPWVD dB representation:** Residual circular cross-term artifacts survive even after dual smoothing.
- **Exponential noise model:** CV = 1.11 versus ideal 1.00. The model is approximate on real biological data (11% deviation).
- **Strongest burst rejected:** Amplifier saturation (44 flat samples) made the highest-power segment unusable.
- These are the method's honest limits, reported transparently rather than hidden or excused.

---

## C.6: The Scope Boundary and Future Work

- This report derives and applies: DFT, windowing, Welch, STFT, autocorrelation, WVD, SPWVD.
- It does **not** derive: FIR/IIR filter design, adaptive filtering, wavelet transforms, independent component analysis.
- The natural next step: use the window functions from Lab 3 to design an FIR bandpass filter via the frequency-sampling method.
- Bandpass pre-filtering before the WVD would dramatically reduce cross-terms, making the SPWVD effective on longer, multi-component segments.
- The scope boundary is a principled choice, not an oversight. The connection to future work is explicit.

---

## Part 4: Closing

The full progression applied to a real neonatal EEG signal, from the simplest tool to the most advanced.

---

## The Full Progression on Real EEG

$$\text{DFT} \xrightarrow{\text{+time}} \text{STFT} \xrightarrow{\text{+sharpness}} \text{WVD} \xrightarrow{\text{+smoothing}} \text{SPWVD}$$

- **DFT (C.2):** Identified delta dominance (91.8%) and the 1/f spectral slope ($-3.18$).
- **STFT (C.3):** Revealed that delta is discontinuous (bursty): 19% burst, 81% quiet, max/median = 17x.
- **WVD (C.5):** Failed on real EEG (49% negative values), confirming the theoretical prediction of cross-term catastrophe.
- **SPWVD (C.5):** Produced the sharpest readable view of individual delta cycles within a burst.
- Each step provided information that the previous tool could not access.

---

## What Worked: Detailed Summary

| Method | Context | Result |
|--------|---------|--------|
| Welch PSD | CZ, full recording | 91.8% delta, 1/f slope $= -3.18$ |
| Band power | All channels | Whole-brain delta synchrony confirmed |
| STFT burst detection | CZ, full recording | 19% burst, median burst duration 8 s |
| Cross-correlation | Auxiliary vs CZ | All $\rho < 0.03$ (no contamination) |
| CV diagnostic | Alpha band, CZ | CV = 1.11 (near-exponential noise floor) |
| SPWVD (linear) | 2 s clean burst | Individual delta cycles resolved |

---

## What Did Not Work: Honest Limits

- **Raw WVD on real EEG:** 49% negative values, cross-term soup, completely unreadable.
- **SPWVD dB representation:** Residual circular cross-term artifacts survive even after dual smoothing.
- **Exponential noise model:** CV = 1.11 versus ideal 1.00; the model is approximate on real biological data.
- **Strongest burst:** Amplifier saturation (44 flat samples) — an instrument artifact, not a physiological event.
- These are the method's honest limits, reported transparently rather than hidden or excused.

---

## The Scope Boundary: Filtering as Future Work

- The SPWVD's remaining limitation on real EEG is multi-component interference within each segment.
- A bandpass filter would isolate the delta band before WVD computation, reducing cross-terms dramatically.
- Filter design requires FIR/IIR theory that was not derived in this report.
- The natural next step: use the window functions from Lab 3 to design an FIR filter via the frequency-sampling method.
- This future work is well-motivated by the current results and directly connected to the existing theory.

---

## Closing Claim

The progression **DFT $\rightarrow$ STFT $\rightarrow$ SPWVD** provides increasingly detailed views of the same neonatal EEG signal.

- The SPWVD achieves **sub-Heisenberg resolution** on selected clean segments.
- Its practical value lies in the **linear-scale representation** of individual burst cycles.
- The dB representation is partially compromised by residual cross-terms in the absence of pre-filtering.
- The raw WVD is theoretically optimal for single components but practically unusable on real multi-component EEG.
- All claims are supported by quantitative measurements verified against Volume A predictions.

---

## Clinical Context

- The burst pattern observed is consistent with **normal discontinuous neonatal EEG activity** (trace discontinue).
- Reference: Lamblin et al. (1999) describe discontinuous activity as a hallmark of neonatal EEG during quiet sleep.
- Reference: Andre et al. (2010) provide normative parameters for burst-interburst patterns in preterm and term neonates.
- **No clinical diagnosis is made or implied.** These are signal-processing observations on a single recording.
- The association between signal features and physiology is cited, not asserted by the authors.

---

## References and Acknowledgments

- Lamblin, M.D. et al. (1999). Electroencephalography of the premature and term newborn.
- Andre, M. et al. (2010). Electroencephalography in premature and full-term infants.
- Oppenheim, A.V. & Schafer, R.W. (2010). Discrete-Time Signal Processing.
- Cohen, L. (1995). Time-Frequency Analysis.
- Boashash, B. (2016). Time-Frequency Signal Analysis and Processing.

All code, figures, and analysis are reproducible from the project repository.

**Thank you.**
