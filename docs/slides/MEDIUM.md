---
marp: true
theme: default
paginate: true
math: mathjax
style: |
  section h2 { margin-top: 0.1em; }
---

# From the DFT to the SPWVD

## Time-Frequency Analysis Applied to Neonatal EEG

**Authors:** Nguyen Duc Hung - 20233960, Bui Phuong Duy - 23233957, Tran Viet Bach - 23233954

Digital Biosignal Processing - Final Report

---

## Three-Volume Structure

- **Volume A** covers the theory: DFT fundamentals, windowing, spectral statistics, the STFT, autocorrelation, the WVD, and the SPWVD.
- **Volume B** contains eight laboratories that validate each theory section on controlled model signals.
- **Volume C** applies the full toolkit to a real neonatal EEG recording (19 minutes, 24 channels).
- Each volume builds on the previous one: theory derives the tools, labs verify their properties, and the application tests their limits on real data.

---

## Why EEG Needs Time-Frequency Analysis

- EEG signals are **non-stationary**: their spectral content changes on the scale of seconds (bursts, spindles, transitions between sleep stages).
- EEG signals are **multi-component**: multiple rhythms (delta, theta, alpha, beta) coexist at any moment, each with independent temporal dynamics.
- EEG signals are **noisy**: the signal-to-noise ratio is low, and physiological artifacts (eye blinks, muscle activity) contaminate the recording.
- A single DFT of the entire recording produces a global average that masks all temporal structure.
- The progression from DFT to STFT to SPWVD addresses each of these challenges incrementally.

---

<!-- Part 1 -->

## The Tool Progression at a Glance

DFT *(no time)* $\rightarrow$ STFT *(Heisenberg-limited)* $\rightarrow$ WVD *(cross-terms)* $\rightarrow$ SPWVD *(controllable)*

- Each tool in this progression addresses one specific limitation of the previous tool.
- The DFT provides a global spectrum but discards all temporal information.
- The STFT restores the time axis but introduces an inescapable resolution tradeoff (Heisenberg uncertainty).
- The WVD bypasses Heisenberg but generates destructive interference patterns (cross-terms) for multi-component signals.
- The SPWVD suppresses cross-terms through two independent smoothing windows, recovering usability.

---

## The DFT: Definition and Frequency Bins

$$X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j2\pi kn/N} \qquad \text{(A.5)}$$

$$\Delta f = \frac{f_s}{N} \qquad \text{(A.6)}$$

- The DFT decomposes a finite-length discrete signal into $N$ frequency components uniformly spaced by $\Delta f$.
- The bin spacing $\Delta f = f_s / N$ is entirely determined by the sampling rate and the number of samples.
- The DFT provides exact magnitudes only when a signal frequency coincides with a bin center (the on-bin condition).
- Increasing $N$ (recording longer) is the only way to reduce $\Delta f$ and improve true frequency resolution.

---

## On-Bin versus Off-Bin: Implications for Real Signals

- An **on-bin** frequency ($f = k \cdot \Delta f$ for integer $k$) concentrates all its energy in a single DFT bin.
- An **off-bin** frequency leaks energy into every bin, with a sinc-like envelope that decays slowly.
- Real-world signals almost never align with the DFT grid, so spectral leakage is the default condition.
- Zero-padding interpolates between existing bins (smoother-looking spectrum) but does **not** reduce leakage or improve resolution.
- The only true remedy for leakage is windowing, which reshapes the leakage pattern at the cost of main-lobe width.

---

## Windowing: The Dirichlet Kernel Concept

- Multiplying a signal by a rectangular window of length $M$ produces the Dirichlet kernel $D_M(f)$ in the frequency domain.
- The Dirichlet kernel has a main lobe of width $4\pi/M$ (in normalized angular frequency) and side lobes that decay at only $-20$ dB per octave.
- Every cosine-sum window (Hann, Hamming, Blackman) is constructed by adding shifted copies of the Dirichlet kernel.
- The shifted copies cancel the side lobes at the expense of widening the main lobe.
- The window choice is the fundamental tradeoff: narrow main lobe (resolution) versus low side lobes (leakage suppression).

---

## Window Comparison: Resolution versus Leakage

![w:900](../../results/graphs/lab3/figure_B_14.png)
**Figure B.14** - Window spectra comparison: Rectangular, Hann, Hamming, Blackman.

- First side-lobe levels: $-13$ dB (Rectangular), $-32$ dB (Hann), $-43$ dB (Hamming), $-58$ dB (Blackman).
- Main-lobe widths (in bins): 2 (Rectangular), 4 (Hann), 4 (Hamming), 6 (Blackman).
- The Hann window offers the best compromise for general time-frequency analysis.

---

## Spectral Statistics: The Exponential Distribution

$$P_{fa} = e^{-\gamma} \qquad \text{(A.37)}$$

- The power in each DFT bin of white Gaussian noise follows an exponential distribution with coefficient of variation CV = 1.0.
- A spectral peak is statistically significant only if it exceeds the noise floor by a threshold $\gamma$ set by the acceptable false-alarm probability $P_{fa}$.
- At $P_{fa} = 0.01$: the threshold is $\gamma = 4.6$ (the peak must be 4.6 times the mean noise power).
- This framework converts spectral peaks from visual impressions into quantitative detections.

---

## Welch Averaging: Trading Resolution for Reliability

- Welch's method divides the signal into $K$ overlapping segments, computes a periodogram for each, and averages.
- Averaging reduces the variance of the spectral estimate by a factor of approximately $1/K$.
- The cost is reduced frequency resolution: each segment is shorter than the full signal, so $\Delta f$ increases.
- For EEG analysis, this tradeoff is essential because single-segment periodograms are too noisy to interpret.
- The overlap (typically 50%) partially recovers the lost data at segment boundaries.

---

## The STFT: Adding a Time Axis

$$\Delta t \cdot \Delta f \geq \beta \qquad \text{(A.49)}$$

- The STFT applies a sliding window of length $M$ to the signal and computes the DFT of each windowed segment.
- The result is a time-frequency representation: power as a function of both time and frequency.
- Time resolution is $\Delta t \approx M/f_s$; frequency resolution is $\Delta f \approx f_s/M$ (for Hann: $\Delta f = 2f_s/M$).
- The Heisenberg uncertainty principle states that the product $\Delta t \cdot \Delta f$ cannot be reduced below a constant $\beta$ that depends on the window shape.

---

## Heisenberg Uncertainty: Practical Implications

- For the Hann window: $\beta = 2$, meaning $\Delta t \cdot \Delta f \geq 2$.
- A window of $M = 256$ samples at $f_s = 250$ Hz gives $\Delta t = 1.02$ s and $\Delta f = 1.95$ Hz.
- Halving the window ($M = 128$) halves $\Delta t$ to 0.51 s but doubles $\Delta f$ to 3.91 Hz.
- There is no window shape that escapes this tradeoff; it is fundamental to the STFT framework.
- The only escape is to abandon the linear time-frequency framework entirely, which motivates the WVD.

---

## The WVD: Definition

$$W_x[n, k] = \sum_{m} R_x[n, m] \, e^{-j2\pi km/N_f} \qquad \text{(A.61)}$$

where $R_x[n, m] = z[n+m] \cdot z^*[n-m]$ is the instantaneous autocorrelation.

- The WVD computes the DFT of the instantaneous autocorrelation at each time index $n$.
- It is a **quadratic** (bilinear) transform: it operates on products of the signal with itself.
- For a single linear chirp, the WVD produces a perfectly sharp diagonal line with zero spreading.
- The analytic signal $z[n]$ (via the Hilbert transform) is used to eliminate the DC self-interference term.

---

## The WVD: Chirp Sharpness and the Cross-Term Problem

- A single-component signal produces a WVD that exactly traces the instantaneous frequency with no Heisenberg blurring.
- For two components, the bilinearity generates a **cross-term** located at the midpoint between them in both time and frequency.
- Cross-terms oscillate at a rate proportional to the distance between the components.
- Cross-terms can be as energetic as the auto-terms (the real components), making visual interpretation impossible.
- For $N$ components, there are $N(N-1)/2$ cross-terms versus only $N$ auto-terms; the problem grows quadratically.

---

## The SPWVD: Definition

$$\text{SPWVD}_x[n, k] = \sum_m h[m] \left(\sum_p g[p] \, z[n\!+\!p\!+\!m] \, z^*[n\!+\!p\!-\!m]\right) e^{-j2\pi km/N_f} \quad \text{(A.72)}$$

- The SPWVD adds two independent smoothing windows to the WVD computation.
- **Lag window** $h[m]$: applied in the lag (frequency) direction, it smooths along the frequency axis.
- **Time window** $g[p]$: applied in the time direction, it smooths along the time axis.
- Each window suppresses one type of cross-term while trading away resolution in its respective axis.

---

## The SPWVD: Two Independent Knobs

- Cross-terms from components separated in **frequency** oscillate along the **time** axis; the time window $g$ suppresses them.
- Cross-terms from components separated in **time** oscillate along the **frequency** axis; the lag window $h$ suppresses them.
- Unlike the STFT, where a single window length controls both resolutions simultaneously, the SPWVD decouples them.
- Increasing $g$ blurs time resolution but suppresses time-oscillating ghosts; increasing $h$ blurs frequency resolution but suppresses frequency-oscillating ghosts.
- The SPWVD is not Heisenberg-free, but it offers a **two-dimensional** tradeoff instead of the STFT's one-dimensional tradeoff.

---

## Tool Progression: Summary Table

| Tool | Adds | Costs | Fundamental Limitation |
| ---- | ---- | ----- | ---------------------- |
| DFT | Global spectrum | No time information | Cannot detect non-stationarity |
| STFT | Time-frequency plane | Heisenberg limit | $\Delta t \cdot \Delta f \geq \beta$ |
| WVD | Sub-Heisenberg sharpness | Cross-terms | Unusable for $\geq 2$ components |
| SPWVD | Controllable smoothing | Some resolution loss | Residual ghosts on complex signals |

Each row addresses a specific limitation of the row above it.

---

## The Story Arc of Volume A

- The DFT provides the foundational spectral decomposition but cannot represent time-varying behavior.
- The STFT solves the time-localization problem but introduces an inescapable resolution tradeoff (Heisenberg).
- The WVD bypasses Heisenberg for single components but introduces cross-terms that scale quadratically.
- The SPWVD controls cross-terms through two independent smoothing windows, recovering usability at the cost of some resolution.
- This is not a progression toward perfection; it is a progression toward **controlled tradeoffs** suited to the analysis of real, multi-component signals like EEG.

---

<!-- Part 2: Volume B Highlights -->

## Lab 1 Purpose: What Happens When a Frequency Falls Between DFT Bins?

- The DFT produces exact magnitudes only when a signal frequency coincides with a bin center (the on-bin condition, $f = k \cdot f_s / N$).
- Lab 1 asks: what happens when a frequency is off-bin, and can zero-padding fix the resulting distortion?
- A sum of sinusoids at 10.0 Hz (on-bin) and 10.5 Hz (off-bin) isolates this variable at $f_s = 250$ Hz, $N = 250$.
- The theoretical prediction is that off-bin frequencies leak energy across all bins following a sinc-like envelope, and that zero-padding interpolates between bins without adding new information.

---

## Lab 1 Result: Off-Bin Leakage and the Limits of Zero-Padding

![w:900](../../results/graphs/lab1/figure_B_02.png)
**Figure B.2** - A 10.5 Hz tone (maximally off-bin at $\Delta f = 1$ Hz) leaks energy across all DFT bins.

- The off-bin tone spreads energy across the entire spectrum, while the on-bin tone concentrates perfectly in a single bin.
- Zero-padding interpolates between existing bins and produces a smoother curve, but no new spectral information is created.
- Two tones separated by 0.5 Hz (below $\Delta f = 1.0$ Hz) cannot be resolved regardless of the zero-padding factor applied.
- True resolution requires acquiring more signal (increasing $N$), not appending zeros.

---

## Lab 1 Result: Zero-Padding Failure Demonstrated

![w:900](../../results/graphs/lab1/figure_B_08.png)
**Figure B.8** - Zero-padding applied to 10 + 10.5 Hz: the spectrum appears resolved but the apparent peaks are interpolation artifacts.

- The smoother curve produced by zero-padding can create a visual illusion of two distinct peaks where only spectral leakage exists.
- This confirms the fundamental distinction between spectral interpolation (cosmetic) and spectral resolution (physical).
- Resolution is determined solely by the observation duration $T = N/f_s$, yielding $\Delta f_{\min} = f_s/N$.

---

## Lab 2 Purpose: How Do We Distinguish a Real Peak from Noise?

- Given a noisy spectrum, a spectral peak might represent a real signal or merely a noise fluctuation.
- Lab 2 asks: what is the statistical distribution of noise power in DFT bins, and how does Welch averaging improve reliability?
- The theoretical prediction is that bin power of white Gaussian noise follows an exponential distribution with CV = 1.0, meaning single-trial spectral estimates have enormous variability.
- Welch averaging should reduce variance by approximately $1/\sqrt{K}$ at the cost of frequency resolution.

---

## Lab 2 Result: The Exponential Distribution Confirmed

![w:900](../../results/graphs/lab2/figure_B_01.png)
**Figure B.1** - Histogram of bin powers for 10,000 white Gaussian noise realizations: exponential distribution confirmed (CV = 1.0).

- The coefficient of variation equals 1.0 exactly, confirming the exponential distribution predicted by theory.
- A single-trial spectral estimate has a standard deviation equal to the mean, making individual peaks unreliable.
- A peak at 2x the mean noise level has a 37% probability of being a noise fluctuation ($P_{fa} = e^{-2}$).
- The detection threshold $\gamma = -\ln(P_{fa})$ provides a principled criterion for distinguishing signal from noise.

---

## Lab 2 Result: Welch Averaging Reduces Variance

![w:900](../../results/graphs/lab2/figure_B_04.png)
**Figure B.4** - Welch periodogram comparison: increasing $K$ segments reduces spectral variance at the cost of frequency resolution.

- With $K = 1$ segment, the spectrum is noisy and unreliable (CV = 1.0).
- With $K = 8$ segments, variance drops by approximately $1/\sqrt{8} \approx 2.8$ times, making peaks statistically meaningful.
- The cost is that each segment is $1/8$ the original length, so frequency resolution degrades by a factor of 8.
- For detection problems, the variance reduction is more valuable than the lost resolution.

---

## Lab 3 Purpose: Why Do Different Windows Have Different Side-Lobe Levels?

- All cosine-sum windows (Hann, Hamming, Blackman) are constructed from the same Dirichlet kernel, yet they exhibit dramatically different side-lobe levels and decay rates.
- Lab 3 asks: what is the mathematical mechanism that explains this cancellation, and can it be derived from first principles?
- The theoretical prediction is that every cosine-sum window shares a common numerator (the Dirichlet kernel) and achieves side-lobe suppression through shifted-kernel cancellation.
- The pure sine form (Equation B.20) makes this cancellation mechanism explicit.

---

## Lab 3 Result: Kernel Anatomy and Cancellation Mechanism

![w:900](../../results/graphs/lab3/figure_B_11.png)
**Figure B.11** - The Dirichlet kernel anatomy: main lobe, side lobes, null spacing ($2\pi/M$), and $-20$ dB/octave decay rate.

- The rectangular window's side lobes decay at only $-20$ dB/octave due to the discontinuous signal edges.
- Adding shifted copies of the Dirichlet kernel cancels the side lobes at the expense of widening the main lobe.
- The pure sine form (Equation B.20) reveals that the cancellation depends on the exact placement of kernel nulls relative to adjacent side-lobe peaks.
- This first-principles derivation explains all cosine-sum windows within a single unified framework.

---

## Lab 3 Result: Window Comparison

![w:900](../../results/graphs/lab3/figure_B_14.png)
**Figure B.14** - Four windows applied to the same signal: the resolution-leakage tradeoff visualized.

- First side-lobe levels: $-13$ dB (Rectangular), $-32$ dB (Hann), $-43$ dB (Hamming), $-58$ dB (Blackman).
- Main-lobe widths in bins: 2 (Rectangular), 4 (Hann), 4 (Hamming), 6 (Blackman).
- The Hann window offers the best compromise for general time-frequency analysis: moderate main-lobe width with $-60$ dB/octave decay.
- For EEG analysis (broad spectral features, wide dynamic range), the Hann window is the standard choice.

---

## Lab 4 Purpose: Can We See Frequency Change over Time?

- The DFT provides a global spectrum but discards all temporal information; it cannot detect non-stationarity.
- Lab 4 asks: can the STFT resolve time-frequency structure, and what fundamental limit constrains its precision?
- A linear chirp (1 to 80 Hz over 1200 seconds) serves as the ideal test signal because any blurring is immediately visible as deviation from the known thin diagonal.
- The theoretical prediction is that Heisenberg uncertainty ($\Delta t \cdot \Delta f \geq \beta$) forces an inescapable tradeoff between time and frequency resolution.

---

## Lab 4 Result: The Heisenberg Tradeoff Made Visible

![w:900](../../results/graphs/lab4/figure_B_18.png)
**Figure B.18** - Same chirp, two window lengths: short ($M = 250$, 1.0 s) versus long ($M = 1250$, 5.0 s).

- The short window tracks the chirp in time with $\Delta t = 1.0$ s but frequency is blurred ($\Delta f = 2.0$ Hz).
- The long window resolves frequency with $\Delta f = 0.4$ Hz but time is blurred ($\Delta t = 5.0$ s).
- The product $\Delta t \cdot \Delta f = 2.0$ in both cases, confirming $\beta = 2$ for the Hann window.
- No single window length can achieve both sharp time and sharp frequency localization simultaneously; this is the Heisenberg bound made visible.

---

## Lab 5 Purpose: At What Frequency Separation Do Two Tones Merge?

- Lab 3 derived the main-lobe width of each window; Lab 4 demonstrated the Heisenberg tradeoff without testing the exact resolution limit.
- Lab 5 asks: at what frequency separation do two stationary tones become indistinguishable on the spectrogram?
- The theoretical prediction is $\Delta f_{\min} = \beta \cdot f_s / M$, derived directly from the window main-lobe width established in Lab 3.
- The test uses two tones at varying separations, with the Hann window ($\beta = 2$, so $\Delta f_{\min} = 2f_s/M$).

---

## Lab 5 Result: Resolved versus Merged

![w:900](../../results/graphs/lab5/figure_B_28.png)
**Figure B.28** - Two stationary tones on the spectrogram: resolved (top, $\Delta f > 2f_s/M$) versus merged (bottom, $\Delta f < 2f_s/M$).

- For $M = 1250$ at $f_s = 250$ Hz: $\Delta f_{\min} = 2 \times 250 / 1250 = 0.4$ Hz.
- Tones at 30.0 Hz and 30.5 Hz ($\Delta f = 0.5$ Hz $> 0.4$ Hz) are resolved as two distinct horizontal lines.
- Tones at 30.0 Hz and 30.3 Hz ($\Delta f = 0.3$ Hz $< 0.4$ Hz) merge into a single blob.
- The transition is sharp and quantitative, confirming the prediction from Lab 3's window theory.

---

## Lab 6 Purpose: Can We Detect Periodicity and Shared Components without Phase?

- The autocorrelation discards phase information by design, retaining only the periodicity structure of a signal.
- Lab 6 asks: what information is preserved (periodicity) and what is lost (phase), and can cross-correlation detect shared frequency content between two independent signals?
- The theoretical prediction is that signals differing only in phase produce identical autocorrelation functions, and that cross-correlation measures shared periodic content as a function of lag.

---

## Lab 6 Result: Phase-Blindness Demonstrated

![w:900](../../results/graphs/lab6/figure_B_33.png)
**Figure B.33** - Two signals with different phase offsets produce identical autocorrelation functions.

- The autocorrelation $R_x[m] = \sum_n x[n] \cdot x[n+m]$ depends only on periodicity, not on the absolute phase of the components.
- This is confirmed by the Wiener-Khinchin theorem: $|X[k]|^2 = \text{DFT}\{R_x[m]\}$, showing that the power spectrum (magnitude squared) contains no phase.
- Phase-blindness is a feature, not a limitation: it extracts what is repeatable and discards what is unique.

---

## Lab 6 Result: Cross-Correlation Detects Shared Content

![w:900](../../results/graphs/lab6/figure_B_37.png)
**Figure B.37** - Cross-correlation between two signals: $\rho = 0.107$ for shared tone, $\rho \approx 0$ for independent signals.

- The cross-correlation $R_{xy}[m]$ quantifies the similarity between two signals as a function of lag.
- When two signals share a common frequency component, the normalized cross-correlation is measurably nonzero ($\rho = 0.107$).
- When signals are fully independent, $\rho \approx 0$, providing a clear criterion for channel independence.
- In Volume C, this technique validates that auxiliary channels (ECG, EMG) do not contaminate the EEG signal.

---

## Lab 7 Purpose: Can the WVD Bypass the Heisenberg Limit?

- The STFT is fundamentally limited by $\Delta t \cdot \Delta f \geq \beta$, meaning no window can achieve simultaneous sharpness in both time and frequency.
- Lab 7 asks: can the WVD produce sub-Heisenberg localization, and what is the cost for multi-component signals?
- The same linear chirp from Lab 4 serves as the test signal, enabling a direct comparison between the STFT and WVD under identical conditions.
- The theoretical prediction is that the WVD achieves perfect localization for a single component but generates cross-terms at the midpoint between any two components.

---

## Lab 7 Result: Razor-Sharp Chirp versus Heisenberg-Blurred STFT

![w:900](../../results/graphs/lab7/figure_B_40.png)
**Figure B.40** - STFT (left) versus WVD (right) on a single linear chirp: the WVD traces the instantaneous frequency exactly.

- The STFT produces a thick, blurred diagonal whose width reflects the Heisenberg tradeoff.
- The WVD produces a razor-sharp line with no time-frequency spreading whatsoever.
- This demonstrates that the Heisenberg limit is a property of the STFT framework, not of time-frequency analysis in general.
- The analytic signal (Hilbert transform) eliminates the DC self-interference term, producing a clean single-component representation.

---

## Lab 7 Result: Cross-Terms Are as Energetic as the Signal

![w:900](../../results/graphs/lab7/figure_B_42.png)
**Figure B.42** - Chirp + constant tone: the WVD generates an oscillating interference pattern at the midpoint between the two components.

- For two components, the bilinearity generates a cross-term located at the midpoint in both time and frequency.
- The cross-term energy equals the auto-term energy; it is not a small perturbation but a fundamental corruption.
- The analytic signal removes DC ghosts but cannot suppress multi-component cross-terms.
- For $N$ components, there are $N(N-1)/2$ cross-terms versus only $N$ auto-terms, making the WVD unusable for real-world signals.

---

## Lab 8 Purpose: Can Cross-Terms Be Suppressed without Losing Sharpness?

- The WVD is unusable for multi-component signals due to cross-terms that are as energetic as the signal itself.
- Lab 8 asks: can the SPWVD suppress cross-terms while preserving the WVD's superior localization, and how do its two independent windows ($h$ for frequency, $g$ for time) achieve this?
- The theoretical prediction is that the PWVD alone (lag window only) fails on time-oscillating ghosts, and the full SPWVD (both windows) is required to suppress both ghost types.
- The duality principle predicts that the ghost type depends on how components are separated in the time-frequency plane.

---

## Lab 8 Result: WVD to PWVD to SPWVD Progression

![w:900](../../results/graphs/lab8/figure_B_46.png)
**Figure B.46** - Step-by-step ghost suppression: WVD (top), PWVD (middle), SPWVD (bottom).

- The WVD shows razor-sharp auto-terms but cross-terms fill the space between the two chirps.
- The PWVD (lag window $h$ only) suppresses frequency-oscillating ghosts but time-oscillating ghosts survive.
- The full SPWVD (both $h$ and $g$) suppresses both ghost types, recovering a readable time-frequency representation.
- This confirms that both windows are independently necessary for a multi-component signal.

---

## Lab 8 Result: Two-Knob Sweep

![w:900](../../results/graphs/lab8/figure_B_47.png)
**Figure B.47** - SPWVD with different window configurations: strong frequency smoothing (left) versus strong time smoothing (right).

- Case 1 ($h = 101$, $g = 5$): strong frequency smoothing preserves time resolution but leaves time-oscillating ghosts visible.
- Case 2 ($h = 25$, $g = 31$): strong time smoothing suppresses all ghosts but blurs the time axis.
- The two knobs are genuinely independent: adjusting one does not affect the other axis.
- This two-dimensional control is the SPWVD's primary advantage over the STFT's single-window design.

---

## Lab 8 Result: The Duality of Ghost Types

![w:900](../../results/graphs/lab8/figure_B_49.png)
**Figure B.49** - Two time-separated impulses: frequency-oscillating ghosts suppressed by the PWVD lag window.

- Components separated in time produce cross-terms that oscillate in the frequency direction; the lag window $h$ suppresses these.
- Components separated in frequency produce cross-terms that oscillate in the time direction; the time window $g$ suppresses these.
- This duality explains why the PWVD alone is insufficient: it addresses only one of the two ghost types.
- The SPWVD is required to suppress both ghost types simultaneously, regardless of how components are arranged.

---

## Volume B Summary: From Theory to Measurement

- **Labs 1-3** establish the DFT framework: bins, leakage, noise statistics, and the unified window theory.
- **Labs 4-5** validate the STFT: Heisenberg uncertainty is measured and the resolution limit is confirmed quantitatively.
- **Lab 6** derives autocorrelation properties (phase-blindness, cross-correlation) that underpin the WVD definition.
- **Labs 7-8** demonstrate the WVD/SPWVD: sub-Heisenberg sharpness, the cross-term problem, and the two-knob solution.
- Every prediction from Volume A is confirmed quantitatively in Volume B, establishing confidence for Volume C.

---

<!-- Part 3 -->

## Volume C: From Model Signals to Real EEG

- Volume B validated every tool on controlled signals where the ground truth is known.
- Volume C applies the same tools to a real neonatal EEG recording where the ground truth is unknown.
- The analysis is **adaptive-directed**: C.1 examines the data first, then each subsequent section follows from what was discovered.
- No filtering is applied. Only tools derived and validated in Volumes A and B are used.
- The goal is to characterize the signal's time-frequency structure, not to make clinical diagnoses.

---

## The Dataset: Overview

- **Subject:** NORB00055 (neonatal), European Data Format (EDF), from the Zenodo neonatal EEG repository.
- **Sampling rate:** 200 Hz. **Duration:** 1140 seconds (19 minutes). **Total channels:** 24.
- **EEG channels:** 19 electrodes in the standard 10-20 system (Fp1, Fp2, F3, F4, C3, C4, P3, P4, O1, O2, F7, F8, T3, T4, T5, T6, FZ, CZ, PZ).
- **Auxiliary channels:** 5 (ECG, EMG, EOG, and reference channels).
- Data loaded via MNE-Python. All amplitudes reported in physical units ($\mu$V).

---

## The Dataset: Channel Selection and Rationale

- **Primary channel: CZ** (vertex electrode), located at the top center of the scalp.
- CZ captures activity from both hemispheres with minimal regional bias, making it the least biased single-channel starting point.
- The 10-20 system provides standardized electrode placement that enables cross-study comparison.
- Auxiliary channels (ECG, EMG, EOG) are used exclusively for artifact verification, never for brain signal analysis.
- **No filtering was applied.** Only tools derived in Volumes A and B are used for analysis.

---

## C.1 Triage: Time-Domain Inspection

![w:900](../../results/graphs/volume_c/c1/figure_C_02.png)
**Figure C.2** - CZ time-domain signal (full 19 minutes): visible burst-suppression pattern.

- The raw signal shows alternating high-amplitude bursts and low-amplitude quiet periods.
- Peak amplitudes reach approximately $\pm 200$ $\mu$V during bursts.
- The burst pattern is irregular: burst durations and inter-burst intervals vary throughout the recording.
- This visual inspection immediately suggests non-stationarity, motivating time-frequency analysis.

---

## C.1 Triage: Band Power Distribution

![w:900](../../results/graphs/volume_c/c1/figure_C_06.png)
**Figure C.6** - Band power distribution computed via Welch's method: 91.8% delta.

- Delta (0.5-4 Hz): **91.8%** of total power dominates the recording.
- Theta (4-8 Hz): 5.7%. Alpha (8-13 Hz): 1.0%. Beta (13-30 Hz): 1.0%.
- The absence of alpha and beta rhythms is consistent with neonatal EEG physiology.
- Triage conclusion: the signal is delta-dominated, bursty, and non-stationary.
- Analysis direction: investigate whether the delta activity is continuous or episodic.

---

## C.1 Triage: Analysis Direction

- The triage findings dictate the direction for all subsequent sections.
- **C.2** investigates whether the delta dominance is a true oscillation or simply the low-frequency tail of 1/f noise.
- **C.3** uses the STFT to determine whether delta power is continuous or arrives in bursts.
- **C.4** validates channel cleanliness before applying the sensitive WVD/SPWVD analysis.
- **C.5** applies the full WVD/SPWVD toolkit to a carefully selected burst segment.

---

## C.2: The 1/f Slope of the Power Spectrum

![w:900](../../results/graphs/volume_c/c2/figure_C_09.png)
**Figure C.9** - Log-log PSD with linear fit: slope = $-3.18$ from 5 to 40 Hz.

- The PSD follows $1/f^{3.18}$ across the 5-40 Hz range, steeper than pink noise ($1/f^1$).
- A steep slope indicates that high-frequency power drops off rapidly, consistent with neonatal EEG.
- The delta band (0.5-4 Hz) contains peaks that sit **below** the 1/f extrapolation line.
- This means the delta peaks represent genuine oscillatory activity, not merely the low-frequency tail of 1/f noise.

---

## C.2: Delta Zoom - Continuous or Bursty?

![w:900](../../results/graphs/volume_c/c2/figure_C_10.png)
**Figure C.10** - Zoomed PSD in the 0-5 Hz range: delta peak structure.

- The DFT reveals that delta power is concentrated around 0.4-0.6 Hz, not spread uniformly across the band.
- However, the stationary DFT cannot determine whether this delta activity is continuous or intermittent.
- A continuous delta oscillation and a series of delta bursts can produce identical DFT spectra.
- Resolving this ambiguity requires the time axis, which motivates the STFT analysis in C.3.

---

## C.3: The Full-Recording Spectrogram

![w:900](../../results/graphs/volume_c/c3/figure_C_12.png)
**Figure C.12** - Full-recording spectrogram (Hann window, $M = 1000$ samples = 5.0 s): delta bursts as vertical stripes.

- Delta power is **not continuous**. It arrives in discrete bursts visible as bright vertical stripes.
- The bursts span multiple frequency bands simultaneously (broadband events, not narrowband oscillations).
- Quiet periods between bursts show near-zero spectral power across all frequencies.
- This answers the C.2 question definitively: **discontinuous delta activity**, not a continuous oscillation.

---

## C.3: Burst Detection Methodology

![w:900](../../results/graphs/volume_c/c3/figure_C_16.png)
**Figure C.16** - Delta power time series overlaid on the CZ signal, with burst markers (threshold: $2 \times$ median).

- Burst threshold: $2 \times$ median delta power = 3479 $\mu$V$^2$.
- **19%** of the recording is classified as burst; **81%** is classified as quiet (inter-burst).
- The maximum-to-median ratio is 17x, indicating highly energetic burst events.
- Delta-theta correlation: $\rho = 0.33$, showing that theta and delta are partially but not fully coupled.

---

## C.3: Burst Statistics

- Median burst duration: approximately 3-5 seconds (consistent with neonatal burst-suppression literature).
- Inter-burst intervals vary from 5 to 30 seconds, reflecting the irregular nature of neonatal discontinuous activity.
- The burst pattern is whole-brain synchronous: all 19 EEG channels show simultaneous activation during bursts.
- This synchrony is expected for neonatal EEG, where cortical activity is less lateralized than in adults.
- The STFT successfully characterizes the temporal structure that the DFT could not resolve.

---

## C.3: What the STFT Answered

- The DFT (C.2) established that delta peaks exist but could not determine their temporal pattern.
- The STFT resolves this ambiguity: delta activity is **discontinuous**, arriving in high-energy bursts separated by quiet periods.
- The burst-suppression pattern (19% burst, 81% quiet) is the dominant non-stationary feature.
- This pattern is consistent with the neonatal EEG literature on discontinuous normal activity (Lamblin et al., 1999).
- The next question: is the CZ channel clean enough for the sensitive WVD/SPWVD analysis?

---

## C.4: Cross-Correlation with Auxiliary Channels

![w:900](../../results/graphs/volume_c/c4/figure_C_20.png)
**Figure C.20** - Cross-correlation between CZ and all auxiliary channels.

- All auxiliary channels show $\rho < 0.03$ versus CZ, indicating negligible linear coupling.
- ECG, EMG, and EOG signals are statistically independent of the CZ brain signal.
- This confirms that the bursts observed in CZ are of cerebral origin, not cardiac or muscular artifacts.
- EEG inter-channel correlations range from $\rho = 0.47$ to $0.78$, reflecting shared brain activity (expected for neonatal recordings).

---

## C.4: Noise Floor Characterization

- The alpha band (8-13 Hz) serves as the noise-floor test band because it contains no physiological signal in neonatal EEG.
- Measured CV = 1.11 on the alpha-band power distribution (theoretical exponential: CV = 1.00).
- The 11% deviation indicates a close but imperfect fit to the exponential noise model.
- Possible causes: residual non-stationarity within each Welch segment, or minor spectral coloring.
- **Verdict:** CZ is sufficiently clean for WVD analysis. No auxiliary contamination was detected, and artifacts were not removed; instead, C.5 selects a clean segment.

---

## C.5: Segment Selection - The Rejected Burst

![w:900](../../results/graphs/volume_c/c5/figure_C_22.png)
**Figure C.22** - REJECTED: the strongest burst (8.7x median, $t = 842.5$ s) exhibits amplifier saturation.

- The strongest burst in the recording has a delta power of 15,597 $\mu$V$^2$ (8.7x median).
- Inspection reveals **44 consecutive flat samples** in the time-domain waveform (amplifier clipping).
- Amplifier saturation is a hardware artifact, not a physiological event.
- This burst must be rejected for WVD analysis because clipping introduces artificial harmonics.

---

## C.5: Segment Selection - The Accepted Burst

![w:900](../../results/graphs/volume_c/c5/figure_C_23.png)
**Figure C.23** - ACCEPTED: the 75th percentile burst ($t = 65.0$ s, 5435 $\mu$V$^2$, 3.0x median) is clean.

- The 75th percentile burst is energetic enough to be representative but not extreme.
- No flat samples, no clipping, no obvious artifacts in the time-domain waveform.
- Duration: approximately 2 seconds, containing 2-3 visible delta cycles.
- This is the segment selected for WVD and SPWVD analysis.
- The transparent reporting of both the rejection and acceptance criteria ensures reproducibility.

---

## C.5: From Lab Validation to Real Data

- Labs 7 and 8 demonstrated WVD cross-terms on controlled model signals with known components.
- The real EEG burst contains approximately 2-3 delta cycles plus broadband background activity.
- The question is whether the SPWVD can produce a readable time-frequency representation of a real burst.
- The parameters from Lab 8 (lag window $h$ and time window $g$) provide the starting point for tuning.
- Success criterion: the linear-scale SPWVD must show the delta burst without being dominated by cross-terms.

---

## C.5: Raw WVD on Real EEG

![w:900](../../results/graphs/volume_c/c5/figure_C_25.png)
**Figure C.25** - Raw WVD of the clean delta burst: cross-term contamination dominates.

- **49% of all WVD values are negative.** The WVD is not a valid power distribution for this signal.
- Oscillating cross-terms fill the entire time-frequency plane, obscuring the auto-terms.
- Even this short (2 s), relatively simple segment (approximately 3 components) generates severe contamination.
- This confirms Section A.7.3 and Lab 7: the raw WVD is unusable on real EEG without smoothing.

---

## C.5: The SPWVD Applied to the Delta Burst

![w:900](../../results/graphs/volume_c/c5/figure_C_26.png)
**Figure C.26** - SPWVD of the same burst segment: cross-terms suppressed, delta structure visible.

- The lag window ($h$) and time window ($g$) suppress the oscillating cross-terms.
- The delta burst is now visible as a concentrated energy blob in the 0.5-3 Hz range.
- Time localization is sharper than the STFT: the burst onset and offset are more precisely delineated.
- Some residual artifacts remain visible in the dB representation (circular patterns from incompletely suppressed cross-terms).

---

## C.5: Three-Way Comparison

![w:900](../../results/graphs/volume_c/c5/figure_C_27.png)
**Figure C.27** - STFT versus WVD versus SPWVD on the same 2-second burst segment.

- **STFT** (top): readable but blurred; the burst occupies a wide time-frequency region due to Heisenberg spreading.
- **WVD** (middle): theoretically sharp but practically unusable due to cross-term contamination.
- **SPWVD** (bottom, linear): the **sharpest readable representation**, combining suppressed cross-terms with sub-Heisenberg localization.
- The SPWVD dB panel still shows residual circular artifacts, indicating incomplete cross-term removal.

---

## C.5: Why Not Filter before the WVD?

- A bandpass filter (0.5-4 Hz) would isolate the delta component before WVD computation, dramatically reducing cross-terms.
- Fewer input components means fewer cross-terms: $N(N-1)/2$ drops rapidly as $N$ decreases.
- **However:** Volumes A and B did not derive filter theory (no FIR design, no passband/stopband specification, no group delay analysis).
- Applying a tool without first deriving and validating it violates the report's foundational methodology.
- **Segment selection** is the only preprocessing available within our derived toolkit, and it produces a usable (though imperfect) result.

---

## C.5: The Scope Boundary

- The natural next step would be: derive FIR filter design using Lab 3's window functions, validate the filter on model signals, then apply bandpass filtering before the WVD on real EEG.
- This would form a complete pipeline: bandpass $\rightarrow$ analytic signal $\rightarrow$ SPWVD $\rightarrow$ clean time-frequency representation.
- This pipeline is well-established in the EEG literature but falls outside the scope of this report.
- The report demonstrates the tools as derived, including their limitations on real data.
- Future work: extend the theoretical foundation to include filter design, completing the signal processing chain.

---

<!-- Part 4 -->

## What Worked

| Tool | Applied To | Key Finding |
| ---- | ---------- | ----------- |
| DFT + Welch | Full recording | Delta dominance (91.8%), 1/f slope = $-3.18$ |
| STFT | Full recording | Burst pattern: 19% burst, 81% quiet |
| Cross-correlation | Auxiliary vs CZ | CZ is clean ($\rho < 0.03$) |
| SPWVD | 2 s burst segment | Sharper burst localization than the STFT |

- Every tool from Volumes A and B found meaningful application on real neonatal EEG.
- Each tool addressed a specific limitation of the previous one, exactly as the theory predicted.

---

## What Did Not Work

- The **raw WVD** was completely unusable on real EEG: 49% negative values and cross-terms obscuring all structure.
- The **SPWVD dB panel** retained circular artifacts from residual cross-terms, limiting the dB representation's utility.
- The **exponential noise model** was approximate on real data (CV = 1.11 versus theoretical 1.00), reducing confidence in detection thresholds.
- The **strongest burst** turned out to be amplifier saturation, not physiology, requiring data-driven segment rejection.
- These are inherent limitations of the methods as applied without filtering, not implementation errors.

---

## What We Learned

- Time-frequency analysis is not a single tool but a **progression** of tools, each with explicit tradeoffs.
- The Heisenberg limit is real and visible on spectrograms; it is not merely a theoretical abstraction.
- Cross-terms are the fundamental cost of bilinear methods (WVD); they cannot be ignored or hand-waved away.
- Real data always contains surprises (amplifier saturation, imperfect noise models) that theory alone does not predict.
- Transparent methodology means reporting failures alongside successes.

---

## The Scope Boundary: Filtering

- The report's methodology requires that every tool be derived theoretically, validated on model signals, then applied to real data.
- Digital filter design (FIR/IIR) was not included in the theoretical framework (Volumes A and B).
- Therefore, bandpass filtering cannot be applied to the EEG signal, even though it would improve WVD/SPWVD results.
- This boundary is a deliberate methodological choice, not an oversight.
- The report demonstrates what time-frequency analysis can achieve **without** filtering, establishing a clear baseline for future work.

---

## Closing Claim

The progression **DFT $\rightarrow$ STFT $\rightarrow$ SPWVD** provides increasingly detailed views of the same neonatal EEG signal, each tool resolving an ambiguity that the previous tool could not address.

- The SPWVD achieves **sub-Heisenberg resolution** on selected clean segments, confirming its theoretical advantage.
- Its practical value lies in the **linear-scale representation**; the dB representation remains partially compromised by residual cross-terms.
- The neonatal EEG is consistent with **normal discontinuous neonatal activity**: delta-dominated, bursty, and whole-brain synchronous (Lamblin et al., 1999; Andre et al., 2010).

**No clinical diagnosis is made or implied.** These are signal-processing observations on a single recording.

---

## References

- Lamblin, M. D. et al. (1999). Electroencephalography of the premature and full-term newborn. *Neurophysiologie Clinique*, 29(3), 167-219.
- Andre, M. et al. (2010). Electroencephalography in premature and full-term infants. *Neurophysiologie Clinique*, 40(2), 59-124.
- Cohen, L. (1995). *Time-Frequency Analysis*. Prentice Hall.
- Boashash, B. (2016). *Time-Frequency Signal Analysis and Processing* (2nd ed.). Academic Press.
- Flandrin, P. (1998). *Time-Frequency/Time-Scale Analysis*. Academic Press.
