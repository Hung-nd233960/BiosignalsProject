---
marp: true
theme: default
paginate: true
math: mathjax
---

# From the DFT to the SPWVD
## Time-Frequency Analysis Applied to Neonatal EEG

**Authors:** Nguyen Duc Hung - 20233960, Bui Phuong Duy - 23233957, Tran Viet Bach - 23233954

Digital Biosignal Processing - Final Report

---

## Three-volume structure

- **Volume A** - Theory: DFT, windowing, STFT, autocorrelation, WVD, SPWVD
- **Volume B** - 8 labs validating each theory section on model signals
- **Volume C** - Application to real neonatal EEG (19 min, 24 channels)

**The arc:** each tool addresses a specific limitation of the previous one.

DFT *(no time)* $\rightarrow$ STFT *(Heisenberg-limited)* $\rightarrow$ WVD *(cross-terms)* $\rightarrow$ SPWVD *(controllable)*

---

<!-- Part 1: Volume A Highlights -->

## The DFT and its frequency bins

$$X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j2\pi kn/N} \qquad \text{(A.5)}$$

- Bin spacing: $\Delta f = f_s / N$ - fixed by the signal length
- **On-bin** ($f$ is a multiple of $\Delta f$): energy captured perfectly in one bin
- **Off-bin**: energy leaks across all bins (spectral leakage)
- Zero-padding interpolates the spectrum but **does not improve resolution**

---

## Windowing and the Dirichlet kernel

![w:900](../results/graphs/lab3/figure_B_14.png)
**Figure B.14** - Window spectra: rectangular, Hann, Hamming, Blackman

- Rectangular window $\rightarrow$ Dirichlet kernel; main lobe width $= 4\pi/M$
- Cosine-sum windows trade wider main lobe for suppressed side lobes
- First side-lobe: $-13$ dB (rect), $-32$ dB (Hann), $-43$ dB (Hamming), $-58$ dB (Blackman)

---

## Spectral statistics: when is a peak real?

$$P_{fa} = e^{-\gamma} \qquad \text{(A.37)}$$

- Bin power of white noise follows an **exponential distribution** (CV = 1.0)
- A peak is real only if it exceeds the noise floor by a threshold set by $P_{fa}$
- **Welch averaging** reduces variance ($\sigma \propto 1/\sqrt{K}$) at the cost of frequency resolution
- Not every spectral peak is a signal - statistics provide the decision rule

---

## The STFT and Heisenberg uncertainty

$$\Delta t \cdot \Delta f \geq \beta \qquad \text{(A.49)}$$

- The STFT adds a **time axis** by windowing the signal before the DFT
- Short window ($M$ small): good time resolution, poor frequency resolution
- Long window ($M$ large): good frequency, poor time
- **Cannot have both simultaneously** - the window length is the single knob
- For Hann: $\beta = 2$, so $\Delta t \cdot \Delta f \geq 2$

---

## The WVD: perfect resolution, fatal flaw

$$W_x[n, k] = \sum_{m} R_x[n, m] \, e^{-j2\pi km/N_f} \qquad \text{(A.61)}$$

where $R_x[n, m] = z[n+m] \cdot z^*[n-m]$ is the instantaneous autocorrelation.

- Single chirp $\rightarrow$ **razor-sharp diagonal** (bypasses Heisenberg)
- Multi-component $\rightarrow$ **cross-terms** at the midpoint between every pair
- Cross-terms oscillate and can be **as energetic as the real components**
- Analytic signal (Hilbert transform) removes the DC self-ghost

---

## The SPWVD: two independent smoothing knobs

$$\text{SPWVD}_x[n, k] = \sum_m h[m] \left(\sum_p g[p] \, z[n+p+m] \, z^*[n+p-m]\right) e^{-j2\pi km/N_f} \qquad \text{(A.72)}$$

- **Lag window** $h[m]$: smooths frequency axis $\rightarrow$ kills frequency-oscillating ghosts
- **Time window** $g[p]$: smooths time axis $\rightarrow$ kills time-oscillating ghosts
- Each knob trades resolution for ghost suppression **independently**
- Unlike the STFT, time and frequency resolution are decoupled

---

## Tool progression summary

| Tool | Adds | Costs | Limitation |
|------|------|-------|------------|
| DFT | Global spectrum | No time info | Cannot see non-stationarity |
| STFT | Time-frequency | Heisenberg limit | $\Delta t \cdot \Delta f \geq \beta$ |
| WVD | Sub-Heisenberg | Cross-terms | Unusable for $\geq 2$ components |
| SPWVD | Controllable smoothing | Some resolution loss | Residual ghosts on complex signals |

**Each tool fixes a specific limitation of the previous one.** This is the spine of the report.

---

<!-- Part 2: Volume B Highlights -->

## Lab 1: DFT bins in action

![w:900](../results/graphs/lab1/figure_B_08.png)
**Figure B.8** - Zero-padding 10 + 10.5 Hz: appears resolved but is not (below $\Delta f_{\min}$)

- 10.5 Hz is maximally off-bin $\rightarrow$ energy smeared across all bins
- Zero-padding adds spectral samples but **cannot create resolution**
- True resolution requires a longer signal, not more zeros

---

## Lab 3: window spectra from first principles

![w:900](../results/graphs/lab3/figure_B_11.png)
**Figure B.11** - Dirichlet kernel anatomy: main lobe, side lobes, null spacing

- All cosine-sum windows share the Dirichlet kernel as their numerator
- The pure sine form (Equation B.20) derives the full spectrum from the kernel
- Decay rate: $-20$ dB/octave (rect, Hamming) vs $-60$ dB/octave (Hann, Blackman)

---

## Lab 4: Heisenberg is visible on the spectrogram

![w:900](../results/graphs/lab4/figure_B_18.png)
**Figure B.18** - Chirp spectrogram: short window (left) vs long window (right)

- Short window tracks the chirp in time but **frequency is blurred**
- Long window resolves frequency but **smears the chirp trajectory**
- The Heisenberg tradeoff is not abstract - it is visible

---

## Lab 5: two-tone resolution confirmed

![w:900](../results/graphs/lab5/figure_B_28.png)
**Figure B.28** - Two tones on spectrogram: resolved (top) vs merged (bottom)

- Hann main-lobe width determines the minimum resolvable frequency separation
- Tones closer than $2f_s/M$ merge into one blob
- Lab 3 predicts the limit; Lab 5 confirms it on the spectrogram

---

## Lab 6: autocorrelation and phase-blindness

![w:900](../results/graphs/lab6/figure_B_33.png)
**Figure B.33** - Two signals with different phases produce identical autocorrelations

- Autocorrelation reveals **periodicity** but discards **phase**
- Wiener-Khinchin theorem: DFT of autocorrelation = power spectrum (Equation A.56)
- The instantaneous autocorrelation is the WVD's building block (A.7)

---

## Lab 7: WVD sharpness on a single chirp

![w:900](../results/graphs/lab7/figure_B_40.png)
**Figure B.40** - STFT (left) vs WVD (right) on a single chirp

- STFT: thick, blurred diagonal (Heisenberg-limited)
- WVD: **razor-sharp line** tracking the instantaneous frequency exactly
- For a single component, the WVD completely bypasses Heisenberg

---

## Lab 7: cross-terms ruin multi-component signals

![w:900](../results/graphs/lab7/figure_B_42.png)
**Figure B.42** - Chirp + tone: clean STFT (left) vs corrupted WVD (right)

- STFT: clean superposition of chirp and tone
- WVD: sharp components **plus oscillating ghost** at the midpoint frequency
- The ghost is as energetic as the real components - the WVD is unusable

---

## Lab 8: WVD to PWVD to SPWVD progression

![w:900](../results/graphs/lab8/figure_B_46.png)
**Figure B.46** - Step-by-step ghost suppression (linear left, dB right)

- **WVD** (top): sharp trajectories, heavy cross-terms
- **PWVD** (middle): lag window $h$ smooths frequency - time ghosts **survive**
- **SPWVD** (bottom): time window $g$ added - **both ghost types suppressed**

---

## Lab 8: the duality of ghost types

![w:900](../results/graphs/lab8/figure_B_49.png)
**Figure B.49** - Two impulses: frequency-oscillating ghosts, PWVD suppresses them

- Components separated in **time** $\rightarrow$ frequency-oscillating ghost $\rightarrow$ PWVD kills it
- Components separated in **frequency** $\rightarrow$ time-oscillating ghost $\rightarrow$ needs SPWVD
- The SPWVD handles both; the PWVD handles only one axis

---

## Lab 8: two-knob sweep

![w:900](../results/graphs/lab8/figure_B_47.png)
**Figure B.47** - SPWVD: strong freq smoothing (top) vs strong time smoothing (bottom)

- Case 1 (h=101, g=5): sharp frequency, time ghosts survive ($T_g < 1/\Delta f$)
- Case 2 (h=25, g=31): clean but frequency blurred
- Unlike STFT, the two axes are **independently controlled**

---

<!-- Part 3: Volume C Highlights -->

## The dataset: neonatal EEG

- **Subject:** NORB00055 (neonatal), European Data Format (EDF)
- **Sampling rate:** 200 Hz, **duration:** 1140 s (19 min), **channels:** 24 (19 EEG + 5 auxiliary)
- **Primary channel:** CZ (vertex) - least biased starting point
- **No filtering applied** - only tools derived in Volumes A and B are used
- Loaded via MNE; amplitude in $\mu$V throughout

---

## C.1 Triage: delta dominance

![w:900](../results/graphs/volume_c/c1/figure_C_06.png)
**Figure C.6** - Band power distribution: 91.8% delta

- Delta (0.5-4 Hz): **91.8%**, theta: 5.7%, alpha: 1.0%, beta: 1.0%
- No alpha or beta rhythms - consistent with neonatal EEG
- Whole-brain synchronous: all channels show the same burst pattern
- Triage decision: delta-dominated, bursty, non-stationary

---

## C.2: is delta rhythmic or 1/f?

![w:900](../results/graphs/volume_c/c2/figure_C_09.png)
**Figure C.9** - Log-log PSD with 1/f fit (slope = -3.18)

- PSD follows $1/f^{3.18}$ from 5-40 Hz - steeper than pink noise ($1/f$)
- Delta peaks at 0.4-0.6 Hz sit **below** the 1/f extrapolation
- Cannot tell from the DFT alone if delta is continuous or bursty
- Need the time axis $\rightarrow$ C.3 (STFT)

---

## C.3: STFT reveals the burst pattern

![w:900](../results/graphs/volume_c/c3/figure_C_12.png)
**Figure C.12** - Full-recording spectrogram: delta bursts visible as vertical stripes

- Delta power is **not continuous** - it comes in bursts
- 19% of the recording is burst, 81% is quiet (threshold: 2x median)
- Max/median ratio: 17x - bursts are highly energetic
- Answers C.2's open question: **discontinuous delta, not continuous oscillation**

---

## C.3: burst detection quantified

![w:900](../results/graphs/volume_c/c3/figure_C_16.png)
**Figure C.16** - Delta power overlaid on CZ time domain with burst markers

- Burst threshold: $2 \times$ median delta power = 3479 $\mu$V$^2$
- Delta-theta correlation: $\rho = 0.33$ (partially independent)
- The burst structure is the dominant non-stationary feature of this EEG

---

## C.4: is CZ clean for WVD analysis?

- **Cross-correlation (Lab 6):** all auxiliary channels $\rho < 0.03$ vs CZ
- **EEG inter-channel:** $\rho = 0.47 - 0.78$ (shared brain activity, expected)
- **Noise floor (Lab 2):** CV = 1.11 on alpha band (close to exponential 1.00)
  - 11% deviation: close but not exact - honest characterization of model limits
- **Verdict:** CZ is clean. No auxiliary contamination detected.
- Artifacts were **not removed** - C.5 selects a clean segment instead

---

## C.5: segment selection - the honest story

![w:900](../results/graphs/volume_c/c5/figure_C_22.png)
**Figure C.22** - REJECTED: strongest burst (8.7x median) shows amplifier saturation

- Strongest burst at t = 842.5 s: 15597 $\mu$V$^2$ - but **44 flat samples** (clipping)
- Fell back to 75th percentile: t = 65.0 s, 5435 $\mu$V$^2$ (3.0x median), clean
- The strongest burst in the recording is an **artifact**, not physiology
- Data-driven selection including the failure: transparent methodology

---

## C.5: raw WVD on real EEG

![w:900](../results/graphs/volume_c/c5/figure_C_25.png)
**Figure C.25** - Raw WVD of clean delta burst: cross-term contamination

- **49% of values are negative** - the WVD is not a true power distribution
- Oscillating cross-terms fill the entire time-frequency plane
- Even a 2-second segment with ~3 components generates severe contamination
- Confirms A.7.3 and Lab 7: the raw WVD is unusable on real EEG

---

## C.5: SPWVD - the payoff and its limit

![w:900](../results/graphs/volume_c/c5/figure_C_27.png)
**Figure C.27** - STFT vs WVD vs SPWVD on the same burst segment

- **STFT** (top): blurred but readable
- **WVD** (middle): sharp but corrupted by cross-terms
- **SPWVD** (bottom, linear): **sharpest readable view** of the burst
- SPWVD dB panel: residual circular artifacts (cross-term remnants without filtering)

---

## Why not filter before WVD?

- Bandpass (0.5-4 Hz) would reduce components $\rightarrow$ fewer cross-terms $\rightarrow$ cleaner SPWVD
- **But:** Volumes A-B did not derive filter theory (no FIR, no IIR design)
- Using a tool without deriving it violates the report's principle
- **Segment selection** is the only preprocessing available to us - and it works
- Future work: derive FIR design using Lab 3's windows, then bandpass before WVD

---

<!-- Part 4: Closing -->

## What worked

| Tool | Applied to | Finding |
|------|-----------|---------|
| DFT + Welch | Full recording | Delta dominance (91.8%), 1/f slope = -3.18 |
| STFT | Full recording | Burst pattern: 19% burst, 81% quiet |
| Cross-correlation | Auxiliary vs CZ | CZ is clean ($\rho < 0.03$) |
| SPWVD | 2 s burst segment | Sharper burst localization than STFT |

Each tool addressed a specific limitation of the previous one.

---

## What did not work

- **Raw WVD** unusable on real EEG (49% negative values, cross-term soup)
- **SPWVD dB panel** still contaminated by residual cross-terms (circular patterns)
- **Exponential noise model** approximate on real data (CV = 1.11 vs ideal 1.00)
- **Strongest burst** was amplifier saturation, not a physiological event

These are the method's honest limits, not failures. Reporting them is the point.

---

## Closing claim

The progression **DFT $\rightarrow$ STFT $\rightarrow$ SPWVD** provides increasingly detailed views of the same neonatal EEG signal.

The SPWVD achieves **sub-Heisenberg resolution** on selected clean segments. Its practical value is in the **linear-scale representation**. The dB representation is partially compromised by residual cross-terms without filtering.

The neonatal EEG is consistent with **normal discontinuous neonatal activity**: delta-dominated, bursty, whole-brain synchronous (Lamblin et al., 1999; Andre et al., 2010).

**No clinical diagnosis is made or implied** - these are signal-processing observations.
