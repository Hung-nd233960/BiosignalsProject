# Slides Content Plan

Detailed description of each slide: what it shows, what point it makes, and what figure/equation it uses. Target: ~30-35 slides.

---

## Part 0: Opening (2 slides)

### Slide 1 — Title

- Report title: "From the DFT to the SPWVD: Time-Frequency Analysis Applied to Neonatal EEG"
- Authors: Nguyen Duc Hung - 20233960, Bui Phuong Duy - 23233957, Tran Viet Bach - 23233954
- Course name, date
- **Point:** first impression, sets the scope.

### Slide 2 — Three-volume structure

- Bullet: Volume A = theory (DFT → STFT → WVD → SPWVD)
- Bullet: Volume B = 8 labs validating the theory on model signals
- Bullet: Volume C = application to real neonatal EEG
- Bullet: Arc: each tool addresses a limitation of the previous one
- **Point:** the reader knows the report structure and the progression logic.

---

## Part 1: Volume A Highlights (7 slides)

### Slide 3 — The DFT and its bins

- Equation: $X[k] = \sum_{n=0}^{N-1} x[n] \, e^{-j2\pi kn/N}$ (Equation A.5)
- Bullets: bin spacing $\Delta f = f_s/N$; on-bin = perfect capture, off-bin = leakage; zero-padding interpolates but does not resolve
- **Point:** the DFT's frequency grid is fixed — you get what the bins give you.

### Slide 4 — Windowing and the Dirichlet kernel

- Figure: Lab 3 kernel anatomy (Figure B.10 or B.11 — side-lobe structure)
- Bullets: rectangular window → Dirichlet kernel; main lobe width = $4\pi/M$; side lobes cause leakage; cosine-sum windows (Hann, Hamming, Blackman) trade main lobe width for side-lobe suppression
- **Point:** windowing is the first design choice — you trade resolution for leakage control.

### Slide 5 — Spectral statistics: when is a peak real?

- Bullets: bin power of noise follows exponential distribution (CV = 1.0); detection threshold from false-alarm probability; Welch averaging reduces variance at the cost of resolution
- Equation: $P_{fa} = e^{-\gamma}$ (Equation A.37)
- **Point:** not every peak is a signal — statistics distinguish signal from noise.

### Slide 6 — The STFT and Heisenberg uncertainty

- Equation: $\Delta t \cdot \Delta f \geq \beta$ (Equation A.49)
- Bullets: short window = good time, poor frequency; long window = opposite; cannot have both; the window length is the Heisenberg knob
- **Point:** the STFT adds a time axis but is fundamentally resolution-limited.

### Slide 7 — The WVD: perfect resolution, fatal flaw

- Equation: WVD definition — Equation (A.61)
- Bullets: DFT of the instantaneous autocorrelation; single chirp → razor-sharp diagonal; multi-component → cross-terms at midpoint; 49% negative values on real EEG
- **Point:** the WVD bypasses Heisenberg but generates artifacts that are as strong as the signal.

### Slide 8 — The SPWVD: two independent smoothing knobs

- Equation: SPWVD definition — Equation (A.72)
- Bullets: lag window $h$ smooths frequency (kills frequency-oscillating ghosts); time window $g$ smooths time (kills time-oscillating ghosts); each knob trades resolution for ghost suppression independently
- **Point:** the SPWVD is the usable form of the WVD — the engineer controls the tradeoff per axis.

### Slide 9 — The tool progression summary

- Table or diagram: DFT → STFT → WVD → SPWVD, what each adds, what each costs
- Bullets: DFT (global spectrum, no time); STFT (time-frequency, Heisenberg-limited); WVD (sharp, cross-terms); SPWVD (sharp, controllable ghosts)
- **Point:** each tool fixes a specific limitation of the previous one. This is the spine of the report.

---

## Part 2: Volume B Highlights (10 slides)

### Slide 10 — Lab 1: DFT bins in action

- Figure: off-bin leakage (Figure B.2b) or zero-padding failure (Figure B.7/B.8)
- Bullets: 10.5 Hz maximally off-bin → energy smeared across all bins; zero-padding on 10 + 10.5 Hz → looks resolved but isn't (below $\Delta f_{\min}$)
- **Point:** the DFT's grid matters — off-bin signals leak, and zero-padding cannot create resolution.

### Slide 11 — Lab 3: the Dirichlet kernel and window spectra

- Figure: window comparison (Figure B.14 — Hann/Hamming/Blackman spectra)
- Bullets: all cosine-sum windows share the Dirichlet numerator; first side-lobe: -13 dB (rect), -32 dB (Hann), -43 dB (Hamming), -58 dB (Blackman); decay: -20 dB/octave (rect/Hamming), -60 dB/octave (Hann/Blackman)
- **Point:** derived from first principles, verified numerically — the window zoo has structure.

### Slide 12 — Lab 4: STFT spectrogram of a chirp

- Figure: chirp spectrogram with Heisenberg comparison (Figure B.17 or B.18)
- Bullets: short window tracks the chirp in time but frequency is blurred; long window resolves frequency but smears the chirp; burst reference lines show real timescale
- **Point:** Heisenberg is not abstract — it's visible on the spectrogram.

### Slide 13 — Lab 5: two-tone resolution on the spectrogram

- Figure: two-tone resolution test (Figure B.27 or B.28)
- Bullets: Hann main-lobe width determines minimum resolvable separation; tones closer than $2f_s/M$ merge; confirmed on spectrogram with known frequencies
- **Point:** Lab 3's theory predicts exactly when two tones merge — Lab 5 confirms it.

### Slide 14 — Lab 6: autocorrelation and phase-blindness

- Figure: phase-blind autocorrelation (Figure B.33) or cross-correlation (Figure B.37)
- Bullets: autocorrelation reveals periodicity, discards phase; Wiener-Khinchin: DFT of autocorrelation = power spectrum; cross-correlation detects shared components between signals
- **Point:** autocorrelation → WVD connection. The instantaneous autocorrelation is the WVD's building block.

### Slide 15 — Lab 7: WVD sharpness on a single chirp

- Figure: STFT vs WVD on single chirp (Figure B.40)
- Bullets: STFT = thick blurred diagonal; WVD = razor-sharp line; the WVD completely bypasses Heisenberg for one component
- **Point:** the payoff is real — for a single component, the WVD is perfect.

### Slide 16 — Lab 7: cross-terms ruin multi-component signals

- Figure: chirp + tone STFT vs WVD (Figure B.42)
- Bullets: STFT shows clean superposition; WVD shows sharp components PLUS oscillating ghost at the midpoint; ghost is as energetic as the real components
- **Point:** the fatal flaw — one extra component makes the WVD unusable.

### Slide 17 — Lab 8: WVD → PWVD → SPWVD progression

- Figure: three-row progression (Figure B.46)
- Bullets: WVD = sharp + ghosts; PWVD = frequency ghosts gone, time ghosts survive; SPWVD = both suppressed
- **Point:** step-by-step ghost suppression — each window targets one axis.

### Slide 18 — Lab 8: the duality of ghost types

- Figure: two impulses — frequency-oscillating ghosts (Figure B.49)
- Bullets: components separated in time → frequency-oscillating ghost → PWVD kills it; components separated in frequency → time-oscillating ghost → needs SPWVD
- Table: duality summary (Table B.28)
- **Point:** the ghost type depends on how the components are separated — the SPWVD handles both.

### Slide 19 — Lab 8: two-knob sweep

- Figure: window sweep (Figure B.47)
- Bullets: strong h, weak g → sharp frequency, time ghosts survive; weak h, strong g → clean but frequency blurred; unlike STFT, the axes are independent
- **Point:** the SPWVD gives engineers a knob per axis — the STFT doesn't.

---

## Part 3: Volume C Highlights (10 slides)

### Slide 20 — The dataset

- Bullets: neonatal EEG, 200 Hz, 19 minutes, 24 channels, EDF format; primary channel: CZ (vertex); no filtering applied — only tools derived in A/B
- **Point:** real data, honest constraints.

### Slide 21 — C.1 Triage: delta dominance

- Figure: band power bar chart (Figure C.5 or C.6)
- Bullets: delta = 91.8%, theta = 5.7%, alpha = 1.0%, beta = 1.0%; no alpha/beta rhythms (neonatal); whole-brain synchronous
- **Point:** this EEG is delta-dominated and bursty — sets the direction for C.2-C.5.

### Slide 22 — C.2: is delta rhythmic or 1/f?

- Figure: log-log PSD with 1/f fit (Figure C.9)
- Bullets: slope = -3.18; delta peaks at 0.4-0.6 Hz sit below the 1/f extrapolation; can't tell from the DFT alone if delta is continuous or bursty
- **Point:** the stationary DFT reaches its limit — we need the time axis.

### Slide 23 — C.3: STFT reveals the burst pattern

- Figure: full-recording spectrogram (Figure C.12) or zoomed 60 s (Figure C.13)
- Bullets: delta comes in bursts (19% of time) separated by quiet intervals (81%); burst threshold = 2x median delta power; answers C.2's open question
- **Point:** the STFT adds the time axis and finds the answer: discontinuous delta.

### Slide 24 — C.3: burst detection and statistics

- Figure: burst overlay (Figure C.16)
- Bullets: 19% burst / 81% quiet; max/median ratio = 17x; delta-theta correlation = 0.33 (partially independent)
- **Point:** the burst structure is the dominant non-stationary feature — now quantified.

### Slide 25 — C.4: is CZ clean?

- Bullets: auxiliary channels (25+, 26+, 27+): all ρ < 0.03 vs CZ; EEG inter-channel: ρ = 0.47-0.78; noise floor test: CV = 1.11 (close to exponential 1.00); verdict: CZ is clean for WVD analysis
- **Point:** cross-correlation (Lab 6) validates the channel before we apply the expensive tool.

### Slide 26 — C.5: segment selection — the honest story

- Figure: rejected burst (Figure C.22) + accepted burst (Figure C.23) — side by side or sequential
- Bullets: strongest burst (8.7x median) → amplifier saturation (44 flat samples) → REJECTED; 75th percentile (3.0x median) → clean delta cycles → ACCEPTED
- **Point:** data-driven selection, including the failure. The strongest burst is an artifact.

### Slide 27 — C.5: raw WVD on real EEG — cross-term soup

- Figure: raw WVD (Figure C.25)
- Bullets: 49% of values negative; oscillating cross-terms fill the plane; multi-component real EEG has far more interactions than Lab 7's two components
- **Point:** the raw WVD is unusable on real EEG — confirms A.7.3.

### Slide 28 — C.5: SPWVD — the payoff

- Figure: three-way comparison (Figure C.27)
- Bullets: STFT = blurred; WVD = sharp but unreadable; SPWVD linear = sharpest readable view; dB panel still has residual artifacts (circular patterns = cross-term remnants)
- **Point:** the SPWVD delivers sharper burst localization than the STFT. The honest limit: dB panel is compromised without filtering.

### Slide 29 — C.5: why not filter?

- Bullets: bandpass would reduce components → fewer cross-terms → cleaner SPWVD; but Volumes A-B did not derive filter theory (no FIR, no IIR); using a tool without deriving it violates the report's principle; segment selection is the only preprocessing available — and it works
- **Point:** scope boundary, not a failure. Future work: derive FIR design using Lab 3's windows, then bandpass before WVD.

---

## Part 4: Closing (2-3 slides)

### Slide 30 — What worked

- Table or bullets: DFT found delta dominance; STFT resolved bursts; cross-correlation confirmed clean channel; SPWVD sharpened burst timing
- **Point:** each tool addressed a specific limitation of the previous one.

### Slide 31 — What did not work

- Bullets: raw WVD unusable on real EEG; SPWVD dB panel contaminated by residual cross-terms; exponential noise model approximate (CV = 1.11); strongest burst was an artifact
- **Point:** honest reporting — these are the method's limits, not failures.

### Slide 32 — The closing claim

- The progression DFT → STFT → SPWVD provides increasingly detailed views of the same signal
- The SPWVD achieves sub-Heisenberg resolution on selected segments; its practical value is in the linear-scale representation
- The neonatal EEG is consistent with normal discontinuous activity (cited: Lamblin et al., 1999; Andre et al., 2010)
- No clinical diagnosis made or implied
- **Point:** defensible, honest, in signal-processing language.
