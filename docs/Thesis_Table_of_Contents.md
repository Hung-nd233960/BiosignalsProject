# Digital Signal Processing: Time-Frequency Analysis from the DFT to the SPWVD, with Application to EEG

**Author:** Nguyễn Đức Hùng — 20233960

---

## How to read this report (½–1 page)

A short framing note stating the gap this report fills. The reader is assumed to know **Signals & Systems** as taught for the masses (circuit theory, control, feedback, modulation) — so they have the DTFT, but are shaky on **bins, Nyquist derived from real (not ideal) signals, windowing, the STFT**, and have never seen the **Wigner-Ville family**. This report rebuilds those specific concepts from the ground up and carries them through to a real EEG signal.

- The structure: **Volume A** (theory + Appendix A), **Volume B** (labs + Appendix B window derivations), **Volume C** (application to a real EEG signal + Appendix C EEG background).
- Each Volume A section opens with a one-line **"gap callout"**: what Signals & Systems gives you, and what this section adds.
- Each lab in Volume B is cross-referenced to the theory section it tests (↔ A.x).
- Appendix A (the six-signal taxonomy) is the reference grid that every lab and every Volume C section points back to.
- A **unified slide deck** accompanies all three volumes.

---

# Volume A — Theoretical Background

> The spine: **sampling → DFT/bins → leakage/windowing → statistics → STFT/uncertainty → autocorrelation → WVD → PWVD/SPWVD.**
> Each section names the gap it fills relative to a standard Signals & Systems course.

## A.1 Signal Theory

*Gap: Signals & Systems teaches sampling via ideal band-limited signals; this rebuilds it from real, finite signals.*

- **A.1.1** What a discrete signal *means* — samples vs. the continuous original; the loss and what is kept
- **A.1.2** Discrete frequency derived: why normalized frequency is **f / fs** (the unit circle picture)
- **A.1.3** Nyquist and aliasing — derived from the f/fs circle, not assumed from an ideal filter
- **A.1.4** Sampling **count** vs. sampling **duration** — the two independent knobs and what each controls

## A.2 The DTFT and the DFT

*Gap: they know the DTFT; they have never been told clearly what a bin is.*

- **A.2.1** From the DTFT (the part they know) to the DFT as its sampled version
- **A.2.2** What a **bin** is: sampling the continuous DTFT; bin spacing = fs / N
- **A.2.3** Resolution vs. bin count — the distinction most people conflate (zero-padding adds bins, not resolution)
- **A.2.4** The DFT as an orthonormal basis; Parseval / energy conservation

## A.3 Leakage and Windowing

*Gap: windowing is never taught; the rectangular cut-off is invisible to them.*

- **A.3.1** The hidden rectangular window in the DFT's summation limits
- **A.3.2** Leakage = the window's own spectrum convolved onto each tone (on-bin vs. off-bin)
- **A.3.3** The cosine-sum family — Rectangular, Hann, Hamming, Blackman — formulas and the design goal behind each
- **A.3.4** Main-lobe vs. side-lobe tradeoff; edge-smoothness → side-lobe rolloff rate (the 1/ω^p rule)
- **A.3.5** The resolution limit: **Δf_min ≈ β · (fs / N)**, with β the main-lobe width factor
- **A.3.6** *(deferred)* the Gaussian window and minimum-uncertainty — flagged here, used in A.8

## A.4 Statistics and the DFT

*Gap: spectral estimation as a statistical act is never framed; they treat the DFT as exact.*

- **A.4.1** Each bin as a random variable; magnitude and phase distributions under noise
- **A.4.2** The noise floor; detecting a tone against it; exponential power statistics
- **A.4.3** Inference from bin statistics — what the distribution of bins tells you about the underlying signal
- **A.4.4** Averaging to reduce variance (the idea behind Welch's method)

## A.5 The STFT and Spectrograms

*Gap: the STFT is never taught; the time-frequency plane is new to them.*

- **A.5.1** The windowed DFT slid in time; the time-frequency plane as the central object
- **A.5.2** **The uncertainty principle: Δt · Δf ≥ 1/4π** — the central tradeoff of the whole report; window length as the slider
- **A.5.3** Parameters beyond the window: hop size, overlap, and the COLA condition for reconstruction
- **A.5.4** Reading a spectrogram; the cost of each parameter choice

## A.6 Autocorrelation

*Gap: introduced as the bridge that makes the Wigner-Ville family understandable.*

- **A.6.1** Definition; the periodicity detector; r[0] = signal energy
- **A.6.2** Wiener–Khinchin: autocorrelation and the power spectrum are a Fourier pair
- **A.6.3** Why autocorrelation keeps frequency but discards phase (and why its transform is the *power* spectrum)
- **A.6.4** Autocorrelation signatures of a tone, an impulse, and noise (sets up the taxonomy)

## A.7 The Wigner-Ville Distribution

*Gap: entirely new; the payoff of the autocorrelation groundwork.*

- **A.7.1** From the static autocorrelation to the **instantaneous autocorrelation**; the WVD definition
- **A.7.2** The chirp → delta: how the WVD beats the STFT's blur for a single component
- **A.7.3** The quadratic nature and **cross-terms**: midpoint location, oscillation at the difference frequency
- **A.7.4** Narrowband vs. broadband ghosts — the duality (tones vs. impulses)
- **A.7.5** Why the **analytic signal** (Hilbert transform) is required — removing the DC self-ghost

## A.8 The PWVD and the SPWVD

*Gap: the endpoint of the report; the usable form of the WVD.*

- **A.8.1** PWVD: the lag window **h** → smoothing in frequency → killing frequency-oscillating ghosts
- **A.8.2** The limitation: PWVD smooths only frequency; time-oscillating ghosts survive
- **A.8.3** SPWVD: the time window **g** → two independent smoothing knobs
- **A.8.4** The resolution-vs-ghost tradeoff, now controllable per axis; choosing windows (Hann / Gaussian) and the minimum smoothing rule
- **A.8.5** *Closing frame:* **Cohen's class** — the spectrogram and the raw WVD as the two extremes, SPWVD in between; wavelets noted as a future direction

## Appendix A — Signal Taxonomy (reference grid)

A table of the **six signal archetypes** × their behaviour under each transform. Referenced from every lab and every Volume C section.

| archetype | time-domain form | DFT | STFT | WVD | SPWVD |
| --- | --- | --- | --- | --- | --- |
| Single tone | sinusoid | single spike | horizontal line | sharp horizontal line (+ DC ghost if real) | clean line |
| Mixed tone | sum of sinusoids | several spikes | parallel lines | lines + midpoint ghosts | lines, ghosts suppressed |
| Chirp (fluctuating tone) | swept-frequency | smeared band | blurred diagonal | razor-sharp diagonal | sharp diagonal |
| Multi-chirp (mixed fluctuating) | sum of chirps | overlapping smears | crossing blurred diagonals | diagonals + moving ghosts | diagonals, ghosts suppressed |
| Transient / pulse | impulse / short burst | flat (broadband) | vertical line | vertical line (+ broadband ghosts) | vertical line, ghosts suppressed |
| Noise | random process | flat, random phase | uniform speckle | speckle (+ self-interference) | speckle, not removed by smoothing |

*Each cell expanded in the appendix with the formula and the signature; the "why" cross-referenced to Volume A.*

---

# Volume B — Application and Derivation (Labs)

> Each lab derives the mathematics, implements it in code, and verifies the theory graphically.
> Labs are paired to the theory sections they test.

## Lab signal constraints (EEG realism)

All model signals used in the labs must satisfy the following constraints so that every demonstration operates in a regime representative of real EEG:

| Constraint | Requirement |
| --- | --- |
| **Maximum frequency** | All signal components below **100 Hz** |
| **Signal duration** | At least **1 200 s** (20 minutes) |
| **Sampling frequency** | **200–300 Hz** (use **250 Hz** unless a specific lab requires otherwise) |

*Rationale:* scalp EEG is band-limited to roughly 0.5–100 Hz; clinical recordings typically last 20+ minutes at 250 Hz sampling. Constraining the labs this way ensures that bin spacing, resolution limits, and window choices encountered in Volume B transfer directly to the real EEG analysis in Volume C.

## Lab report template

Every lab follows a fixed six-section structure:

1. **Introduction** — Gloss the corresponding Volume A theory (cite section numbers). State what this lab tests and what we expect to confirm.
2. **Setup** — The model signal's mathematical formula first (the specification), then the code that implements it.
3. **Parameters** — Table listing every value used (fs, N, frequencies, amplitudes, window type, SNR, seed, etc.).
4. **Results** — Figures (time-domain first, then spectral/dual-stack). Observations on what the plots show.
5. **Verification** — Compare measured results against Volume A predictions quantitatively (e.g. "Equation (A.6) predicts Δf = 0.05 Hz; measured: 0.05 Hz").
6. **Conclusion** — What was confirmed, any surprises, bridge to the next lab.

---

## B.1 — Lab 1: The DFT of Basic Signals  *(↔ A.1, A.2)*

- **Introduction:** the DFT as a sampled DTFT (A.2.1); bins, bin spacing, resolution (A.2.2–A.2.3)
- **Setup:** single tone and a dual-tone "chord" — math and code; 1-second segment for on-bin vs. maximally off-bin comparison
- **Parameters:** fs, N, f0 (on-bin 10.0 Hz and off-bin 10.5 Hz), amplitudes, duration
- **Results:** on-bin → clean spike, off-bin → leakage across entire band; dual-tone resolution; zero-padding adds bins not resolution
- **Verification:** bin spacing matches Equation (A.6); on-bin = zero leakage; off-bin = Dirichlet kernel side lobes visible
- **Conclusion:** what the DFT reveals and what it misses; leakage motivates windowing (Lab 3)

## B.2 — Lab 2: Statistics on a Noisy Signal  *(↔ A.4)*

- **Introduction:** DFT bins under noise are random variables (A.4.1); the noise floor and spectral detection (A.4.2); the periodogram variance problem (A.4.3)
- **Setup:** known tone buried in noise at controlled SNR — math and code
- **Parameters:** fs, N, f0, A, σ (noise), SNR, seed, detection threshold γ
- **Results:** magnitude/phase distributions of bins; noise floor estimation; detection via Equation (A.28); Welch averaging progression
- **Verification:** bin power follows exponential distribution (A.25); detection threshold γ matches false-alarm probability (Table A.3); Welch reduces variance by 1/L (A.32)
- **Conclusion:** spectral detection works where time-domain thresholding fails; Welch trades resolution for reliability

## B.3 — Lab 3: Windowing and the Dirichlet Kernel  *(↔ A.3)*

- **Introduction:** Lab 1 showed leakage from off-bin tones under the rectangular window. This lab derives *why* leakage occurs (the Dirichlet kernel) and *how* the cosine-sum windows suppress it.
- **Setup:** derive the DFT of the rectangular window from the geometric series → Dirichlet kernel; expand Hann, Hamming, Blackman into pure sine form (shifted Dirichlet kernels); all at $M = 256$ (typical EEG epoch). Code: `src/appendix_b/appendix_b.py`.
  - AB.1: Dirichlet kernel derivation (geometric series → sin/sin closed form)
  - AB.2: Graphical anatomy — main lobe, side lobes, nulls, skew of maxima
  - AB.3.1: First side-lobe strength (visual + $k = 1.5$ approximation)
  - AB.3.2: Decay rate (analytical $1/\omega$ + regression on 6 maxima)
  - AB.3.3: Cosine-sum identity → pure sine forms for Hann, Hamming, Blackman; side-lobe cancellation mechanism; Hamming rolloff explained
  - AB.3.4: Comparative graph of all window spectra (linear scale)
- **Parameters:** M = 256, fs = 250, zero-pad = 2048×, windows: Rectangular, Hann, Hamming, Blackman
- **Results:** Figures AB.1–AB.6; side-lobe data table; regression (slope = −0.985, R² = 0.99997, −5.9 dB/oct)
- **Verification:** first side-lobe −13.3 dB (textbook −13.0); decay 5.9 dB/oct (theoretical 6.0); Hann/Blackman side-lobe cancellation confirmed from the pure sine form
- **Conclusion:** all Table A.1 properties are derived from first principles; the pure sine form is the analytical tool for understanding any cosine-sum window

## B.4 — Lab 4: The STFT of a Fluctuating Signal  *(↔ A.5)*

- **Introduction:** the STFT as a windowed DFT slid in time (A.5.1); the uncertainty principle (A.5.2); hop size and COLA (A.5.3)
- **Setup:** a chirp (or rhythm that drifts) — math and code
- **Parameters:** fs, N, f0, μ (chirp rate), window type, window length M (swept), hop size H, overlap
- **Results:** spectrograms at multiple window lengths; the Heisenberg tradeoff visualized; effect of overlap on sample coverage
- **Verification:** Δt·Δf product is constant across window lengths (A.39); COLA condition confirmed for Hann at 50% overlap (A.42)
- **Conclusion:** the STFT forces a choice; no single window captures multi-scale features; the spectrogram is now the reader's primary tool → used in Lab 5

## B.5 — Lab 5: Two-Tone Resolution on the Spectrogram  *(↔ A.3, A.5, Lab 3)*

- **Introduction:** Lab 3 derived the resolution limit $\Delta f_{\min} \approx \beta \cdot f_s/M$ from the Dirichlet kernel. Lab 4 introduced the spectrogram. This lab confirms the resolution limit visually: two stationary tones on a spectrogram — too close = one line, far enough = two lines.
- **Setup:** two stationary tones with swept separation, STFT with Rectangular, Hann, Blackman windows — code
- **Parameters:** fs, M (window length), f1, f2 (separation sweep), window type, hop size, overlap
- **Results:** spectrograms at several separations per window; one horizontal line splitting into two at the predicted $\beta \cdot f_s/M$
- **Verification:** Rectangular splits at $\beta = 1$ limit, Hann at $\beta = 2$, Blackman at $\beta = 3$ — matching Lab 3 derivations and Table A.1
- **Conclusion:** the resolution-leakage tradeoff confirmed visually on spectrograms; window choice for EEG band separation justified

## B.6 — Lab 6: Autocorrelation of a Noisy Signal  *(↔ A.6)*

- **Introduction:** autocorrelation as periodicity detector (A.6.1); Wiener-Khinchin theorem (A.6.2); phase-blindness (A.6.3)
- **Setup:** tone buried in noise — math and code; autocorrelation computed from definition
- **Parameters:** fs, N, f0, A, σ (noise), SNR, seed
- **Results:** autocorrelation of tone-in-noise: lag-0 spike vs. periodic peaks; Wiener-Khinchin verified numerically (DFT of autocorrelation vs. power spectrum)
- **Verification:** r[0] equals total energy (A.45); DFT of r[l] matches |X[k]|² (A.47); periodicity visible at lags where noise contribution vanishes
- **Conclusion:** autocorrelation exposes periodicity that noise hides; phase is lost → bridge to WVD (instantaneous autocorrelation)

## B.7 — Lab 7: The WVD and its Tradeoffs  *(↔ A.7)*

- **Introduction:** instantaneous autocorrelation and the WVD definition (A.7.1); single-component sharpness (A.7.2); cross-terms (A.7.3)
- **Setup:** single chirp; two-component signal (chirp + tone) — math and code
- **Parameters:** fs, N, f0, μ (chirp rate), component frequencies, analytic signal flag
- **Results:** WVD of chirp (sharp diagonal) vs. STFT (blurred diagonal); WVD of two components → cross-terms at midpoint with difference-frequency oscillation
- **Verification:** chirp trajectory matches instantaneous frequency (AA.4); cross-term location matches midpoint prediction (A.7.3); analytic signal removes DC ghost (A.7.5)
- **Conclusion:** the WVD's super-resolution is real but cross-terms make it unusable for multi-component signals → motivation for SPWVD

## B.8 — Lab 8: The SPWVD and its Tradeoffs  *(↔ A.8)*

- **Introduction:** PWVD lag window (A.8.1); SPWVD two-knob smoothing (A.8.3); resolution-vs-ghost tradeoff (A.8.4)
- **Setup:** multi-component signal from Lab 7 — same signal, now processed with PWVD and SPWVD
- **Parameters:** fs, N, signal from Lab 7, lag window h (type, length), time window g (type, length)
- **Results:** WVD → PWVD → SPWVD progression; two-knob (g, h) sweep; each ghost suppressed in the axis it oscillates
- **Verification:** PWVD suppresses frequency-oscillating ghosts but not time-oscillating ones (A.8.2); SPWVD suppresses both; resolution traded per axis independently (A.8.4)
- **Conclusion:** the SPWVD is the usable form of the WVD; the two knobs give explicit control over the resolution-ghost tradeoff → ready for real EEG (Volume C)

---

# Volume C — Application to a Real EEG Signal

> **Approach: adaptive-directed.** The instructor's brief is "show what you find," so the *question* is supplied by the analysis itself. Volume C starts by characterizing the data, decides which of the six archetypes the EEG resembles, then deploys the right tool for what was found. The unifying thesis: **EEG can be read as a superposition of the six signal archetypes, and the time-frequency theory of Volumes A–B predicts which tool reveals each feature.** Medical associations are *cited and descriptive*, never asserted as diagnoses.

## Volume C standards

All Volume B rules apply (code format, graphs, numbering, dual-stack, six-section template). Additional EEG-specific standards:

- **Physical units mandatory.** Amplitude in µV, PSD in µV²/Hz, time in s. No arbitrary or normalized units.
- **Data-first, not math-first.** Each section starts with the data (channel, time range, raw waveform), then applies the tool.
- **Adaptive-directed.** C.1 (triage) decides the direction. Each subsequent section justifies its existence from C.1's findings.
- **Reproducibility.** Dataset source, file format, sampling rate, channel selection, preprocessing — all stated explicitly. MNE for all EEG I/O.
- **Artifact handling stated per section.** Removed or not, how, and what effect.
- **Parameter justification.** Welch segment length, window choice, overlap — justified from Labs 2–3, not default.
- **WVD/SPWVD on selected segments only.** Selection criteria stated. Full-record WVD is not attempted.
- **Clinical language.** Associations cited from literature, never asserted as diagnoses.

## C.1 First Look and Triage  *(decide the direction from the data)*

- **C.1.1** The EEG signal: source, sampling rate, duration, channel(s) — stated concretely
- **C.1.2** **Time-domain plot first** (the report's own principle) — what is visible by eye
- **C.1.3** A first broad DFT — which standard bands dominate (δ 0.5–4, θ 4–8, α 8–13, β 13–30, γ >30 Hz)
- **C.1.4** **Triage decision:** which of the six archetypes does this EEG most resemble? → sets the centerpiece for C.2–C.5

## C.2 Stationary Characterization — DFT and Band Power  *(tone / mixed-tone archetypes)*

- Dominant rhythm(s); band-power estimate using the window chosen for the *actual* band spacing present (justified from Lab 3)
- The resolution limit applied to real band separation
- **Finding:** the resting spectral signature of the signal

## C.3 Time-Varying Characterization — STFT Spectrogram  *(chirp / non-stationary archetypes)*

- Does the rhythm persist, come and go, or drift in frequency over time?
- The Heisenberg choice made explicitly for *this* signal (resolve the band vs. catch the timing)
- **Finding:** the temporal structure the stationary DFT could not show

## C.4 Events and Artifacts — Statistics and Transient Detection  *(transient / noise archetypes)*

- Time-localized broadband events flagged against the noise floor (blinks, movement, spikes)
- Statistical detection from Volume A.4 used to separate signal from artifact
- **Finding:** what is genuine signal vs. what is artifact

## C.5 High-Resolution Time-Frequency — WVD / SPWVD  *(the payoff and the honest test)*

- Applied to a **selected, cleaned segment** (analytic signal, one or two components) — the correct way to use the WVD family on EEG
- Why the raw WVD fails on the full multi-component noisy record (cross-term soup), and how the SPWVD's two knobs (tuned per Lab 7) recover a readable picture
- **Finding:** the sharpest time-frequency view of one chosen feature, with the tradeoff characterized honestly

## C.6 Synthesis — What We Found

- Which archetype–tool pairings held up on real data; which EEG features needed which tool
- The honest limits: where the WVD family helped and where it did not
- The defensible closing claim — in signal-processing language — with any clinical association supported by a citation, not asserted

## Appendix C — EEG Background (domain knowledge for Volume C)

> *Purpose: give the reader (and the author) the minimum, cited EEG background needed to interpret Volume C honestly. This appendix is descriptive domain knowledge, kept separate from the DSP contribution of Volumes A–B. All clinical statements are referenced; nothing here is an original medical claim.*

**C.1 What EEG is**

- C.1.1 The physical signal: scalp potentials from cortical activity; microvolt scale
- C.1.2 The 10–20 electrode system; what the channel names (Fp, F, C, P, O, T) mean and where they sit
- C.1.3 Sampling, filtering, and referencing as they affect the recorded signal (ties to Volume A.1)
- C.1.4 Non-brain channels in a recording: ECG, EMG, EOG — what they are and why they appear

**C.2 The EEG frequency bands**

- C.2.1 Delta (0.5–4 Hz), theta (4–8), alpha (8–13), beta (13–30), gamma (>30) — the standard bands
- C.2.2 What each band is associated with (state, arousal, sleep) — adult baseline
- C.2.3 Band power as the quantity of interest; why this connects to the DFT (Volume A.2, A.4)

**C.3 Infant and neonatal EEG (the relevant case)**

- C.3.1 How the immature brain's EEG differs from the adult: delta dominance, lower frequencies, larger amplitudes
- C.3.2 **Continuity vs. discontinuity** — what "bursty" / discontinuous activity means
- C.3.3 Tracé discontinu and tracé alternant: the normal discontinuous patterns of prematurity and term sleep
- C.3.4 Why discontinuity is a **maturity marker** — it decreases as the brain develops
- C.3.5 The clinical caution: what these patterns do and do not indicate (and the limits of what DSP alone can claim)

**C.4 Common EEG features and artifacts**

- C.4.1 Rhythms vs. transients: sharp waves, spikes, K-complexes, sleep spindles (which exist at which age)
- C.4.2 Artifacts: movement, electrode pop, ECG contamination, EMG, eye movement — how each looks
- C.4.3 Mapping EEG features to the six signal archetypes of Appendix A (the bridge to the DSP analysis)

**C.5 The dataset used in this report**

- C.5.1 Source, subject (infant), montage, sampling rate, duration — the concrete facts
- C.5.2 First characterization: delta-dominated, discontinuous, low-passed to ~28 Hz; ECG/EMG on auxiliary channels
- C.5.3 What to expect, and which DSP tool each expected feature calls for (the plan that Volume C executes)

**References** — EEG textbooks and infant-EEG literature for every clinical statement made.

---

## Formatting conventions

- **Equations** are numbered sequentially per volume: `(A.1)`, `(A.2)`, … / `(B.1)`, … / `(C.1)`, …
- **Figures** are numbered sequentially per volume: `Figure A.1`, `Figure A.2`, … / `Figure B.1`, … / `Figure C.1`, …
- **Tables** are numbered sequentially per volume: `Table A.1`, `Table A.2`, … / `Table B.1`, … / `Table C.1`, …
- All equations, figures, and tables are **referenced by number** in the text — never by position ("above", "below", "the following").

---

## Notes on scope and structure (for our discussion)

1. **Three volumes:** Volume A (theory + Appendix A), Volume B (8 labs + Appendix C EEG background in Volume C). One unified slide deck covers all three.
2. **A↔B pairing is explicit** — every theory section has a matching lab. Lab 3 (windowing) is the derivation-heavy lab; Lab 5 (two-tone spectrogram) is its experimental confirmation. The spine: Lab 1 (DFT) → Lab 2 (stats) → Lab 3 (windows) → Lab 4 (STFT) → Lab 5 (resolution) → Lab 6 (autocorrelation) → Lab 7 (WVD) → Lab 8 (SPWVD).
3. **Heisenberg (A.5.2) is the conceptual hinge** — it closes the STFT half and motivates the entire WVD half.
4. **Cohen's class (A.8.5) is the unifying frame** — stated once, it ties spectrogram/WVD/PWVD/SPWVD into one family.
5. **Volume C is adaptive-directed** — the centerpiece (C.2 vs C.3 vs C.4 vs C.5) is decided by the C.1 triage of the actual EEG, not assumed in advance.
