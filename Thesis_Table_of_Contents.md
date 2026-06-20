# Digital Signal Processing: Time-Frequency Analysis from the DFT to the SPWVD, with Application to EEG

**Author:** Nguyễn Đức Hùng — 20233960

---

## How to read this report (½–1 page)

A short framing note stating the gap this report fills. The reader is assumed to know **Signals & Systems** as taught for the masses (circuit theory, control, feedback, modulation) — so they have the DTFT, but are shaky on **bins, Nyquist derived from real (not ideal) signals, windowing, the STFT**, and have never seen the **Wigner-Ville family**. This report rebuilds those specific concepts from the ground up and carries them through to a real EEG signal.

- The structure: **Volume A** (theory + appendices), **Volume B** (labs that derive and verify the theory), **Volume C** (application to a real EEG signal + EEG background appendix).
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

## Appendix B — Hand-written Derivations

- **B.1** The Dirichlet kernel from the finite geometric series (rectangular-window spectrum)
- **B.2** The Hann window coefficients from the edge conditions (value and slope vanish), and the M vs M−1 true-zeros subtlety
- **B.3** *(optional)* The uncertainty bound Δt·Δf ≥ 1/4π via Cauchy–Schwarz (continuous case; note that it does not transfer cleanly to the discrete/periodic domain)

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

## B.1 — Lab 1: The DFT of Basic Signals  *(↔ A.1, A.2)*

- Single tone and a dual-tone "chord": math and code
- On-bin vs. off-bin placement; reading bins, resolution, and the negative-frequency twin
- Verifies: bin spacing, resolution vs. bin count

## B.2 — Lab 2: Statistics on a Noisy Signal  *(↔ A.4)*

- Add noise to a known signal; test whether statistical detection recovers it
- The noise floor, magnitude/phase distributions, averaging to reduce variance
- Verifies: a tone can be detected below the apparent noise by statistics

## B.3 — Lab 3: The DFT of Window Functions  *(↔ A.3)*

- Rectangular, Hann, Blackman: math and code; verify their spectra graphically
- Derive (by heuristic / mathematical approximation): side-lobe level, rolloff rate, main-lobe width, and the **minimum distance between two tones to be distinguishable**
- Verifies: the resolution limit Δf_min ≈ β·(fs/N) and the leakage-vs-resolution tradeoff

## B.4 — Lab 4: The STFT of a Fluctuating Signal  *(↔ A.5)*

- A chirp (or rhythm that drifts): build the spectrogram
- Demonstrate the **Heisenberg tradeoff** by sweeping window length; the effect of windowing and overlap
- Verifies: Δt·Δf is fixed; window length slides you along the curve

## B.5 — Lab 5: Autocorrelation of a Noisy Signal  *(↔ A.6)*

- Autocorrelation of a tone buried in noise; the lag-0 spike vs. the periodic structure
- Wiener–Khinchin checked numerically (autocorrelation ↔ power spectrum)
- Verifies: autocorrelation exposes periodicity that noise hides

## B.6 — Lab 6: The WVD and its Tradeoffs  *(↔ A.7)*

- WVD of a chirp (sharp) vs. the STFT (blurred) — the direct comparison
- WVD of two components → the cross-terms appear; midpoint and oscillation verified
- Verifies: the WVD's super-resolution and its cross-term cost

## B.7 — Lab 7: The SPWVD and its Tradeoffs  *(↔ A.8)*

- WVD → PWVD → SPWVD progression; the two-knob (g, h) sweep
- Smoothing each ghost in the axis it oscillates; resolution traded per axis
- Verifies: the two independent knobs and the cross-term suppression

---

# Volume C — Application to a Real EEG Signal

> **Approach: adaptive-directed.** The instructor's brief is "show what you find," so the *question* is supplied by the analysis itself. Volume C starts by characterizing the data, decides which of the six archetypes the EEG resembles, then deploys the right tool for what was found. The unifying thesis: **EEG can be read as a superposition of the six signal archetypes, and the time-frequency theory of Volumes A–B predicts which tool reveals each feature.** Medical associations are *cited and descriptive*, never asserted as diagnoses.

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

## Notes on scope and structure (for our discussion)

1. **Three volumes:** Volume A (theory + Appendices A–B), Volume B (labs), Volume C (EEG application + Appendix C). One unified slide deck covers all three.
2. **A↔B pairing is explicit** — every theory section has a matching lab, cross-referenced. The spine A.1→A.2→A.3→A.5→A.6→A.7→A.8 is the irreducible core; A.4/Lab 2 (statistics) is the most detachable if time is short.
3. **Heisenberg (A.5.2) is the conceptual hinge** — it closes the STFT half and motivates the entire WVD half.
4. **Cohen's class (A.8.5) is the unifying frame** — stated once, it ties spectrogram/WVD/PWVD/SPWVD into one family.
5. **Volume C is adaptive-directed** — the centerpiece (C.2 vs C.3 vs C.4 vs C.5) is decided by the C.1 triage of the actual EEG, not assumed in advance.
6. **Already built in "Part 3":** A.3, A.5, A.7, A.8 and Labs 3, 4, 6, 7 are largely done — roughly the windowing-through-SPWVD half of the thesis.
