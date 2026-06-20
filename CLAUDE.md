# Project Context

We are undergraduate Biomedical Engineering students taking a Digital Biosignal Processing class. This is the final report on signal processing. Our team chose to study the basics of the DFT and STFT, then expand into the Wigner-Ville Distribution (WVD) and Smoothed Pseudo Wigner-Ville Distribution (SPWVD) as the advanced topic to learn and apply. The application target is real EEG data.

The report is split into three volumes (see `Thesis_Table_of_Contents.md`):

- **Volume A** — Theory (A.1–A.8) + Appendices A–B
- **Volume B** — Labs (B.1–B.7) with EEG-realism signal constraints
- **Volume C** — Application to real EEG + Appendix C (EEG background)

A unified slide deck accompanies all three volumes.

---

## Tone and Writing Style

- **Academic but honest.** Serious enough for a university report, but written as students explaining concepts we learned — not as textbook authors. No false authority; if something is approximate or heuristic, say so.
- **Neutral.** No hype, no hedging. State what the method does, what it costs, and what we observed.

---

## Code Standards

- **Show your tools.** Every code block must explicitly display the libraries/imports being used. If a calculation is built from scratch (no library call), show the implementation explicitly.
- **Comment every line.** Code explanation must include inline comments stating what each line does.
- **Compact and modular.** Wrap reusable operations in functions. Each function does one thing. Code should be easy to copy into another notebook or script and run.
- **Reproducibility.** Any signal involving randomness (noise, stochastic processes) must use a fixed random seed. State the seed explicitly.

---

## Model Signals (Volume B Labs)

All model signals must satisfy the lab constraints (see Volume B header in `Thesis_Table_of_Contents.md`):

- All frequency components below **100 Hz**
- Signal duration at least **1 200 s** (20 minutes)
- Sampling frequency **200–300 Hz** (default to **250 Hz** unless a lab specifically requires otherwise)

Additional requirements:

- **Well-behaved and continuous.** No discontinuities in the signal. If a signal is piecewise, ensure smooth transitions.
- **Math first.** Before any code, show the explicit mathematical formula for the model signal (e.g. the sum-of-sines expression). The formula is the specification; the code implements it.
- **Time-domain first.** Every signal must have its time-domain waveform plotted and inspected before any spectral or time-frequency analysis. This is the control step — no skipping it.

---

## Graphs and Figures

### Quality

- Render all figures at **300+ DPI**.
- Export as **PNG at 300 DPI** for raster figures, or **PDF/SVG** for vector figures.
- All axes must have **labels with physical units** (Hz, s, µV, dB, µV²/Hz, etc.) — not just variable names.
- Use **perceptually uniform colormaps** (`viridis`, `inferno`, `plasma`, `magma`). **Never use `jet` or `rainbow`** — they distort perception and are not accessible.
- Colors should be high-contrast and distinguishable in grayscale print.

### The Dual-Stack Rule

Every spectral or time-frequency plot must appear as a **dual stack: linear scale first, then log/dB scale below it.**

- The **linear-scale plot is primary** — it shows the physical quantity in its real units and is always presented first.
- The **log/dB-scale plot is secondary** — it reveals dynamic range and low-level structure, but is always shown *after* the linear representation.

For real EEG signals (Volume C), this is non-negotiable:

- Amplitude in **µV** (linear first).
- Power spectral density in **µV²/Hz** (linear first), then dB scale.
- Time-frequency representations in physical units (linear first), then log/dB.
- The log/dB plot exists to support interpretation, not to replace the physical representation.

For model signals (Volume B), the same dual-stack applies. Linear scale grounds the reader in the actual signal; log/dB scale exposes features that linear compresses.

---

## Interpretation and Analysis

- **Physical interpretation governs.** When working with real signals (EEG), all numbers, axes, and discussion must respect physical units and physiological meaning. A "peak at 10 Hz" is an alpha-band observation, not just a bin number.
- **Honest reporting.** If a method fails or produces artifacts (e.g. WVD cross-terms on multi-component EEG), say so. The point is to characterize the tools, not to sell them.
- **Clinical claims are cited, not asserted.** Any association between a signal feature and a clinical condition must reference published literature. We are engineering students, not clinicians.
