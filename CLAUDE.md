# Project Context

We are undergraduate Biomedical Engineering students taking a Digital Biosignal Processing class. This is the final report on signal processing. Our team chose to study the basics of the DFT and STFT, then expand into the Wigner-Ville Distribution (WVD) and Smoothed Pseudo Wigner-Ville Distribution (SPWVD) as the advanced topic to learn and apply. The application target is real EEG data.

The report is split into three volumes (see `TABLE_OF_CONTENTS.md`):

- **Volume A** — Theory (A.1–A.8) + Appendix A (signal taxonomy)
- **Volume B** — Labs (B.1–B.8) with EEG-realism signal constraints
- **Volume C** — Application to real EEG + Appendix C (EEG background)

A unified slide deck accompanies all three volumes.

---

## Tone and Writing Style

- **Academic but honest.** Serious enough for a university report, but written as students explaining concepts we learned — not as textbook authors. No false authority; if something is approximate or heuristic, say so.
- **Neutral.** No hype, no hedging. State what the method does, what it costs, and what we observed.
- **Discrete realm.** All mathematical formulas, derivations, and interpretations are explicitly discrete-time and discrete-frequency. Use summations (not integrals), sample indices (not continuous time), and normalized/discrete frequency. The continuous-time form may be mentioned for motivation, but the working formula is always the discrete version.
- **No "gap callouts" in the reports.** The gap framing ("Signals & Systems gives you X, this section adds Y") is internal planning language in `TABLE_OF_CONTENTS.md`. Do NOT include it in the actual volumes (volume_A.md, volume_B.md, volume_C.md) or slides. The reports teach the material directly without referencing what the reader is assumed to be missing.
- **No `>` blockquotes in submitted reports.** Blockquotes (`> text`) are internal editorial notes and planning reminders. They must be removed before submission. The `prettier.py` script strips them automatically.

---

## Environment

- **Conda (miniforge3).** All code runs inside the `biosignals` conda environment defined in `environment.yml` at the project root.
- **Setup from scratch:** `conda env create -f environment.yml` then activate as above.
- **Adding dependencies:** add to `environment.yml`, then `conda env update -f environment.yml --prune`. Never `pip install` directly — keep the environment file as the single source of truth.
- **Core stack:** Python (≥3.11), NumPy, SciPy, Matplotlib, JupyterLab, MNE (EEG analysis for Volume C). R and rpy2 are included for collaborators who prefer R.

---

## Reference Files

- `TABLE_OF_CONTENTS.md` — master plan with all section outlines, lab templates, and internal planning notes (gap callouts). This is the source of truth for structure.
- `CLAUDE.md` — this file. Project standards, conventions, and rules for all contributors (human and AI).
- `src/README.md` — code contributor guide: directory layout, import rules, naming conventions.
- `environment.yml` — conda environment definition. Single source of truth for dependencies.
- `template/reference.docx` — Word template for pandoc: Roboto body, Consolas code, margins (top 2cm, left 2cm, right 1cm, bottom 1cm), page numbers bottom center.
- `docs/METHODOLOGY.md` — records every rule in this file, why it was added, and what problem it solved. **When adding or changing a rule in CLAUDE.md, update METHODOLOGY.md with the reason.**

---

## Building Documents

**Touch up markdown before converting** (removes horizontal rules, replaces fancy dashes):

```bash
python build_docs/prettier.py docs/volume_B.md
```

**Markdown to Word (.docx)** via pandoc with the reference template:

```bash
cd docs && ~/miniforge3/bin/pandoc volume_B.md \
    -o ../output/volume_B.docx \
    --from=markdown+tex_math_dollars+pipe_tables \
    --mathml \
    --reference-doc=../template/reference.docx
```

Or use the build script: `bash build_docs/build_docx.sh`

**Markdown to LaTeX/PDF** via pandoc + tectonic:

```bash
cd docs && ~/miniforge3/bin/pandoc volume_B.md \
    -o ../output/volume_B.tex \
    --from=markdown+tex_math_dollars --standalone \
    -V geometry:margin=2cm
cd ../output && ~/miniforge3/bin/tectonic volume_B.tex
```

**Slides** via Marp (markdown to HTML):

```bash
export PATH=~/miniforge3/bin:$PATH
cd docs && marp slides.md -o ../output/slides.html --allow-local-files
```

For PDF/PPTX slides, install Chromium and use `marp slides.md -o ../output/slides.pdf`. Or print-to-PDF from the HTML in a browser.

---

## Code Standards

- **Show your tools.** Every code block must explicitly display the libraries/imports being used. If a calculation is built from scratch (no library call), show the implementation explicitly. When using external functions (e.g. `scipy.stats.linregress`, `scipy.signal.argrelextrema`), state briefly what they do and what method they use — not just "fit a line" or "find maxima."
- **Comment every line.** Code explanation must include inline comments stating what each line does.
- **Compact and modular.** Wrap reusable operations in functions. Each function does one thing. Code should be easy to copy into another notebook or script and run.
- **Reproducibility.** Any signal involving randomness (noise, stochastic processes) must use a fixed random seed. State the seed explicitly.
- **Import from `src/common/`.** Constants (`FS`, `DURATION`, `SEED`, `DPI`), signal generators, window functions, and plotting utilities all live in `src/common/`. Never hardcode these values in lab code. See `src/README.md` for the full contributor guide.
- **Every figure must have its code shown in the report.** The instructor does not go to GitHub. Before each figure, show the code that produces it — the function call, the parameters, the data flow. Pattern: **code block → figure → interpretation.** No figure should appear without the reader seeing exactly how it was generated.
- **Numbers that drive decisions must be traceable to code.** If a number determines the analysis direction (e.g. "91.8% delta" drives the triage, "1/f slope = -3.18" determines signal vs. noise), the code that computes it must appear immediately before, or the raw output must be shown. Descriptive statistics (amplitude range, std) are fine as prose. The critical chain is: **code → number → decision.**

---

## Model Signals (Volume B Labs)

All model signals must satisfy the lab constraints (see Volume B header in `TABLE_OF_CONTENTS.md`):

- All frequency components below **100 Hz**
- Signal duration at least **1 200 s** (20 minutes)
- Sampling frequency **200–300 Hz** (default to **250 Hz** unless a lab specifically requires otherwise)

Additional requirements:

- **Well-behaved and continuous.** No discontinuities in the signal. If a signal is piecewise, ensure smooth transitions.
- **Math first.** Before any code, show the explicit mathematical formula for the model signal (e.g. the sum-of-sines expression). The formula is the specification; the code implements it.
- **Time-domain first.** Every signal must have its time-domain waveform plotted and inspected before any spectral or time-frequency analysis. This is the control step — no skipping it.

---

## Lab Report Structure (Volume B)

Every lab follows a fixed six-section template:

1. **Introduction** — Gloss the corresponding Volume A theory (cite section numbers). State what this lab tests and what we expect to confirm.
2. **Setup** — The model signal's mathematical formula (the specification), then the code that implements it. Explain the construction.
3. **Parameters** — A table listing every value used: `fs`, `N`, frequencies, amplitudes, window type, SNR, seed, etc. One glance = full reproducibility.
4. **Results** — Figures (time-domain first, then spectral/dual-stack per the Dual-Stack Rule). Observations on what the plots show.
5. **Verification** — Compare measured results against Volume A predictions quantitatively. E.g. "Equation (A.6) predicts Δf = 0.05 Hz; measured: 0.05 Hz." This is where the "Verifies: ..." claim from the ToC is fulfilled.
6. **Conclusion** — What was confirmed, any surprises or deviations, bridge to the next lab.

---

## Numbering

All tables, figures, and equations must be **explicitly numbered** using the volume prefix. The three types have **independent counters** — Equation (B.1), Figure B.1, and Table B.1 are three different objects.

- **Equations:** always in parentheses → `(A.1)`, `(B.1)`, `(C.1)`, … Referenced as `Equation (B.3)`.
- **Figures:** always with "Figure" prefix → `Figure A.1`, `Figure B.1`, `Figure C.1`, … Referenced as `Figure B.3`.
- **Tables:** always with "Table" prefix → `Table A.1`, `Table B.1`, `Table C.1`, … Referenced as `Table B.3`.

Each counter is sequential within its type within each volume. **Never write a bare number** — always include "Equation", "Figure", or "Table" so the type is unambiguous. Cross-volume references use the full label (e.g. "Equation (A.3)" when referenced from Volume B).

---

## Graphs and Figures

### Quality

- Render all figures at **300+ DPI**.
- Export **all figures in both PNG (300 DPI) and SVG**. PNG for raster embedding; SVG for scalable/print quality. Both are generated automatically by `save_figure()` and the appendix save loops.
- All axes must have **labels with physical units** (Hz, s, µV, dB, µV²/Hz, etc.) — not just variable names. Use proper Unicode symbols in axis labels: **µ** (not u), **²** (not ^2), **³** (not ^3). Write `"µV²/Hz"` not `"uV^2/Hz"`.
- Use **perceptually uniform colormaps** (`viridis`, `inferno`, `plasma`, `magma`). **Never use `jet` or `rainbow`** — they distort perception and are not accessible.
- Colors should be high-contrast and distinguishable in grayscale print.
- **Every spectrogram and heatmap must have a colorbar** with labeled units (e.g. "Power (dB)", "Power (linear)"). The scale choice (linear or dB) must be **justified in the report text** — e.g. "dB is used because the burst and chirp differ in amplitude" or "linear is used because both tones have equal amplitude." When multiple spectrograms appear in one figure, use **one shared colorbar** with a shared color range — not one per panel.
- **Window lengths in samples must always state the equivalent duration** — e.g. "$M = 256$ (1.024 s)" or "$M = 1250$ (5.0 s)". Never a bare sample count.

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

---

## Volume C Standards (Real EEG)

All Volume B rules (code format, graph quality, numbering, dual-stack, lab template) apply to Volume C. The following additional standards are specific to real EEG data.

### Use `src/common/` — no exceptions

Volume C code **must** import from `src/common/`, the same infrastructure used in Volume B. Do not redefine constants, rewrite plotting functions, or hardcode parameters.

- **`config.py`** — `FS`, `DPI`, `COLORMAP`, `FIGURE_FORMATS`, `SEED`, and EEG-specific constants (`EEG_BANDS`, `DATASET_PATH`).
- **`plotting.py`** — `plot_time_domain`, `plot_dual_stack_spectrum`, `plot_spectrogram`, `save_figure`. Change axis labels to physical units (µV, µV²/Hz) via function arguments, not by rewriting the functions.
- **`windows.py`** — same windows, same periodic convention, justified by Labs 3 and 5.
- **`eeg.py`** — EEG loading (MNE wrapper), band-power computation (Welch with justified parameters). All EEG I/O goes through this module.

If Volume C needs a new utility, add it to `src/common/`, not to the lab code.

### Primary Channel

- **CZ (vertex)** is the default channel for all single-channel analysis. It sits at the top center of the scalp and picks up activity from both hemispheres with minimal regional bias — the least biased starting point for triage.
- Use other channels only when the analysis specifically requires it (e.g., comparing hemispheres, investigating a regional feature, or examining auxiliary channels for artifact identification).
- When switching from CZ, state which channel and why.
- For auxiliary/non-brain channels (ECG, EMG, EOG), label axes with both the channel name and the suspected signal type — e.g. "25+ (suspected ECG)" — so the reader never has to guess what a channel contains.

### Physical Units

- Amplitude: **µV** — never arbitrary units or "normalized."
- Power spectral density: **µV²/Hz** (linear), then dB.
- Time-frequency: physical units on both axes (s, Hz) and colorbar (µV²/Hz or dB).
- Every axis, every colorbar, every number in the text carries its unit.

### Data-First Structure

- Volume B starts each lab with a formula (math-first). Volume C starts each section with the **data**: what channel, what time range, what the raw signal looks like.
- The tool (DFT, STFT, WVD) comes after — it is applied to the data, not the specification of the data.
- **Appendix C (EEG background)** must be read before C.1. It provides the band definitions, electrode system, and dataset description.

### Adaptive-Directed Analysis

- C.1 (triage) examines the data and decides the direction for C.2–C.5.
- Each subsequent section must explain **why** it was chosen based on what C.1 found — not presented as a predetermined sequence.

### Reproducibility

- State the **dataset source**, file format, sampling rate, channel selection, montage, and any preprocessing (filtering, rereferencing).
- Use **MNE** for all EEG loading, channel selection, and basic preprocessing. Show the MNE calls explicitly in the report.
- Someone with the same dataset and code must be able to reproduce every figure.

### Artifact Handling

- Every analysis section must state whether artifacts were removed, how, and what effect this has.
- If artifacts are NOT removed (e.g., to demonstrate detection in C.4), say so explicitly.

### Parameter Justification

- Band-power analysis (Welch) must justify segment length, window, and overlap with explicit reference to the resolution-variance tradeoff (Lab 2) and the resolution limit (Lab 3).
- No default parameters without justification — state why each choice is appropriate for **this** EEG signal.

### WVD/SPWVD on Selected Segments

- The raw WVD cannot handle the full multi-component noisy EEG record (cross-term soup).
- Select **clean, short segments** for WVD/SPWVD analysis. State the selection criteria (which channel, which time range, why).
- This is the correct use of the tool, not a limitation to apologize for.

### Clinical Language

- "The 10 Hz peak corresponds to alpha rhythm" — fine (established EEG knowledge, cite if needed).
- "This patient has abnormal alpha" — **NOT allowed** (clinical diagnosis, outside our scope).
- Every association between a signal feature and a clinical condition must reference published literature with a citation.
