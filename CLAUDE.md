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
- `docs/DOCUMENT_STANDARD.md` — LaTeX/PDF formatting standard: page layout, figure handling (no floats), equation tags, numbering, build pipeline, and known issues. **All LaTeX output must conform to this file.**

---

## Building Documents

**The primary output format is LaTeX** (for Overleaf). Word (.docx) is secondary.

**Step 1: Touch up markdown** (removes horizontal rules, blockquotes, fancy dashes):

```bash
python build_docs/prettier.py docs/volume_B.md
```

**Step 2a: Markdown to LaTeX** (primary) -- with embedded images for Overleaf:

```bash
~/miniforge3/bin/pandoc docs/volume_B.md \
    -o output/volume_B.tex \
    --from=markdown+tex_math_dollars+pipe_tables --standalone \
    -V geometry:"top=2cm, bottom=2cm, left=2cm, right=1cm" \
    --extract-media=output/media_B \
    --resource-path=docs

# Fix image paths and remove floating figures (prevents double-numbering)
python3 build_docs/fix_tex.py output/volume_B.tex
```

Upload to Overleaf: `volume_B.tex` + the `media_B/` folder (side by side). Images are extracted as hashed PNGs -- no external paths needed. The `fix_tex.py` script strips the `output/` prefix from media paths and replaces `\begin{figure}` floats with inline `\begin{center}` blocks so our manual "Figure B.X" labels don't conflict with LaTeX's auto-numbering.

**Step 2b: Markdown to Word (.docx)** (secondary):

```bash
~/miniforge3/bin/pandoc docs/volume_B.md \
    -o output/volume_B.docx \
    --from=markdown+tex_math_dollars+pipe_tables \
    --mathml \
    --reference-doc=template/reference.docx \
    --resource-path=docs
```

**Step 2c: LaTeX to PDF** locally via tectonic:

```bash
cd output && ~/miniforge3/bin/tectonic volume_B.tex

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
