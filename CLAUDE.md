# Project Context

We are undergraduate Biomedical Engineering students taking a Digital Biosignal Processing class. This is the final report on signal processing. Our team chose to study the basics of the DFT and STFT, then expand into the Wigner-Ville Distribution (WVD) and Smoothed Pseudo Wigner-Ville Distribution (SPWVD) as the advanced topic to learn and apply. The application target is real EEG data.

The report is split into three volumes (see `docs/standards/TABLE_OF_CONTENTS.md`):

- **Volume A** — Theory (A.1–A.8) + Appendix A (signal taxonomy)
- **Volume B** — Labs (B.1–B.8) with EEG-realism signal constraints
- **Volume C** — Application to real EEG + Appendix C (EEG background)

A unified slide deck accompanies all three volumes.

---

## Tone and Writing Style

- **Academic but honest.** Serious enough for a university report, but written as students explaining concepts we learned — not as textbook authors. No false authority; if something is approximate or heuristic, say so.
- **Neutral.** No hype, no hedging. State what the method does, what it costs, and what we observed.
- **Discrete realm.** All mathematical formulas, derivations, and interpretations are explicitly discrete-time and discrete-frequency. Use summations (not integrals), sample indices (not continuous time), and normalized/discrete frequency. The continuous-time form may be mentioned for motivation, but the working formula is always the discrete version.
- **No "gap callouts" in the reports.** The gap framing ("Signals & Systems gives you X, this section adds Y") is internal planning language in `docs/standards/TABLE_OF_CONTENTS.md`. Do NOT include it in the actual volumes (`docs/reports/volume_A.md`, `volume_B.md`, `volume_C.md`) or slides. The reports teach the material directly without referencing what the reader is assumed to be missing.
- **No `>` blockquotes in submitted reports.** Blockquotes (`> text`) are internal editorial notes and planning reminders. They must be removed before submission. The `prettier.py` script strips them automatically.

---

## Environment

- **Conda (miniforge3).** All code runs inside the `biosignals` conda environment defined in `environment.yml` at the project root.
- **Setup from scratch:** `conda env create -f environment.yml` then activate as above.
- **Adding dependencies:** add to `environment.yml`, then `conda env update -f environment.yml --prune`. Never `pip install` directly — keep the environment file as the single source of truth.
- **Core stack:** Python (≥3.11), NumPy, SciPy, Matplotlib, JupyterLab, MNE (EEG analysis for Volume C). R and rpy2 are included for collaborators who prefer R.

---

## Reference Files

- `docs/standards/TABLE_OF_CONTENTS.md` — master plan with all section outlines, lab templates, and internal planning notes (gap callouts). This is the source of truth for structure.
- `CLAUDE.md` — this file. Project standards, conventions, and rules for all contributors (human and AI).
- `src/README.md` — code contributor guide: directory layout, import rules, naming conventions.
- `environment.yml` — conda environment definition. Single source of truth for dependencies.
- `template/reference.docx` — Word template for pandoc: Roboto body, Consolas code, margins (top 2cm, left 2cm, right 1cm, bottom 1cm), page numbers bottom center.
- `docs/METHODOLOGY.md` — records every rule in this file, why it was added, and what problem it solved. **When adding or changing a rule in CLAUDE.md, update METHODOLOGY.md with the reason.**
- `docs/standards/DOCUMENT_STANDARD.md` — LaTeX/PDF formatting standard: page layout, figure handling (no floats), equation tags, numbering, build pipeline, and known issues. **All LaTeX output must conform to this file.**

---

## Building Documents

**The primary output format is LaTeX** (for Overleaf). Word (.docx) is secondary.

**Step 1: Touch up markdown** (removes horizontal rules, blockquotes, fancy dashes):

```bash
python build_docs/prettier.py docs/reports/volume_B.md
```

**Step 2a: Markdown to LaTeX** (primary) -- with embedded images for Overleaf:

```bash
~/miniforge3/bin/pandoc docs/reports/volume_B.md \
    -o output/volume_B.tex \
    --from=markdown+tex_math_dollars+pipe_tables --standalone \
    -V geometry:"top=2cm, bottom=2cm, left=2cm, right=1cm" \
    --extract-media=output/media_B \
    --resource-path=docs/reports

# Fix image paths and remove floating figures (prevents double-numbering)
python3 build_docs/fix_tex.py output/volume_B.tex
```

Upload to Overleaf: `volume_B.tex` + the `media_B/` folder (side by side). Images are extracted as hashed PNGs -- no external paths needed. The `fix_tex.py` script strips the `output/` prefix from media paths and replaces `\begin{figure}` floats with inline `\begin{center}` blocks so our manual "Figure B.X" labels don't conflict with LaTeX's auto-numbering.

**Step 2b: Markdown to Word (.docx)** (secondary):

```bash
~/miniforge3/bin/pandoc docs/reports/volume_B.md \
    -o output/volume_B.docx \
    --from=markdown+tex_math_dollars+pipe_tables \
    --mathml \
    --reference-doc=template/reference.docx \
    --resource-path=docs/reports
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
- **Time-domain first — no exceptions.** Every signal must have its time-domain waveform plotted and inspected before any spectral or time-frequency analysis. This applies to all signals: tones, chirps, noise, transients, real EEG, and even pure white noise. A noise signal that looks clipped, saturated, or has unexpected structure would corrupt every downstream result. The time-domain plot is the sanity check — if the signal looks wrong here, nothing after it matters. No signal enters the DFT, STFT, or WVD without being seen in the time domain first.

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

### Logarithmic Scale Labeling

**"dB" versus "10 log₁₀(·)" — case by case.** The label "dB" implies a defined physical reference. Use it only when the reference exists and is stated. Otherwise, label the axis as the mathematical transformation it is.

| Plot type | Reference? | Axis label |
| --- | --- | --- |
| PSD of real EEG (µV²/Hz) | Yes (1 µV²/Hz) | **dB re 1 µV²/Hz** |
| Window spectrum (Lab 3, normalized to peak) | Yes (peak = 0 dB) | **dB (relative to peak)** |
| PSD of model signal (no physical units) | No | **10 log₁₀(power)** |
| WVD / SPWVD (can be negative) | No | **10 log₁₀(\|WVD\|)** |

"dB re 1 µV²/Hz" is a legitimate, defensible unit — same convention as dB SPL (re 20 µPa) in acoustics. The "re" defines the reference; the number converts back to physical units unambiguously.

**Dimensionless signals (Volume B model signals) should not use "dB."** The log transformation is still useful — it compresses dynamic range so a -60 dB side lobe is visible next to a 0 dB main lobe. The information is the same. But the label "dB" implies a physical reference that dimensionless signals do not have. Use **10 log₁₀(·)** instead. Same plot, same numbers, honest label.

### The Dual-Stack Rule

Both linear and log-scale plots are produced for every spectral or time-frequency figure. They are rendered as **separate figures**, not stacked in one panel.

- The **linear-scale plot is primary** — it shows the physical quantity in its real units and is always presented first. It is always numbered and included in the report.
- The **log-scale plot is secondary** — it reveals dynamic range and low-level structure. It is numbered and included **only when it reveals something the linear plot hides** (e.g. burst structure in C.3, side-lobe decay in Lab 3). Otherwise it stays in `results/graphs/` as an unnumbered draft.

**When to include the log-scale plot:**

- C.3 full spectrogram — burst structure invisible in linear, obvious in log. **Include.**
- Lab 3 side-lobe analysis — decay rate measured in dB/octave. **Include.**
- C.2 log-log PSD — 1/f slope is a log-log measurement. **Include.**
- Lab 1 on-bin/off-bin — linear shows everything. **Draft only.**
- WVD/SPWVD — log panel shows misleading circular artifacts. **Draft only, or include with explicit caveat.**

For real EEG signals (Volume C):

- Amplitude in **µV** (linear).
- Power spectral density in **µV²/Hz** (linear), then **dB re 1 µV²/Hz** when justified.
- Time-frequency representations in physical units (linear), then log scale when justified.

For model signals (Volume B), the same principle applies. Linear scale grounds the reader in the actual signal; log scale is included only when it adds information.

---

## PSD Density versus Band Power

**PSD (µV²/Hz) is a density, not energy.** The PSD plot shows how power is distributed across frequency — its shape reveals peaks, slopes, and spectral structure. But comparing the curve height at two frequencies does not tell you which band has more total energy, because bands have different widths.

**To compare energy across bands, integrate:**

$$\text{band power} = \sum_{k \in \text{band}} \text{PSD}[k] \cdot \Delta f \qquad \text{(µV²/Hz × Hz = µV²)}$$

The result is in µV² — actual energy, not density. This is what `compute_band_power()` does. The bar chart (C.1 Figure C.6) shows integrated band power, not PSD heights.

**The chain in the report:**

1. **PSD plot** (µV²/Hz) → for looking at spectral shape (where are peaks, what is the slope)
2. **Integrate over band** → µV² (energy per band)
3. **Divide by total** → % (relative energy distribution: "91.8% delta")

Step 1 is for shape. Steps 2-3 are for energy comparison. Both are needed; neither replaces the other. When discussing energy distribution in the report text, always reference the integrated values, not the PSD curve height.

**scipy's `spectrogram` uses PSD scaling by default** (`scaling='density'`). This matches the Welch convention. The alternative (`scaling='spectrum'`, µV² per bin) shows the identical shape with a different normalization. The choice only matters when reading numbers off the plot — density requires integration, spectrum requires summation. We use density throughout for consistency with `compute_band_power()`.

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
