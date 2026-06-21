# Methodology - How This Report Was Built

This document records the standards and conventions developed during the making of this project. What started as a few basic rules grew into a comprehensive framework as we encountered problems, received instructor feedback, and discovered what worked and what didn't. Each rule exists because something went wrong without it.

## Writing Standards

**Academic but honest tone.** We are students explaining what we learned, not textbook authors. If something is approximate, we say so. If a method fails, we report it.

**Discrete realm only.** All formulas use summations, sample indices, and discrete frequency. Continuous-time forms appear only for motivation - the working formula is always discrete. This decision was made early and enforced throughout.

**No gap callouts in reports.** The table of contents uses "Gap: Signals & Systems gives you X, this section adds Y" as internal planning language. The actual reports teach the material directly. The gap framing was useful for planning but disrespectful to the reader when left in.

**No blockquotes in submitted reports.** The `>` markdown syntax is used for internal editorial notes during drafting. The `prettier.py` script strips them before conversion to Word/PDF.

## Code Standards

**Show your tools.** Every code block displays its imports. If a library function is used (e.g. `scipy.stats.linregress`, `scipy.signal.argrelextrema`), the report states what it does and what method it uses. "Fit a line" is not acceptable - "ordinary least squares via `linregress`" is. This rule was added after recognizing that vague descriptions like "fit a line" or "find maxima" hide the actual computation from the reader.

**Comment every line.** Inline comments explain what each line does. This was a founding rule.

**Single source of truth (`src/common/`).** Constants, signal generators, window functions, plotting utilities, and EEG loading all live in `src/common/`. Labs and Volume C import from there - never redefine `fs = 250` or rewrite the dual-stack plotting function. This was established before writing any lab code and prevented parameter drift across 8 labs.

**Fixed random seed.** Every noise signal uses `SEED = 42` from `config.py`. Results are reproducible across machines and runs.

## Signal Constraints (EEG Realism)

All model signals in Volume B satisfy constraints that make them representative of real EEG:
- Frequencies below 100 Hz
- Duration at least 1200 s (20 minutes)
- Sampling at 250 Hz (default)

These constraints were added at the start to ensure that bin spacing, resolution limits, and window choices encountered in Volume B transfer directly to the real EEG analysis in Volume C. Without them, the labs would demonstrate concepts at unrealistic scales that don't apply to the application.

## Numbering System

**Three independent counters per volume.** Equation (B.1), Figure B.1, and Table B.1 are three different objects. Each counter is sequential within its type within each volume. The rule: never write a bare number - always include "Equation", "Figure", or "Table".

This system was adopted after an early design where equations and figures shared the same counter, causing confusion. The independent counters were modeled on standard academic conventions (IEEE, APA).

**Every table must be numbered.** This rule was enforced retroactively after discovering 15 unlabeled tables in Volume B. A renumbering script was written to fix them.

## Graph Quality

**300+ DPI, perceptually uniform colormaps.** `viridis`, `inferno`, `plasma`, `magma` only. `jet` and `rainbow` are explicitly forbidden - they distort perception and are inaccessible to colorblind readers.

**Proper Unicode in axis labels.** Write µV²/Hz, not uV^2/Hz. This was added after the instructor noted crude labels on early Volume C figures.

**The Dual-Stack Rule.** Every spectral or time-frequency plot appears as a dual stack: linear scale first (physical units), dB second (dynamic range). Linear is primary; dB is supplementary. This was a founding rule, strengthened for Volume C where physical units (µV, µV²/Hz) are non-negotiable.

**Scale choice must be justified.** "dB is used because the burst and chirp differ in amplitude" or "linear is used because both tones have equal amplitude." Added after Lab 5 where we found linear scale was cleaner for equal-amplitude tones.

**Colorbars mandatory on spectrograms.** With labeled units. When multiple spectrograms appear in one figure, one shared colorbar with a shared color range - not one per panel.

**Window lengths always include duration.** Write "$M = 256$ (1.024 s)", never just "$M = 256$".

**SVG export for line plots, PNG-only for spectrograms.** SVG of dense pcolormesh data produces massive files that crash renderers. This was discovered when Lab 4's spectrogram SVG export crashed the system.

## Lab Report Structure

Every Volume B lab follows a fixed six-section template:
1. Introduction - gloss theory, state what the lab tests
2. Setup - formula first, then code
3. Parameters - one table, full reproducibility at a glance
4. Results - time-domain first, then spectral/dual-stack, with figures embedded inline
5. Verification - predicted vs. measured, quantitative
6. Conclusion - what was confirmed, bridge to next lab

This template was designed early and every lab follows it. The verification section is the key differentiator - it forces quantitative comparison against Volume A predictions, not just "the plot looks right."

**Code must be explicit in reports.** The instructor does not go to GitHub. All important functions, parameters, and calls appear inline in the markdown report. This rule was added after the first draft showed only function names without implementations.

## The M vs M-1 Convention

The cosine-sum windows (Hann, Hamming, Blackman) can use either $M$ (periodic) or $M-1$ (symmetric) in the denominator. This distinction consumed significant effort:

- **Volume A** documents both conventions with a library comparison table (numpy, scipy, MATLAB defaults)
- **Volume B Lab 3** derives everything using the periodic convention ($M$) because it enables the shared-numerator factorization in the geometric series
- **Appendix B** provides the symmetric ($M-1$) derivation and proves both converge as $M \to \infty$
- **The practical difference at $M = 256$ is 0.11 dB** - invisible in any measurement

The lesson: when two valid conventions exist, pick one, document why, show they agree, and move on.

## Volume C Standards

Volume C applies the same tools to real EEG data. Additional standards:

- **Data-first, not math-first.** Each section starts with the data, then applies the tool.
- **Primary channel: CZ.** Vertex electrode, least regional bias. Other channels only when justified.
- **Physical units mandatory.** µV, µV²/Hz - never arbitrary or normalized.
- **Parameter justification.** Every Welch segment length, window choice, and overlap must reference specific Volume B labs.
- **Adaptive-directed.** C.1 triage decides the analysis direction. Subsequent sections justify themselves from C.1's findings.
- **Clinical claims cited, not asserted.** We are engineering students, not clinicians.
- **WVD/SPWVD on selected segments only.** Full-record WVD is cross-term soup. Selecting clean segments is the correct use of the tool.

## Code-Figure Traceability

**Every figure must have its code in the report.** This rule was added when Volume C's early drafts showed figures without the code that produced them. The instructor's concern: "how do I know this is what you actually computed?" Without the code, the figure is an unsupported claim. The pattern enforced is: **code block → figure → interpretation.** No figure appears without the reader seeing exactly how it was generated. This applies to both Volume B and Volume C.

## Decision-Driving Numbers

**Numbers that drive analysis decisions must be traceable to code.** Not every number needs a code block — descriptive statistics like "amplitude range ±230 µV" are fine as prose. But numbers that determine what happens next — "91.8% delta" (drives the triage), "1/f slope = -3.18" (determines if delta is signal or noise), "excess ratio = 0.2x" (the verdict) — must have the code that produced them shown immediately before. The critical chain is: **code → number → decision.** This was added to distinguish between aesthetic figures (which need code for traceability) and quantitative claims (which need code for auditability). Both need code, but for different reasons.

## Self-Documenting Standards

**METHODOLOGY.md tracks every rule change.** When a rule is added or changed in CLAUDE.md, the reason is recorded in this file. This was added to prevent the standards from becoming an opaque list of commands - the "why" behind each rule is preserved so future contributors (and future sessions) understand the reasoning, not just the conclusion.

## Build Pipeline

**LaTeX is the primary output format.** The project started with Word (.docx) via pandoc + MathML. However, pandoc's MathML renderer choked on complex equations with `\tag{}` inside `\text{}` blocks, producing 13+ warnings on Volume A alone. LaTeX handles `\tag{}` natively with no issues. Additionally, Overleaf provides collaborative editing, proper typesetting, and PDF export without local TeX installation. The `--extract-media` flag embeds images as hashed PNGs in a `media_X/` folder, making Overleaf upload simple: one .tex file + one image folder.

**Word remains secondary** for quick previews and instructor compatibility (some prefer .docx). The `--resource-path=docs` flag resolves relative image paths from the markdown's perspective.

## fix_tex.py Post-Processing

Pandoc's LaTeX output requires several fixes before it compiles cleanly:

1. **Figure floats removed.** Pandoc wraps images in `\begin{figure}` which auto-numbers them ("Figure 1:"), conflicting with our manual "Figure B.1" labels. `fix_tex.py` replaces all `\begin{figure}` environments with `\begin{center}` blocks.
2. **Media paths stripped.** Pandoc's `--extract-media` produces paths like `output/media_B/hash.png`. The script strips the `output/` prefix so Overleaf sees `media_B/hash.png`.
3. **Clearpage before sections.** Each lab and appendix starts on a new page. The script inserts `\clearpage` before `\subsection{A.\d`, `\subsection{B.\d`, `\subsection{C.\d`, `\section{Appendix`, and `\subsection{Shared Infrastructure`.
4. **Title shrunk.** The report title is too long for `\section` font size. Replaced with `\Large\bfseries` centered block.
5. **Table A.2 column widths.** The 4-column library defaults table was cramped at 25/25/25/25. Changed to 35/25/20/20 for the wider first column.

Each fix was added when a specific problem was encountered. The script is idempotent when run on fresh pandoc output but should NOT be run twice on the same file (clearpage would double).

## WVD Colorbar Labels

The WVD is a quasi-distribution that can produce **negative values** (unlike the STFT power spectrogram which is always non-negative). Colorbars on WVD/SPWVD plots say "WVD value" (linear) and "|WVD| (dB)" rather than "Power" to avoid implying the values are physical power. STFT colorbars correctly say "Power (linear)" and "Power (dB)".

## Segment Selection for WVD on Real EEG

The WVD/SPWVD cannot handle the full multi-component noisy EEG record. The prescribed approach: select clean, short segments using C.3's burst detection (2x median delta power threshold). The selection process is data-driven and includes honest failure reporting: the strongest burst was amplifier saturation (44 flat samples) and was rejected; the 75th percentile burst was clean and was accepted. This methodology is shown transparently in C.5.1, not hidden.

No filtering was applied because Volumes A-B did not derive filter theory. This is a scope boundary, not a failure.

## Slide Standards

Three slide deck versions (SHORT ~35, MEDIUM ~70, LONG ~117) in `docs/slides/`. Standards in `docs/standards/SLIDES_STANDARD.md`:

- **Title Case** on all slide titles (e.g. "Lab 3: Window Spectra from First Principles").
- **Formal language.** Complete sentences in all bullets, no slang, no casual phrasing. Slides are the public face of the report.
- **Purpose → Result structure for labs.** Each lab's slides lead with the question being asked, then show the answer (the key figure). The figure is the star, not the setup.
- **Content density:** 3-5 bullets, 1 figure max, 1-2 equations max, minimum ~20 words per slide.

## Directory Structure

Reorganized `docs/` into three subdirectories:

- `docs/reports/` - Volume A, B, C markdown source
- `docs/slides/` - SHORT.md, MEDIUM.md, LONG.md
- `docs/standards/` - TABLE_OF_CONTENTS, DOCUMENT_STANDARD, SLIDES_STANDARD, SLIDES_CONTENT

This was done because the flat `docs/` directory had 10+ files mixing reports, slides, and standards. The subdirectories make the purpose of each file obvious. Image paths in reports changed from `../results/` to `../../results/` accordingly.

## What We Would Do Differently

- **Number tables from the start.** Retroactive numbering of 15 tables was painful and error-prone.
- **Settle the M vs M-1 convention earlier.** The derivation rework cost significant time.
- **Write `prettier.py` before writing any reports.** The horizontal rules and fancy dashes accumulated and required batch cleanup.
- **Show code for every figure from the first draft.** Adding code retroactively to Volume C required reworking sections that were already "done."
