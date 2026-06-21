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

## Build Pipeline

Reports are written in markdown, then converted:
- **prettier.py** - strips horizontal rules, blockquotes, and fancy dashes before conversion
- **pandoc** - markdown to Word (.docx) with a reference template (Roboto body, Consolas code, margins 2/2/1/1 cm, page numbers bottom center)
- **pandoc + tectonic** - markdown to LaTeX to PDF
- **Marp** - markdown to HTML slides

The build pipeline was set up mid-project when we needed to deliver formatted documents. The key insight: write in markdown (diffable, versionable), convert for submission.

## Code-Figure Traceability

**Every figure must have its code in the report.** This rule was added when Volume C's early drafts showed figures without the code that produced them. The instructor's concern: "how do I know this is what you actually computed?" Without the code, the figure is an unsupported claim. The pattern enforced is: **code block → figure → interpretation.** No figure appears without the reader seeing exactly how it was generated. This applies to both Volume B and Volume C.

## Decision-Driving Numbers

**Numbers that drive analysis decisions must be traceable to code.** Not every number needs a code block — descriptive statistics like "amplitude range ±230 µV" are fine as prose. But numbers that determine what happens next — "91.8% delta" (drives the triage), "1/f slope = -3.18" (determines if delta is signal or noise), "excess ratio = 0.2x" (the verdict) — must have the code that produced them shown immediately before. The critical chain is: **code → number → decision.** This was added to distinguish between aesthetic figures (which need code for traceability) and quantitative claims (which need code for auditability). Both need code, but for different reasons.

## Self-Documenting Standards

**METHODOLOGY.md tracks every rule change.** When a rule is added or changed in CLAUDE.md, the reason is recorded in this file. This was added to prevent the standards from becoming an opaque list of commands - the "why" behind each rule is preserved so future contributors (and future sessions) understand the reasoning, not just the conclusion.

## What We Would Do Differently

- **Number tables from the start.** Retroactive numbering of 15 tables was painful and error-prone.
- **Settle the M vs M-1 convention earlier.** The derivation rework cost significant time.
- **Write `prettier.py` before writing any reports.** The horizontal rules and fancy dashes accumulated and required batch cleanup.
- **Show code for every figure from the first draft.** Adding code retroactively to Volume C required reworking sections that were already "done."
