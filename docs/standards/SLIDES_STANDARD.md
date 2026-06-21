# Slides Standard

Rules for the presentation slide deck accompanying all three volumes.

## Format

- **Tool:** Marp (Markdown to HTML/PDF/PPTX)
- **Theme:** default or uncover
- **Math:** LaTeX via `$...$` and `$$...$$`
- **Images:** embedded from `results/graphs/` relative paths
- **Export:** HTML (primary), PDF (print), PPTX (fallback)
- **Versions:** three lengths in `docs/slides/`:
  - `SHORT.md` (~35 slides) - highlights only, one slide per concept/lab
  - `MEDIUM.md` (~70 slides) - expanded, two slides per lab (setup + result)
  - `LONG.md` (~100-120 slides) - comprehensive teaching deck, full walkthrough

## Content Density

### Text
- **3-5 bullet points per slide.** Each bullet is one sentence. No paragraphs.
- **Minimum ~20 words per slide.** A slide with only a title or 3-5 words is a waste to print. Combine short headers with their content.
- **No wall of text.** If a slide has more than 6 bullets or more than 50 words of body text, split it.

### Figures
- **1 figure per slide maximum.** The figure occupies ~60% of the slide area; bullets occupy ~40% below.
- **Every figure slide must have 2-3 observation bullets** — what the figure shows, not just a caption.
- **No slides that are just a figure with no text** — the audience needs to know what to look at.

### Equations
- **1-2 equations per slide maximum.** Show the final result, not the derivation.
- **Equations are centered, large font.** No inline clutter around them.
- **Cite the source** — e.g. "Equation (A.61)" — so the reader can find the derivation in the report.
- **Derivations live in the report, not on slides.** If a derivation is essential for the talk, spread it across 2-3 slides with one step per slide. Never show a full derivation on one slide.

### Code
- **No code on slides.** Code is in the report (Volume B, Volume C). The slides reference it but do not reproduce it.
- **Exception:** a 2-3 line pseudocode snippet is acceptable if it clarifies the algorithm (e.g. the WVD inner loop), but only if it fits naturally.

## Slide Titles

- **Titles state what the slide proves or shows**, not just a section label.
  - Good: "SPWVD suppresses cross-terms independently per axis"
  - Bad: "Lab 8"
  - Bad: "Results"
- **Title Case after section or lab labels.** Capitalize each major word in the title.
  - Good: "Lab 3: Window Spectra from First Principles"
  - Bad: "Lab 3: Window spectra from first principles"
  - Good: "C.5: Segment Selection - The Honest Story"
- **Titles are one line.** If the title wraps, shorten it.

## Language and Tone

- **Extra formal.** Slides are the public face of the report. Language must be precise and professional - no casual phrasing, no abbreviations without definition, no sentence fragments.
  - Good: "The SPWVD suppresses both time-oscillating and frequency-oscillating cross-terms."
  - Bad: "SPWVD kills both ghost types."
- **Complete sentences in bullets.** Every bullet is a grammatically complete sentence, not a note or keyword.
- **No slang, no hedging, no hype.** State facts. "The raw WVD produces 49% negative values" not "The WVD completely fails."

## Title Positioning

- **Titles are positioned at the top of the slide** with extra spacing above. In Marp, use a custom style or `<!-- _class: lead -->` for emphasis slides. For regular slides, the default `## Title` heading provides adequate top positioning. If the title appears too low, add a Marp directive:
  ```
  <style>section h2 { margin-top: 0; }</style>
  ```

## Narrative Rules

- **Slides are not the report again.** The reader has Volumes A-C. Slides show highlights and the narrative thread.
- **Each slide must advance the story.** If removing a slide doesn't break the narrative, the slide is unnecessary.
- **The story arc:** problem (EEG is non-stationary) → tools (DFT → STFT → WVD → SPWVD) → application (real neonatal EEG) → honest result (what worked, what didn't).
- **No gap callouts.** Same rule as the reports — teach the material, don't reference what the reader is assumed to be missing.
- **No blockquotes.** Same as reports.

## Figure References

- Figures on slides are the same figures from the report. Reference them by number: "Figure B.40", "Figure C.27".
- Do not create new figures for slides. If a slide needs a visual that doesn't exist, create it in the lab code first (reproducibility), then reference it.

## Printing

- Slides will be printed. Design for readability on paper:
  - High contrast (dark text on light background)
  - No gradient backgrounds
  - Figures must be legible at slide-on-paper size (~17 cm wide)
  - Avoid color-only distinctions (some printers are grayscale)
