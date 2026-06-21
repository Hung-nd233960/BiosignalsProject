# Document Standard — LaTeX Output

This file defines how the final PDF should look. All volumes (A, B, C) follow the same standard. The pipeline is: Markdown → pandoc → fix_tex.py → tectonic → PDF.

## Page Layout

- **Paper:** A4
- **Margins:** top 2 cm, bottom 2 cm, left 2 cm, right 1 cm
- **Page numbers:** bottom center (automatic via LaTeX)
- **Font:** LaTeX default (Computer Modern). Pandoc standalone template handles this.

## Title Block

Every volume starts with:

- Report title: "From the DFT to the SPWVD: Time-Frequency Analysis Applied to Neonatal EEG"
- Volume subtitle: e.g. "Volume B — Application and Derivation (Labs)"
- Author: Nguyen Duc Hung — 20233960

## Figures

- **No floating figures.** Pandoc wraps images in `\begin{figure}` which auto-numbers them ("Figure 1:", "Figure 2:"). This conflicts with our manual numbering ("Figure B.1", "Figure B.39"). The `fix_tex.py` script replaces all `\begin{figure}` environments with `\begin{center}` blocks.
- **Captions are bold text** below the image, using our manual labels: `\textbf{Figure B.1 - Description}`.
- **Images are inline** — they appear exactly where the markdown places them, not repositioned by LaTeX's float algorithm.
- **Image format:** PNG (300 DPI), extracted by pandoc's `--extract-media` into a `media_X/` folder (where X = A, B, or C).
- **No SVG in LaTeX.** Pandoc extracts PNGs only. SVGs are kept in `results/graphs/` for other uses.

## Equations

- **Numbered with `\tag{}`** — e.g. `\tag{A.15}`, `\tag{B.3}`. These are our manual tags, not LaTeX's `\label/\ref` system.
- **No inline units inside equations.** Patterns like `\quad \text{(Hz)} \tag{A.15}` break pandoc's MathML parser (for docx) and clutter the equation. State units in the surrounding prose instead.
- **Display math only for tagged equations.** Inline math (`$...$`) for expressions in prose; display math (`$$...$$`) for numbered equations.
- **Fallback for stubborn tags:** if pandoc chokes on `\tag{X.Y}` in a specific equation, use the workaround: `\qquad \text{(X.Y)}` instead of `\tag{X.Y}`.

## Tables

- **Pipe tables in markdown** — pandoc converts them to `\begin{longtable}`. These do NOT conflict with our numbering (no auto-numbering for tables in the pandoc template).
- **Table labels** are bold text above the table: `**Table B.25 — Lab 7 parameters**`.
- **No floating tables.** Pandoc's default longtable is inline, which is correct.

## Code Blocks

- **Fenced code blocks** (triple backticks with `python` language tag) → pandoc converts to `\begin{Shaded}\begin{Highlighting}` with syntax coloring.
- **Every figure has its code block immediately before it.** Pattern: code → figure → interpretation.
- Code blocks should not break across pages in ways that separate the code from its figure. If this happens, it's a manual fix in the tex.

## Numbering (Three Independent Counters)

Each volume has three independent sequential counters:

- **Equations:** `(A.1)`, `(B.1)`, `(C.1)` — in parentheses, referenced as "Equation (B.3)"
- **Figures:** `Figure A.1`, `Figure B.1`, `Figure C.1` — referenced as "Figure B.3"
- **Tables:** `Table A.1`, `Table B.1`, `Table C.1` — referenced as "Table B.3"

These are all manual — we do NOT use LaTeX's `\label`/`\ref` system because pandoc-generated tex doesn't support cross-references cleanly.

## Build Pipeline

```bash
# Step 1: Clean markdown
python3 build_docs/prettier.py docs/reports/volume_B.md

# Step 2: Markdown → LaTeX (with embedded images)
~/miniforge3/bin/pandoc docs/reports/volume_B.md \
    -o output/volume_B.tex \
    --from=markdown+tex_math_dollars+pipe_tables --standalone \
    -V geometry:"top=2cm, bottom=2cm, left=2cm, right=1cm" \
    --extract-media=output/media_B \
    --resource-path=docs/reports

# Step 3: Fix figures + image paths
python3 build_docs/fix_tex.py output/volume_B.tex

# Step 4: LaTeX → PDF
cd output && ~/miniforge3/bin/tectonic volume_B.tex
```

For Volume A: substitute `volume_A`, `media_A`. For Volume C: substitute `volume_C`, `media_C`.

## Known Issues

- **`\pandocbounded` macro:** pandoc 3.x emits this wrapper around `\includegraphics`. The pandoc standalone template defines it. If using a custom template, ensure it includes the definition.
- **Long tables:** tables wider than the text width will overflow. Keep parameter tables to 3-4 columns.
- **Page breaks around figures:** LaTeX may insert page breaks between a code block and its figure. This is cosmetic and can be fixed with manual `\newpage` or `\clearpage` in the tex if needed.
- **Overleaf compile timeout:** 54+ images can exceed Overleaf's free-tier compile time. Use local tectonic instead, or Overleaf's paid tier.
