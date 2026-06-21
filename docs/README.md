# docs/ Directory Guide

## Structure

```
docs/
  reports/       — The three report volumes (markdown source)
  slides/        — Presentation slide decks (Marp markdown)
  standards/     — Formatting rules, content plans, and the master table of contents
  METHODOLOGY.md — Why each rule exists (history of standards evolution)
  OPINIONS.md    — Resolved design decisions
```

## Reports (`reports/`)

| File | Content |
|------|---------|
| `volume_A.md` | Theory: DFT, windowing, STFT, autocorrelation, WVD, SPWVD + Appendix A (signal taxonomy) |
| `volume_B.md` | 8 labs + Appendix B1 (M vs M-1 convention) + Appendix B2 (CV of archetypes) |
| `volume_C.md` | Application to neonatal EEG: triage, DFT, STFT, artifacts, WVD/SPWVD, synthesis |

## Slides (`slides/`)

| File | Content |
|------|---------|
| `slides.md` | Marp slide deck covering all three volumes (~32 slides) |

## Standards (`standards/`)

| File | Purpose |
|------|---------|
| `TABLE_OF_CONTENTS.md` | Master plan: all section outlines, lab templates, internal planning notes |
| `DOCUMENT_STANDARD.md` | LaTeX/PDF formatting: page layout, figure handling, equation tags, build pipeline |
| `SLIDES_STANDARD.md` | Slide rules: content density, formatting, narrative structure |
| `SLIDES_CONTENT.md` | Slide-by-slide content plan: what each slide shows and why |

## Building

All builds run from the **project root** (`Biosignals_Project/`), not from `docs/`.

```bash
# Reports: markdown → LaTeX → PDF
python3 build_docs/prettier.py docs/reports/volume_B.md
~/miniforge3/bin/pandoc docs/reports/volume_B.md \
    -o output/volume_B.tex \
    --from=markdown+tex_math_dollars+pipe_tables --standalone \
    -V geometry:"top=2cm, bottom=2cm, left=2cm, right=1cm" \
    --extract-media=output/media_B \
    --resource-path=docs/reports
python3 build_docs/fix_tex.py output/volume_B.tex
cd output && ~/miniforge3/bin/tectonic volume_B.tex

# Slides: markdown → HTML
npx @marp-team/marp-cli docs/slides/slides.md -o output/slides.html --allow-local-files
```

Substitute `volume_A`/`media_A` or `volume_C`/`media_C` for the other volumes.

See `CLAUDE.md` (project root) for the full standards and `standards/DOCUMENT_STANDARD.md` for LaTeX-specific rules.
