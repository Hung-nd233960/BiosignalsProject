# TODO - Priority Ranking

## Completed

- [x] **Volume A** - Theory (A.1-A.8 + Appendix A). Exported to LaTeX, zero warnings.
- [x] **Volume B** - Labs 1-8 + Appendix B1 (M vs M-1) + Appendix B2 (CV archetypes). Code and reports complete. Exported to LaTeX.
- [x] **Volume C** - C.1 triage, C.2 stationary DFT, C.3 STFT bursts, C.4 artifacts, C.5 WVD/SPWVD, C.6 synthesis. Code and reports complete. Exported to LaTeX.
- [x] **Labs 7-8** - WVD and SPWVD code (`src/lab7_wvd/`, `src/lab8_spwvd/`), dual-stack figures, reports in Volume B.
- [x] **WVD infrastructure** - `src/common/wvd.py` with `wigner_ville()` and `smoothed_pseudo_wigner_ville()`.
- [x] **Build pipeline** - prettier.py → pandoc → fix_tex.py → tectonic. LaTeX is primary, docx secondary.
- [x] **fix_tex.py** - strips figure floats, fixes media paths, inserts clearpage, shrinks title, fixes Table A.2 widths.
- [x] **LaTeX Workshop** - configured in VS Code with tectonic backend.
- [x] **Directory rework** - `docs/reports/`, `docs/slides/`, `docs/standards/`.
- [x] **Slides** - SHORT.md (~35), MEDIUM.md (~70), LONG.md (~117). Standards and content plan in `docs/standards/`.

## Critical (remaining)

- [ ] **Appendix C** - EEG background (being written by collaborators - integrate when ready)
- [ ] **Final equation/table/figure audit** - renumber across all three volumes after all content is finalized. Check for collisions or gaps.
- [ ] **Final LaTeX build** - rebuild all three volumes, verify PDF output (page breaks, figures, equations).

## Important (strengthens the report)

- [ ] **Proof of orthonormality** of DFT basis vectors (add to A.2.4 or Appendix)
- [ ] **Gaussian window deep dive** - expand A.3.6: minimum uncertainty, connection to SPWVD smoothing (A.8)
- [ ] **Figures for Volume A** - illustrative graphs for A.1 (sampling, aliasing), A.2 (bin spacing visual), A.5 (spectrogram reading guide)

## Nice to Have (scope extensions)

- [ ] **Z-Transformation and Filtering** - would enable bandpass before WVD (the scope boundary identified in C.6)
- [ ] **Wavelets** - brief comparison to STFT/WVD in A.8.5 (Cohen's class). Not a full treatment.
- [ ] **Hilbert-Huang Transform** - empirical mode decomposition as a data-driven alternative.
