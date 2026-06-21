# TODO - Priority Ranking

## Critical (must complete for a valid report)

- [ ] **Appendix C** - EEG background (being written by someone else - integrate when ready)
- [ ] **C.1 Triage** - load dataset, time-domain plot, broad DFT, band power, triage decision
- [ ] **C.2 Stationary characterization** - DFT + Welch band power with justified parameters
- [ ] **C.3 Time-varying characterization** - STFT spectrogram, Heisenberg choice for this signal
- [ ] **C.4 Events and artifacts** - transient detection, artifact identification (channels 25+, 26+, 27+)
- [ ] **C.5 WVD/SPWVD** - selected clean segment, honest comparison to STFT
- [ ] **C.6 Synthesis** - which archetype-tool pairings held up, honest limits
- [ ] **Volume B Labs 7-8** - WVD and SPWVD code + reports (needed before C.5)
- [ ] **Volume A sections A.7-A.8** - WVD and SPWVD theory (needed before Labs 7-8)
- [ ] **Unified slides** - extend slides.md to cover Volumes B and C

## Important (strengthens the report)

- [ ] **Proof of orthonormality** of DFT basis vectors (add to A.2.4 or Appendix)
- [ ] **Gaussian window deep dive** - expand A.3.6 into a full treatment: minimum uncertainty, connection to SPWVD smoothing (A.8)
- [ ] **Figures for Volume A** - illustrative graphs for A.1 (sampling, aliasing, unit circle), A.2 (bin spacing visual), A.5 (spectrogram reading guide)
- [ ] **Lab 3 confirmation** - numerically verify M vs M-1 convergence (referenced in Appendix B but not yet computed)

## Nice to have (scope extensions)

- [ ] **Wavelets** - brief comparison to STFT/WVD in A.8.5 (Cohen's class closing frame). Not a full treatment - just note as a future direction.
- [ ] **Hilbert-Huang Transform** - empirical mode decomposition as a data-driven alternative. Mention in C.6 synthesis if relevant to EEG findings.
