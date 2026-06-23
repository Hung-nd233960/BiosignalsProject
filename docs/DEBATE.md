# DEBATE.md — Open Questions and Unsettled Arguments

Things to try, test, or resolve before finalizing. Each entry has the question, the arguments for each side, and what experiment or evidence would settle it.

---

## 1. Spectrogram Visual Smoothness

**The argument:** "Your sampling rate is high, why is the spectrogram not smooth enough?"

**Our position:** The spectrogram's coarseness reflects the actual time-frequency resolution of the STFT at the chosen window length. The visual grid is the Heisenberg tradeoff — not a rendering bug. Smoothing the display would hide the real resolution.

**What we already do:**
- `shading="gouraud"` (interpolated, not flat blocks)
- 300 DPI rendering
- 50% overlap (COLA for Hann)

**What we could try:**
- [ ] Increase overlap to 75% or 87.5% — more STFT frames, smoother time axis, but does not improve actual resolution
- [ ] Increase `nfft` (zero-pad the FFT) — smoother frequency axis, interpolation only
- [ ] Side-by-side comparison: 50% overlap vs 87.5% overlap on the same signal, show that the information content is identical

**The risk:** over-smoothing makes the spectrogram look good but misrepresents the actual resolution. The reader thinks the tool is sharper than it is.

**To settle:** produce both versions on one figure (Lab 4 chirp), label one "actual resolution" and the other "interpolated display." If the instructor prefers the smooth one, use it but add a note that visual smoothness is interpolation, not resolution.

---

## 2. dB Labeling on Dimensionless Signals

**The argument:** "dB doesn't mean anything physically if you just scale it."

**Our position (settled):** Use "dB re 1 µV²/Hz" for real EEG PSD (physical reference). Use "10 log₁₀(·)" for dimensionless model signals and WVD values. Same plot, honest label.

**Status:** Resolved. Rule added to CLAUDE.md.

---

## 3. Dual-Stack: Always Show Both, or Only When Justified?

**The argument:** "Having both linear and dB for every figure is excessive."

**Our position (settled):** Produce both as separate figures. Include the log-scale plot in the report only when it reveals something the linear plot hides. Otherwise keep it as an unnumbered draft in `results/graphs/`.

**Status:** Resolved. Rule updated in CLAUDE.md.

---

## 4. Burst Threshold: Why 2x Median?

**The argument:** No formal justification for choosing 2x median as the burst detection threshold.

**Our position:** It's a pragmatic choice. The bursts are so energetic (17x median peak) that any threshold between 1.5x and 3x gives roughly the same result.

**What we could try:**
- [ ] Sweep the threshold from 1.5x to 4x, plot burst percentage vs threshold — if the curve is flat in the 1.5-3x range, the choice is robust
- [ ] Connect to Lab 2's detection framework: 2x median ≈ γ = 1.4, P_fa ≈ 25%. Not rigorous, but shows the connection.
- [ ] Cite neonatal EEG literature if a similar threshold is used

**To settle:** run the sweep experiment, add one sentence to C.3 showing robustness.

---

## 5. Code Shown in Report — How Much Is Too Much?

**The argument:** The instructor wants to see the code, but showing every function call makes the report very long.

**Our position:** Show the code that produces each figure (the critical chain: code → figure → interpretation). Infrastructure code (config.py, plotting.py) is shown once in Volume B's "Basis Functions" section, not repeated.

**What we could try:**
- [ ] Collapse long code blocks into shorter "key lines only" versions with comments referencing the full source
- [ ] Move some code to an appendix

**To settle:** instructor preference. Ask directly.

---

## 6. WVD dB Panel Artifacts (Circular Patterns)

**The argument:** The dB panel of the SPWVD on real EEG shows circular/elliptical patterns that could be mistaken for brain activity.

**Our position:** These are residual cross-term artifacts, not physiology. The linear panel is the reliable deliverable. The dB panel is included with an explicit caveat, or excluded entirely.

**What we could try:**
- [ ] Label the dB panel explicitly: "Visualization only — circular patterns are residual cross-term artifacts, not neural activity"
- [ ] Show both labeled and unlabeled to the instructor, ask which is clearer
- [ ] Drop the WVD/SPWVD dB panel entirely from the report, keep in drafts

**To settle:** instructor feedback on whether the caveat is sufficient or the panel should be removed.
