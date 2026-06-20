# src/ - Code Standards for Contributors

## Environment setup

This project uses **Conda** (miniforge3). The environment is defined in `environment.yml` at the project root.

```bash
# First time setup
conda env create -f environment.yml

# After adding a dependency to environment.yml
conda env update -f environment.yml --prune
```

Never `pip install` directly into the environment. All dependencies go in `environment.yml` so the environment is reproducible across machines.

## Directory layout

```
src/
  common/              shared utilities - the single source of truth
    config.py          constants: FS, DURATION, SEED, DPI, colormaps, EEG_BANDS, paths
    signals.py         model signal generators (one per Appendix A archetype)
    windows.py         window functions (cosine-sum family + Gaussian)
    plotting.py        dual-stack plotting, spectrogram, time-domain plots
    eeg.py             EEG loading (MNE), band-power (Welch), channel selection
  lab1_dft/            Lab 1: DFT basics
  lab2_statistics/     Lab 2: statistics on noisy signal
  lab3_windows/        Lab 3: windowing and the Dirichlet kernel
  lab4_stft/           Lab 4: STFT of a fluctuating signal
  lab5_resolution/     Lab 5: two-tone resolution on the spectrogram
  lab6_autocorrelation/ Lab 6: autocorrelation of a noisy signal
  lab7_wvd/            Lab 7: WVD and its tradeoffs
  lab8_spwvd/          Lab 8: SPWVD and its tradeoffs

build_docs/            document build scripts
  prettier.py          touch up markdown before conversion
  make_template.py     generate reference.docx template
  build_docx.sh        pandoc build script

docs/                  all markdown content
  volume_A.md          Volume A - theory
  volume_B.md          Volume B - labs + Appendix B
  slides.md            unified slide deck (Marp)
  Thesis_Table_of_Contents.md  master plan (internal)

template/              reference.docx for pandoc (Roboto, Consolas, margins)
data/                  raw EEG data (Volume C)
notebooks/             .ipynb mirrors of each lab
results/graphs/labN/   exported figures (PNG + SVG)
output/                generated documents (docx, pdf, html) - gitignored
```

## Rules

1. **Import from `common/`, never redefine constants.** Use `from src.common import FS, DURATION, SEED` or `from src.common import make_tone, plot_dual_stack_spectrum`. Never hardcode `fs = 250` or `dpi = 300` in lab code.

2. **All signals go through `signals.py`.** The generators enforce the Volume B constraints (frequencies < 100 Hz, fs = 250 Hz, duration >= 1200 s). If you need a new signal type, add it to `signals.py`, not to your lab file.

3. **All plots go through `plotting.py`.** This enforces DPI, axis labels with units, perceptually uniform colormaps, and the dual-stack rule (linear first, dB second). Call `plot_dual_stack_spectrum()` or `plot_spectrogram()` - don't configure matplotlib yourself.

4. **All EEG I/O goes through `eeg.py`.** Use `load_eeg()` for loading (returns data in uV), `compute_band_power()` for Welch PSD, `get_channel_data()` for channel selection. Never call MNE directly in lab code.

5. **Fixed random seed.** Any noise uses `SEED` from `config.py`. Pass it explicitly: `make_noise(seed=SEED)`. Results must be reproducible across machines.

6. **Notebooks mirror, don't own.** Lab logic lives in `src/labN_*/`. Notebooks import from `src/` and add markdown narrative. The notebook is the presentation; `src/` is the testable, diffable logic.

7. **Figure naming.** Export via `save_figure(fig, lab_number=1, fig_id="01")`. Figures go to `results/graphs/labN/figure_B_XX.png` (+ SVG for line plots). Use `raster_only=True` for spectrograms (SVG of dense heatmaps crashes renderers).

8. **Comment every line.** Per CLAUDE.md: inline comments on what each line does. Keep functions compact and single-purpose.

9. **Three independent counters.** Equations (B.1), Figures B.1, Tables B.1 are different objects. Always prefix with "Equation", "Figure", or "Table" - never a bare number.

10. **Building documents.** Run `python build_docs/prettier.py docs/volume_B.md` before converting. Then `bash build_docs/build_docx.sh` for Word, or `marp docs/slides.md -o output/slides.html` for slides.
