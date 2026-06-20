# src/ — Code Standards for Contributors

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
  common/              shared utilities — the single source of truth
    config.py          constants: FS, DURATION, SEED, DPI, colormaps, paths
    signals.py         model signal generators (one per Appendix A archetype)
    windows.py         window functions (cosine-sum family + Gaussian)
    plotting.py        dual-stack plotting, spectrogram, time-domain plots
  lab1_dft/            one folder per lab
  lab2_statistics/
  lab3_windows/
  lab4_stft/
  lab5_autocorrelation/
  lab6_wvd/
  lab7_spwvd/

notebooks/             .ipynb mirrors of each lab (narrative + code + inline figures)

results/graphs/labN/   exported figures, named by lab and figure number
```

## Rules

1. **Import from `common/`, never redefine constants.** Use `from src.common import FS, DURATION, SEED` or `from src.common import make_tone, plot_dual_stack_spectrum`. Never hardcode `fs = 250` or `dpi = 300` in lab code.

2. **All signals go through `signals.py`.** The generators enforce the Volume B constraints (frequencies < 100 Hz, fs = 250 Hz, duration >= 1200 s). If you need a new signal type, add it to `signals.py`, not to your lab file.

3. **All plots go through `plotting.py`.** This enforces DPI, axis labels with units, perceptually uniform colormaps, and the dual-stack rule (linear first, dB second). Call `plot_dual_stack_spectrum()` or `plot_spectrogram()` — don't configure matplotlib yourself.

4. **Fixed random seed.** Any noise uses `SEED` from `config.py`. Pass it explicitly: `make_noise(seed=SEED)`. Results must be reproducible across machines.

5. **Notebooks mirror, don't own.** Lab logic lives in `src/labN_*/`. Notebooks import from `src/` and add markdown narrative. The notebook is the presentation; `src/` is the testable, diffable logic.

6. **Figure naming.** Export via `save_figure(fig, lab_number=1, fig_id="01")` → `results/graphs/lab1/figure_B_01.png`. The figure ID matches the Volume B numbering in the report.

7. **Comment every line.** Per CLAUDE.md: inline comments on what each line does. Keep functions compact and single-purpose.
