"""
Project-wide constants. Import from here — never hardcode these values in lab code.
"""

import os

# --- Sampling and signal parameters (Volume B lab constraints) ---
FS = 250  # sampling frequency (Hz) — default for all labs
DURATION = 1200  # minimum signal duration (s) — 20 minutes
N_SAMPLES = FS * DURATION  # total samples at default fs and duration
F_MAX = 100  # maximum allowed signal frequency (Hz)

# --- Reproducibility ---
SEED = 42  # fixed random seed for all noise generation

# --- Figure export ---
DPI = 300  # minimum render DPI
FIGURE_FORMATS = ["png", "svg"]  # export in both raster and vector
COLORMAP = "viridis"  # default perceptually uniform colormap
COLORMAP_DIVERGING = "inferno"  # alternative for diverging data
FORBIDDEN_COLORMAPS = {"jet", "rainbow", "hsv"}  # never use these

# --- EEG frequency bands (Volume C) ---
EEG_BANDS = {
    "delta": (0.5, 4),     # δ: 0.5–4 Hz
    "theta": (4, 8),       # θ: 4–8 Hz
    "alpha": (8, 13),      # α: 8–13 Hz
    "beta": (13, 30),      # β: 13–30 Hz
    "gamma": (30, 100),    # γ: 30–100 Hz (capped at F_MAX)
}

# --- Paths (relative to project root) ---
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GRAPH_DIR = os.path.join(PROJECT_ROOT, "results", "graphs")
DATA_DIR = os.path.join(PROJECT_ROOT, "data")             # raw EEG data (Volume C)
