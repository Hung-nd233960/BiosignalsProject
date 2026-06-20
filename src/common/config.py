"""
Project-wide constants. Import from here — never hardcode these values in lab code.
"""
import os

# --- Sampling and signal parameters (Volume B lab constraints) ---
FS = 250              # sampling frequency (Hz) — default for all labs
DURATION = 1200       # minimum signal duration (s) — 20 minutes
N_SAMPLES = FS * DURATION  # total samples at default fs and duration
F_MAX = 100           # maximum allowed signal frequency (Hz)

# --- Reproducibility ---
SEED = 42             # fixed random seed for all noise generation

# --- Figure export ---
DPI = 300                         # minimum render DPI
FIGURE_FORMATS = ["png", "svg"]   # export in both raster and vector
COLORMAP = "viridis"              # default perceptually uniform colormap
COLORMAP_DIVERGING = "inferno"    # alternative for diverging data
FORBIDDEN_COLORMAPS = {"jet", "rainbow", "hsv"}  # never use these

# --- Paths (relative to project root) ---
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GRAPH_DIR = os.path.join(PROJECT_ROOT, "results", "graphs")
