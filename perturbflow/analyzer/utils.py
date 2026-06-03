# =============================================================================
# perturbflow/analyzer/utils.py
#
# Shared utility helpers used across all PerturbScope modules.
# Contains directory management, numeric normalization, and safe DataFrame
# column access.  No biological logic lives here.
# =============================================================================

from __future__ import annotations

from pathlib import Path
import numpy as np


def ensure_dir(path: str | Path) -> Path:
    # Create the directory (and all parents) if it does not already exist.
    # Returns the resolved Path object so callers can chain it.
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def minmax_01(values: np.ndarray) -> np.ndarray:
    # Rescale a 1-D array to the [0, 1] range using min-max normalisation.
    # NaN values are ignored when computing the range.
    # If all values are equal the function returns an all-zero array.
    arr = np.asarray(values, dtype=float)
    lo = np.nanmin(arr)
    hi = np.nanmax(arr)
    if np.isclose(hi, lo):
        return np.zeros_like(arr, dtype=float)
    return (arr - lo) / (hi - lo)


def safe_col(obs, key: str, default):
    # Return obs[key] if it exists, otherwise initialise it to `default` first.
    # Useful for lazily creating optional metadata columns without overwriting
    # values that may already be present in the AnnData .obs table.
    if key not in obs:
        obs[key] = default
    return obs[key]
