# =============================================================================
# perturbflow/analyzer/io.py
#
# Data ingestion and standardisation module.
#
# Responsibility: load a Perturb-seq dataset from an .h5ad file and produce a
# uniformly annotated AnnData object that all downstream modules can rely on.
#
# Guaranteed .obs columns after standardise_adata():
#   perturbation        -- gene name being targeted, or "control"
#   guide_ids           -- list of sgRNA IDs detected in the cell
#   is_control          -- True for non-targeting / control cells
#   n_guides            -- number of guides detected per cell
#   guide_confidence    -- initial confidence score (0 for controls, 1 otherwise;
#                          refined later by the scoring module)
#
# Guaranteed .uns keys:
#   perturbflow_version        -- version string
#   perturbation_targets       -- sorted list of all perturbation labels
#   n_cells_per_perturbation   -- dict mapping label -> cell count
# =============================================================================

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from . import __version__


def _control_like(value: str) -> bool:
    # Return True if the value looks like a non-targeting / control label.
    # Checks for common naming conventions: control, ctrl, NTC, scramble, etc.
    v = str(value).strip().lower()
    tags = ["control", "ctrl", "ntc", "non", "safe", "scramble", "wt"]
    return any(tag in v for tag in tags)


def _split_guides(value) -> list[str]:
    # Parse a guide-ID field into a list of individual guide strings.
    # Handles None, pre-split lists, or delimited strings (;  ,  |).
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    raw = str(value).strip()
    if not raw or raw.lower() == "nan":
        return []
    for sep in [";", ",", "|"]:
        if sep in raw:
            return [p.strip() for p in raw.split(sep) if p.strip()]
    return [raw]


def read_data(path: str, perturbation_col: Optional[str] = None):
    # Load an .h5ad file from disk and immediately standardise its metadata.
    # `perturbation_col` can name a specific .obs column holding perturbation
    # labels; if None, the function auto-detects from a list of common names.
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("scanpy is required to read AnnData files") from exc

    p = Path(path)
    if p.suffix.lower() != ".h5ad":
        raise ValueError("Currently only .h5ad inputs are supported")
    adata = sc.read_h5ad(path)
    return standardize_adata(adata, perturbation_col=perturbation_col)


def standardize_adata(adata, perturbation_col: Optional[str] = None):
    # Normalise an already-loaded AnnData so that all downstream modules see
    # the same column names regardless of the original dataset conventions.
    #
    # Step 1: resolve which .obs column holds perturbation labels, trying a
    #         priority list of common names; fall back to "control" for all cells.
    candidates = [
        perturbation_col,
        "perturbation",
        "target_gene",
        "gene_target",
        "perturb",
        "sgRNA_target",
        "guide_target",
        "guide",
    ]
    chosen = None
    for key in candidates:
        if key and key in adata.obs.columns:
            chosen = key
            break

    if chosen is None:
        adata.obs["perturbation"] = "control"
    else:
        adata.obs["perturbation"] = adata.obs[chosen].astype(str).fillna("control")

    # Step 2: locate guide-ID column and parse it into a list-per-cell.
    guide_col = None
    for key in ["guide_ids", "guide", "guides", "sgRNA", "sgRNA_id"]:
        if key in adata.obs.columns:
            guide_col = key
            break

    if guide_col is None:
        # Infer guide IDs from the perturbation label when no dedicated column exists.
        adata.obs["guide_ids"] = adata.obs["perturbation"].apply(
            lambda x: [] if _control_like(x) else [str(x)]
        )
    else:
        adata.obs["guide_ids"] = adata.obs[guide_col].apply(_split_guides)

    # Step 3: add derived boolean / integer columns used throughout the pipeline.
    adata.obs["is_control"] = adata.obs["perturbation"].map(_control_like).astype(bool)
    adata.obs["n_guides"] = adata.obs["guide_ids"].apply(len).astype(int)

    # Step 4: add a placeholder guide_confidence if not already present.
    #         The scoring module will overwrite this with a calibrated score.
    if "guide_confidence" not in adata.obs.columns:
        adata.obs["guide_confidence"] = np.where(adata.obs["is_control"], 0.0, 1.0)

    # Step 5: store dataset-level metadata in .uns for traceability.
    targets = sorted(pd.Series(adata.obs["perturbation"]).astype(str).unique().tolist())
    adata.uns["perturbflow_version"] = __version__
    adata.uns["perturbation_targets"] = targets
    adata.uns["n_cells_per_perturbation"] = (
        adata.obs["perturbation"].value_counts().to_dict()
    )
    return adata
