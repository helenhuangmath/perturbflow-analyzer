# =============================================================================
# perturbflow/analyzer/state.py
#
# Dual-level effect decomposition — core novelty of PerturbFlow Analyzer.
#
# Standard pipelines report a single "perturbation effect".  This module
# separates that effect into two orthogonal components:
#
#   transcriptional_score  -- how much gene expression changes *within* each
#                             cell state (computed as mean |log2FC| across states)
#   state_shift_score      -- how much the *proportion* of cells in each state
#                             changes relative to controls (total variation
#                             distance between composition vectors)
#
# A perturbation that forces cells to differentiate might have a large
# state_shift_score but a small transcriptional_score; a gene that regulates
# a pathway without changing identity shows the opposite pattern.
#
# Outputs written to adata.uns["effect_scores"][perturbation]:
#   transcriptional  -- float score
#   state_shift      -- float score
#
# Returns (adata, effect_df) where effect_df is a tidy DataFrame with one
# row per perturbation and a dominant_effect_type label.
# =============================================================================

from __future__ import annotations

from collections import defaultdict
import numpy as np
import pandas as pd


def _to_dense(x):
    # Convert a sparse matrix slice to a dense numpy array.
    return x.A if hasattr(x, "A") else np.asarray(x)


def compute_effect_decomposition(adata, min_cells=5):
    # Decompose the effect of each perturbation into transcriptional and
    # state-shift components.
    #
    # Args:
    #   min_cells -- minimum cells required in a (perturbation, state) group
    #               to include that state in the transcriptional calculation.
    #
    # Returns (adata, effect_df).

    # Ensure cell_state exists; fall back to a single dummy state.
    if "cell_state" not in adata.obs.columns:
        adata.obs["cell_state"] = "state0"

    controls = adata.obs["is_control"].values
    if controls.sum() == 0:
        raise ValueError("At least one control cell is required for effect decomposition")

    scores = defaultdict(dict)
    rows = []
    states = adata.obs["cell_state"].astype(str).unique().tolist()

    for perturb in sorted(adata.obs["perturbation"].astype(str).unique()):
        if perturb.lower() == "control":
            continue
        p_mask = adata.obs["perturbation"].values == perturb
        if p_mask.sum() < min_cells:
            continue

        # -- Within-state transcriptional effect --
        # For each cell state, compute the mean per-gene |log2FC| between
        # perturbed and control cells.  Average across states.
        state_effects = []
        for st in states:
            s_mask = adata.obs["cell_state"].astype(str).values == st
            p_idx = p_mask & s_mask
            c_idx = controls & s_mask
            if p_idx.sum() < min_cells or c_idx.sum() < min_cells:
                continue
            xp = _to_dense(adata.X[p_idx, :]).mean(axis=0)
            xc = _to_dense(adata.X[c_idx, :]).mean(axis=0)
            lfc = np.log2((xp + 1e-3) / (xc + 1e-3))
            state_effects.append(np.mean(np.abs(lfc)))
        transcriptional = float(np.mean(state_effects)) if state_effects else 0.0

        # -- Between-state composition shift --
        # Compare the cell-type composition vectors of perturbed vs control
        # using total variation distance (= half the L1 distance).
        comp_p = (
            adata.obs.loc[p_mask, "cell_state"].astype(str).value_counts(normalize=True)
        )
        comp_c = (
            adata.obs.loc[controls, "cell_state"].astype(str).value_counts(normalize=True)
        )
        all_states = sorted(set(comp_p.index).union(comp_c.index))
        vp = np.array([float(comp_p.get(s, 0.0)) for s in all_states])
        vc = np.array([float(comp_c.get(s, 0.0)) for s in all_states])
        state_shift = float(np.sum(np.abs(vp - vc)) / 2.0)

        scores[perturb]["transcriptional"] = transcriptional
        scores[perturb]["state_shift"] = state_shift

        # Classify the dominant effect type using a 20% margin.
        if transcriptional > state_shift * 1.2:
            dominant = "transcriptional"
        elif state_shift > transcriptional * 1.2:
            dominant = "state_shift"
        elif transcriptional == 0.0 and state_shift == 0.0:
            dominant = "neither"
        else:
            dominant = "both"

        rows.append(
            {
                "perturbation": perturb,
                "transcriptional_score": transcriptional,
                "state_shift_score": state_shift,
                "dominant_effect_type": dominant,
            }
        )

    adata.uns.setdefault("effect_scores", {})
    adata.uns["effect_scores"].update({k: dict(v) for k, v in scores.items()})
    return adata, pd.DataFrame(rows)
