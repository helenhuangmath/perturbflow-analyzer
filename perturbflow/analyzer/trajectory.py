# =============================================================================
# perturbflow/analyzer/trajectory.py
#
# Trajectory and cell-fate analysis module.
#
# Maps how each perturbation shifts cells along a developmental or activation
# trajectory, capturing both the direction (earlier vs later pseudotime) and
# the degree of fate redistribution.
#
# Pseudotime inference:
#   If adata.obs["pseudotime"] is already present (e.g. from diffusion
#   pseudotime or Monocle), it is used directly.  Otherwise the first PCA
#   component is rescaled to [0,1] as a surrogate.
#
# Outputs per perturbation (returned as DataFrame and stored in adata.uns):
#   pseudotime_shift         -- mean pseudotime of perturbed cells minus mean
#                               pseudotime of control cells (signed)
#   branch_probability_change -- maximum absolute change in cell-state fraction
#   fate_bias_score          -- Jensen-Shannon divergence between perturbed and
#                               control cell-state composition vectors
#   commitment_shift_index   -- absolute pseudotime shift (unsigned magnitude)
# =============================================================================

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


def compute_trajectory_effects(adata):
    # Compute per-perturbation trajectory metrics and store results.
    #
    # Returns (adata, trajectory_df) where trajectory_df has one row per
    # perturbation and the four metrics described above.
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca is required for trajectory effects")

    # Use PC1 as a pseudotime surrogate if no real pseudotime is available.
    if "pseudotime" not in adata.obs.columns:
        pc1 = adata.obsm["X_pca"][:, 0]
        mn, mx = float(pc1.min()), float(pc1.max())
        adata.obs["pseudotime"] = (pc1 - mn) / (mx - mn + 1e-8)

    if "cell_state" not in adata.obs.columns:
        adata.obs["cell_state"] = "state0"

    ctrl = adata.obs["is_control"].values
    # Baseline pseudotime: mean pseudotime across all control cells.
    ctrl_pt = adata.obs.loc[ctrl, "pseudotime"].mean()

    rows = []
    for perturb in sorted(adata.obs["perturbation"].astype(str).unique()):
        if perturb.lower() == "control":
            continue
        pm = adata.obs["perturbation"].values == perturb
        if pm.sum() < 5:
            continue

        # Signed pseudotime shift: positive = shifted toward later / more
        # differentiated states; negative = shifted toward earlier states.
        pt_shift = float(adata.obs.loc[pm, "pseudotime"].mean() - ctrl_pt)

        # Cell-state composition vectors for perturbed vs control cells.
        comp_p = adata.obs.loc[pm, "cell_state"].astype(str).value_counts(normalize=True)
        comp_c = adata.obs.loc[ctrl, "cell_state"].astype(str).value_counts(normalize=True)
        all_states = sorted(set(comp_p.index).union(comp_c.index))
        vp = np.array([float(comp_p.get(s, 0.0)) for s in all_states])
        vc = np.array([float(comp_c.get(s, 0.0)) for s in all_states])
        # Jensen-Shannon divergence as a symmetric, bounded fate-bias measure.
        js = float(jensenshannon(vp + 1e-8, vc + 1e-8))

        rows.append(
            {
                "perturbation": perturb,
                "pseudotime_shift": pt_shift,
                "branch_probability_change": float(np.max(np.abs(vp - vc))) if len(vp) else 0.0,
                "fate_bias_score": js,
                "commitment_shift_index": abs(pt_shift),
            }
        )

    df = pd.DataFrame(rows)
    adata.uns["trajectory_effects"] = df.to_dict(orient="records")
    return adata, df
