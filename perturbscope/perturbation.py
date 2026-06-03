# =============================================================================
# perturbscope/perturbation.py
#
# Effective perturbation scoring module.
#
# Many cells in a Perturb-seq experiment do not actually get perturbed despite
# carrying a guide RNA — they "escape".  This module detects and quantifies
# that distinction using a three-feature Gaussian Mixture Model (GMM) fit
# separately for each perturbation target.
#
# Features used per cell:
#   1. target_expr_reduction  -- how much the target gene is knocked down
#   2. perturbation_signature -- distance from control centroid in PCA space
#   3. 1 - neighbor_control_fraction -- how isolated from controls the cell is
#
# Outputs written to adata.obs:
#   neighbor_control_fraction -- fraction of k nearest neighbours that are controls
#   perturbation_signature    -- normalised distance from control centroid
#   perturbation_score        -- P(truly perturbed) from GMM posterior
#   escape_probability        -- 1 - perturbation_score
#   perturb_class             -- categorical: "perturbed" / "escaped" / "ambiguous"
#
# Outputs written to adata.uns:
#   perturbation_stats -- per-perturbation cell count summary
# =============================================================================

from __future__ import annotations

import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

from .utils import minmax_01


def score_effective_perturbation(adata, n_neighbors=20):
    # Score each cell on how effectively it was perturbed.
    #
    # Args:
    #   n_neighbors -- number of nearest neighbours to use when computing
    #                  neighbor_control_fraction.
    #
    # Returns the modified AnnData.
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca is required. Run preprocessing first.")

    x = adata.obsm["X_pca"]

    # -- Step 1: neighbour control fraction --
    # For each cell, find its k nearest neighbours in PCA space and record
    # what fraction are control cells.  A highly perturbed cell should sit far
    # from controls, so this fraction will be low.
    nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, adata.n_obs))
    nn.fit(x)
    neigh = nn.kneighbors(return_distance=False)

    is_ctrl = adata.obs["is_control"].values.astype(bool)
    frac = np.zeros(adata.n_obs, dtype=float)
    for i in range(adata.n_obs):
        idx = neigh[i][1:]  # exclude the cell itself (index 0)
        frac[i] = is_ctrl[idx].mean() if len(idx) else 0.0
    adata.obs["neighbor_control_fraction"] = frac

    # -- Step 2: perturbation signature score --
    # Use the Euclidean distance from the mean control centroid in PCA space
    # as a proxy for how far a cell has moved away from the unperturbed state.
    # Distances are min-max normalised to [0,1] across all cells.
    ctrl_mask = adata.obs["is_control"].values
    ctrl_centroid = x[ctrl_mask].mean(axis=0) if ctrl_mask.sum() else np.zeros(x.shape[1])
    dist = np.linalg.norm(x - ctrl_centroid.reshape(1, -1), axis=1)
    adata.obs["perturbation_signature"] = minmax_01(dist)

    # -- Step 3: retrieve earlier QC metrics for the feature matrix --
    ter = adata.obs.get("target_expr_reduction", 0.0)
    if not hasattr(ter, "values"):
        ter = np.zeros(adata.n_obs)
    ter = np.nan_to_num(np.asarray(ter, dtype=float), nan=0.0)
    sig = np.asarray(adata.obs["perturbation_signature"].values, dtype=float)
    neigh_eff = 1.0 - np.asarray(adata.obs["neighbor_control_fraction"].values, dtype=float)

    # Initialise output columns with conservative defaults.
    adata.obs["perturbation_score"] = 0.0
    adata.obs["escape_probability"] = 1.0
    adata.obs["perturb_class"] = "ambiguous"

    # -- Step 4: per-perturbation GMM classification --
    # Fit a 2-component GMM for each perturbation separately.
    # The component with the higher mean feature values is labelled "perturbed".
    stats = {}
    for perturb in sorted(adata.obs["perturbation"].astype(str).unique()):
        if perturb.lower() == "control":
            continue
        mask = adata.obs["perturbation"].values == perturb
        if mask.sum() < 5:  # too few cells to fit a GMM reliably
            continue
        feat = np.column_stack([ter[mask], sig[mask], neigh_eff[mask]])
        if np.allclose(feat.std(axis=0), 0):
            # All features are identical: cannot separate modes.
            score = np.zeros(mask.sum())
        else:
            gmm = GaussianMixture(n_components=2, random_state=0)
            gmm.fit(feat)
            probs = gmm.predict_proba(feat)
            means = gmm.means_.mean(axis=1)  # average feature value per component
            perturbed_comp = int(np.argmax(means))  # component with higher mean = perturbed
            score = probs[:, perturbed_comp]

        adata.obs.loc[mask, "perturbation_score"] = score
        adata.obs.loc[mask, "escape_probability"] = 1.0 - score
        # Threshold-based hard classification: >0.66 perturbed, <0.33 escaped.
        cls = np.where(score > 0.66, "perturbed", np.where(score < 0.33, "escaped", "ambiguous"))
        adata.obs.loc[mask, "perturb_class"] = cls
        stats[perturb] = {
            "n_cells": int(mask.sum()),
            "n_perturbed": int((cls == "perturbed").sum()),
            "n_escaped": int((cls == "escaped").sum()),
            "n_ambiguous": int((cls == "ambiguous").sum()),
            "mean_perturbation_score": float(np.mean(score)),
        }

    adata.uns["perturbation_stats"] = stats
    return adata
