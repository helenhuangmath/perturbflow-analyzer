# =============================================================================
# perturbscope/preprocessing.py
#
# Normalisation, dimensionality reduction, and cell-state clustering.
#
# This module transforms raw counts into two complementary representations:
#
#   X_bio           -- standard PCA embedding capturing the biological manifold,
#                      used for clustering and trajectory analysis.
#   X_perturb_resid -- perturbation residual space: PCA coords minus the mean
#                      control centroid, highlighting deviation from baseline.
#
# Additional outputs:
#   X_umap   (in adata.obsm) -- 2-D UMAP for visualisation (scanpy path only)
#   cell_state (in adata.obs) -- cluster / leiden community label
#
# Two execution paths are supported:
#   Primary   -- full scanpy pipeline (normalise → HVG → scale → PCA →
#                neighbours → UMAP → Leiden clustering)
#   Fallback  -- sklearn-only pipeline (log-normalise → PCA → KMeans)
#                used when scanpy is not available or fails
# =============================================================================

from __future__ import annotations

import numpy as np


def normalize_and_embed(adata, random_state=0, leiden_resolution=0.5):
    # Normalise expression, compute PCA, UMAP, and cell-state clusters.
    #
    # Args:
    #   random_state       -- seed for PCA, UMAP, KMeans; ensures reproducibility.
    #   leiden_resolution  -- Leiden community-detection resolution (default 0.5).
    #
    # Returns the modified AnnData with obsm keys X_pca, X_bio, X_perturb_resid
    # (and X_umap when scanpy is used) plus obs column cell_state.
    try:
        import scanpy as sc

        # Normalise total counts to 10 000 per cell, then log-transform.
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # Select the top highly variable genes to reduce dimensionality and noise.
        sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars), inplace=True)
        ad = adata[:, adata.var["highly_variable"]].copy() if "highly_variable" in adata.var else adata
        sc.pp.scale(ad, max_value=10)
        sc.tl.pca(ad, n_comps=min(50, ad.n_vars - 1, ad.n_obs - 1), random_state=random_state)
        adata.obsm["X_pca"] = ad.obsm["X_pca"]
        # Build kNN graph and embed in 2-D UMAP for visualisation.
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
        sc.tl.umap(adata, random_state=random_state)
        # Leiden community detection defines cell states used in effect decomposition.
        sc.tl.leiden(adata, key_added="cell_state", resolution=leiden_resolution)
    except Exception:
        # Fallback: sklearn-only path (no scanpy dependency).
        from sklearn.decomposition import PCA
        from sklearn.cluster import KMeans

        x = adata.X.A if hasattr(adata.X, "A") else np.asarray(adata.X)
        x = np.log1p(x)
        x = x - x.mean(axis=0, keepdims=True)  # centre genes
        pca = PCA(n_components=min(30, x.shape[1] - 1, x.shape[0] - 1), random_state=random_state)
        adata.obsm["X_pca"] = pca.fit_transform(x)
        # Use KMeans as a substitute for Leiden clustering.
        km = KMeans(n_clusters=min(10, max(2, adata.n_obs // 50)), random_state=random_state, n_init=10)
        adata.obs["cell_state"] = km.fit_predict(adata.obsm["X_pca"]).astype(str)

    # X_bio is a copy of the biological PCA embedding kept separate from X_perturb_resid.
    adata.obsm["X_bio"] = adata.obsm["X_pca"].copy()

    # X_perturb_resid centres the embedding on the mean control cell, so that
    # distances in this space reflect deviation from the unperturbed baseline.
    ctrl_mask = adata.obs["is_control"].values
    if ctrl_mask.sum() > 0:
        centroid = adata.obsm["X_pca"][ctrl_mask].mean(axis=0, keepdims=True)
        adata.obsm["X_perturb_resid"] = adata.obsm["X_pca"] - centroid
    else:
        adata.obsm["X_perturb_resid"] = adata.obsm["X_pca"].copy()
    return adata
