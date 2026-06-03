# =============================================================================
# perturbscope/qc.py
#
# Quality-control module for Perturb-seq data.
#
# Beyond standard scRNA-seq QC (gene counts, mitochondrial fraction), this
# module computes three Perturb-seq-specific metrics that are novel to
# PerturbScope:
#
#   perturbation_burden       -- proxy for how heavily a cell is perturbed,
#                                derived from guide count × total UMIs,
#                                min-max normalised per perturbation group.
#   target_expr_reduction     -- fraction by which the target gene's expression
#                                is reduced relative to control cells; closer to
#                                1 means strong knock-down (TER formula).
#   guide_confidence_score    -- composite score [0,1] combining perturbation
#                                burden and target expression reduction.
#
# Outputs written to adata.obs:
#   n_genes_by_counts, total_counts, pct_counts_mt  (standard)
#   perturbation_burden, target_expr_reduction, guide_confidence_score  (novel)
#   qc_pass  -- boolean flag: True if the cell passes all QC thresholds
# =============================================================================

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from .utils import minmax_01

# ---------------------------------------------------------------------------
# Global font constants (mirror eda.py — increase here to resize QC plots)
# ---------------------------------------------------------------------------
_TITLE_FS  = 15
_LABEL_FS  = 13
_TICK_FS   = 11
_LEGEND_FS = 11


def _to_dense(x):
    # Convert a sparse matrix slice to a dense numpy array.
    # Works transparently on both scipy sparse (.A) and dense arrays.
    return x.A if hasattr(x, "A") else np.asarray(x)


def _target_expr_reduction(adata) -> np.ndarray:
    # Compute Target Expression Reduction (TER) for every perturbed cell.
    #
    # Formula (per target gene g, per perturbed cell i):
    #   TER_i = clip(1 - expr_i(g) / median_ctrl(g), 0, 1)
    #
    # For combinatorial perturbations (multi-gene targets) the TER values
    # across all target genes are averaged.
    # Control cells receive NaN (they have no target gene to reduce).
    out = np.full(adata.n_obs, np.nan, dtype=float)
    if "perturbation" not in adata.obs.columns:
        return out

    controls = adata.obs["is_control"].values
    if controls.sum() == 0:
        return out

    for perturb in sorted(adata.obs["perturbation"].astype(str).unique()):
        if perturb.lower() == "control":
            continue
        # Split combinatorial perturbations (e.g. "GENE1+GENE2") into parts.
        genes = [g.strip() for g in str(perturb).replace(",", "+").split("+") if g.strip()]
        pert_mask = adata.obs["perturbation"].values == perturb
        vals = []
        for g in genes:
            if g not in adata.var_names:
                continue
            idx = int(np.where(adata.var_names == g)[0][0])
            ctrl_expr = _to_dense(adata.X[controls, idx]).reshape(-1)
            cell_expr = _to_dense(adata.X[pert_mask, idx]).reshape(-1)
            baseline = float(np.median(ctrl_expr)) + 1e-6  # pseudocount avoids div/0
            ter = np.clip(1.0 - (cell_expr / baseline), 0.0, 1.0)
            vals.append(ter)
        if vals:
            # Average TER across all target genes for this perturbation.
            out[pert_mask] = np.nanmean(np.vstack(vals), axis=0)
    return out


def run_qc(adata, min_genes=200, max_pct_mt=20.0, min_cells_per_perturbation=20):
    # Run the full QC pipeline on an AnnData object and annotate each cell.
    #
    # Args:
    #   min_genes                   -- minimum detected genes to keep a cell
    #   max_pct_mt                  -- maximum mitochondrial % to keep a cell
    #   min_cells_per_perturbation  -- perturbations with fewer cells are excluded
    #
    # Returns the same AnnData with new columns added to .obs.

    # -- Standard QC metrics via scanpy (falls back to manual if unavailable) --
    try:
        import scanpy as sc

        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    except Exception:
        # Fallback: compute minimal metrics without scanpy.
        x = _to_dense(adata.X)
        adata.obs["n_genes_by_counts"] = (x > 0).sum(axis=1)
        adata.obs["total_counts"] = x.sum(axis=1)
        adata.obs["pct_counts_mt"] = 0.0

    # -- Novel metric 1: perturbation_burden --
    # Raw burden = number of guides in cell × total UMI count.
    # This captures how "loaded" a cell is with perturbation reagent.
    adata.obs["perturbation_burden_raw"] = adata.obs["n_guides"].values * adata.obs[
        "total_counts"
    ].astype(float).values

    # Min-max normalise burden within each perturbation group so values are
    # comparable across groups with different guide efficiencies.
    adata.obs["perturbation_burden"] = 0.0
    for perturb, idx in adata.obs.groupby("perturbation").groups.items():
        arr = adata.obs.loc[idx, "perturbation_burden_raw"].values
        adata.obs.loc[idx, "perturbation_burden"] = minmax_01(arr)

    # -- Novel metric 2: target_expr_reduction (TER) --
    adata.obs["target_expr_reduction"] = _target_expr_reduction(adata)

    # -- Novel metric 3: guide_confidence_score --
    # Equal-weight average of burden and TER, clamped to [0,1].
    ter = adata.obs["target_expr_reduction"].fillna(0.0).values
    burden = adata.obs["perturbation_burden"].fillna(0.0).values
    adata.obs["guide_confidence_score"] = np.clip(0.5 * burden + 0.5 * ter, 0.0, 1.0)

    # -- QC pass/fail flag --
    # A cell passes if it meets gene count, MT%, and perturbation group size thresholds.
    counts = adata.obs["perturbation"].value_counts()
    keep_perturb = counts[counts >= min_cells_per_perturbation].index
    qc_pass = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mt)
        & (adata.obs["perturbation"].isin(keep_perturb))
    )
    adata.obs["qc_pass"] = qc_pass
    return adata


def plot_qc_summary(adata, out_png: str | None = None):
    # Generate a 2×2 violin-plot panel showing the four key QC metrics grouped
    # by perturbation label (top-10 by cell count, always including controls).
    # Saves to `out_png` if a path is provided; always returns the Figure.
    metrics = [
        "n_genes_by_counts",
        "total_counts",
        "target_expr_reduction",
        "guide_confidence_score",
    ]
    # Select the top-10 most common perturbations plus the control group.
    top = adata.obs["perturbation"].value_counts().head(10).index.tolist()
    if "control" in adata.obs["perturbation"].values and "control" not in top:
        top.append("control")
    sub = adata.obs[adata.obs["perturbation"].isin(top)].copy()

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.ravel()
    for i, m in enumerate(metrics):
        ax = axes[i]
        sns.violinplot(
            data=sub, x="perturbation", y=m, ax=ax,
            inner=None, cut=0,
        )
        # Overlay mean (black diamond) and median (red hlines) per group.
        for j, p in enumerate(top):
            vals = sub.loc[sub["perturbation"] == p, m].dropna().values
            if len(vals) == 0:
                continue
            ax.scatter(j, float(vals.mean()), marker="D", color="black",
                       s=40, zorder=5, label="mean" if j == 0 else "_")
            half = 0.18
            ax.hlines(float(np.median(vals)), j - half, j + half,
                      colors="#e63946", linewidths=2, zorder=5,
                      label="median" if j == 0 else "_")
        ax.tick_params(axis="x", rotation=45, labelsize=_TICK_FS)
        ax.tick_params(axis="y", labelsize=_TICK_FS)
        ax.set_title(m, fontsize=_TITLE_FS)
        ax.set_xlabel("Perturbation", fontsize=_LABEL_FS)
        ax.set_ylabel(m, fontsize=_LABEL_FS)
        # Show mean/median legend only on the first panel to avoid clutter.
        if i == 0:
            ax.legend(loc="upper right", fontsize=_LEGEND_FS)

    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=180, bbox_inches="tight")
    return fig
