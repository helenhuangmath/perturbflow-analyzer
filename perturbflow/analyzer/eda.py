# =============================================================================
# perturbflow/analyzer/eda.py
#
# Exploratory data analysis (EDA) module for Perturb-seq data.
#
# Produces a suite of diagnostic plots and statistics AFTER QC and
# normalisation but BEFORE the more involved scoring/effects steps, so that
# analysts can inspect the data before committing to the full pipeline.
#
# Outputs (all written to <output_dir>/plots/ and <output_dir>/csv/):
#
#   plots/
#     eda_cells_per_perturbation.png             -- bar chart: n cells per group
#     eda_cluster_proportions.png                -- stacked bar: cluster composition
#                                                   per perturbation (top 30 groups)
#     eda_gene_by_cell_heatmap.png               -- top HVGs x sampled cells
#     eda_gene_by_perturbation_heatmap.png       -- pseudobulk mean per perturbation
#     eda_gene_correlation_heatmap.png           -- top HVG x HVG Pearson correlation
#     umap_cell_state.png                        -- UMAP coloured by cell_state cluster
#     eda_umap_perturbation.png                  -- UMAP coloured by perturbation label
#     eda_clustered_gene_by_pert_heatmap.png     -- top HVGs x perturbations with
#                                                   gene/pert cluster colour sidebars
#     eda_gene_pert_cluster_summary_heatmap.png  -- gene cluster x pert cluster
#                                                   mean z-score summary heatmap
#   csv/
#     eda_cells_per_perturbation.csv        -- raw cell-count table
#     eda_cluster_proportions.csv           -- per-perturbation cluster fractions
#     eda_gene_clusters.csv                 -- gene → cluster_id assignments
#     eda_pert_clusters.csv                 -- perturbation → cluster_id assignments
#
# All plots use seaborn / matplotlib only; no scanpy plotting is called so
# the module works on the sklearn-only fallback preprocessing path too.
# =============================================================================

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from .utils import ensure_dir

# ---------------------------------------------------------------------------
# Global font / style constants — increase here to resize all EDA plots
# ---------------------------------------------------------------------------
_TITLE_FS  = 18   # plot title
_LABEL_FS  = 15   # axis labels (xlabel / ylabel)
_TICK_FS   = 13   # axis tick labels
_LEGEND_FS = 13   # legend text
_CBAR_FS   = 13   # colorbar label
_ANNOT_FS  = 12   # in-plot annotations

_FONT_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "axes.titlesize": _TITLE_FS,
    "axes.labelsize": _LABEL_FS,
    "xtick.labelsize": _TICK_FS,
    "ytick.labelsize": _TICK_FS,
    "legend.fontsize": _LEGEND_FS,
    "legend.title_fontsize": _LEGEND_FS,
}
sns.set_theme(style="white", font="Arial", rc=_FONT_RC)
plt.rcParams.update(_FONT_RC)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _to_dense(x):
    return x.A if hasattr(x, "A") else np.asarray(x)


def _apply_plot_style() -> None:
    sns.set_theme(style="white", font="Arial", rc=_FONT_RC)
    plt.rcParams.update(_FONT_RC)


def _top_hvg_indices(adata, n: int) -> np.ndarray:
    # Return the indices of the top-n most variable genes.
    # Prefers scanpy's highly_variable flag when available; falls back to
    # computing per-gene variance on a cell subsample.
    if "highly_variable" in adata.var.columns:
        idx = np.where(adata.var["highly_variable"].values)[0]
        if len(idx) >= n:
            return idx[:n]
        if len(idx) > 0:
            return idx  # fewer HVGs than requested — use them all
    n_sample = min(5000, adata.n_obs)
    np.random.seed(0)
    sample_idx = np.random.choice(adata.n_obs, size=n_sample, replace=False)
    x = _to_dense(adata.X[np.sort(sample_idx), :])
    gene_var = x.var(axis=0)
    return np.argsort(gene_var)[::-1][:n]


def _eda_gene_indices(adata, n_top_genes: int) -> np.ndarray:
    return _top_hvg_indices(adata, min(n_top_genes, adata.n_vars))


_CONTROL_COLOR = "#9aa0a6"  # neutral grey for control/non-targeting groups


_CONTROL_KEYS = {"control", "ctrl", "nontargeting", "non-targeting", "nt",
                 "scramble", "safe-targeting", "safe_targeting"}


def _spectral_palette(labels) -> dict:
    """Spectral colormap palette for UMAP-style categorical plots.

    Control labels map to ``_CONTROL_COLOR``; non-control labels are spread
    across the Spectral colormap so neighbouring categories get visually
    distinct hues across the red→blue spectrum.
    """
    non_ctrl = sorted({str(l) for l in labels if str(l).lower() not in _CONTROL_KEYS})
    cols = sns.color_palette("Spectral", max(2, len(non_ctrl)))
    palette: dict = {}
    for i, p in enumerate(non_ctrl):
        palette[p] = cols[i % len(cols)]
    for l in labels:
        if str(l).lower() in _CONTROL_KEYS:
            palette[l] = _CONTROL_COLOR
    return palette


def _perturbation_palette(labels) -> dict:
    uniq = sorted({str(x) for x in labels})
    base = (
        sns.color_palette("tab20", 20)
        + sns.color_palette("tab20b", 20)
        + sns.color_palette("tab20c", 20)
        + sns.color_palette("Set3", 12)
    )
    palette = {p: base[i % len(base)] for i, p in enumerate(uniq)}
    for p in uniq:
        if p.lower() in {
            "control", "ctrl", "nontargeting", "non-targeting", "nt",
            "scramble", "safe-targeting", "safe_targeting",
        }:
            palette[p] = _CONTROL_COLOR
    return palette


def _pseudobulk_mean(adata, gene_idx: np.ndarray) -> pd.DataFrame:
    # Compute per-perturbation mean expression for the given gene indices.
    # Returns a DataFrame shaped (perturbations, genes).
    gene_names = adata.var_names[gene_idx].tolist()
    groups = sorted(adata.obs["perturbation"].astype(str).unique())
    rows = {}
    for g in groups:
        mask = adata.obs["perturbation"].values == g
        if mask.sum() == 0:
            continue
        x = _to_dense(adata.X[np.where(mask)[0], :][:, gene_idx])
        rows[g] = x.mean(axis=0)
    return pd.DataFrame(rows, index=gene_names).T  # (perturbation, gene)


def _zscore_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Column-wise z-score with clipping for heatmaps/similarity plots."""
    with np.errstate(invalid="ignore"):
        z = (df - df.mean(axis=0)) / (df.std(axis=0) + 1e-6)
    return z.replace([np.inf, -np.inf], 0).fillna(0).clip(-3, 3)


def _corr_for_display(corr: np.ndarray) -> np.ndarray:
    """Keep correlation sign while expanding weak-but-real colour contrast."""
    return np.sign(corr) * np.sqrt(np.abs(corr))


# ---------------------------------------------------------------------------
# Individual plot functions  (each returns silently on error/missing data)
# ---------------------------------------------------------------------------

def _plot_cells_per_perturbation(adata, plots_dir: Path, tables_dir: Path) -> pd.DataFrame:
    # Bar chart: number of cells per perturbation group.
    counts = (
        adata.obs["perturbation"]
        .astype(str)
        .value_counts()
        .rename_axis("perturbation")
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    counts.to_csv(tables_dir / "eda_cells_per_perturbation.csv", index=False)

    n_groups = len(counts)
    fig_width = max(10, min(0.35 * n_groups, 60))
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    pert_palette = _perturbation_palette(counts["perturbation"])
    colors = [pert_palette[str(p)] for p in counts["perturbation"]]
    ax.bar(range(n_groups), counts["n_cells"], color=colors, edgecolor="none", width=0.8)
    ax.set_xticks(range(n_groups))
    tick_fs = max(6, min(_TICK_FS, 140 // max(n_groups, 1)))
    ax.set_xticklabels(counts["perturbation"], rotation=90, fontsize=tick_fs)
    ax.tick_params(axis="y", labelsize=_TICK_FS)
    ax.set_ylabel("Number of cells", fontsize=_LABEL_FS)
    ax.set_title(
        f"Cells per perturbation group  "
        f"({n_groups} groups, {adata.n_obs:,} total cells)",
        fontsize=_TITLE_FS,
    )
    n_legend_cols = max(1, min(4, int(np.ceil(n_groups / 18))))
    leg_fs = max(6, min(_LEGEND_FS, 130 // max(n_groups, 1)))
    handles = [mpatches.Patch(color=pert_palette[p], label=p) for p in counts["perturbation"]]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        fontsize=leg_fs,
        title="Perturbation",
        title_fontsize=leg_fs + 1,
        ncol=n_legend_cols,
        frameon=True,
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "eda_cells_per_perturbation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return counts


def _plot_cluster_proportions(
    adata, plots_dir: Path, tables_dir: Path
) -> Optional[pd.DataFrame]:
    # Stacked bar of cell-state proportions per perturbation (top 30 groups).
    if "cell_state" not in adata.obs.columns:
        return None
    ct = pd.crosstab(
        adata.obs["perturbation"].astype(str),
        adata.obs["cell_state"].astype(str),
        normalize="index",
    )
    ct.to_csv(tables_dir / "eda_cluster_proportions.csv")

    # Limit to the 30 most-abundant perturbations for readability.
    top_perts = (
        adata.obs["perturbation"].astype(str).value_counts().head(30).index.tolist()
    )
    ct_plot = ct.loc[ct.index.isin(top_perts)]
    if ct_plot.empty:
        return ct

    states = [str(col) for col in ct_plot.columns]
    state_palette = dict(zip(states, sns.color_palette("Spectral", max(2, len(states)))))

    n_perts = len(ct_plot)
    fig_width = max(10, min(0.45 * n_perts, 50))
    fig, ax = plt.subplots(figsize=(fig_width, 5))
    bottom = np.zeros(n_perts)
    for i, col in enumerate(ct_plot.columns):
        state = str(col)
        ax.bar(
            range(n_perts),
            ct_plot[col].values,
            bottom=bottom,
            label=state,
            color=state_palette[state],
            edgecolor="none",
        )
        bottom += ct_plot[col].values

    ax.set_xticks(range(n_perts))
    ax.set_xticklabels(ct_plot.index.tolist(), rotation=90, fontsize=_TICK_FS)
    ax.tick_params(axis="y", labelsize=_TICK_FS)
    ax.set_ylabel("Proportion", fontsize=_LABEL_FS)
    ax.set_title("Cell-state cluster proportions per perturbation (top 30 groups)", fontsize=_TITLE_FS)
    ax.legend(
        title="state",
        title_fontsize=_LEGEND_FS,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=_LEGEND_FS,
    )
    fig.tight_layout()
    fig.savefig(
        plots_dir / "eda_cluster_proportions.png", dpi=150, bbox_inches="tight"
    )
    fig.savefig(
        plots_dir / "eda_cluster_proportions_states_spectral.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)
    return ct


def _plot_gene_by_cell_heatmap(
    adata, plots_dir: Path, n_top_genes: int, max_cells: int
) -> None:
    # Heatmap: top HVGs (rows) x all cells (columns).
    # Columns are annotated by perturbation group via a colour sidebar.
    gene_idx = _eda_gene_indices(adata, n_top_genes)
    gene_names = adata.var_names[gene_idx].tolist()
    if not gene_names:
        return
    n_cells_total = adata.n_obs
    if n_cells_total > max_cells:
        rng = np.random.default_rng(0)
        cell_idx = np.sort(rng.choice(n_cells_total, size=max_cells, replace=False))
    else:
        cell_idx = np.arange(adata.n_obs)
    n_cells = len(cell_idx)

    x = _to_dense(adata.X[cell_idx, :][:, gene_idx])
    # z-score each gene across sampled cells for visual comparability.
    with np.errstate(invalid="ignore"):
        mean = x.mean(axis=0, keepdims=True)
        std = x.std(axis=0, keepdims=True) + 1e-6
        x_z = np.clip((x - mean) / std, -3, 3)

    cell_names = [str(adata.obs_names[i]) for i in cell_idx]
    cell_df = pd.DataFrame(x_z.T, index=gene_names, columns=cell_names)

    # Colour sidebar: one colour per perturbation.
    perts = adata.obs["perturbation"].astype(str).values[cell_idx]
    uniq = sorted(set(perts))
    n_uniq = len(uniq)
    pert_palette = _perturbation_palette(uniq)
    col_colors = pd.Series(perts, index=cell_names, name="perturbation").map(pert_palette)

    fig_w = 10.0
    fig_h = 13.0
    try:
        g = sns.clustermap(
            cell_df,
            col_colors=col_colors,
            col_cluster=True,
            row_cluster=True,
            cmap="viridis",
            vmin=-3,
            vmax=3,
            xticklabels=False,
            yticklabels=False,
            figsize=(fig_w, fig_h),
            cbar_kws={"label": "z-score", "shrink": 0.5},
        )
        g.ax_heatmap.set_xlabel("")
        g.ax_heatmap.set_ylabel(f"Top {len(gene_names)} HVGs", fontsize=_LABEL_FS + 4)
        g.ax_heatmap.tick_params(axis="x", length=0)
        g.ax_heatmap.tick_params(axis="y", labelsize=_TICK_FS + 4)
        g.ax_col_dendrogram.set_title("Gene × Cell Expression", fontsize=_TITLE_FS + 6, pad=10)
        # Resize colorbar label
        g.cax.yaxis.label.set_size(_CBAR_FS + 4)
        g.cax.tick_params(labelsize=_TICK_FS + 3)
        # Legend for perturbation row-colour sidebar
        n_legend_cols = max(1, int(np.ceil(n_uniq / 12)))
        leg_fs = max(10, min(_LEGEND_FS + 3, 180 // max(n_uniq, 1)))
        legend_handles = [mpatches.Patch(color=pert_palette[p], label=p) for p in uniq]
        g.fig.legend(
            handles=legend_handles,
            title="Perturbation",
            title_fontsize=leg_fs + 2,
            loc="lower left",
            fontsize=leg_fs,
            ncol=n_legend_cols,
            frameon=True,
        )
        g.fig.savefig(
            plots_dir / "eda_gene_by_cell_heatmap.png", dpi=120, bbox_inches="tight"
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")


def _plot_gene_by_perturbation_heatmap(
    adata, plots_dir: Path, n_top_genes: int
) -> None:
    # Clustered heatmap: top HVGs (rows) x perturbations (columns).
    # Values are z-scored pseudobulk means so colours are comparable.
    gene_idx = _eda_gene_indices(adata, n_top_genes)
    label_genes = set(adata.var_names[_top_hvg_indices(adata, min(20, adata.n_vars))].tolist())
    pb = _pseudobulk_mean(adata, gene_idx)
    if pb.empty or pb.shape[0] < 2:
        return

    # z-score each gene (column) across perturbations.
    with np.errstate(invalid="ignore"):
        pb_z = (pb - pb.mean(axis=0)) / (pb.std(axis=0) + 1e-6)
    pb_z = pb_z.clip(-3, 3)

    mat = pb_z.T
    n_rows, n_cols = mat.shape
    try:
        g = sns.clustermap(
            mat,
            col_cluster=True,
            row_cluster=True,
            cmap="viridis",
            vmin=-3,
            vmax=3,
            xticklabels=True,
            yticklabels=False,
            figsize=(max(10, min(28, n_cols * 0.6 + 5)), max(10, min(24, n_rows / 220))),
            cbar_kws={"label": "z-score (across perturbations)", "shrink": 0.5},
        )
        g.ax_heatmap.set_xlabel("Perturbation", fontsize=_LABEL_FS)
        g.ax_heatmap.set_ylabel(f"Top {n_rows} HVGs", fontsize=_LABEL_FS)
        g.ax_heatmap.tick_params(axis="x", labelsize=_TICK_FS)
        g.ax_heatmap.tick_params(axis="y", labelsize=max(_TICK_FS, 13))
        g.ax_col_dendrogram.set_title("Gene × Perturbation pseudobulk expression heatmap", fontsize=_TITLE_FS, pad=8)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        g.fig.savefig(
            plots_dir / "eda_gene_by_perturbation_heatmap.png",
            dpi=120,
            bbox_inches="tight",
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")


def _plot_gene_correlation_heatmap(
    adata, plots_dir: Path, n_top_genes: int,
    gene_cluster_s: pd.Series | None = None,
    n_modules: int = 8,
) -> None:
    # Top-HVG gene-gene Pearson correlation across all cells. Module side bars
    # use the gene clusters from `_plot_clustered_gene_pert_heatmap` (passed in
    # via ``gene_cluster_s``) so the colour groupings on this heatmap match
    # eda_gene_clusters.csv exactly. Falls back to ad-hoc clustering on
    # 1-|r| distance only when gene_cluster_s is unavailable.
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    gene_idx = _eda_gene_indices(adata, n_top_genes)
    gene_names = adata.var_names[gene_idx].tolist()
    if len(gene_names) < 2:
        return

    x = _to_dense(adata.X[:, :][:, gene_idx]).astype(np.float32, copy=False)
    corr = np.corrcoef(x.T)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)

    n_g = len(gene_names)

    # ---- Use the EDA-shared gene_cluster assignments if provided ----
    labels: np.ndarray | None = None
    if gene_cluster_s is not None and not gene_cluster_s.empty:
        try:
            mapped = gene_cluster_s.reindex(gene_names)
            if mapped.notna().all():
                labels = mapped.astype(int).to_numpy()
        except Exception:
            labels = None

    if labels is None:
        # Fallback: ad-hoc hierarchical clustering on 1 − |r|.
        try:
            dist = np.clip(1.0 - np.abs(corr), 0.0, None)
            np.fill_diagonal(dist, 0.0)
            condensed = squareform(dist, checks=False)
            Z = linkage(condensed, method="average")
            n_clust = max(2, min(n_modules, n_g - 1))
            labels = fcluster(Z, t=n_clust, criterion="maxclust")
        except Exception:
            labels = np.ones(n_g, dtype=int)
    order = np.argsort(labels, kind="stable")

    corr_sorted = corr[np.ix_(order, order)]
    labels_sorted = labels[order]
    n_modules_real = int(labels_sorted.max())
    palette = sns.color_palette("tab10", min(n_modules_real, 10))
    module_colors = np.array(
        [palette[(int(c) - 1) % len(palette)] for c in labels_sorted]
    )

    # ---- Square layout: square heatmap with thin top + left module strips
    #      and a dedicated legend row at the bottom (kept inside the canvas).
    try:
        # Bigger fonts than the EDA defaults so the heatmap reads well at
        # native size and when downscaled in the report's heatmaps tab.
        _TITLE_FS_HM   = 24
        _LABEL_FS_HM   = 20
        _CBAR_FS_HM    = 18
        _LEGEND_FS_HM  = 18

        fig = plt.figure(figsize=(15, 15))
        gs = fig.add_gridspec(
            3, 3,
            # Thin top + left module strips so the side bars don't dominate.
            width_ratios=[0.035, 1.0, 0.05],
            height_ratios=[0.035, 1.0, 0.18],
            wspace=0.025, hspace=0.04,
        )
        ax_top   = fig.add_subplot(gs[0, 1])   # top module strip
        ax_left  = fig.add_subplot(gs[1, 0])   # left module strip
        ax_main  = fig.add_subplot(gs[1, 1])   # main square heatmap
        ax_cbar  = fig.add_subplot(gs[1, 2])   # colour-bar slot
        ax_leg   = fig.add_subplot(gs[2, :])   # legend row (full width)
        ax_leg.axis("off")

        im = ax_main.imshow(
            _corr_for_display(corr_sorted),
            cmap="RdBu_r",
            vmin=-0.6,
            vmax=0.6,
            # aspect="auto" so the heatmap fills its gridspec slot exactly
            # → top + left module strips line up with the heatmap edges.
            aspect="auto",
            interpolation="nearest",
        )
        ax_main.set_xticks([])
        ax_main.set_yticks([])
        ax_main.set_xlabel("Genes (sorted by module)", fontsize=_LABEL_FS_HM)
        ax_main.set_ylabel("Genes (sorted by module)", fontsize=_LABEL_FS_HM)
        boundaries = np.where(labels_sorted[1:] != labels_sorted[:-1])[0] + 0.5
        for b in boundaries:
            ax_main.axhline(b, color="#202020", lw=0.7, alpha=0.7)
            ax_main.axvline(b, color="#202020", lw=0.7, alpha=0.7)

        ax_top.imshow(module_colors[np.newaxis, :, :], aspect="auto")
        ax_top.set_xticks([]); ax_top.set_yticks([])
        ax_top.set_title(
            f"Gene–Gene Correlation (all cells, all perturbations) — "
            f"top {n_g} HVGs, {n_modules_real} modules",
            fontsize=_TITLE_FS_HM, pad=12,
        )

        ax_left.imshow(module_colors[:, np.newaxis, :], aspect="auto")
        ax_left.set_xticks([]); ax_left.set_yticks([])
        # Row title intentionally omitted; module identity is shown by the
        # legend row at the bottom.

        cbar = fig.colorbar(im, cax=ax_cbar, label="Pearson r (contrast-enhanced)")
        cbar.ax.tick_params(labelsize=_CBAR_FS_HM)
        cbar.ax.yaxis.label.set_size(_CBAR_FS_HM)

        # Module legend laid out in 2–3 rows so it stays compact and is
        # placed inside the dedicated bottom-row axis (never clipped).
        handles = [
            mpatches.Patch(color=palette[(c - 1) % len(palette)],
                           label=f"Module {c}")
            for c in range(1, n_modules_real + 1)
        ]
        # Pick ncol so the legend uses 2 or 3 rows: ceil(n/3) ≤ rows ≤ ceil(n/2).
        ncol = max(1, (n_modules_real + 1) // 2)        # → 2 rows for 4–8 modules
        if n_modules_real >= 7:
            ncol = max(1, (n_modules_real + 2) // 3)    # → 3 rows for 7+ modules
        ax_leg.legend(
            handles=handles, title="Gene module",
            loc="upper center", ncol=ncol,
            fontsize=_LEGEND_FS_HM, title_fontsize=_LEGEND_FS_HM + 1,
            frameon=True, handlelength=2.0, handletextpad=0.6,
            columnspacing=1.4, labelspacing=0.5,
        )

        fig.savefig(
            plots_dir / "eda_gene_correlation_heatmap.png",
            dpi=120,
        )
        plt.close(fig)
    except Exception:
        plt.close("all")


def _plot_umap_by_cell_state(adata, plots_dir: Path) -> None:
    # UMAP scatter coloured by cell_state cluster label.
    if "X_umap" not in adata.obsm or "cell_state" not in adata.obs.columns:
        return
    um = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    um["cell_state"] = adata.obs["cell_state"].astype(str).values
    states = sorted(um["cell_state"].unique())
    n_s = len(states)
    pal = sns.color_palette("Spectral", max(2, n_s))
    palette = dict(zip(states, pal))

    _LBL_UMAP_FS = 36
    _TICK_UMAP_FS = 30
    _TITLE_UMAP_FS = 34
    # Fixed legend font so every UMAP (eda + report) renders identically.
    _LEG_UMAP_FS = 24

    fig, ax = plt.subplots(figsize=(11, 11))
    for s in states:
        mask = um["cell_state"] == s
        ax.scatter(
            um.loc[mask, "UMAP1"],
            um.loc[mask, "UMAP2"],
            s=3,
            alpha=0.65,
            color=palette[s],
            label=s,
            edgecolors="none",
            rasterized=True,
        )
    # Annotate each cluster with its label at the centroid position.
    for s in states:
        mask = um["cell_state"] == s
        cx = um.loc[mask, "UMAP1"].mean()
        cy = um.loc[mask, "UMAP2"].mean()
        ax.text(
            cx, cy, s,
            fontsize=14,
            fontweight="bold",
            ha="center", va="center",
            color="white",
            bbox=dict(boxstyle="round,pad=0.25", fc=palette[s], ec="none", alpha=0.9),
        )
    leg = ax.legend(
        title="cluster",
        title_fontsize=_LEG_UMAP_FS + 1,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=_LEG_UMAP_FS,
        markerscale=3,
        frameon=True,
        framealpha=0.9,
    )
    leg.get_frame().set_edgecolor("#888888")
    ax.set_xlabel("UMAP 1", fontsize=_LBL_UMAP_FS)
    ax.set_ylabel("UMAP 2", fontsize=_LBL_UMAP_FS)
    ax.tick_params(labelsize=_TICK_UMAP_FS)
    ax.set_aspect("equal", adjustable="box")  # square plot box
    ax.set_title("UMAP colored by cell_state (cluster)", fontsize=_TITLE_UMAP_FS)
    fig.tight_layout(rect=[0, 0, 0.78, 1])
    fig.savefig(plots_dir / "umap_cell_state.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_umap_by_perturbation(adata, plots_dir: Path) -> None:
    # UMAP scatter coloured by perturbation label. Control rendered with the
    # neutral grey and plotted first so non-control points sit on top.
    if "X_umap" not in adata.obsm or "perturbation" not in adata.obs.columns:
        return
    um = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    um["perturbation"] = adata.obs["perturbation"].astype(str).values
    perts = sorted(um["perturbation"].unique())
    # Plot control first (so non-control sits on top); palette already maps
    # control → _CONTROL_COLOR.
    ctrl_keys = {"control", "ctrl", "nontargeting", "non-targeting", "nt",
                 "scramble", "safe-targeting", "safe_targeting"}
    perts_ordered = [p for p in perts if p.lower() in ctrl_keys] + \
                    [p for p in perts if p.lower() not in ctrl_keys]
    n_p = len(perts)
    palette = _spectral_palette(perts)

    # Legend lives outside the axes so it never covers the UMAP cloud.
    _LBL_UMAP_FS = 36
    _TICK_UMAP_FS = 30
    _TITLE_UMAP_FS = 34
    # Fixed legend font so every UMAP (eda + report) renders identically.
    _LEG_UMAP_FS = 24
    n_leg_cols = max(1, min(2, n_p // 25 + 1))

    fig, ax = plt.subplots(figsize=(11, 11))
    for p in perts_ordered:
        mask = um["perturbation"] == p
        ax.scatter(
            um.loc[mask, "UMAP1"],
            um.loc[mask, "UMAP2"],
            s=3,
            alpha=0.65,
            color=palette[p],
            label=p,
            edgecolors="none",
            rasterized=True,
        )
    leg = ax.legend(
        title="perturbation",
        title_fontsize=_LEG_UMAP_FS + 1,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=_LEG_UMAP_FS,
        markerscale=3,
        ncol=n_leg_cols,
        frameon=True,
        framealpha=0.9,
    )
    leg.get_frame().set_edgecolor("#888888")
    ax.set_xlabel("UMAP 1", fontsize=_LBL_UMAP_FS)
    ax.set_ylabel("UMAP 2", fontsize=_LBL_UMAP_FS)
    ax.tick_params(labelsize=_TICK_UMAP_FS)
    ax.set_aspect("equal", adjustable="box")  # square plot box
    ax.set_title("UMAP colored by perturbation", fontsize=_TITLE_UMAP_FS)
    fig.tight_layout(rect=[0, 0, 0.76, 1])
    fig.savefig(plots_dir / "eda_umap_perturbation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_perturbation_similarity_heatmap(
    adata, plots_dir: Path, n_top_genes: int
) -> None:
    """Clustered heatmap of pairwise perturbation cosine similarity.

    Computed from pseudobulk mean expression profiles of the top HVGs.
    Perturbations that cluster together in this heatmap share similar
    transcriptional effects (i.e. are functionally related).
    """
    gene_idx = _eda_gene_indices(adata, n_top_genes)
    pb = _pseudobulk_mean(adata, gene_idx)
    if pb.empty or pb.shape[0] < 3:
        return

    pb_z = _zscore_columns(pb)
    mat = pb_z.values  # (n_perturbations, genes z-scored across perturbations)
    norms = np.linalg.norm(mat, axis=1, keepdims=True) + 1e-12
    mat_norm = mat / norms
    sim = mat_norm @ mat_norm.T          # (n_perturbations, n_perturbations)
    np.fill_diagonal(sim, 1.0)
    sim_df = pd.DataFrame(sim, index=pb.index, columns=pb.index)

    n_p = len(sim_df)
    # More compact square footprint and tighter cell size.
    cell_size = max(0.3, min(0.55, 18 / max(n_p, 1)))
    fig_side = max(6, min(9, n_p * cell_size))
    try:
        g = sns.clustermap(
            sim_df,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            xticklabels=(n_p <= 80),
            yticklabels=(n_p <= 80),
            figsize=(fig_side, fig_side),
            cbar_kws={"label": "Cosine similarity", "shrink": 0.5},
        )
        # Larger floor and ceiling so labels stay readable even though the
        # canvas got smaller.
        tick_fs = max(20, min(28, 460 // max(n_p, 1)))
        g.ax_col_dendrogram.set_title(
            f"Perturbation–perturbation cosine similarity  (z-scored pseudobulk, top {len(gene_idx)} HVGs)",
            fontsize=_TITLE_FS,
            pad=8,
        )
        g.ax_heatmap.set_xlabel("Perturbation", fontsize=_LABEL_FS)
        g.ax_heatmap.set_ylabel("Perturbation", fontsize=_LABEL_FS)
        g.ax_heatmap.tick_params(axis="x", labelsize=tick_fs, rotation=45)
        for lab in g.ax_heatmap.get_xticklabels():
            lab.set_ha("right")
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_fs)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        g.fig.savefig(
            plots_dir / "eda_perturbation_similarity.png",
            dpi=120,
            bbox_inches="tight",
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")


# ---------------------------------------------------------------------------
# Per-perturbation gene–gene correlation heatmaps
# ---------------------------------------------------------------------------

def _plot_gene_corr_per_perturbation(
    adata, plots_dir: Path, n_top_genes: int
) -> None:
    """For each perturbation, save a 3-panel PNG:
      left  = control gene–gene Pearson correlation (same gene order)
      centre = perturbation gene–gene Pearson correlation
      right  = difference (pert − ctrl)
    Output: eda_gene_corr_vs_<pert>.png
    """
    gene_idx = _eda_gene_indices(adata, n_top_genes)
    gene_names = adata.var_names[gene_idx].tolist()
    if len(gene_names) < 2:
        return

    perts_all = adata.obs["perturbation"].astype(str).values

    # Identify control group (most common label matching known keywords)
    _ctrl_kw = {"control", "ctrl", "nontargeting", "non-targeting", "nt",
                "scramble", "safe-targeting", "safe_targeting"}
    ctrl_mask = np.array([p.lower() in _ctrl_kw for p in perts_all])
    if ctrl_mask.sum() < 3:
        groups, counts = np.unique(perts_all, return_counts=True)
        ctrl_label = groups[np.argmax(counts)]
        ctrl_mask = perts_all == ctrl_label

    ctrl_x = _to_dense(adata.X[np.where(ctrl_mask)[0], :][:, gene_idx])
    ctrl_corr = np.corrcoef(ctrl_x.T)
    ctrl_corr = np.nan_to_num(ctrl_corr, nan=0.0)
    np.fill_diagonal(ctrl_corr, 1.0)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        dist = np.clip(1.0 - np.abs(ctrl_corr), 0.0, None)
        np.fill_diagonal(dist, 0.0)
        order = leaves_list(linkage(squareform(dist, checks=False), method="average"))
    except Exception:
        order = np.arange(len(gene_names))

    n_g = len(gene_names)

    unique_perts = [p for p in sorted(set(perts_all)) if not (p.lower() in _ctrl_kw)]
    for pert in unique_perts:
        pert_mask = perts_all == pert
        if pert_mask.sum() < 3:
            continue

        pert_x = _to_dense(adata.X[np.where(pert_mask)[0], :][:, gene_idx])
        pert_corr = np.corrcoef(pert_x.T)
        pert_corr = np.nan_to_num(pert_corr, nan=0.0)
        np.fill_diagonal(pert_corr, 1.0)
        ctrl_ord = ctrl_corr[np.ix_(order, order)]
        pert_ord = pert_corr[np.ix_(order, order)]
        diff_ord = pert_ord - ctrl_ord
        diff_lim = float(np.nanpercentile(np.abs(diff_ord), 98))
        diff_lim = max(0.15, min(0.75, diff_lim))

        try:
            fig_h = 6
            fig, axes = plt.subplots(1, 3, figsize=(fig_h * 3, fig_h))
            for ax, mat, title, cmap, vmin, vmax, cbar_lbl in [
                (axes[0], _corr_for_display(ctrl_ord), "Control",       "RdBu_r",  -1, 1, "Pearson r (contrast)"),
                (axes[1], _corr_for_display(pert_ord), pert,             "RdBu_r",  -1, 1, "Pearson r (contrast)"),
                (axes[2], diff_ord, f"\u0394 ({pert}\u2212Ctrl)", "coolwarm", -diff_lim, diff_lim, "\u0394r"),
            ]:
                im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax,
                               aspect="equal", interpolation="nearest")
                ax.set_title(title, fontsize=_LABEL_FS + 1, pad=8)
                ax.set_xticks([])
                ax.set_yticks([])
                fig.colorbar(im, ax=ax, shrink=0.6, label=cbar_lbl)

            fig.suptitle(
                f"Gene\u2013Gene Correlation: Control vs {pert}  (top {n_g} HVGs)",
                fontsize=_TITLE_FS, y=1.01,
            )
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            safe = pert.replace("/", "_").replace(" ", "_").replace("+", "_")
            fig.savefig(
                plots_dir / f"eda_gene_corr_vs_{safe}.png",
                dpi=120, bbox_inches="tight",
            )
        except Exception:
            pass
        finally:
            plt.close("all")


# ---------------------------------------------------------------------------
# Clustered gene × perturbation heatmap with cluster annotations
# ---------------------------------------------------------------------------

def _plot_clustered_gene_pert_heatmap(
    adata,
    plots_dir: Path,
    tables_dir: Path,
    n_top_genes: int,
    n_gene_clusters: int = 10,
    n_pert_clusters: int = 10,
):
    """Top-HVG × perturbations clustermap with coloured cluster sidebars.

    Rows = top HVGs, Columns = perturbations. Both axes are clustered with
    Ward hierarchical clustering and then flat-cut into N clusters; the cluster
    membership is encoded as coloured sidebars (left = gene clusters,
    top = perturbation clusters).

    Saves
    -----
    plots/eda_clustered_gene_by_pert_heatmap.png
    csv/eda_gene_clusters.csv      -- one row per gene: gene, gene_cluster
    csv/eda_pert_clusters.csv      -- one row per pert: perturbation, pert_cluster

    Returns
    -------
    (gene_cluster_series, pert_cluster_series) — used by the summary heatmap.
    None, None on failure.
    """
    from scipy.cluster.hierarchy import linkage, fcluster

    gene_idx = _eda_gene_indices(adata, n_top_genes)
    pb = _pseudobulk_mean(adata, gene_idx)          # (perturbations × genes)
    if pb.empty or pb.shape[0] < 2 or pb.shape[1] < 2:
        return None, None

    # z-score each gene (column) across perturbations
    with np.errstate(invalid="ignore"):
        pb_z = (pb - pb.mean(axis=0)) / (pb.std(axis=0) + 1e-6)
    pb_z = pb_z.clip(-3, 3)

    # Transpose so rows = genes, cols = perturbations  (Figure A orientation)
    mat = pb_z.T   # (genes × perturbations)

    n_gc = min(n_gene_clusters, mat.shape[0])
    n_pc = min(n_pert_clusters, mat.shape[1])

    # Hierarchical clustering for cluster-label assignment
    try:
        Z_genes = linkage(mat.values, method="ward")
        gene_labels = fcluster(Z_genes, t=n_gc, criterion="maxclust")
    except Exception:
        gene_labels = np.ones(mat.shape[0], dtype=int)

    try:
        Z_perts = linkage(mat.values.T, method="ward")
        pert_labels = fcluster(Z_perts, t=n_pc, criterion="maxclust")
    except Exception:
        pert_labels = np.ones(mat.shape[1], dtype=int)
    control_mask = np.array([str(p).lower() == "control" for p in mat.columns])
    if control_mask.any():
        pert_labels = np.asarray(pert_labels, dtype=int) + 1
        pert_labels[control_mask] = 1

    gene_cluster_s = pd.Series(gene_labels, index=mat.index, name="gene_cluster")
    pert_cluster_s = pd.Series(pert_labels, index=mat.columns, name="pert_cluster")

    # Save cluster assignments to CSV
    (gene_cluster_s
     .rename_axis("gene")
     .reset_index()
     .to_csv(tables_dir / "eda_gene_clusters.csv", index=False))
    (pert_cluster_s
     .rename_axis("perturbation")
     .reset_index()
     .to_csv(tables_dir / "eda_pert_clusters.csv", index=False))

    # Build cluster colour palettes
    n_gc_uniq = gene_cluster_s.nunique()
    n_pc_uniq = pert_cluster_s.nunique()
    gc_pal = dict(zip(
        sorted(gene_cluster_s.unique()),
        (sns.color_palette("tab20", min(n_gc_uniq, 20)) * 2)[:n_gc_uniq],
    ))
    pc_pal = dict(zip(
        sorted(pert_cluster_s.unique()),
        (sns.color_palette("Set2", 8) + sns.color_palette("tab20", 20))[:n_pc_uniq],
    ))
    if control_mask.any() and 1 in pc_pal:
        pc_pal[1] = _CONTROL_COLOR

    row_colors_s = gene_cluster_s.map(lambda c: gc_pal.get(c, "#aaaaaa"))
    row_colors_s.name = "Gene cluster"
    col_colors_s = pert_cluster_s.map(lambda c: pc_pal.get(c, "#aaaaaa"))
    col_colors_s.name = "Pert cluster"

    n_rows, n_cols = mat.shape
    try:
        tick_col = n_cols <= 80
        tick_fs = max(11, min(18, 1200 // max(n_rows, 1)))
        fig_w = 12.5
        fig_h = 18.0  # extra vertical room so legends never overlap x-ticks
        g = sns.clustermap(
            mat,
            row_colors=row_colors_s,
            col_colors=col_colors_s,
            col_cluster=True,
            row_cluster=True,
            cmap="RdBu_r",
            center=0,
            vmin=-3,
            vmax=3,
            xticklabels=tick_col,
            yticklabels=False,
            figsize=(fig_w, fig_h),
            # Narrower z-score colour-scale bar.
            cbar_pos=(0.02, 0.84, 0.012, 0.14),
            cbar_kws={"label": "z-score"},
        )
        g.ax_heatmap.set_xlabel(f"{n_cols} perturbations", fontsize=_LABEL_FS + 10)
        g.ax_heatmap.set_ylabel(f"Top {n_rows} HVGs", fontsize=_LABEL_FS + 10)
        g.ax_heatmap.tick_params(axis="x", labelsize=max(tick_fs + 8, _TICK_FS + 9), rotation=90)
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_fs + 6)
        g.ax_col_dendrogram.set_title(
            "All Genes \u00d7 Perturbations Clusters",
            fontsize=_TITLE_FS + 12, pad=14,
        )
        g.cax.yaxis.label.set_size(_CBAR_FS + 8)
        g.cax.tick_params(labelsize=_TICK_FS + 6)
        # Bigger annotation text on the row + column module side bars.
        if hasattr(g, "ax_row_colors") and g.ax_row_colors is not None:
            g.ax_row_colors.tick_params(labelsize=_LABEL_FS + 6)
            for lbl in g.ax_row_colors.get_xticklabels():
                lbl.set_rotation(0); lbl.set_ha("right")
        if hasattr(g, "ax_col_colors") and g.ax_col_colors is not None:
            g.ax_col_colors.tick_params(labelsize=_LABEL_FS + 6)

        # Cluster colour legends live below the clustermap so they never cover
        # labels, dendrograms, or the heatmap body.
        gene_handles = [
            mpatches.Patch(color=gc_pal[c], label=f"GC{c}")
            for c in sorted(gc_pal)
        ]
        pert_handles = [
            mpatches.Patch(color=pc_pal[c], label=("control" if c == 1 and control_mask.any() else f"PC{c}"))
            for c in sorted(pc_pal)
        ]
        # Reserve a deeper bottom strip so the rotated x-tick labels
        # (perturbations) and the two legends never overlap.
        g.fig.subplots_adjust(bottom=0.32)
        leg1 = g.fig.legend(
            handles=gene_handles, title="Gene cluster",
            loc="lower center", bbox_to_anchor=(0.5, 0.085),
            ncol=min(8, n_gc_uniq), fontsize=19, title_fontsize=21, frameon=True,
        )
        g.fig.legend(
            handles=pert_handles, title="Pert cluster",
            loc="lower center", bbox_to_anchor=(0.5, 0.005),
            ncol=min(8, n_pc_uniq), fontsize=19, title_fontsize=21, frameon=True,
        )
        g.fig.add_artist(leg1)

        g.fig.savefig(
            plots_dir / "eda_clustered_gene_by_pert_heatmap.png",
            dpi=120, bbox_inches="tight",
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")

    return gene_cluster_s, pert_cluster_s


def _plot_gene_pert_cluster_summary_heatmap(
    adata,
    plots_dir: Path,
    n_top_genes: int,
    gene_cluster_s,
    pert_cluster_s,
) -> None:
    """Gene-expression cluster × perturbation cluster summary heatmap.

    Each cell shows the mean z-score of all genes in gene-cluster GCi and all
    perturbations in perturbation-cluster PCj (Figure B orientation).
    Values are annotated inside each cell.

    Saves: plots/eda_gene_pert_cluster_summary_heatmap.png
    """
    if gene_cluster_s is None or pert_cluster_s is None:
        return
    if not isinstance(gene_cluster_s, pd.Series) or gene_cluster_s.empty:
        return
    if not isinstance(pert_cluster_s, pd.Series) or pert_cluster_s.empty:
        return

    gene_idx = _eda_gene_indices(adata, n_top_genes)
    pb = _pseudobulk_mean(adata, gene_idx)          # (perturbations × genes)
    if pb.empty:
        return

    with np.errstate(invalid="ignore"):
        pb_z = (pb - pb.mean(axis=0)) / (pb.std(axis=0) + 1e-6)
    pb_z = pb_z.clip(-3, 3)

    mat = pb_z.T   # (genes × perturbations)

    genes_common = mat.index.intersection(gene_cluster_s.index)
    perts_common = mat.columns.intersection(pert_cluster_s.index)
    if genes_common.empty or perts_common.empty:
        return

    mat    = mat.loc[genes_common, perts_common]
    gc     = gene_cluster_s.loc[genes_common]
    pc     = pert_cluster_s.loc[perts_common]

    gene_clust_ids = sorted(gc.unique())
    pert_clust_ids = sorted(pc.unique())
    n_gc = len(gene_clust_ids)
    n_pc = len(pert_clust_ids)

    # Build summary matrix: mean z-score per (gene_cluster, pert_cluster)
    summary = np.zeros((n_gc, n_pc))
    for i, gc_id in enumerate(gene_clust_ids):
        gene_mask = gc == gc_id
        for j, pc_id in enumerate(pert_clust_ids):
            pert_mask = pc == pc_id
            vals = mat.loc[gene_mask, pert_mask].values
            summary[i, j] = float(vals.mean()) if vals.size > 0 else 0.0

    vabs = max(0.1, float(np.abs(summary).max()))

    fig_w = max(8, n_pc * 0.75 + 2)
    fig_h = max(6, n_gc * 0.6 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(summary, cmap="RdBu_r", aspect="auto",
                   vmin=-vabs, vmax=vabs, interpolation="nearest")

    ax.set_xticks(range(n_pc))
    ax.set_xticklabels([f"PC{c}" for c in pert_clust_ids],
                       rotation=45, ha="right", fontsize=_TICK_FS)
    ax.set_yticks(range(n_gc))
    ax.set_yticklabels([f"GC{c}" for c in gene_clust_ids], fontsize=_TICK_FS)
    ax.set_xlabel("Perturbation Clusters", fontsize=_LABEL_FS)
    ax.set_ylabel("Gene Expression Clusters", fontsize=_LABEL_FS)
    ax.set_title(
        f"Gene Expression Cluster \u00d7 Perturbation Cluster Summary\n"
        f"(mean z-score | {n_gc} gene clusters \u00d7 {n_pc} perturbation clusters)",
        fontsize=14,
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, label="mean z-score")
    cbar.ax.tick_params(labelsize=_CBAR_FS)

    # Annotate each cell with its mean z-score value
    annot_fs = max(7, min(_ANNOT_FS, 70 // max(n_gc, n_pc, 1)))
    thresh = vabs * 0.55
    for i in range(n_gc):
        for j in range(n_pc):
            val = summary[i, j]
            color = "white" if abs(val) > thresh else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=annot_fs, color=color)

    fig.tight_layout()
    fig.savefig(
        plots_dir / "eda_gene_pert_cluster_summary_heatmap.png",
        dpi=150, bbox_inches="tight",
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_eda(
    adata,
    output_dir: str,
    n_top_genes: int = 500,
    max_cells_heatmap: int = 500,
    n_gene_clusters: int = 10,
    n_pert_clusters: int = 10,
) -> pd.DataFrame:
    """Run the full EDA suite and write all plots and tables to output_dir.

    Args:
        adata              -- AnnData after QC + normalisation (at minimum).
                             If cell_state and X_umap are present (after
                             preprocess step) the cluster plots are included.
        output_dir         -- root directory; plots/ and csv/ are created.
        n_top_genes        -- number of top HVGs used in all three heatmaps.
        max_cells_heatmap  -- cap on cells sampled for the gene×cell heatmap.
        n_gene_clusters    -- number of gene clusters for clustered heatmap.
        n_pert_clusters    -- number of perturbation clusters for clustered heatmap.

    Returns:
        counts DataFrame with columns [perturbation, n_cells].
    """
    _apply_plot_style()
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    counts = _plot_cells_per_perturbation(adata, plots, tables)
    _plot_cluster_proportions(adata, plots, tables)
    _plot_gene_by_cell_heatmap(adata, plots, n_top_genes, max_cells_heatmap)
    _plot_gene_by_perturbation_heatmap(adata, plots, n_top_genes)
    # Compute clustered gene × perturbation modules first so the gene
    # correlation heatmap can reuse the same `gene_cluster` assignments.
    gene_cluster_s, pert_cluster_s = _plot_clustered_gene_pert_heatmap(
        adata, plots, tables, n_top_genes, n_gene_clusters, n_pert_clusters,
    )
    _plot_gene_correlation_heatmap(adata, plots, n_top_genes,
                                   gene_cluster_s=gene_cluster_s)
    _plot_perturbation_similarity_heatmap(adata, plots, n_top_genes)
    _plot_umap_by_cell_state(adata, plots)
    _plot_umap_by_perturbation(adata, plots)
    _plot_gene_corr_per_perturbation(adata, plots, n_top_genes)
    _plot_gene_pert_cluster_summary_heatmap(
        adata, plots, n_top_genes, gene_cluster_s, pert_cluster_s,
    )

    return counts
