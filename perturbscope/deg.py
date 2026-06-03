# =============================================================================
# perturbscope/deg.py
#
# Differential expression analysis module for Perturb-seq data.
#
# Identifies the top perturbations by transcriptional impact and runs DEG
# analysis (each perturbation vs all control cells) using Welch's t-test on
# log-normalised counts (already in adata.X after preprocessing).
# All statistics are computed in pure NumPy/Python — no scipy required.
#
# Outputs written to <output_dir>/plots/ and <output_dir>/csv/:
#
#   plots/
#     deg_volcano_<perturbation>.png        -- volcano plot (one per perturbation)
#     deg_top_perturbations_heatmap.png     -- log2FC heatmap across top perturbations
#   csv/
#     deg_<perturbation>.csv                -- full gene-level DEG table per perturbation
#                                              columns: gene, log2fc, mean_perturbed,
#                                              mean_control, pval, padj,
#                                              neg_log10_padj, significant
#     deg_summary.csv                       -- one row per perturbation:
#                                              n_cells, n_de_up, n_de_down,
#                                              n_de_total, top_gene
#
# Usage in pipeline:
#   from .deg import identify_top_perturbations, run_deg_analysis
#   top_perts = identify_top_perturbations(adata, effect_df=effect_df, n_top=5)
#   deg_results = run_deg_analysis(adata, output_dir, perturbations=top_perts)
# =============================================================================

from __future__ import annotations

from math import erf, sqrt
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from .utils import ensure_dir
from .pathways import annotate_deg, run_enrichment

# ---------------------------------------------------------------------------
# Global font constants (mirror eda.py — increase here to resize DEG plots)
# ---------------------------------------------------------------------------
_TITLE_FS  = 15
_LABEL_FS  = 13
_TICK_FS   = 11
_LEGEND_FS = 11
_CBAR_FS   = 11
_ANNOT_FS  = 10


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _to_dense(x):
    return x.A if hasattr(x, "A") else np.asarray(x)


def _two_sided_pval(z: np.ndarray) -> np.ndarray:
    # Normal-approximation two-sided p-value (avoids scipy dependency).
    # Accurate enough for gene ranking; should not be used for small t.
    out = np.empty(z.shape, dtype=float)
    for i, zi in enumerate(z.ravel()):
        out.ravel()[i] = 2.0 * (1.0 - 0.5 * (1.0 + erf(abs(zi) / sqrt(2.0))))
    return out


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    # Benjamini-Hochberg FDR correction (pure NumPy; no scipy).
    n = len(pvals)
    if n == 0:
        return pvals.copy()
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    bh_raw = np.minimum(1.0, pvals * n / ranks)
    # Enforce monotonicity via reverse cumulative minimum.
    bh_sorted = bh_raw[order]
    cummin = np.minimum.accumulate(bh_sorted[::-1])[::-1]
    result = np.empty(n)
    result[order] = cummin
    return result


def _compute_deg(adata, perturbation: str) -> pd.DataFrame:
    # Compute genome-wide DEG stats for one perturbation vs all controls.
    #
    # Method: Welch's t-test on log-normalised counts (already in adata.X
    # after the preprocessing step).  p-values are adjusted with BH-FDR.
    #
    # Returns a DataFrame sorted by |log2fc| descending with columns:
    #   gene, log2fc, mean_perturbed, mean_control, pval, padj,
    #   neg_log10_padj, significant
    p_mask = adata.obs["perturbation"].values == perturbation
    c_mask = adata.obs["is_control"].values
    if p_mask.sum() < 3 or c_mask.sum() < 3:
        return pd.DataFrame()

    xp = _to_dense(adata.X[p_mask, :])
    xc = _to_dense(adata.X[c_mask, :])

    mp = xp.mean(axis=0)
    mc = xc.mean(axis=0)
    # adata.X already holds log1p-normalised values (from sc.pp.log1p in
    # preprocessing).  The correct log2 fold-change is therefore the
    # difference in log-space converted from natural-log to log2.
    # Using np.log2((mp+eps)/(mc+eps)) would treat log-space means as raw
    # counts and produce spuriously extreme values (up to ±9) when one group
    # mean is 0.  The formula below gives biologically plausible ±2–4 FC.
    log2fc = (mp - mc) / np.log(2)

    # Welch variance: var_p/n_p + var_c/n_c
    vp = xp.var(axis=0, ddof=1)
    vc = xc.var(axis=0, ddof=1)
    se = np.sqrt(vp / max(int(p_mask.sum()), 1) + vc / max(int(c_mask.sum()), 1)) + 1e-12
    t_stat = (mp - mc) / se

    pvals = _two_sided_pval(t_stat)
    padj = _benjamini_hochberg(pvals)

    df = pd.DataFrame(
        {
            "gene": adata.var_names.tolist(),
            "log2fc": log2fc,
            "mean_perturbed": mp,
            "mean_control": mc,
            "pval": pvals,
            "padj": padj,
            "neg_log10_padj": -np.log10(np.maximum(padj, 1e-300)),
        }
    )
    df["significant"] = (df["padj"] < 0.05) & (df["log2fc"].abs() >= 0.5)
    df = df.sort_values("log2fc", key=abs, ascending=False).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Top-perturbation selection
# ---------------------------------------------------------------------------

def identify_top_perturbations(
    adata,
    effect_df: Optional[pd.DataFrame] = None,
    n_top: int = 5,
) -> List[str]:
    """Return the n_top most impactful non-control perturbation names.

    Ranking priority:
      1. If effect_df (from compute_effect_decomposition) is provided, rank by
         |transcriptional_score| + |state_shift_score| (combined effect size).
      2. Otherwise fall back to ranking by cell count (most abundant groups).

    Controls are always excluded.
    """
    non_ctrl = [
        p
        for p in adata.obs["perturbation"].astype(str).unique()
        if str(p).lower() != "control"
    ]
    if not non_ctrl:
        return []

    if (
        effect_df is not None
        and not effect_df.empty
        and "perturbation" in effect_df.columns
    ):
        sub = effect_df[effect_df["perturbation"].isin(non_ctrl)].copy()
        for col in ("transcriptional_score", "state_shift_score"):
            if col not in sub.columns:
                sub[col] = 0.0
        sub["_rank"] = sub["transcriptional_score"].abs() + sub["state_shift_score"].abs()
        return sub.sort_values("_rank", ascending=False)["perturbation"].head(n_top).tolist()

    # Fallback: most cells.
    counts = (
        adata.obs[adata.obs["perturbation"].isin(non_ctrl)]["perturbation"]
        .value_counts()
    )
    return counts.head(n_top).index.tolist()


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

def _plot_volcano(
    deg_df: pd.DataFrame,
    perturbation: str,
    plots_dir: Path,
    logfc_threshold: float,
    pval_threshold: float,
) -> None:
    # Draw a volcano plot for one perturbation vs control.
    if deg_df.empty:
        return

    def _colour(row):
        if row["padj"] < pval_threshold and row["log2fc"] >= logfc_threshold:
            return "#e63946"  # up-regulated
        if row["padj"] < pval_threshold and row["log2fc"] <= -logfc_threshold:
            return "#457b9d"  # down-regulated
        return "#adb5bd"      # not significant

    colors = deg_df.apply(_colour, axis=1)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(
        deg_df["log2fc"],
        deg_df["neg_log10_padj"],
        c=colors,
        s=10,
        alpha=0.7,
        rasterized=True,
    )

    # Label top 10 significant genes by |log2fc|.
    top = deg_df[deg_df["significant"]].head(10)
    for _, row in top.iterrows():
        ax.annotate(
            row["gene"],
            xy=(row["log2fc"], row["neg_log10_padj"]),
            fontsize=_ANNOT_FS,
            xytext=(4, 0),
            textcoords="offset points",
            color="#333333",
        )

    # Show mean and median of log2FC as vertical reference lines.
    fc_mean   = float(deg_df["log2fc"].mean())
    fc_median = float(deg_df["log2fc"].median())
    ax.axvline(fc_mean,   color="#555555", lw=1.2, ls="-",  label=f"mean FC = {fc_mean:.2f}")
    ax.axvline(fc_median, color="#e76f51", lw=1.2, ls="--", label=f"median FC = {fc_median:.2f}")

    sig = deg_df["padj"] < pval_threshold
    n_up = int((sig & (deg_df["log2fc"] >= logfc_threshold)).sum())
    n_dn = int((sig & (deg_df["log2fc"] <= -logfc_threshold)).sum())
    handles = [
        mpatches.Patch(color="#e63946", label=f"Up ({n_up})"),
        mpatches.Patch(color="#457b9d", label=f"Down ({n_dn})"),
        mpatches.Patch(color="#adb5bd", label="NS"),
    ]
    legend1 = ax.legend(handles=handles, loc="upper left", fontsize=_LEGEND_FS)
    ax.add_artist(legend1)  # keep patch legend
    ax.legend(loc="upper right", fontsize=_LEGEND_FS)  # mean/median lines legend
    ax.set_xlabel("log₂ fold-change  (perturbed / control)", fontsize=_LABEL_FS)
    ax.set_ylabel("−log₁₀(padj)", fontsize=_LABEL_FS)
    ax.tick_params(labelsize=_TICK_FS)
    ax.set_title(f"Volcano: {perturbation}  vs  control", fontsize=_TITLE_FS)
    fig.tight_layout()

    safe = perturbation.replace("/", "_").replace(" ", "_").replace("+", "_")
    fig.savefig(
        plots_dir / f"deg_volcano_{safe}.png", dpi=150, bbox_inches="tight"
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# Top-DEG summary heatmap
# ---------------------------------------------------------------------------

def _plot_deg_heatmap(
    adata,
    deg_results: Dict[str, pd.DataFrame],
    plots_dir: Path,
    n_top_genes: int,
) -> None:
    # Clustered heatmap of log2FC values: rows = perturbations,
    # columns = union of the top-n DEGs from each perturbation.
    if not deg_results:
        return

    gene_set: List[str] = []
    for df in deg_results.values():
        if df.empty:
            continue
        for g in df.head(n_top_genes)["gene"].tolist():
            if g not in gene_set:
                gene_set.append(g)

    if not gene_set:
        return

    records = {}
    for pert, df in deg_results.items():
        if df.empty:
            records[pert] = pd.Series(np.nan, index=gene_set)
        else:
            records[pert] = df.set_index("gene")["log2fc"].reindex(gene_set)

    mat = pd.DataFrame(records).T  # (perturbation × gene)
    if mat.empty:
        return

    n_r, n_c = mat.shape
    try:
        g = sns.clustermap(
            mat.fillna(0),
            cmap="RdBu_r",
            center=0,
            vmin=-2,
            vmax=2,
            xticklabels=False,
            yticklabels=True,
            figsize=(max(12, n_c // 3), max(6, n_r)),
            cbar_kws={"label": "log₂FC"},
        )
        g.ax_heatmap.set_xlabel("Genes (top DEGs)", fontsize=_LABEL_FS)
        g.ax_heatmap.set_ylabel("Perturbation", fontsize=_LABEL_FS)
        g.ax_heatmap.tick_params(axis="x", labelsize=_TICK_FS)
        g.ax_heatmap.tick_params(axis="y", labelsize=_TICK_FS)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        g.fig.suptitle("DEG log₂FC: top perturbations vs control", y=1.01, fontsize=_TITLE_FS)
        g.fig.savefig(
            plots_dir / "deg_top_perturbations_heatmap.png",
            dpi=120,
            bbox_inches="tight",
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_deg_analysis(
    adata,
    output_dir: str,
    perturbations: Optional[List[str]] = None,
    effect_df: Optional[pd.DataFrame] = None,
    n_top_perturbations: int = 5,
    logfc_threshold: float = 0.5,
    pval_threshold: float = 0.05,
    n_top_deg_heatmap: int = 20,
) -> Dict[str, pd.DataFrame]:
    """Run Welch t-test DEG analysis for the top perturbations vs control.

    Args:
        adata                -- AnnData after QC + normalisation.
        output_dir           -- root output directory.
        perturbations        -- explicit list of perturbations to analyse.
                                When None, identify_top_perturbations() is
                                called with n_top=n_top_perturbations.
        effect_df            -- DataFrame from compute_effect_decomposition();
                                used for smarter ranking (optional).
        n_top_perturbations  -- how many perturbations to auto-select.
        logfc_threshold      -- |log2fc| cutoff for "significant" label.
        pval_threshold       -- BH-adjusted p-value cutoff.
        n_top_deg_heatmap    -- top DEGs per perturbation included in the
                                summary heatmap.

    Returns:
        dict mapping perturbation name -> DEG DataFrame (or empty DataFrame
        if too few cells were available).
    """
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    if perturbations is None:
        perturbations = identify_top_perturbations(
            adata, effect_df=effect_df, n_top=n_top_perturbations
        )

    deg_results: Dict[str, pd.DataFrame] = {}
    summary_rows = []

    for pert in perturbations:
        df = _compute_deg(adata, pert)
        deg_results[pert] = df
        if df.empty:
            continue

        safe = pert.replace("/", "_").replace(" ", "_").replace("+", "_")
        df = annotate_deg(df, adata.var_names.tolist())
        df.to_csv(tables / f"deg_{safe}.csv", index=False)
        sig_genes = df.loc[df["significant"], "gene"].tolist()
        # If fewer than 5 significant DEGs, fall back to top-200 by p-value
        if len(sig_genes) >= 5:
            query_genes = sig_genes
            query_type = "sig_deg"
        else:
            query_genes = (
                df.sort_values("padj").head(200)["gene"].tolist()
            )
            query_type = "top200_pval"
        enrich_df = run_enrichment(
            query_genes,
            adata.var_names.tolist(),
            top_n=50,
        )
        if not enrich_df.empty:
            enrich_df["query_type"] = query_type
            enrich_df.to_csv(tables / f"deg_enrichment_{safe}.csv", index=False)
        _plot_volcano(df, pert, plots, logfc_threshold, pval_threshold)

        sig = df["padj"] < pval_threshold
        n_up = int((sig & (df["log2fc"] >= logfc_threshold)).sum())
        n_dn = int((sig & (df["log2fc"] <= -logfc_threshold)).sum())
        top_gene = df[df["significant"]].iloc[0]["gene"] if df["significant"].any() else ""
        summary_rows.append(
            {
                "perturbation": pert,
                "n_cells": int((adata.obs["perturbation"].values == pert).sum()),
                "n_de_up": n_up,
                "n_de_down": n_dn,
                "n_de_total": n_up + n_dn,
                "top_gene": top_gene,
            }
        )

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(tables / "deg_summary.csv", index=False)

    _plot_deg_heatmap(adata, deg_results, plots, n_top_deg_heatmap)

    return deg_results
