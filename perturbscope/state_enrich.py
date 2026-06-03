# =============================================================================
# perturbscope/state_enrich.py
#
# Cell state enrichment / depletion analysis for Perturb-seq data.
#
# For each (perturbation, cell_state) pair, tests whether the proportion of
# cells from that perturbation in the state is significantly enriched or
# depleted relative to the control population, using a chi-square test on a
# 2×2 contingency table (with Yates correction).
#
# Sign convention (matching published figures):
#   Positive value → perturbation is ENRICHED in that cell state  (red)
#   Negative value → perturbation is DEPLETED in that cell state  (blue)
#
# Outputs written to <output_dir>/plots/ and <output_dir>/csv/:
#
#   csv/
#     state_enrich_results.csv   -- long-form table (pert × state, pval, qval, …)
#     state_enrich_matrix.csv    -- wide pert × state matrix of signed -log10(q)
#   plots/
#     state_enrich_heatmap.png   -- clustered heatmap (≤50 perts shown)
#
# Usage in pipeline (wired by pipeline.py):
#   from .state_enrich import run_state_enrichment
#   run_state_enrichment(adata, output_dir)
# =============================================================================

from __future__ import annotations

from math import erfc, log10, sqrt
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from .utils import ensure_dir

# ---------------------------------------------------------------------------
# Font constants (mirror other modules)
# ---------------------------------------------------------------------------
_TITLE_FS = 15
_LABEL_FS = 13
_TICK_FS = 10
_CBAR_FS = 11


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _chi2_2x2_pval(a: int, b: int, c: int, d: int) -> float:
    """Yates-corrected chi-square p-value for a 2×2 contingency table.

    Table layout:
        pert in state | pert not in state   →  a  b
        ctrl in state | ctrl not in state   →  c  d
    """
    N = a + b + c + d
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    if N == 0 or row1 == 0 or row2 == 0 or col1 == 0 or col2 == 0:
        return 1.0
    # Yates-corrected chi2
    chi2 = N * max(0.0, abs(a * d - b * c) - N / 2.0) ** 2 / (row1 * row2 * col1 * col2)
    # chi2(df=1) survival: erfc(sqrt(chi2/2))
    return min(float(erfc(sqrt(chi2 / 2.0))), 1.0)


def _bh(pvals: list) -> list:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    result = [1.0] * n
    cummin = 1.0
    for rank, idx in enumerate(reversed(order)):
        raw = pvals[idx] * n / (n - rank)
        cummin = min(cummin, raw)
        result[idx] = cummin
    return result


def _signed_log10q(q: float, sign: float) -> float:
    return round(-log10(max(q, 1e-300)) * sign, 4)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _plot_state_enrich_heatmap(mat: pd.DataFrame, plots_dir: Path) -> None:
    """Clustered heatmap of cell state enrichment / depletion.

    Rows = perturbations, columns = cell states.
    Color = signed -log10(q-value); red = enriched, blue = depleted.
    """
    if mat.empty or mat.shape[0] < 2:
        return

    n_r, n_c = mat.shape
    figw = max(5 + n_c * 0.55, 6)
    figh = max(4 + n_r * 0.18, 5)
    tick_fs = max(5, min(_TICK_FS, 180 // max(n_r, 1)))

    try:
        vmax = min(max(abs(mat.values).max(), 1.0), 6.0)
        g = sns.clustermap(
            mat,
            col_cluster=(n_c > 2),
            row_cluster=(n_r > 2),
            cmap="RdBu_r",
            center=0,
            vmin=-vmax,
            vmax=vmax,
            figsize=(figw, figh),
            xticklabels=True,
            yticklabels=(n_r <= 80),
            cbar_kws={"label": "signed −log₁₀(q)", "shrink": 0.5},
            linewidths=0.5 if n_r < 30 else 0,
        )
        g.ax_heatmap.set_title(
            "Cell state enrichment or depletion", fontsize=_TITLE_FS, pad=12
        )
        g.ax_heatmap.set_xlabel("Cell state", fontsize=_LABEL_FS)
        g.ax_heatmap.set_ylabel("Perturbed genes, pairs and triplets", fontsize=_LABEL_FS)
        g.ax_heatmap.tick_params(axis="x", labelsize=_TICK_FS, rotation=0)
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_fs)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        g.fig.savefig(
            plots_dir / "state_enrich_heatmap.png", dpi=120, bbox_inches="tight"
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_state_enrichment(
    adata,
    output_dir: str,
    min_cells: int = 5,
    fdr_threshold: float = 0.05,
    max_perts_plot: int = 60,
) -> pd.DataFrame:
    """Cell state enrichment / depletion analysis for each perturbation.

    Args:
        adata           -- AnnData after QC + cell_state assignment.
        output_dir      -- Root output directory.
        min_cells       -- Minimum cells in perturbation group to test.
        fdr_threshold   -- BH-FDR threshold for significance flag.
        max_perts_plot  -- Maximum perturbations shown in the heatmap
                           (top by max |signal| are kept).

    Returns:
        DataFrame: perturbation × cell_state of signed -log10(q-value).
    """
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    if "cell_state" not in adata.obs.columns:
        return pd.DataFrame()

    states = sorted(adata.obs["cell_state"].astype(str).unique())
    ctrl_mask = adata.obs["is_control"].values.astype(bool)
    n_ctrl = int(ctrl_mask.sum())
    if n_ctrl < min_cells:
        return pd.DataFrame()

    ctrl_obs = adata.obs.loc[ctrl_mask, "cell_state"].astype(str)
    ctrl_counts = {s: int((ctrl_obs == s).sum()) for s in states}

    perts = [
        p for p in sorted(adata.obs["perturbation"].astype(str).unique())
        if str(p).lower() not in ("control", "nan")
    ]

    records = []
    for pert in perts:
        p_mask = adata.obs["perturbation"].values == pert
        n_pert = int(p_mask.sum())
        if n_pert < min_cells:
            continue
        pert_obs = adata.obs.loc[p_mask, "cell_state"].astype(str)
        pert_counts = {s: int((pert_obs == s).sum()) for s in states}

        for s in states:
            a = pert_counts.get(s, 0)
            b = n_pert - a
            c = ctrl_counts.get(s, 0)
            d = n_ctrl - c
            pval = _chi2_2x2_pval(a, b, c, d)
            obs_frac = a / n_pert if n_pert > 0 else 0.0
            exp_frac = c / n_ctrl if n_ctrl > 0 else 0.0
            sign = 1.0 if obs_frac >= exp_frac else -1.0
            records.append({
                "perturbation": pert,
                "cell_state": s,
                "n_pert_in_state": a,
                "n_pert": n_pert,
                "n_ctrl_in_state": c,
                "n_ctrl": n_ctrl,
                "obs_frac": round(obs_frac, 4),
                "exp_frac": round(exp_frac, 4),
                "pval": pval,
                "sign": sign,
            })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    qvals = _bh(df["pval"].tolist())
    df["qval"] = qvals
    df["significant"] = df["qval"] < fdr_threshold
    df["signed_log10q"] = [
        _signed_log10q(q, s) for q, s in zip(df["qval"], df["sign"])
    ]
    df.to_csv(tables / "state_enrich_results.csv", index=False)

    # Pivot to (pert × state) matrix
    mat = (
        df.pivot(index="perturbation", columns="cell_state", values="signed_log10q")
        .fillna(0)
    )
    mat.columns.name = None
    mat.to_csv(tables / "state_enrich_matrix.csv")

    # Limit heatmap rows to the top perturbations by max |signal|
    mat_plot = mat.reindex(
        mat.abs().max(axis=1).sort_values(ascending=False).index
    ).head(max_perts_plot)

    _plot_state_enrich_heatmap(mat_plot, plots)

    return mat
