# =============================================================================
# perturbflow/analyzer/cscore.py
#
# Connectivity Score (C-score): a metric quantifying the structural rewiring
# of gene co-expression networks induced by each perturbation.
#
# The C-score decomposes into three complementary components:
#
#   C_gain(p)   = |E_p \ E_c| / (|E_c| + 1)
#                 Fraction of edges gained (new edges in perturbation not in ctrl)
#
#   C_loss(p)   = |E_c \ E_p| / (|E_c| + 1)
#                 Fraction of edges lost (ctrl edges absent in perturbation)
#
#   C_shift(p)  = (1 / |E_c ∩ E_p|) * Σ |r_ij^p − r_ij^c|   for shared edges
#                 Mean absolute weight shift on preserved edges
#
#   C_total(p)  = C_gain + C_loss + C_shift
#
# Additionally, a gene-level "hub rewiring" score is computed:
#   hub_change(g) = degree(g, G_pert) − degree(g, G_ctrl)
#
# Outputs (written to <output_dir>/csv/ and <output_dir>/plots/):
#
#   csv/
#     cscore_summary.csv         -- per-perturbation C-score table
#     cscore_edge_<safe>.csv     -- per-edge status (gained/lost/shifted/shared)
#     cscore_gene_hub.csv        -- per-gene-per-perturbation hub degree change
#
#   plots/
#     cscore_ranked_bar.png      -- perturbations ranked by C_total, stacked
#     cscore_decomposition.png   -- scatter: C_gain vs C_loss, size=C_shift
#     cscore_vs_deg.png          -- C_total vs n_DE_genes scatter
#     cscore_gene_hub_heatmap.png-- gene × pert hub-degree-change heatmap
#     cscore_module_rewiring.png -- per gene-module fraction gained/lost
#
# Usage:
#   from .cscore import run_cscore
#   run_cscore(adata, output_dir, perturbations=top_perts, corr_threshold=0.4)
# =============================================================================

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from .utils import ensure_dir

# ---------------------------------------------------------------------------
# Global style: Arial font for all figures
# ---------------------------------------------------------------------------
import matplotlib as _mpl
_mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})

# ---------------------------------------------------------------------------
# Global font size constants (mirrors eda.py / genenet.py)
# ---------------------------------------------------------------------------
_TITLE_FS  = 18
_LABEL_FS  = 15
_TICK_FS   = 13
_LEGEND_FS = 13
_CBAR_FS   = 13

# C-score summary panels are embedded in the interactive report, so keep the
# canvas compact while making axis text large enough to read after scaling.
_CS_AXIS_FS = 17
_CS_TICK_FS = 15
_CS_LEGEND_FS = 14

# Publication-quality color palette (colorblind-friendly, muted)
_COL_GAIN  = "#C4616B"   # muted rose   — new co-expression edges
_COL_LOSS  = "#4392B4"   # steel blue   — lost edges
_COL_SHIFT = "#E8A838"   # warm amber   — weight-shift on preserved edges


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _to_dense(x):
    return x.A if hasattr(x, "A") else np.asarray(x)


def _pearson_corrmat(x: np.ndarray) -> np.ndarray:
    corr = np.corrcoef(x.T)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def _edge_sets(corr: np.ndarray, gene_names: List[str], threshold: float):
    """Return {(i,j): r_ij} for |r| >= threshold (i < j always)."""
    edges = {}
    n = len(gene_names)
    for i in range(n):
        for j in range(i + 1, n):
            r = float(corr[i, j])
            if abs(r) >= threshold:
                edges[(i, j)] = r
    return edges


def _compute_cscore(
    ctrl_edges: dict, pert_edges: dict
) -> dict:
    """Compute C-score components from two edge dicts {(i,j): r}."""
    ec = set(ctrl_edges)
    ep = set(pert_edges)

    n_ctrl = len(ec)
    n_pert = len(ep)
    gained = ep - ec
    lost   = ec - ep
    shared = ec & ep

    c_gain  = len(gained) / (n_ctrl + 1)
    c_loss  = len(lost)   / (n_ctrl + 1)

    if shared:
        c_shift = float(np.mean([abs(pert_edges[e] - ctrl_edges[e]) for e in shared]))
    else:
        c_shift = 0.0

    c_total = c_gain + c_loss + c_shift

    return {
        "c_total":    round(c_total, 5),
        "c_gain":     round(c_gain, 5),
        "c_loss":     round(c_loss, 5),
        "c_shift":    round(c_shift, 5),
        "n_edges_ctrl":  n_ctrl,
        "n_edges_pert":  n_pert,
        "n_gained":   len(gained),
        "n_lost":     len(lost),
        "n_shared":   len(shared),
    }


def _hub_degrees(edges: dict, n_genes: int) -> np.ndarray:
    """Node degree array from an edge dict."""
    deg = np.zeros(n_genes, dtype=int)
    for i, j in edges:
        deg[i] += 1
        deg[j] += 1
    return deg


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _plot_ranked_bar(summary_df: pd.DataFrame, plots_dir: Path) -> None:
    """Stacked bar: perturbations ranked by C_total."""
    df = summary_df.sort_values("c_total", ascending=True)
    n = len(df)
    if n == 0:
        return

    fig_h = max(5.0, n * 0.34 + 1.8)
    fig, ax = plt.subplots(figsize=(7.8, fig_h))

    # Stacked horizontal bar
    bars = [
        ("c_gain",  _COL_GAIN,  "Gained"),
        ("c_loss",  _COL_LOSS,  "Lost"),
        ("c_shift", _COL_SHIFT, "Weight-shift"),
    ]
    left = np.zeros(n)
    ys = np.arange(n)
    for col, color, label in bars:
        vals = df[col].values
        ax.barh(ys, vals, left=left, color=color, label=label, height=0.42, alpha=0.88)
        left += vals

    ax.set_yticks(ys)
    ax.set_yticklabels(df["perturbation"].tolist(), fontsize=max(12, _CS_TICK_FS - max(0, n - 20)))
    ax.set_xlabel("C-score (gain + loss + shift)", fontsize=_CS_AXIS_FS)
    ax.set_title("Connectivity Score — perturbations ranked by C_total", fontsize=_TITLE_FS, pad=10)
    ax.legend(loc="lower right", fontsize=_CS_LEGEND_FS, frameon=True)
    ax.tick_params(axis="x", labelsize=_CS_TICK_FS)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(plots_dir / "cscore_ranked_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_decomposition(summary_df: pd.DataFrame, plots_dir: Path) -> None:
    """Scatter: C_gain vs C_loss, dot size ~ C_shift, labelled."""
    df = summary_df.copy()
    if df.empty:
        return

    shift_norm = (df["c_shift"] - df["c_shift"].min())
    shift_range = df["c_shift"].max() - df["c_shift"].min()
    if shift_range > 0:
        shift_norm /= shift_range
    sizes = (shift_norm * 300 + 40).values

    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    sc = ax.scatter(
        df["c_gain"], df["c_loss"],
        s=sizes, c=df["c_total"], cmap="YlOrRd",
        alpha=0.85, edgecolors="#333", linewidths=0.5, zorder=3,
    )
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7)
    cbar.set_label("C_total", fontsize=_CS_AXIS_FS)
    cbar.ax.tick_params(labelsize=_CS_TICK_FS)

    # Label only the top-10 perturbations by C_total to keep the figure
    # readable when there are many points.
    top10 = df.nlargest(10, "c_total")
    for _, row in top10.iterrows():
        ax.text(row["c_gain"] + 0.001, row["c_loss"] + 0.001, row["perturbation"],
                fontsize=max(11, _CS_TICK_FS - 2), ha="left", va="bottom", color="#333")

    # Diagonal guide (gain == loss → symmetric rewiring)
    lim = max(df[["c_gain", "c_loss"]].values.max() * 1.1, 0.01)
    ax.plot([0, lim], [0, lim], "--", color="#aaa", linewidth=1, label="C_gain = C_loss")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("C_gain (new edges)", fontsize=_CS_AXIS_FS)
    ax.set_ylabel("C_loss (lost edges)", fontsize=_CS_AXIS_FS)
    ax.set_title("C-score Decomposition\n(dot size = C_shift, colour = C_total)",
                 fontsize=_TITLE_FS, pad=10)
    ax.legend(fontsize=_CS_LEGEND_FS)
    ax.tick_params(labelsize=_CS_TICK_FS)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(plots_dir / "cscore_decomposition.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_cscore_vs_deg(
    summary_df: pd.DataFrame,
    deg_csv: Optional[Path],
    plots_dir: Path,
) -> None:
    """Scatter: C_total vs n_DE_genes, labelled."""
    df = summary_df.copy()
    if deg_csv and deg_csv.exists():
        deg_df = pd.read_csv(deg_csv)
        df = df.merge(deg_df[["perturbation", "n_de_total"]], on="perturbation", how="left")
    else:
        df["n_de_total"] = np.nan

    if df["n_de_total"].isna().all():
        return

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(df["n_de_total"], df["c_total"],
               s=60, color="#4c78a8", alpha=0.8, edgecolors="#333", linewidths=0.5)
    # Label only the top-10 perturbations by C_total. Tiny offsets keep the
    # text close to its dot (no adjustText available in this env, so the
    # smaller label set is the simple anti-overlap strategy).
    valid = df.dropna(subset=["n_de_total"])
    top10 = valid.nlargest(10, "c_total")
    for _, row in top10.iterrows():
        ax.text(row["n_de_total"] + 0.3, row["c_total"] + 0.001, row["perturbation"],
                fontsize=max(7, _TICK_FS - 2), color="#333")

    ax.set_xlabel("Number of DE genes", fontsize=_LABEL_FS)
    ax.set_ylabel("C_total", fontsize=_LABEL_FS)
    ax.set_title("Network Rewiring vs Differential Expression",
                 fontsize=_TITLE_FS - 4, pad=8)
    ax.tick_params(labelsize=_TICK_FS)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(plots_dir / "cscore_vs_deg.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_gene_hub_heatmap(
    hub_df: pd.DataFrame,
    plots_dir: Path,
) -> None:
    """Gene × perturbation heatmap of hub-degree change."""
    if hub_df.empty:
        return

    # Pivot: rows=gene, cols=perturbation, values=hub_change
    pivot = hub_df.pivot_table(
        index="gene", columns="perturbation", values="hub_change", fill_value=0
    )
    if pivot.empty:
        return

    # Sort rows by total absolute change
    pivot = pivot.loc[pivot.abs().sum(axis=1).sort_values(ascending=False).index]
    n_genes, n_perts = pivot.shape

    # Show top 25 genes for readability
    pivot = pivot.head(25)
    n_genes = len(pivot)

    tick_fs_r = max(18, min(_TICK_FS + 8, 400 // max(n_genes, 1)))
    tick_fs_c = max(16, min(_TICK_FS + 6, 280 // max(n_perts, 1)))
    vmax = float(pivot.abs().max().max()) or 1.0

    try:
        # Taller figure so the rotated column labels + horizontal colorbar
        # both fit below the heatmap without colliding with each other.
        g = sns.clustermap(
            pivot,
            cmap="RdBu_r", center=0,
            vmin=-vmax, vmax=vmax,
            xticklabels=True, yticklabels=False,
            figsize=(max(7, n_perts * 1.4 + 2), max(6.5, n_genes * 0.18 + 2.5)),
            cbar_pos=None,           # colorbar added manually below
            dendrogram_ratio=(0.10, 0.10),
        )
        # Vertical column labels — they occupy a fixed bottom strip and never
        # bleed sideways into the colorbar.
        g.ax_heatmap.tick_params(axis="x", labelsize=tick_fs_c, rotation=90)
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_fs_r)
        for lbl in g.ax_heatmap.get_xticklabels():
            lbl.set_rotation(90)
            lbl.set_ha("center")
            lbl.set_va("top")
        g.ax_col_dendrogram.set_title(
            "Gene Hub Rewiring (Δ degree per perturbation, top 25 genes)",
            fontsize=_TITLE_FS + 6, pad=10,
        )

        # Reserve a generous bottom strip: top of strip = column labels,
        # bottom of strip = horizontal colorbar. Heatmap is shifted up so
        # labels live at y∈[0.10, 0.22] and colorbar at y=[0.04, 0.06].
        g.fig.subplots_adjust(bottom=0.26)
        ax_cbar = g.fig.add_axes([0.20, 0.04, 0.55, 0.022])
        norm = plt.Normalize(vmin=-vmax, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
        sm.set_array([])
        cb = g.fig.colorbar(sm, cax=ax_cbar, orientation="horizontal")
        cb.set_label("Δ degree (pert − ctrl)", fontsize=_CBAR_FS + 6)
        cb.ax.tick_params(labelsize=_TICK_FS + 4)
        g.fig.savefig(plots_dir / "cscore_gene_hub_heatmap.png", dpi=120,
                      bbox_inches="tight")
        plt.close(g.fig)
    except Exception:
        plt.close("all")


def _plot_module_rewiring(
    summary_df: pd.DataFrame,
    module_df: pd.DataFrame,
    plots_dir: Path,
) -> None:
    """Edge gain / loss horizontal bars matching the ranked-bar layout.

    Mirrors `_plot_ranked_bar` so the two figures sit side-by-side at the
    same orientation and tick / label / title sizes. Perturbations on the
    y-axis (sorted by C_gain ascending so the largest is at the top of the
    plot), C_gain and C_loss as grouped horizontal bars.
    """
    if module_df.empty or summary_df.empty:
        return
    df = summary_df[["perturbation", "c_gain", "c_loss"]].sort_values(
        "c_gain", ascending=True
    )
    n = len(df)
    if n == 0:
        return

    fig_h = max(5.0, n * 0.34 + 1.8)
    fig, ax = plt.subplots(figsize=(7.8, fig_h))

    ys = np.arange(n)
    bar_h = 0.4
    ax.barh(ys - bar_h / 2, df["c_gain"].values, height=bar_h,
            color=_COL_GAIN, label="C_gain", alpha=0.88)
    ax.barh(ys + bar_h / 2, df["c_loss"].values, height=bar_h,
            color=_COL_LOSS, label="C_loss", alpha=0.88)

    ax.set_yticks(ys)
    ax.set_yticklabels(df["perturbation"].tolist(),
                       fontsize=max(12, _CS_TICK_FS - max(0, n - 20)))
    ax.set_xlabel("C-score component", fontsize=_CS_AXIS_FS)
    ax.set_title("Edge Gain vs Loss per Perturbation", fontsize=_TITLE_FS, pad=10)
    ax.legend(loc="lower right", fontsize=_CS_LEGEND_FS, frameon=True)
    ax.tick_params(axis="x", labelsize=_CS_TICK_FS)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(plots_dir / "cscore_module_rewiring.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_cscore(
    adata,
    output_dir: str,
    perturbations: Optional[List[str]] = None,
    n_top_genes: int = 1000,
    corr_threshold: float = 0.4,
) -> pd.DataFrame:
    """Compute the Connectivity Score for each perturbation vs control.

    Args:
        adata           -- AnnData after QC + normalisation.
        output_dir      -- root output directory.
        perturbations   -- list of perturbation labels to analyse.
                           Defaults to all non-control perturbations.
        n_top_genes     -- number of top HVGs to build networks from.
        corr_threshold  -- |Pearson r| threshold for an edge to exist.

    Returns:
        DataFrame with columns [perturbation, c_total, c_gain, c_loss, c_shift, ...]
        (also written to csv/cscore_summary.csv).
    """
    out    = ensure_dir(output_dir)
    plots  = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    # ── Identify control cells ────────────────────────────────────────────
    perts_arr = adata.obs["perturbation"].astype(str).values
    _ctrl_kw = {"control", "ctrl", "nontargeting", "non-targeting", "nt",
                "scramble", "safe-targeting", "safe_targeting"}
    ctrl_mask = np.array([p.lower() in _ctrl_kw for p in perts_arr])
    if ctrl_mask.sum() < 3:
        # Fall back to the most-abundant label
        groups, counts = np.unique(perts_arr, return_counts=True)
        ctrl_label = groups[np.argmax(counts)]
        ctrl_mask = perts_arr == ctrl_label

    # ── Select top HVGs ───────────────────────────────────────────────────
    if "highly_variable" in adata.var.columns:
        hv_idx = np.where(adata.var["highly_variable"].values)[0][:n_top_genes]
    else:
        xc_full = _to_dense(adata.X[ctrl_mask, :])
        hv_idx = np.argsort(xc_full.var(axis=0))[::-1][:n_top_genes]
    gene_names: List[str] = adata.var_names[hv_idx].tolist()
    n_g = len(gene_names)
    if n_g < 4:
        return pd.DataFrame()

    # ── Control network ───────────────────────────────────────────────────
    ctrl_expr = _to_dense(adata.X[ctrl_mask, :][:, hv_idx])
    ctrl_corr = _pearson_corrmat(ctrl_expr)
    ctrl_edges = _edge_sets(ctrl_corr, gene_names, corr_threshold)
    ctrl_deg   = _hub_degrees(ctrl_edges, n_g)

    # ── Perturbation list ─────────────────────────────────────────────────
    all_perts = [p for p in np.unique(perts_arr) if p.lower() not in _ctrl_kw]
    if perturbations is None:
        perturbations = all_perts
    else:
        perturbations = [p for p in perturbations if p in set(all_perts)]

    # ── Load gene-module assignments (from genenet step if available) ──────
    module_csv = tables / "genenet_gene_clusters.csv"
    module_df  = pd.read_csv(module_csv) if module_csv.exists() else pd.DataFrame()

    # ── Per-perturbation loop ─────────────────────────────────────────────
    summary_rows  = []
    edge_dfs      = []
    hub_rows      = []

    for pert in perturbations:
        pert_mask = perts_arr == pert
        if int(pert_mask.sum()) < 3:
            continue

        pert_expr  = _to_dense(adata.X[pert_mask, :][:, hv_idx])
        pert_corr  = _pearson_corrmat(pert_expr)
        pert_edges = _edge_sets(pert_corr, gene_names, corr_threshold)

        # C-score
        cs = _compute_cscore(ctrl_edges, pert_edges)
        cs["perturbation"] = pert
        summary_rows.append(cs)

        # Edge-level table
        safe = pert.replace("/", "_").replace(" ", "_").replace("+", "_")
        union = set(ctrl_edges) | set(pert_edges)
        erows = []
        for i, j in union:
            rc = ctrl_edges.get((i, j), 0.0)
            rp = pert_edges.get((i, j), 0.0)
            if (i, j) in ctrl_edges and (i, j) in pert_edges:
                status = "shared"
            elif (i, j) in ctrl_edges:
                status = "lost"
            else:
                status = "gained"
            erows.append({
                "gene_A": gene_names[i], "gene_B": gene_names[j],
                "status": status, "ctrl_r": round(rc, 4),
                "pert_r": round(rp, 4), "delta_r": round(rp - rc, 4),
            })
        edf = pd.DataFrame(erows)
        edf.to_csv(tables / f"cscore_edges_{safe}.csv", index=False)
        edge_dfs.append(edf)

        # Hub degrees
        pert_deg = _hub_degrees(pert_edges, n_g)
        for gi, gname in enumerate(gene_names):
            hub_rows.append({
                "perturbation": pert,
                "gene": gname,
                "degree_ctrl": int(ctrl_deg[gi]),
                "degree_pert": int(pert_deg[gi]),
                "hub_change":  int(pert_deg[gi]) - int(ctrl_deg[gi]),
            })

    if not summary_rows:
        return pd.DataFrame()

    summary_df = pd.DataFrame(summary_rows)
    cols_order = ["perturbation", "c_total", "c_gain", "c_loss", "c_shift",
                  "n_edges_ctrl", "n_edges_pert", "n_gained", "n_lost", "n_shared"]
    summary_df = summary_df[[c for c in cols_order if c in summary_df.columns]]
    summary_df.to_csv(tables / "cscore_summary.csv", index=False)

    hub_df = pd.DataFrame(hub_rows)
    hub_df.to_csv(tables / "cscore_gene_hub.csv", index=False)

    # ── Figures ───────────────────────────────────────────────────────────
    _plot_ranked_bar(summary_df, plots)
    _plot_decomposition(summary_df, plots)
    _plot_cscore_vs_deg(summary_df, tables / "deg_summary.csv", plots)
    _plot_gene_hub_heatmap(hub_df, plots)
    _plot_module_rewiring(summary_df, module_df, plots)

    return summary_df
