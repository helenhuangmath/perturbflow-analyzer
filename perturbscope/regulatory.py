# =============================================================================
# perturbscope/regulatory.py
#
# Perturbation effect correlation and TF-TF regulatory network analysis.
#
# Three analyses, all derived from the DEG log2FC / padj tables:
#
# 1. Perturbation effect correlation (before / after cell-state correction)
#    Pearson correlation between per-perturbation log2FC effect vectors.
#    "Before" = raw log2FC across all cells.
#    "After"  = cell-state-weighted log2FC (removes composition confound).
#    Guides targeting functionally related TFs cluster together.
#
# 2. TF-TF regulatory heatmap
#    For each perturbed gene (guide, rows) and each other perturbed gene
#    that appears in the DEG table (columns), extract the signed -log10(padj):
#      Positive (red)  = KO increases target → TF normally inhibits target.
#      Negative (blue) = KO decreases target → TF normally activates target.
#
# 3. TF regulatory network (directed graph)
#    Nodes = perturbed genes; edges = significant regulatory relationships.
#    Edge colour:
#      Blue = activating  (KO → target ↓, so TF→target is activating).
#      Red  = inhibiting  (KO → target ↑, so TF→target is inhibiting).
#    Spring-layout positions are stored in the JSON so the interactive
#    report can render a Plotly network without re-running networkx.
#
# Outputs written to <output_dir>/plots/ and <output_dir>/csv/:
#
#   csv/
#     pert_effect_corr_before.csv      -- pert×pert Pearson corr (raw)
#     pert_effect_corr_after.csv       -- pert×pert Pearson corr (state-adj.)
#     tf_regulatory_matrix.csv         -- signed -log10(q) guide×gene matrix
#     tf_regulatory_network.json       -- nodes, positions, edges (for report)
#   plots/
#     pert_effect_corr_before.png
#     pert_effect_corr_after.png
#     tf_regulatory_heatmap.png
#     tf_regulatory_network.png
#
# Usage in pipeline (wired by pipeline.py):
#   from .regulatory import run_regulatory_analysis
#   run_regulatory_analysis(adata, output_dir)
# =============================================================================

from __future__ import annotations

import json
from math import log10
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from .utils import ensure_dir

# ---------------------------------------------------------------------------
# Font constants
# ---------------------------------------------------------------------------
_TITLE_FS = 20
_LABEL_FS = 17
_TICK_FS = 13
_CBAR_FS = 15
_LEGEND_FS = 15
_PLOT_DPI = 90
_PNG_KWARGS = {"optimize": True, "compress_level": 9}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _to_dense(x):
    return x.A if hasattr(x, "A") else np.asarray(x)


def _save_png(fig, out_path: Path, tight: bool = True) -> None:
    """Write compact PNGs while keeping text large enough for reports."""
    bbox = "tight" if tight else None
    try:
        fig.savefig(
            out_path,
            dpi=_PLOT_DPI,
            bbox_inches=bbox,
            pil_kwargs=_PNG_KWARGS,
        )
    except TypeError:
        fig.savefig(out_path, dpi=_PLOT_DPI, bbox_inches=bbox)


def _pearson_rows(M: np.ndarray) -> np.ndarray:
    """Pearson correlation between rows of M (NaN-safe)."""
    M_c = M - M.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(M_c, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1.0, norms)
    M_n = M_c / norms
    return np.clip(M_n @ M_n.T, -1.0, 1.0)


def _bh(pvals: list) -> list:
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


# ---------------------------------------------------------------------------
# Load DEG matrices
# ---------------------------------------------------------------------------

def _load_deg_matrix(csv_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load all per-perturbation DEG CSVs and build:
      - log2fc_df : (pert × gene) log2FC matrix
      - padj_df   : (pert × gene) adjusted-p matrix
    """
    lfc_dict: Dict[str, pd.Series] = {}
    padj_dict: Dict[str, pd.Series] = {}

    for f in sorted(csv_dir.glob("deg_*.csv")):
        if "summary" in f.stem or "enrichment" in f.stem:
            continue
        pert_name = f.stem[4:]  # strip "deg_" prefix
        try:
            df = pd.read_csv(f)
            if "gene" not in df.columns:
                continue
            df = df.set_index("gene")
            if "log2fc" in df.columns:
                lfc_dict[pert_name] = df["log2fc"]
            if "padj" in df.columns:
                padj_dict[pert_name] = df["padj"]
        except Exception:
            continue

    if not lfc_dict:
        return pd.DataFrame(), pd.DataFrame()

    log2fc_df = pd.DataFrame(lfc_dict).T.fillna(0.0)
    padj_df = pd.DataFrame(padj_dict).T.fillna(1.0)
    return log2fc_df, padj_df


# ---------------------------------------------------------------------------
# 1. Perturbation effect correlation
# ---------------------------------------------------------------------------

def _compute_pert_corr(log2fc_df: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlation between perturbation effect vectors (rows = perts)."""
    if log2fc_df.shape[0] < 2 or log2fc_df.shape[1] < 2:
        return pd.DataFrame()
    corr = _pearson_rows(log2fc_df.values.astype(float))
    return pd.DataFrame(corr, index=log2fc_df.index, columns=log2fc_df.index)


def _compute_pert_corr_adjusted(
    adata,
    perts: List[str],
    genes: List[str],
) -> pd.DataFrame:
    """Cell-state-adjusted perturbation effect correlation.

    For each perturbation computes a cell-state-weighted log2FC:
        log2fc_adj[p][g] = Σ_s (w_s · log2fc(p,g | state s))
    where w_s = proportion of control cells in state s.
    This removes the state-composition confound from the correlation.
    """
    if "cell_state" not in adata.obs.columns:
        return pd.DataFrame()

    states = adata.obs["cell_state"].astype(str).unique().tolist()
    ctrl_mask = adata.obs["is_control"].values.astype(bool)
    n_ctrl = int(ctrl_mask.sum())
    if n_ctrl == 0:
        return pd.DataFrame()

    ctrl_states = adata.obs.loc[ctrl_mask, "cell_state"].astype(str)
    ctrl_state_fracs = {s: float((ctrl_states == s).sum()) / n_ctrl for s in states}

    gene_idx_map = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    valid_genes = [g for g in genes if g in gene_idx_map]
    if not valid_genes:
        return pd.DataFrame()
    gidx = np.array([gene_idx_map[g] for g in valid_genes])

    adj_lfc: Dict[str, np.ndarray] = {}
    for pert in perts:
        p_mask = adata.obs["perturbation"].values == pert
        if p_mask.sum() < 5:
            continue
        weighted = np.zeros(len(valid_genes))
        total_w = 0.0
        for s in states:
            s_mask = adata.obs["cell_state"].astype(str).values == s
            ps = p_mask & s_mask
            cs = ctrl_mask & s_mask
            if ps.sum() < 3 or cs.sum() < 3:
                continue
            xp = _to_dense(adata.X[ps, :][:, gidx]).mean(axis=0)
            xc = _to_dense(adata.X[cs, :][:, gidx]).mean(axis=0)
            lfc = np.log2((xp + 1e-3) / (xc + 1e-3))
            w = ctrl_state_fracs.get(s, 0.0)
            weighted += w * lfc
            total_w += w
        if total_w > 0:
            adj_lfc[pert] = weighted / total_w

    if not adj_lfc:
        return pd.DataFrame()

    adj_df = pd.DataFrame(adj_lfc, index=valid_genes).T
    return _compute_pert_corr(adj_df)


# ---------------------------------------------------------------------------
# 2. TF-TF regulatory heatmap
# ---------------------------------------------------------------------------

def _build_tf_regulatory_matrix(
    log2fc_df: pd.DataFrame,
    padj_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build signed -log10(padj) guide × TF-gene regulatory matrix.

    Only genes that also appear in the perturbation set are used as columns
    (they are both perturbed TFs and observed targets).

    Sign convention:
      +value (red)  = KO increases target expression (TF normally inhibits).
      -value (blue) = KO decreases target expression (TF normally activates).
    """
    perts_set = set(log2fc_df.index)
    genes_set = set(log2fc_df.columns)
    tf_genes = sorted(perts_set & genes_set)

    if not tf_genes:
        return pd.DataFrame()

    rows = []
    row_names = []
    for pert in log2fc_df.index:
        row = []
        for gene in tf_genes:
            lfc = float(log2fc_df.loc[pert, gene])
            padj = float(padj_df.loc[pert, gene]) if gene in padj_df.columns else 1.0
            sign = 1.0 if lfc >= 0 else -1.0
            val = round(-log10(max(padj, 1e-300)) * sign, 4)
            row.append(val)
        rows.append(row)
        row_names.append(pert)

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows, index=row_names, columns=tf_genes)


# ---------------------------------------------------------------------------
# 3. TF regulatory network
# ---------------------------------------------------------------------------

def _build_tf_network(
    log2fc_df: pd.DataFrame,
    padj_df: pd.DataFrame,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.3,
) -> dict:
    """Build directed TF regulatory network.

    Edge convention (matching biology):
      KO of TF_A → TF_B goes DOWN → A activates B → blue "activating" edge.
      KO of TF_A → TF_B goes UP   → A inhibits  B → red  "inhibiting" edge.
    """
    perts_set = set(log2fc_df.index)
    tf_genes = sorted(perts_set & set(log2fc_df.columns))
    if len(tf_genes) < 2:
        return {"nodes": [], "positions": {}, "edges_activating": [], "edges_inhibiting": []}

    edges_act: list = []
    edges_inh: list = []
    for pert in log2fc_df.index:
        for gene in tf_genes:
            if gene == pert:
                continue
            lfc = float(log2fc_df.loc[pert, gene])
            padj = float(padj_df.loc[pert, gene]) if gene in padj_df.columns else 1.0
            if padj >= fdr_threshold or abs(lfc) < lfc_threshold:
                continue
            edge = {
                "source": pert,
                "target": gene,
                "weight": round(abs(lfc), 3),
                "lfc": round(lfc, 3),
            }
            if lfc < 0:
                edges_act.append(edge)  # KO→↓ = TF activates target
            else:
                edges_inh.append(edge)  # KO→↑ = TF inhibits target

    # Compute spring-layout positions with networkx
    positions: dict = {}
    try:
        import networkx as nx
        G = nx.DiGraph()
        G.add_nodes_from(tf_genes)
        for e in edges_act + edges_inh:
            G.add_edge(e["source"], e["target"], weight=e["weight"])
        raw_pos = nx.spring_layout(G, seed=0, k=2.0)
        positions = {
            g: [round(float(xy[0]), 4), round(float(xy[1]), 4)]
            for g, xy in raw_pos.items()
        }
    except Exception:
        pass

    return {
        "nodes": tf_genes,
        "positions": positions,
        "edges_activating": edges_act,
        "edges_inhibiting": edges_inh,
    }


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _plot_corr_heatmap(corr_df: pd.DataFrame, title: str, out_path: Path) -> None:
    """Clustered heatmap of perturbation-effect Pearson correlation."""
    if corr_df.empty:
        return
    n = len(corr_df)
    figsize = (max(6, n * 0.4 + 2), max(5, n * 0.35 + 2))
    tick_fs = max(5, min(_TICK_FS, 120 // max(n, 1)))
    try:
        g = sns.clustermap(
            corr_df,
            cmap="RdBu_r",
            center=0,
            vmin=-1.0,
            vmax=1.0,
            figsize=figsize,
            xticklabels=(n <= 60),
            yticklabels=(n <= 60),
            cbar_kws={"label": "Pearson correlation of β between guides", "shrink": 0.5},
            linewidths=0.5 if n < 20 else 0,
        )
        g.ax_heatmap.set_title(title, fontsize=_TITLE_FS, pad=12)
        g.ax_heatmap.set_xlabel("Guides", fontsize=_LABEL_FS)
        g.ax_heatmap.set_ylabel(
            "Effect of guide on its target expression", fontsize=_LABEL_FS
        )
        g.ax_heatmap.tick_params(axis="x", labelsize=tick_fs, rotation=90)
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_fs)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        _save_png(g.fig, out_path)
        plt.close(g.fig)
    except Exception:
        plt.close("all")


def _plot_tf_regulatory_heatmap(tf_mat: pd.DataFrame, out_path: Path) -> None:
    """Clustered heatmap: guides × TF genes, signed -log10(q)."""
    if tf_mat.empty:
        return
    n_r, n_c = tf_mat.shape
    figw = max(4.5 + n_c * 0.28, 5.5)
    figh = max(5.5 + n_r * 0.85, 7.0)
    tick_r = max(5, min(_TICK_FS, 160 // max(n_r, 1)))
    tick_c = max(5, min(_TICK_FS, 160 // max(n_c, 1)))
    vmax = min(max(abs(tf_mat.values).max(), 1.0), 8.0)

    try:
        g = sns.clustermap(
            tf_mat,
            row_cluster=(n_r > 2),
            col_cluster=(n_c > 2),
            cmap="RdBu_r",
            center=0,
            vmin=-vmax,
            vmax=vmax,
            figsize=(figw, figh),
            xticklabels=False,
            yticklabels=(n_r <= 80),
            cbar_kws={"label": "signed −log₁₀(q-value)", "shrink": 0.5},
            linewidths=0.5 if n_r < 30 else 0,
        )
        g.ax_heatmap.set_title(
            "TF regulatory relationships", fontsize=_TITLE_FS, pad=12
        )
        g.ax_heatmap.set_xlabel("Genes", fontsize=_LABEL_FS)
        g.ax_heatmap.set_ylabel("Guides", fontsize=_LABEL_FS)
        g.ax_heatmap.tick_params(axis="x", labelsize=tick_c, rotation=90)
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_r)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        handles = [
            mpatches.Patch(color="#1a9850", label="Activating (KO→↓)"),
            mpatches.Patch(color="#d73027", label="Inhibiting (KO→↑)"),
        ]
        g.fig.legend(
            handles=handles, loc="lower left", fontsize=_LEGEND_FS,
            frameon=True, title="Effect of KO"
        )
        _save_png(g.fig, out_path)
        plt.close(g.fig)
    except Exception:
        plt.close("all")


def _plot_tf_network(network: dict, out_path: Path) -> None:
    """Directed TF regulatory network.

    Edge colours follow the publication-grade ColorBrewer RdYlGn pair:
      green = activating (KO → target ↓ → TF normally activates target)
      red   = inhibiting (KO → target ↑ → TF normally inhibits target)
    Always writes a PNG. If the regulatory map is empty, the PNG is a
    text panel explaining why (so the report has a concrete artefact to
    show instead of a missing image).
    """
    import networkx as nx

    _C_ACT = "#1a9850"   # green — activating (matches GRN palette)
    _C_INH = "#d73027"   # red   — inhibiting

    nodes = network.get("nodes", [])

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    for e in network.get("edges_activating", []):
        G.add_edge(e["source"], e["target"], weight=e["weight"], kind="act")
    for e in network.get("edges_inhibiting", []):
        G.add_edge(e["source"], e["target"], weight=e["weight"], kind="inh")

    # ---- Empty-network placeholder PNG ----
    if not nodes or not G.edges():
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.text(
            0.5, 0.55, "No significant TF→target edges",
            fontsize=22, fontweight="bold",
            ha="center", va="center", color="#444444",
        )
        ax.text(
            0.5, 0.40,
            "Either no perturbed gene is in HumanTFs, or no edges met the\n"
            "FDR / |log₂FC| thresholds. Lower the regulatory_fdr_threshold\n"
            "or regulatory_lfc_threshold in the config to see more edges.",
            fontsize=14, ha="center", va="center", color="#666666",
        )
        ax.set_title("TF regulatory network", fontsize=18, pad=10)
        ax.set_xlabel("Layout X", fontsize=12)
        ax.set_ylabel("Layout Y", fontsize=12)
        ax.grid(True, color="#e7edf2", linewidth=0.8)
        fig.tight_layout()
        _save_png(fig, out_path, tight=False)
        plt.close(fig)
        return

    pos = network.get("positions") or {}
    if not pos:
        try:
            pos = nx.spring_layout(G, seed=0, k=2.0)
        except Exception:
            pos = {n: [float(i), 0.0] for i, n in enumerate(nodes)}

    in_deg = dict(G.in_degree())
    out_deg = dict(G.out_degree())
    node_sizes = [max(800, (in_deg.get(n, 0) + out_deg.get(n, 0)) * 200 + 800)
                  for n in nodes]

    fig, ax = plt.subplots(figsize=(9, 9))
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_size=node_sizes,
                           node_color="#f4a261", alpha=0.9,
                           edgecolors="#1d3557", linewidths=1.5, ax=ax)
    if len(nodes) <= 80:
        nx.draw_networkx_labels(
            G, pos, font_size=18, font_weight="bold", ax=ax,
        )

    act_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("kind") == "act"]
    inh_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("kind") == "inh"]

    def _w(u, v):
        return max(0.8, G[u][v].get("weight", 0.5) * 1.6)

    if act_edges:
        nx.draw_networkx_edges(
            G, pos, edgelist=act_edges,
            edge_color=_C_ACT, width=[_w(u, v) for u, v in act_edges],
            alpha=0.8, arrows=True, arrowsize=12,
            connectionstyle="arc3,rad=0.12", ax=ax,
        )
    if inh_edges:
        nx.draw_networkx_edges(
            G, pos, edgelist=inh_edges,
            edge_color=_C_INH, width=[_w(u, v) for u, v in inh_edges],
            alpha=0.8, arrows=True, arrowsize=12,
            connectionstyle="arc3,rad=-0.12", ax=ax,
        )

    handles = [
        mpatches.Patch(color=_C_ACT, label="Activating (KO → target ↓)"),
        mpatches.Patch(color=_C_INH, label="Inhibiting (KO → target ↑)"),
    ]
    ax.legend(handles=handles, loc="upper left",
              fontsize=10, frameon=True)
    ax.set_title("TF regulatory network", fontsize=16)
    ax.set_xlabel("Layout X", fontsize=12)
    ax.set_ylabel("Layout Y", fontsize=12)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, color="#e7edf2", linewidth=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    _save_png(fig, out_path, tight=False)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_regulatory_analysis(
    adata,
    output_dir: str,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.3,
) -> dict:
    """Run perturbation effect correlation and TF regulatory network analysis.

    Requires that the DEG step has already been run (reads deg_*.csv files).

    Args:
        adata           -- AnnData after QC + normalisation + DEG step.
        output_dir      -- Root output directory (same as pipeline output).
        fdr_threshold   -- FDR cutoff for TF regulatory edges.
        lfc_threshold   -- |log2FC| cutoff for TF regulatory edges.

    Returns:
        dict with keys: corr_before, corr_after, tf_regulatory_matrix, tf_network.
    """
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    log2fc_df, padj_df = _load_deg_matrix(tables)
    if log2fc_df.empty:
        return {}

    perts = log2fc_df.index.tolist()
    genes = log2fc_df.columns.tolist()
    results: dict = {}

    # 1. Perturbation effect correlation — before cell-state correction
    corr_before = _compute_pert_corr(log2fc_df)
    if not corr_before.empty:
        corr_before.to_csv(tables / "pert_effect_corr_before.csv")
        _plot_corr_heatmap(
            corr_before,
            "Without modeling cell states",
            plots / "pert_effect_corr_before.png",
        )
        results["corr_before"] = corr_before

    # 2. Perturbation effect correlation — after cell-state correction
    corr_after = _compute_pert_corr_adjusted(adata, perts, genes)
    if not corr_after.empty:
        corr_after.to_csv(tables / "pert_effect_corr_after.csv")
        _plot_corr_heatmap(
            corr_after,
            "After modeling cell states",
            plots / "pert_effect_corr_after.png",
        )
        results["corr_after"] = corr_after

    # 3. TF-TF regulatory heatmap
    tf_mat = _build_tf_regulatory_matrix(log2fc_df, padj_df)
    if not tf_mat.empty:
        tf_mat.to_csv(tables / "tf_regulatory_matrix.csv")
        _plot_tf_regulatory_heatmap(tf_mat, plots / "tf_regulatory_heatmap.png")
        results["tf_regulatory_matrix"] = tf_mat

    # 4. TF regulatory network (JSON for interactive report + static PNG)
    network = _build_tf_network(log2fc_df, padj_df, fdr_threshold, lfc_threshold)
    if network["nodes"]:
        (tables / "tf_regulatory_network.json").write_text(
            json.dumps(network, ensure_ascii=False, separators=(",", ":"))
        )
        _plot_tf_network(network, plots / "tf_regulatory_network.png")
        results["tf_network"] = network

    return results
