# =============================================================================
# perturbflow/analyzer/genenet.py
#
# Gene co-expression network module for Perturb-seq data.
#
# For each top perturbation, this module:
#   1. Clusters genes (based on control-cell expression) via hierarchical
#      clustering of Pearson correlation → identifies gene modules.
#   2. Plots a side-by-side heatmap: control vs perturbed (same gene order,
#      sorted by module) so expression changes are immediately visible.
#   3. Builds a co-expression network (nodes = genes, edges = |r| > threshold)
#      for both conditions using NetworkX, then draws:
#        • Control network — edges coloured by correlation sign.
#        • Perturbed network — same layout, edges reflect perturbed correlations.
#        • Differential overlay — gained edges (red), lost edges (green),
#          unchanged edges (grey), and nodes scaled by |mean expression change|.
#   4. Writes a CSV of gene-cluster assignments.
#
# Outputs written to <output_dir>/plots/ and <output_dir>/csv/:
#
#   plots/
#     genenet_control_heatmap.png           -- gene cluster heatmap (control)
#     genenet_<pert>_heatmap_comparison.png -- control vs perturbed side-by-side
#     genenet_control_network.png           -- network graph (control)
#     genenet_<pert>_network.png            -- network graph (perturbed)
#     genenet_<pert>_diff_network.png       -- differential network
#   csv/
#     genenet_gene_clusters.csv             -- gene, cluster_id, mean_ctrl,
#                                             mean_pert, log2fc columns
#
# Usage in pipeline (wired by pipeline.py):
#   from .genenet import run_gene_network
#   run_gene_network(adata, output_dir, perturbations=top_perts)
# =============================================================================

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib import font_manager
import seaborn as sns

from .utils import ensure_dir

# ---------------------------------------------------------------------------
# Global font constants
# ---------------------------------------------------------------------------
_TITLE_FS  = 18
_LABEL_FS  = 15
_TICK_FS   = 13
_LEGEND_FS = 13
_CBAR_FS   = 13


def _preferred_network_font() -> str:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for family in ("Arial", "Liberation Sans", "DejaVu Sans"):
        if family in available:
            return family
    return "sans-serif"

# ---------------------------------------------------------------------------
# Publication-grade colour palette
# ---------------------------------------------------------------------------
# Static co-expression panels use red for positive correlations and blue for
# negative correlations. Differential panels keep red/green for gained/lost
# edges because those are edge-status categories, not correlation signs.
_C_RED   = "#ef8a8a"   # light red (positive correlation / gained edge)
_C_BLUE  = "#67a9cf"   # light blue (negative correlation)
_C_GREEN = "#91cf60"   # light green (lost edge)
_C_GREY  = "#cccccc"

# Sign of Pearson r in static co-expression panels:
#   positive r  → red,
#   negative r  → blue.
_C_POS_R = _C_RED
_C_NEG_R = _C_BLUE
# Differential edge status:
#   gained edge → red   ("added" in perturbation),
#   lost   edge → green ("removed" in perturbation),
#   unchanged   → grey.
_C_DIFF_GAINED = _C_RED
_C_DIFF_LOST   = _C_GREEN
_C_DIFF_SHARED = _C_GREY
# Node outlines used to keep nodes legible against light backgrounds.
_C_NODE_OUTLINE = "#1d3557"
_C_TF_OUTLINE   = "#000000"

# ---------------------------------------------------------------------------
# Node filter: drop genes that are isolated in BOTH the control and the
# perturbation networks. A node is kept whenever it has at least one edge
# (|r| >= corr_threshold) in EITHER network. The filter uses the same
# threshold as edge construction so it stays consistent — the legacy
# _NODE_FILTER_THRESHOLD constant is no longer used.
# ---------------------------------------------------------------------------
_NODE_FILTER_THRESHOLD = 0.2  # retained only for backward compatibility


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _to_dense(x):
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def _subset_expr(adata, mask: np.ndarray, gene_idx: np.ndarray) -> np.ndarray:
    """Return dense (cells, genes) array for a boolean cell mask."""
    return _to_dense(adata.X[np.where(mask)[0], :][:, gene_idx])


def _pearson_corrmat(x: np.ndarray) -> np.ndarray:
    """Fast Pearson correlation matrix for a (cells × genes) array."""
    corr = np.corrcoef(x.T)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def _cluster_genes(corr: np.ndarray, n_clusters: int) -> np.ndarray:
    """Hierarchical clustering on 1-|corr| distance; returns cluster labels."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    dist = np.clip(1.0 - np.abs(corr), 0, None)
    np.fill_diagonal(dist, 0.0)
    # squareform needs a condensed distance vector.
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    return labels  # 1-indexed


def _sorted_order(labels: np.ndarray) -> np.ndarray:
    """Argsort so that genes within the same cluster are contiguous."""
    return np.argsort(labels, kind="stable")


def _critical_r(n_cells: int, alpha: float = 0.05, n_tests: int = 1) -> float:
    """Two-sided critical |Pearson r| at significance ``α/n_tests`` for ``n_cells``
    samples (Bonferroni-corrected when ``n_tests > 1``).

    Uses ``t = r·√((n-2)/(1-r²))`` ↔ Student-t with df = n − 2. Returns 1.0
    (nothing passes) when ``n_cells ≤ 3`` or any input is invalid.
    """
    if n_cells is None or n_cells <= 3:
        return 1.0
    from scipy.stats import t as student_t
    df = int(n_cells) - 2
    alpha_eff = float(alpha) / max(1, int(n_tests))
    if alpha_eff <= 0 or alpha_eff >= 1:
        return 1.0
    t_crit = float(student_t.ppf(1.0 - alpha_eff / 2.0, df))
    return float(np.sqrt(t_crit * t_crit / (t_crit * t_crit + df)))


def _build_network(
    corr: np.ndarray,
    gene_names: List[str],
    threshold: float,
    n_cells: Optional[int] = None,
    alpha: float = 0.05,
    multiple_testing: str = "none",
):
    """Build a NetworkX graph from a correlation matrix.

    By default an edge is added whenever ``|r| >= threshold`` — the user's
    threshold is honoured exactly, regardless of ``n_cells``. To opt-in to
    sample-size-aware Bonferroni correction (effective threshold raised to
    ``max(threshold, r_crit(n_cells, α/n_tests))``), pass
    ``multiple_testing="bonferroni"`` and an integer ``n_cells``.

    Edge attributes: ``weight`` (signed r), ``positive`` (bool sign).
    """
    import networkx as nx

    n_genes = len(gene_names)
    n_tests = max(1, n_genes * (n_genes - 1) // 2)

    if multiple_testing == "bonferroni" and n_cells is not None:
        rc = _critical_r(int(n_cells), alpha, n_tests=n_tests)
        eff = max(float(threshold), rc)
    else:
        rc = None
        eff = float(threshold)

    G = nx.Graph()
    G.add_nodes_from(gene_names)
    G.graph["threshold"] = float(threshold)
    G.graph["effective_threshold"] = float(eff)
    G.graph["n_cells"] = None if n_cells is None else int(n_cells)
    G.graph["n_tests"] = n_tests
    G.graph["multiple_testing"] = multiple_testing
    G.graph["alpha"] = float(alpha)
    G.graph["r_crit"] = rc
    n = n_genes
    for i in range(n):
        for j in range(i + 1, n):
            r = float(corr[i, j])
            if abs(r) >= eff:
                G.add_edge(gene_names[i], gene_names[j],
                           weight=r, positive=(r >= 0))
    return G


def _spring_layout(G, seed: int = 0):
    """Return a dict of node positions using spring layout."""
    import networkx as nx
    return nx.spring_layout(G, weight="weight", seed=seed, k=2.4, scale=1.8, iterations=100)


def _safe_name(text: str) -> str:
    return str(text).replace("/", "_").replace(" ", "_").replace("+", "_")


def _deg_prioritized_gene_indices(
    adata,
    tables_dir: Path,
    perturbation: str,
    hv_idx: np.ndarray,
    n_top_genes: int,
    deg_fraction: float = 0.75,
) -> tuple[np.ndarray, dict]:
    """Select network genes using top DEGs first, then HVGs as background.

    Global HVGs are excellent for unsupervised structure, but they can make
    a perturbation-specific network look biologically generic. For a focused
    network, start with genes whose expression actually changes under this
    perturbation, then fill the panel with high-variance genes so there is
    enough co-expression structure to draw.
    """
    var_lookup = {str(g): i for i, g in enumerate(adata.var_names)}
    selected: list[int] = []
    source: dict[str, str] = {}

    def _add_gene(gene: str, tag: str) -> None:
        idx = var_lookup.get(str(gene))
        if idx is None or idx in selected:
            return
        selected.append(idx)
        source[str(gene)] = tag

    safe = _safe_name(perturbation)
    deg_path = tables_dir / f"deg_{safe}.csv"
    n_deg_target = max(10, int(round(n_top_genes * deg_fraction)))
    if deg_path.exists():
        try:
            deg = pd.read_csv(deg_path)
            if "gene" in deg.columns:
                deg["_abs_log2fc"] = pd.to_numeric(
                    deg.get("log2fc", 0), errors="coerce"
                ).abs().fillna(0)
                deg["_padj"] = pd.to_numeric(
                    deg.get("padj", 1), errors="coerce"
                ).fillna(1)
                if "significant" in deg.columns:
                    deg["_sig"] = deg["significant"].astype(str).str.lower().isin(
                        {"true", "1", "yes"}
                    ).astype(int)
                else:
                    deg["_sig"] = (deg["_padj"] <= 0.05).astype(int)
                deg = deg.sort_values(
                    ["_sig", "_padj", "_abs_log2fc"],
                    ascending=[False, True, False],
                )
                for gene in deg["gene"].astype(str).head(n_deg_target):
                    _add_gene(gene, "DEG")
        except Exception:
            pass

    # Always include the perturbed gene if it is measured.
    _add_gene(perturbation, "perturbed_gene")

    for idx in hv_idx:
        if len(selected) >= n_top_genes:
            break
        gene = str(adata.var_names[int(idx)])
        _add_gene(gene, "HVG")

    return np.asarray(selected[:n_top_genes], dtype=int), source


def _draw_network(
    G,
    pos: dict,
    ax,
    node_sizes: dict,
    title: str,
    show_labels: bool = True,
    default_node_size: float = 900,
    node_color="#2a9d8f",
    label_font_size: int = 10,
) -> None:
    """Draw a co-expression network on the given axis.

    ``node_sizes`` is a per-gene dict; values are passed straight to
    NetworkX (matplotlib points squared). The default values used by callers
    here are intentionally large because each network is rendered as its own
    square panel and the genes are typically a small set (≤ 50).
    """
    import networkx as nx

    node_list = list(G.nodes())
    sizes = [node_sizes.get(n, default_node_size) for n in node_list]

    # Edge colours: positive = coral, negative = steel-blue.
    edge_colors = [
        _C_POS_R if G[u][v].get("positive", True) else _C_NEG_R
        for u, v in G.edges()
    ]
    edge_widths = [
        max(0.8, abs(G[u][v].get("weight", 0)) * 3.0)
        for u, v in G.edges()
    ]

    nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_size=sizes,
                           node_color=node_color, alpha=0.9,
                           edgecolors="#1d3557", linewidths=1.0, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors,
                           width=edge_widths, alpha=0.65, ax=ax)
    if show_labels and len(node_list) <= 80:
        nx.draw_networkx_labels(G, pos, font_size=label_font_size,
                                font_weight="normal", ax=ax)
    ax.set_title(title, fontsize=_TITLE_FS)
    ax.axis("off")


# ---------------------------------------------------------------------------
# Plot 1 — control gene cluster heatmap (one-time, all top-pert genes)
# ---------------------------------------------------------------------------

def _plot_control_heatmap(
    ctrl_corr: np.ndarray,
    gene_names: List[str],
    sorted_idx: np.ndarray,
    cluster_labels: np.ndarray,
    plots_dir: Path,
) -> None:
    """Clustered heatmap of gene-gene correlation in control cells,
    with gene-cluster sidebar annotation."""
    n_g = len(gene_names)
    if n_g < 2:
        return

    ordered_genes = [gene_names[i] for i in sorted_idx]
    ordered_corr = ctrl_corr[np.ix_(sorted_idx, sorted_idx)]
    df = pd.DataFrame(ordered_corr, index=ordered_genes, columns=ordered_genes)

    n_clusters = int(cluster_labels.max())
    palette = sns.color_palette("tab10", min(n_clusters, 10))
    cluster_pal = {str(c): palette[(c - 1) % len(palette)] for c in range(1, n_clusters + 1)}
    row_colors = pd.Series(
        [str(cluster_labels[i]) for i in sorted_idx],
        index=ordered_genes,
        name="Gene module",
    ).map(cluster_pal)

    tick_fs = max(5, min(_TICK_FS, 120 // max(n_g, 1)))
    try:
        g = sns.clustermap(
            df,
            row_cluster=False, col_cluster=False,   # already sorted by module
            row_colors=row_colors,
            col_colors=row_colors,
            cmap="RdBu_r", center=0, vmin=-1, vmax=1,
            xticklabels=False, yticklabels=False,
            figsize=(max(12, n_g // 3), max(10, n_g // 3)),
            cbar_kws={"label": "Pearson r", "shrink": 0.5},
        )
        g.ax_heatmap.set_title("")
        g.fig.suptitle(
            "Control Gene Co-expression Heatmap (modules highlighted)",
            fontsize=_TITLE_FS, y=0.995,
        )
        g.ax_heatmap.tick_params(axis="x", labelsize=tick_fs, rotation=90)
        g.ax_heatmap.tick_params(axis="y", labelsize=tick_fs)
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        g.fig.subplots_adjust(top=0.92)
        # Legend for gene modules.
        handles = [mpatches.Patch(color=cluster_pal[str(c)], label=f"Module {c}")
                   for c in range(1, n_clusters + 1)]
        g.fig.legend(handles=handles, title="Gene module", title_fontsize=_LEGEND_FS,
                     loc="lower left", fontsize=_LEGEND_FS, frameon=True)
        g.fig.savefig(plots_dir / "genenet_control_heatmap.png", dpi=120, bbox_inches="tight")
        plt.close(g.fig)
    except Exception:
        plt.close("all")


# ---------------------------------------------------------------------------
# Plot 2 — side-by-side expression heatmap: control vs perturbed
# ---------------------------------------------------------------------------

def _plot_comparison_heatmap(
    ctrl_expr: np.ndarray,
    pert_expr: np.ndarray,
    gene_names: List[str],
    sorted_idx: np.ndarray,
    cluster_labels: np.ndarray,
    perturbation: str,
    plots_dir: Path,
) -> None:
    """Side-by-side mean-expression heatmaps (control | perturbed) with
    genes sorted by control-derived module and a log2FC panel."""
    ordered_genes = [gene_names[i] for i in sorted_idx]
    n_g = len(ordered_genes)
    if n_g == 0:
        return

    # Per-gene means (z-scored for visual comparability).
    def _row_mean_z(expr, idx):
        m = expr[:, idx].mean(axis=0)
        z = (m - m.mean()) / (m.std() + 1e-8)
        return np.clip(z, -3, 3)

    ctrl_z  = _row_mean_z(ctrl_expr,  sorted_idx)
    pert_z  = _row_mean_z(pert_expr,  sorted_idx)
    log2fc  = (
        (pert_expr[:, sorted_idx].mean(axis=0) - ctrl_expr[:, sorted_idx].mean(axis=0))
        / np.log(2)
    )

    # Build a (3-panel) DataFrame: control z, perturbed z, log2fc.
    panel = pd.DataFrame(
        {"Control (z)": ctrl_z, f"{perturbation} (z)": pert_z, "log₂FC": log2fc},
        index=ordered_genes,
    )

    n_clusters = int(cluster_labels.max())
    palette = sns.color_palette("tab10", min(n_clusters, 10))
    row_colors = pd.Series(
        [str(cluster_labels[i]) for i in sorted_idx],
        index=ordered_genes,
        name="Module",
    ).map({str(c): palette[(c - 1) % len(palette)] for c in range(1, n_clusters + 1)})

    tick_fs = max(5, min(_TICK_FS, 120 // max(n_g, 1)))
    try:
        g = sns.clustermap(
            panel,
            row_cluster=False, col_cluster=False,
            row_colors=row_colors,
            # RdBu_r divergent palette + tighter range so small differences
            # between control and perturbation contrast clearly.
            cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5,
            xticklabels=True, yticklabels=False,
            figsize=(7.5, max(7, n_g // 3)),
            cbar_kws={"label": "z-score / log₂FC", "shrink": 0.5},
        )
        g.ax_heatmap.set_title(
            f"Gene expression — Control vs {perturbation}",
            fontsize=_TITLE_FS, pad=12,
        )
        g.ax_heatmap.tick_params(axis="x", labelsize=_LABEL_FS + 2)
        g.ax_heatmap.tick_params(axis="y", labelsize=max(tick_fs, _TICK_FS))
        g.cax.yaxis.label.set_size(_CBAR_FS)
        g.cax.tick_params(labelsize=_TICK_FS)
        safe = perturbation.replace("/", "_").replace(" ", "_").replace("+", "_")
        g.fig.savefig(
            plots_dir / f"genenet_{safe}_heatmap_comparison.png",
            dpi=120, bbox_inches="tight",
        )
        plt.close(g.fig)
    except Exception:
        plt.close("all")


# ---------------------------------------------------------------------------
# Plot 3 — network graphs: control, perturbed, differential
# ---------------------------------------------------------------------------

def _plot_networks(
    ctrl_corr: np.ndarray,
    pert_corr: np.ndarray,
    gene_names: List[str],
    cluster_labels: np.ndarray,
    ctrl_means: np.ndarray,
    pert_means: np.ndarray,
    perturbation: str,
    plots_dir: Path,
    corr_threshold: float,
    n_ctrl_cells: Optional[int] = None,
    n_pert_cells: Optional[int] = None,
) -> None:
    """Draw control network, perturbed network, and differential network.

    The user-supplied ``corr_threshold`` is combined with a sample-size
    aware critical value (Pearson significance at α = 0.05) so that small
    per-condition samples don't produce near-clique networks driven by
    sampling noise. The control panel is drawn at the control's effective
    threshold, the perturbation panel at its own (typically higher when
    n_pert_cells is small), and the differential panel uses the **larger**
    of the two so a "gained / lost" classification is symmetric.
    """
    import networkx as nx

    # ------------------------------------------------------------------
    # Node pre-filter: keep genes that have at least one edge with
    # abs(r) >= corr_threshold in EITHER the control or the perturbation
    # network. Genes that are isolated in both conditions are dropped before
    # any drawing so the panels stay clean and readable.
    # ------------------------------------------------------------------
    G_filt_ctrl = _build_network(ctrl_corr, gene_names, corr_threshold)
    G_filt_pert = _build_network(pert_corr, gene_names, corr_threshold)
    active_nodes = [
        g for g in gene_names
        if G_filt_ctrl.degree(g) > 0 or G_filt_pert.degree(g) > 0
    ]
    if len(active_nodes) < 2:
        return  # nothing to plot after filtering

    if len(active_nodes) < len(gene_names):
        active_idx = np.array([gene_names.index(g) for g in active_nodes])
        ctrl_corr      = ctrl_corr[np.ix_(active_idx, active_idx)]
        pert_corr      = pert_corr[np.ix_(active_idx, active_idx)]
        cluster_labels = cluster_labels[active_idx]
        ctrl_means     = ctrl_means[active_idx]
        pert_means     = pert_means[active_idx]
        gene_names     = active_nodes
    # ------------------------------------------------------------------

    G_ctrl = _build_network(ctrl_corr, gene_names, corr_threshold,
                            n_cells=n_ctrl_cells)
    G_pert = _build_network(pert_corr, gene_names, corr_threshold,
                            n_cells=n_pert_cells)
    eff_ctrl = G_ctrl.graph.get("effective_threshold", corr_threshold)
    eff_pert = G_pert.graph.get("effective_threshold", corr_threshold)

    if len(G_ctrl.edges()) == 0 and len(G_pert.edges()) == 0:
        return

    # Use the union graph to compute a shared layout.
    G_union = nx.Graph()
    G_union.add_nodes_from(gene_names)
    G_union.add_edges_from(G_ctrl.edges())
    G_union.add_edges_from(G_pert.edges())
    pos = _spring_layout(G_union)

    # ---- Differential network (built first so its degree feeds node sizes) ----
    # Each diff edge carries the relevant Pearson r in `weight` so the panel's
    # line width can reflect correlation magnitude:
    #   shared → max(|r_ctrl|, |r_pert|)
    #   lost   → |r_ctrl|
    #   gained → |r_pert|
    G_diff = nx.Graph()
    G_diff.add_nodes_from(gene_names)
    ctrl_edges = set(
        (u, v) if u < v else (v, u) for u, v in G_ctrl.edges()
    )
    pert_edges = set(
        (u, v) if u < v else (v, u) for u, v in G_pert.edges()
    )
    for u, v in ctrl_edges | pert_edges:
        rc = abs(G_ctrl[u][v]["weight"]) if (u, v) in ctrl_edges else 0.0
        rp = abs(G_pert[u][v]["weight"]) if (u, v) in pert_edges else 0.0
        if (u, v) in ctrl_edges and (u, v) in pert_edges:
            status = "shared"; w = max(rc, rp)
        elif (u, v) in ctrl_edges:
            status = "lost";   w = rc
        else:
            status = "gained"; w = rp
        G_diff.add_edge(u, v, status=status, weight=float(w))

    # Node size = number of significant connections in the relevant graph
    # (Pearson partners surviving the Bonferroni-corrected |r| threshold).
    # Hubs naturally pop out; isolated genes stay small but still visible.
    _STATIC_NODE_BASE = 150
    _STATIC_NODE_SCALE = 45
    _STATIC_NODE_MAX = 950
    _DIFF_NODE_BASE = 130
    _DIFF_NODE_SCALE = 38
    _DIFF_NODE_MAX = 850
    ctrl_deg = dict(G_ctrl.degree())
    pert_deg = dict(G_pert.degree())
    diff_deg = dict(G_diff.degree())
    ctrl_size = {
        g: min(_STATIC_NODE_MAX, _STATIC_NODE_BASE + ctrl_deg.get(g, 0) * _STATIC_NODE_SCALE)
        for g in gene_names
    }
    pert_size = {
        g: min(_STATIC_NODE_MAX, _STATIC_NODE_BASE + pert_deg.get(g, 0) * _STATIC_NODE_SCALE)
        for g in gene_names
    }
    change_size = {
        g: min(_DIFF_NODE_MAX, _DIFF_NODE_BASE + diff_deg.get(g, 0) * _DIFF_NODE_SCALE)
        for g in gene_names
    }

    safe = perturbation.replace("/", "_").replace(" ", "_").replace("+", "_")
    show_lab = len(gene_names) <= 80

    # All three panels are rendered as square (12×12) single-axis figures so
    # they line up cleanly when the report places them in a 3-column grid.
    # Title / legend / label fonts are sized so that the labels remain
    # readable when the PNG is shown ~⅓ of report width.
    _PANEL_SIZE = 10.5
    _FONT_NET = _preferred_network_font()
    _LABEL_FS_NET = 11
    _LABEL_FS_DIFF = 10
    _TITLE_FS_NET = 18
    _LEGEND_FS_NET = 14

    # Per-gene module colours used for both the control and perturbed plots so
    # the same gene is the same colour across the row. Set2 keeps nodes light
    # enough that overlaid labels and edges remain legible while still being
    # distinguishable in print.
    n_clusters = int(cluster_labels.max())
    tab_pal = sns.color_palette("Set2", min(n_clusters, 8))
    name_to_color = {
        g: tab_pal[(cluster_labels[i] - 1) % len(tab_pal)]
        for i, g in enumerate(gene_names)
    }

    edge_legend = [
        mpatches.Patch(color=_C_POS_R, label="Positive correlation"),
        mpatches.Patch(color=_C_NEG_R, label="Negative correlation"),
    ]
    module_handles = [
        mpatches.Patch(color=tab_pal[(c - 1) % len(tab_pal)], label=f"Module {c}")
        for c in range(1, n_clusters + 1)
    ]

    # ---- Figure: control-only network ----
    fig_c, ax_c = plt.subplots(figsize=(_PANEL_SIZE, _PANEL_SIZE))
    ctrl_node_colors = [name_to_color[n] for n in G_ctrl.nodes()]
    ctrl_node_sizes_list = [ctrl_size.get(n, _STATIC_NODE_BASE) for n in G_ctrl.nodes()]
    edge_colors_c = [_C_POS_R if G_ctrl[u][v].get("positive", True) else _C_NEG_R
                     for u, v in G_ctrl.edges()]
    edge_widths_c = [max(0.3, abs(G_ctrl[u][v].get("weight", 0)) * 2.5)
                     for u, v in G_ctrl.edges()]
    nx.draw_networkx_nodes(G_ctrl, pos, node_color=ctrl_node_colors,
                           node_size=ctrl_node_sizes_list, alpha=0.7,
                           edgecolors=_C_NODE_OUTLINE, linewidths=0.8, ax=ax_c)
    nx.draw_networkx_edges(G_ctrl, pos, edge_color=edge_colors_c,
                           width=edge_widths_c, alpha=0.4, ax=ax_c)
    if show_lab:
        nx.draw_networkx_labels(G_ctrl, pos, font_size=_LABEL_FS_NET,
                                font_weight="normal", font_family=_FONT_NET, ax=ax_c)
    ax_c.set_title(
        f"Control — co-expression network\n"
        f"(|r| ≥ {eff_ctrl:.2f}, n={n_ctrl_cells} cells)",
        fontsize=_TITLE_FS_NET,
        fontfamily=_FONT_NET,
    )
    leg_c1 = ax_c.legend(handles=edge_legend, loc="lower right", ncol=1,
                         prop={"family": _FONT_NET, "size": _LEGEND_FS_NET}, frameon=True)
    ax_c.legend(handles=module_handles, title="Node color", loc="lower left",
                ncol=1, prop={"family": _FONT_NET, "size": _LEGEND_FS_NET - 2},
                title_fontproperties={"family": _FONT_NET, "size": _LEGEND_FS_NET - 1},
                frameon=True)
    ax_c.add_artist(leg_c1)
    ax_c.margins(0.18)
    ax_c.axis("off")
    fig_c.tight_layout()
    fig_c.savefig(plots_dir / f"genenet_{safe}_ctrl_network.png", dpi=120)
    plt.close(fig_c)

    # ---- Figure: perturbation-only network ----
    fig_p, ax_p = plt.subplots(figsize=(_PANEL_SIZE, _PANEL_SIZE))
    pert_node_colors = [name_to_color[n] for n in G_pert.nodes()]
    pert_node_sizes_list = [pert_size.get(n, _STATIC_NODE_BASE) for n in G_pert.nodes()]
    edge_colors_p = [_C_POS_R if G_pert[u][v].get("positive", True) else _C_NEG_R
                     for u, v in G_pert.edges()]
    edge_widths_p = [max(0.3, abs(G_pert[u][v].get("weight", 0)) * 2.5)
                     for u, v in G_pert.edges()]
    nx.draw_networkx_nodes(G_pert, pos, node_color=pert_node_colors,
                           node_size=pert_node_sizes_list, alpha=0.7,
                           edgecolors=_C_NODE_OUTLINE, linewidths=0.8, ax=ax_p)
    nx.draw_networkx_edges(G_pert, pos, edge_color=edge_colors_p,
                           width=edge_widths_p, alpha=0.4, ax=ax_p)
    if show_lab:
        nx.draw_networkx_labels(G_pert, pos, font_size=_LABEL_FS_NET,
                                font_weight="normal", font_family=_FONT_NET, ax=ax_p)
    ax_p.set_title(
        f"{perturbation} — co-expression network\n"
        f"(|r| ≥ {eff_pert:.2f}, n={n_pert_cells} cells)",
        fontsize=_TITLE_FS_NET,
        fontfamily=_FONT_NET,
    )
    leg_p1 = ax_p.legend(handles=edge_legend, loc="lower right", ncol=1,
                         prop={"family": _FONT_NET, "size": _LEGEND_FS_NET}, frameon=True)
    ax_p.legend(handles=module_handles, title="Node color", loc="lower left",
                ncol=1, prop={"family": _FONT_NET, "size": _LEGEND_FS_NET - 2},
                title_fontproperties={"family": _FONT_NET, "size": _LEGEND_FS_NET - 1},
                frameon=True)
    ax_p.add_artist(leg_p1)
    ax_p.margins(0.18)
    ax_p.axis("off")
    fig_p.tight_layout()
    fig_p.savefig(plots_dir / f"genenet_{safe}_pert_network.png", dpi=120)
    plt.close(fig_p)

    # ---- Figure: differential network (square, same panel size, smaller nodes) ----
    fig2, ax2 = plt.subplots(figsize=(_PANEL_SIZE, _PANEL_SIZE))
    diff_edge_colors = {
        "shared": _C_DIFF_SHARED,
        "gained": _C_DIFF_GAINED,
        "lost":   _C_DIFF_LOST,
    }
    for status, color in diff_edge_colors.items():
        edges = [(u, v) for u, v, d in G_diff.edges(data=True) if d.get("status") == status]
        if edges:
            widths = [max(0.3, abs(G_diff[u][v].get("weight", 0)) * 2.5)
                      for u, v in edges]
            nx.draw_networkx_edges(G_diff, pos, edgelist=edges,
                                   edge_color=color, width=widths, alpha=0.5, ax=ax2)

    diff_node_colors = [name_to_color[n] for n in G_diff.nodes()]
    diff_node_sizes_list = [change_size.get(n, _DIFF_NODE_BASE) for n in G_diff.nodes()]
    nx.draw_networkx_nodes(G_diff, pos, node_color=diff_node_colors,
                           node_size=diff_node_sizes_list, alpha=0.7,
                           edgecolors=_C_NODE_OUTLINE, linewidths=0.8, ax=ax2)
    if show_lab:
        nx.draw_networkx_labels(G_diff, pos, font_size=_LABEL_FS_DIFF,
                                font_weight="normal", font_family=_FONT_NET, ax=ax2)

    diff_handles = [
        mpatches.Patch(color=_C_DIFF_GAINED, label="Gained edge"),
        mpatches.Patch(color=_C_DIFF_LOST, label="Lost edge"),
        mpatches.Patch(color=_C_DIFF_SHARED, label="Unchanged edge"),
    ]
    leg_d1 = ax2.legend(handles=diff_handles, loc="lower right", ncol=1,
                        prop={"family": _FONT_NET, "size": _LEGEND_FS_NET}, frameon=True)
    ax2.legend(handles=module_handles, title="Node color", loc="lower left",
               ncol=1, prop={"family": _FONT_NET, "size": _LEGEND_FS_NET - 2},
               title_fontproperties={"family": _FONT_NET, "size": _LEGEND_FS_NET - 1},
               frameon=True)
    ax2.add_artist(leg_d1)
    ax2.set_title(
        f"Differential co-expression network\n{perturbation} vs control",
        fontsize=_TITLE_FS_NET,
        fontfamily=_FONT_NET,
    )
    ax2.axis("off")
    ax2.margins(0.18)
    fig2.tight_layout()
    fig2.savefig(plots_dir / f"genenet_{safe}_diff_network.png", dpi=120)
    plt.close(fig2)


# ---------------------------------------------------------------------------
# JSON export for interactive report
# ---------------------------------------------------------------------------

def _save_network_json(
    ctrl_corr: np.ndarray,
    pert_corr: np.ndarray,
    gene_names: List[str],
    cluster_labels: np.ndarray,
    ctrl_means: np.ndarray,
    pert_means: np.ndarray,
    perturbation: str,
    tables_dir: Path,
    corr_threshold: float,
    n_ctrl_cells: Optional[int] = None,
    n_pert_cells: Optional[int] = None,
) -> None:
    """Save network topology (nodes, edges, spring-layout positions) as JSON.

    Produces ``genenet_<safe_pert>_network.json`` in *tables_dir*.
    The interactive HTML report loads this file to render a Plotly network.
    Edges are filtered with the same sample-size aware threshold used by
    :func:`_plot_networks` so the JSON matches the rendered PNG.
    """
    import networkx as nx  # local import keeps module-level deps minimal

    G_ctrl = _build_network(ctrl_corr, gene_names, corr_threshold,
                            n_cells=n_ctrl_cells)
    G_pert = _build_network(pert_corr, gene_names, corr_threshold,
                            n_cells=n_pert_cells)

    G_union = nx.Graph()
    G_union.add_nodes_from(gene_names)
    G_union.add_edges_from(G_ctrl.edges())
    G_union.add_edges_from(G_pert.edges())
    pos = _spring_layout(G_union)

    log2fc = ((pert_means - ctrl_means) / np.log(2)).tolist()

    ctrl_edges = [
        {"s": u, "t": v,
         "r": round(float(G_ctrl[u][v]["weight"]), 4),
         "pos": bool(G_ctrl[u][v]["positive"])}
        for u, v in G_ctrl.edges()
    ]
    pert_edges = [
        {"s": u, "t": v,
         "r": round(float(G_pert[u][v]["weight"]), 4),
         "pos": bool(G_pert[u][v]["positive"])}
        for u, v in G_pert.edges()
    ]

    payload = {
        "perturbation": perturbation,
        "genes": gene_names,
        "cluster_ids": cluster_labels.tolist(),
        "ctrl_means": [round(float(m), 4) for m in ctrl_means],
        "pert_means": [round(float(m), 4) for m in pert_means],
        "log2fc": [round(float(f), 4) for f in log2fc],
        "positions": {
            g: [round(float(xy[0]), 4), round(float(xy[1]), 4)]
            for g, xy in pos.items()
        },
        "ctrl_edges": ctrl_edges,
        "pert_edges": pert_edges,
        "threshold": corr_threshold,
    }

    safe = perturbation.replace("/", "_").replace(" ", "_").replace("+", "_")
    out_path = tables_dir / f"genenet_{safe}_network.json"
    out_path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_gene_network(
    adata,
    output_dir: str,
    perturbations: Optional[List[str]] = None,
    n_top_genes: int = 50,
    n_gene_clusters: int = 5,
    corr_threshold: float = 0.5,
) -> pd.DataFrame:
    """Run gene co-expression network analysis for each perturbation vs control.

    Steps per perturbation:
      1. Select a perturbation-focused panel: top DEGs for that perturbation,
         filled with HVGs when needed.
      2. Compute Pearson gene–gene correlations in control cells and cluster
         them into ``n_gene_clusters`` modules via hierarchical clustering.
      3. Plot the control gene-cluster heatmap (once, shared across all perts).
      4. For each perturbation: side-by-side heatmap (control vs perturbed),
         co-expression network graphs (control + perturbed), and a differential
         network highlighting gained / lost / unchanged edges.
      5. Save gene cluster assignments to CSV.

    Args:
        adata            -- AnnData after QC + normalisation.
        output_dir       -- root output directory.
        perturbations    -- list of perturbation names to analyse.
                            If None, runs on all non-control perturbations
                            (use carefully for large datasets).
        n_top_genes      -- number of top HVGs to include in the networks.
        n_gene_clusters  -- number of gene modules identified in control cells.
        corr_threshold   -- |Pearson r| threshold for drawing an edge.

    Returns:
        DataFrame with columns [gene, cluster_id, mean_ctrl] indexed by gene.
    """
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    ctrl_mask = adata.obs["is_control"].values.astype(bool)
    if ctrl_mask.sum() < 5:
        return pd.DataFrame()

    # --- Baseline gene selection: top HVGs for the shared control overview ---
    if "highly_variable" in adata.var.columns:
        hv_idx = np.where(adata.var["highly_variable"].values)[0]
        hv_idx = hv_idx[: max(n_top_genes, len(hv_idx))]
    else:
        xc_full = _to_dense(adata.X[ctrl_mask, :])
        hv_idx = np.argsort(xc_full.var(axis=0))[::-1]

    base_hv_idx = hv_idx[:n_top_genes]
    gene_names: List[str] = adata.var_names[base_hv_idx].tolist()
    n_g = len(gene_names)
    if n_g < 4:
        return pd.DataFrame()

    # --- Control correlation matrix and gene clustering ---
    ctrl_expr = _to_dense(adata.X[ctrl_mask, :][:, base_hv_idx])
    ctrl_corr = _pearson_corrmat(ctrl_expr)
    n_clust = min(n_gene_clusters, n_g - 1)
    cluster_labels = _cluster_genes(ctrl_corr, n_clust)   # 1-indexed, length n_g
    sorted_idx = _sorted_order(cluster_labels)

    # --- Save gene cluster CSV ---
    ctrl_means = ctrl_expr.mean(axis=0)
    cluster_df = pd.DataFrame(
        {"gene": gene_names, "cluster_id": cluster_labels, "mean_ctrl": ctrl_means}
    )

    # --- Control heatmap (once) ---
    _plot_control_heatmap(ctrl_corr, gene_names, sorted_idx, cluster_labels, plots)

    # --- Determine which perturbations to analyse ---
    all_perts = [
        p for p in adata.obs["perturbation"].astype(str).unique()
        if str(p).lower() != "control"
    ]
    if perturbations is None:
        perturbations = all_perts
    else:
        perturbations = [p for p in perturbations if p in all_perts]

    if not perturbations:
        cluster_df.to_csv(tables / "genenet_gene_clusters.csv", index=False)
        return cluster_df

    n_ctrl_cells = int(ctrl_mask.sum())
    focused_rows = []
    for pert in perturbations:
        pert_mask = adata.obs["perturbation"].values == pert
        n_pert_cells = int(pert_mask.sum())
        if n_pert_cells < 5:
            continue

        gene_idx, selection_source = _deg_prioritized_gene_indices(
            adata, tables, pert, hv_idx, n_top_genes
        )
        if len(gene_idx) < 4:
            continue

        gene_names_p: List[str] = adata.var_names[gene_idx].tolist()
        ctrl_expr_p = _to_dense(adata.X[ctrl_mask, :][:, gene_idx])
        ctrl_corr_p = _pearson_corrmat(ctrl_expr_p)
        n_clust_p = min(n_gene_clusters, len(gene_names_p) - 1)
        cluster_labels_p = _cluster_genes(ctrl_corr_p, n_clust_p)
        sorted_idx_p = _sorted_order(cluster_labels_p)
        ctrl_means_p = ctrl_expr_p.mean(axis=0)

        pert_expr = _to_dense(adata.X[pert_mask, :][:, gene_idx])
        pert_corr = _pearson_corrmat(pert_expr)
        pert_means = pert_expr.mean(axis=0)

        _plot_comparison_heatmap(
            ctrl_expr_p, pert_expr,
            gene_names_p, sorted_idx_p, cluster_labels_p,
            pert, plots,
        )
        _plot_networks(
            ctrl_corr_p, pert_corr,
            gene_names_p, cluster_labels_p,
            ctrl_means_p, pert_means,
            pert, plots, corr_threshold,
            n_ctrl_cells=n_ctrl_cells,
            n_pert_cells=n_pert_cells,
        )

        _save_network_json(
            ctrl_corr_p, pert_corr,
            gene_names_p, cluster_labels_p,
            ctrl_means_p, pert_means,
            pert, tables, corr_threshold,
            n_ctrl_cells=n_ctrl_cells,
            n_pert_cells=n_pert_cells,
        )

        log2fc = (pert_means - ctrl_means_p) / np.log(2)
        for gene, cl, mc, mp, lf in zip(
            gene_names_p, cluster_labels_p, ctrl_means_p, pert_means, log2fc
        ):
            focused_rows.append({
                "perturbation": pert,
                "gene": gene,
                "cluster_id": int(cl),
                "mean_ctrl": float(mc),
                "mean_pert": float(mp),
                "log2fc": float(lf),
                "selection_source": selection_source.get(gene, "HVG"),
            })

    cluster_df.to_csv(tables / "genenet_gene_clusters.csv", index=False)
    if focused_rows:
        pd.DataFrame(focused_rows).to_csv(
            tables / "genenet_focused_gene_clusters.csv", index=False
        )
    return cluster_df


# =============================================================================
# TF-centric gene regulatory network
# =============================================================================
# This is a sibling of run_gene_network that builds the panel from
#   1. the top `n_top_tfs` highly-variable transcription factors (variance
#      computed over control cells, intersected with a HumanTFs symbol list),
#   2. plus the `n_partners_per_tf` most-correlated genes per TF in control
#      cells (|Pearson r| ≥ corr_threshold), drawn from a candidate pool of
#      top-variance genes.
# Every other piece (correlation, edge predicate, layout, differential edges)
# is identical to run_gene_network. TFs are rendered as squares, partners as
# circles; module colour is shared with the generic gene-network panels.
# =============================================================================

# Default source for the TF symbol list (HumanTFs v1.01 master list).
_DEFAULT_TF_LIST_PATH = (
    "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/"
    "PerturbVerse_v4/biological_priors/grn/humantfs_v1.01.csv"
)


def _load_tf_set(path: Optional[str]) -> set:
    """Read a one-column CSV (`tf_symbol`) or one-symbol-per-line text file."""
    p = Path(path) if path else Path(_DEFAULT_TF_LIST_PATH)
    if not p.exists():
        return set()
    try:
        if p.suffix.lower() in {".csv", ".tsv"}:
            df = pd.read_csv(p)
            col = "tf_symbol" if "tf_symbol" in df.columns else df.columns[0]
            return set(df[col].astype(str).str.strip().tolist())
        return set(line.strip() for line in p.read_text().splitlines() if line.strip())
    except Exception:
        return set()


def _zscored_columns(x: np.ndarray) -> np.ndarray:
    """Return per-column z-score (cells × genes); zeros where std == 0."""
    mu = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, keepdims=True)
    safe = np.where(sd < 1e-8, 1.0, sd)
    return ((x - mu) / safe).astype(np.float32, copy=False)


def _select_partners(
    z_matrices: List[np.ndarray],
    tf_local_idx: np.ndarray,
    candidate_local_idx: np.ndarray,
    n_partners: int,
    corr_threshold: float,
) -> dict:
    """Per-TF, return the top-``n_partners`` candidate genes ranked by the
    **maximum** ``|r|`` over a list of conditions.

    A partner is kept if it passes the threshold in *any* of the provided
    z-scored matrices (e.g. control or one of the analysed perturbations).
    Self-correlation is excluded. ``z_matrices`` is a list of per-column
    z-scored expression matrices ``(n_cells × n_pool_genes)`` — one per
    condition — all sharing the same gene-axis ordering as
    ``candidate_local_idx`` and ``tf_local_idx``.

    Returns ``{tf_local_idx: [partner_local_idx, ...]}`` sorted by
    ``max(|r|)`` descending, capped at ``n_partners`` per TF.
    """
    if not z_matrices:
        return {int(t): [] for t in tf_local_idx}

    cand_local_arr = np.asarray(candidate_local_idx)
    K_tf = len(tf_local_idx)
    K_cand = len(cand_local_arr)
    abs_max = np.zeros((K_cand, K_tf), dtype=np.float32)

    for z in z_matrices:
        if z is None or z.shape[0] < 2:
            continue
        tf_z = z[:, tf_local_idx]               # (n_cells, K_tf)
        cand_z = z[:, cand_local_arr]           # (n_cells, K_cand)
        corr = (cand_z.T @ tf_z) / float(z.shape[0])   # (K_cand, K_tf)
        np.maximum(abs_max, np.abs(corr), out=abs_max)

    out = {}
    for k, tf_loc in enumerate(tf_local_idx):
        rs = abs_max[:, k]
        keep = cand_local_arr != tf_loc        # exclude self
        rs_eff = np.where(keep, rs, 0.0)
        order = np.argsort(-rs_eff)
        chosen = []
        for j in order:
            if not keep[j]:
                continue
            if rs_eff[j] < corr_threshold:
                break
            chosen.append(int(cand_local_arr[j]))
            if len(chosen) >= n_partners:
                break
        out[int(tf_loc)] = chosen
    return out


def _draw_tf_network(
    G,
    pos: dict,
    ax,
    tf_set: set,
    name_to_color: dict,
    node_sizes: dict,
    title: str,
    show_labels: bool = True,
    label_font_size: int = 22,
    default_node_size: float = 800,
    title_font_size: int = 24,
    legend_font_size: int = 18,
    module_handles: Optional[List] = None,
) -> None:
    """Draw a TF-centric network. TFs are squares, partners are circles."""
    import networkx as nx

    edge_colors = [
        _C_POS_R if G[u][v].get("positive", True) else _C_NEG_R
        for u, v in G.edges()
    ]
    edge_widths = [
        max(0.3, abs(G[u][v].get("weight", 0)) * 2.5) for u, v in G.edges()
    ]
    nx.draw_networkx_edges(
        G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.4, ax=ax
    )

    tf_nodes = [n for n in G.nodes() if n in tf_set]
    pt_nodes = [n for n in G.nodes() if n not in tf_set]

    if pt_nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=pt_nodes,
            node_color=[name_to_color.get(n, "#9aa0a6") for n in pt_nodes],
            node_size=[node_sizes.get(n, default_node_size) for n in pt_nodes],
            node_shape="o", alpha=0.7,
            edgecolors="#1d3557", linewidths=0.8, ax=ax,
        )
    if tf_nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=tf_nodes,
            node_color=[name_to_color.get(n, "#9aa0a6") for n in tf_nodes],
            node_size=[node_sizes.get(n, default_node_size) * 1.4 for n in tf_nodes],
            node_shape="s", alpha=0.85,
            edgecolors="#000000", linewidths=1.4, ax=ax,
        )

    if show_labels and len(G.nodes()) <= 120:
        # Bigger / bolder labels for TFs.
        if pt_nodes:
            nx.draw_networkx_labels(
                G, pos,
                labels={n: n for n in pt_nodes},
                font_size=label_font_size,
                font_color="#1d3557",
                ax=ax,
            )
        if tf_nodes:
            nx.draw_networkx_labels(
                G, pos,
                labels={n: n for n in tf_nodes},
                font_size=label_font_size + 2,
                font_color="#000000",
                font_weight="bold",
                ax=ax,
            )

    edge_legend = [
        mpatches.Patch(color=_C_POS_R, label="Positive correlation"),
        mpatches.Patch(color=_C_NEG_R, label="Negative correlation"),
    ]
    shape_legend = [
        plt.Line2D([0], [0], marker="s", color="w", label="TF (square)",
                   markerfacecolor="#9aa0a6", markeredgecolor="#000000",
                   markersize=12, markeredgewidth=2),
        plt.Line2D([0], [0], marker="o", color="w", label="Partner (circle)",
                   markerfacecolor="#9aa0a6", markeredgecolor="#1d3557",
                   markersize=12, markeredgewidth=1),
    ]
    leg1 = ax.legend(
        handles=edge_legend, loc="lower left",
        fontsize=legend_font_size, frameon=True,
    )
    ax.add_artist(leg1)
    leg2 = ax.legend(
        handles=shape_legend, loc="upper right",
        fontsize=legend_font_size, frameon=True,
    )
    ax.add_artist(leg2)
    if module_handles:
        ax.legend(
            handles=module_handles, title="Node color (gene module)",
            loc="lower right",
            fontsize=max(10, legend_font_size - 4),
            title_fontsize=max(11, legend_font_size - 3),
            frameon=True, ncol=1,
        )
    ax.set_title(title, fontsize=title_font_size)
    ax.axis("off")


def _draw_tf_diff_network(
    G_diff, pos, ax, tf_set, name_to_color, node_sizes, title,
    show_labels: bool = True, label_font_size: int = 20,
    default_node_size: float = 500,
    title_font_size: int = 24,
    legend_font_size: int = 18,
    module_handles: Optional[List] = None,
) -> None:
    """Differential TF-centric network. Same node encoding as _draw_tf_network
    but edge colour reflects gained / lost / shared status."""
    import networkx as nx

    diff_edge_colors = {
        "shared": _C_DIFF_SHARED,
        "gained": _C_DIFF_GAINED,
        "lost":   _C_DIFF_LOST,
    }
    for status, color in diff_edge_colors.items():
        edges = [(u, v) for u, v, d in G_diff.edges(data=True)
                 if d.get("status") == status]
        if edges:
            widths = [max(0.3, abs(G_diff[u][v].get("weight", 0)) * 2.5)
                      for u, v in edges]
            nx.draw_networkx_edges(G_diff, pos, edgelist=edges,
                                   edge_color=color, width=widths, alpha=0.5, ax=ax)

    tf_nodes = [n for n in G_diff.nodes() if n in tf_set]
    pt_nodes = [n for n in G_diff.nodes() if n not in tf_set]
    if pt_nodes:
        nx.draw_networkx_nodes(
            G_diff, pos, nodelist=pt_nodes,
            node_color=[name_to_color.get(n, "#9aa0a6") for n in pt_nodes],
            node_size=[node_sizes.get(n, default_node_size) for n in pt_nodes],
            node_shape="o", alpha=0.7,
            edgecolors="#1d3557", linewidths=0.8, ax=ax,
        )
    if tf_nodes:
        nx.draw_networkx_nodes(
            G_diff, pos, nodelist=tf_nodes,
            node_color=[name_to_color.get(n, "#9aa0a6") for n in tf_nodes],
            node_size=[node_sizes.get(n, default_node_size) * 1.4 for n in tf_nodes],
            node_shape="s", alpha=0.85,
            edgecolors="#000000", linewidths=1.4, ax=ax,
        )

    if show_labels and len(G_diff.nodes()) <= 120:
        if pt_nodes:
            nx.draw_networkx_labels(
                G_diff, pos, labels={n: n for n in pt_nodes},
                font_size=label_font_size, font_color="#1d3557", ax=ax,
            )
        if tf_nodes:
            nx.draw_networkx_labels(
                G_diff, pos, labels={n: n for n in tf_nodes},
                font_size=label_font_size + 2, font_color="#000000",
                font_weight="bold", ax=ax,
            )

    diff_handles = [
        mpatches.Patch(color=_C_DIFF_GAINED, label="Gained edge"),
        mpatches.Patch(color=_C_DIFF_LOST, label="Lost edge"),
        mpatches.Patch(color=_C_DIFF_SHARED, label="Unchanged edge"),
    ]
    leg_e = ax.legend(handles=diff_handles, loc="lower left", ncol=1,
                      fontsize=legend_font_size, frameon=True)
    ax.add_artist(leg_e)
    if module_handles:
        ax.legend(
            handles=module_handles, title="Node color (gene module)",
            loc="lower right",
            fontsize=max(10, legend_font_size - 4),
            title_fontsize=max(11, legend_font_size - 3),
            frameon=True, ncol=1,
        )
    ax.set_title(title, fontsize=title_font_size)
    ax.axis("off")


def run_tf_gene_network(
    adata,
    output_dir: str,
    perturbations: Optional[List[str]] = None,
    n_top_tfs: int = 10,
    n_partners_per_tf: int = 8,
    corr_threshold: float = 0.5,
    tf_list_path: Optional[str] = None,
    candidate_gene_pool_size: int = 2000,
    n_gene_modules: int = 5,
    max_ctrl_cells: int = 20000,
) -> pd.DataFrame:
    """Build a TF-anchored co-expression network for control + each perturbation.

    Pipeline:
      1. Load HumanTFs symbol list; intersect with ``adata.var_names`` →
         available TFs in this dataset.
      2. Rank available TFs by their variance over control cells; take the
         top ``n_top_tfs``.
      3. For each TF, find the ``n_partners_per_tf`` genes with the
         highest ``|Pearson r|`` to the TF in control cells (drawn from a
         candidate pool = top ``candidate_gene_pool_size`` highly-variable
         genes), keeping only those with ``|r| ≥ corr_threshold``.
      4. Panel = union(top TFs, all selected partners).
      5. Build a thresholded gene–gene correlation graph for the control pool
         and for each perturbation (same hard ``|r| ≥ corr_threshold`` rule
         used by ``run_gene_network``).
      6. Spring layout on the union graph keeps node positions consistent
         across panels.
      7. Render one shared control PNG and, per perturbation, a perturbation
         PNG and a differential PNG. TFs are drawn as squares, partners as
         circles; node colour reflects the gene module computed from the
         control |r| matrix on the panel.

    Outputs (under ``<output_dir>/``):
      plots/
        tfnet_control_network.png             square 10×10 control network
        tfnet_<pert>_pert_network.png         perturbation network (same panel)
        tfnet_<pert>_diff_network.png         differential network
      csv/
        tfnet_seed_tfs.csv                    seed TFs + variance + mean expr
        tfnet_partners.csv                    partner gene → top-TF and |r|
        tfnet_<pert>_network.json             topology + positions for the
                                              interactive report

    Returns the seed-TF table (or an empty DataFrame on failure).
    """
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    if "is_control" not in adata.obs.columns:
        return pd.DataFrame()
    ctrl_mask = adata.obs["is_control"].values.astype(bool)
    if int(ctrl_mask.sum()) < 5:
        return pd.DataFrame()

    tf_set = _load_tf_set(tf_list_path)
    if not tf_set:
        return pd.DataFrame()

    var_names = np.asarray(adata.var_names)
    tf_in_data = np.array([str(g) in tf_set for g in var_names], dtype=bool)
    if not tf_in_data.any():
        return pd.DataFrame()

    # ---- Sub-sample control cells if there are very many (Pearson r is
    #      already accurate at ~10k cells; this caps memory/time).
    ctrl_idx_full = np.where(ctrl_mask)[0]
    if ctrl_idx_full.size > max_ctrl_cells:
        rng = np.random.default_rng(0)
        ctrl_idx_full = np.sort(
            rng.choice(ctrl_idx_full, size=max_ctrl_cells, replace=False)
        )

    # ---- Step 1: candidate gene pool — top-variance genes in control. ----
    if "highly_variable" in adata.var.columns:
        hv_pool = np.where(adata.var["highly_variable"].values)[0]
        if hv_pool.size < 50:                       # too few HVGs annotated
            hv_pool = None
    else:
        hv_pool = None
    if hv_pool is None:
        sub_x = _to_dense(adata.X[ctrl_idx_full, :]).astype(np.float32, copy=False)
        gene_var = sub_x.var(axis=0)
        hv_pool = np.argsort(gene_var)[::-1][:max(candidate_gene_pool_size, 200)]
        ctrl_pool_x = sub_x[:, hv_pool]
    else:
        if hv_pool.size > candidate_gene_pool_size:
            sub_x_full = _to_dense(adata.X[ctrl_idx_full, :][:, hv_pool])\
                .astype(np.float32, copy=False)
            local_var = sub_x_full.var(axis=0)
            keep = np.argsort(local_var)[::-1][:candidate_gene_pool_size]
            hv_pool = hv_pool[keep]
            ctrl_pool_x = sub_x_full[:, keep]
        else:
            ctrl_pool_x = _to_dense(adata.X[ctrl_idx_full, :][:, hv_pool])\
                .astype(np.float32, copy=False)
    pool_gene_names = var_names[hv_pool].tolist()

    # ---- Step 2: top variable TFs, restricted to TFs that survived the pool. ----
    pool_is_tf = np.array([g in tf_set for g in pool_gene_names], dtype=bool)
    if not pool_is_tf.any():
        # If no top-variance TF made it into the candidate pool, force the
        # top variable TFs (anywhere in adata) into the pool.
        tf_global_idx = np.where(tf_in_data)[0]
        sub_x_tfs = _to_dense(adata.X[ctrl_idx_full, :][:, tf_global_idx])\
            .astype(np.float32, copy=False)
        tf_var_global = sub_x_tfs.var(axis=0)
        keep = np.argsort(tf_var_global)[::-1][:n_top_tfs]
        forced_tf_idx = tf_global_idx[keep]
        merged = np.unique(np.concatenate([hv_pool, forced_tf_idx]))
        ctrl_pool_x = _to_dense(adata.X[ctrl_idx_full, :][:, merged])\
            .astype(np.float32, copy=False)
        hv_pool = merged
        pool_gene_names = var_names[hv_pool].tolist()
        pool_is_tf = np.array([g in tf_set for g in pool_gene_names], dtype=bool)

    pool_tf_local = np.where(pool_is_tf)[0]
    pool_tf_var = ctrl_pool_x[:, pool_tf_local].var(axis=0)
    take = min(n_top_tfs, pool_tf_local.size)
    top_tf_local_idx = pool_tf_local[np.argsort(pool_tf_var)[::-1][:take]]
    top_tf_names = [pool_gene_names[i] for i in top_tf_local_idx]

    # ---- Step 3: per-TF top correlated partners. A partner is kept if it
    #             passes |r| >= corr_threshold in ANY analysed condition
    #             (control or one of the perturbations under study), so genes
    #             that only co-vary with a TF inside a perturbation are not
    #             missed. Pre-screen perturbations to those that will be
    #             rendered (≥ 5 cells), to mirror the panel that the rest of
    #             the function builds.
    pool_z_ctrl = _zscored_columns(ctrl_pool_x)
    cand_local_idx = np.arange(len(pool_gene_names))

    _all_perts = [
        p for p in adata.obs["perturbation"].astype(str).unique()
        if str(p).lower() != "control"
    ]
    if perturbations is None:
        _perts_for_partner_scan = _all_perts
    else:
        _perts_for_partner_scan = [p for p in perturbations if p in _all_perts]

    z_matrices = [pool_z_ctrl]
    for _p in _perts_for_partner_scan:
        _pmask = adata.obs["perturbation"].values == _p
        if int(_pmask.sum()) < 5:
            continue
        _ppool = _to_dense(adata.X[_pmask, :][:, hv_pool])\
            .astype(np.float32, copy=False)
        z_matrices.append(_zscored_columns(_ppool))

    partners_by_tf = _select_partners(
        z_matrices, top_tf_local_idx, cand_local_idx,
        n_partners=n_partners_per_tf, corr_threshold=corr_threshold,
    )

    # ---- Step 4: panel = TFs ∪ partners. ----
    panel_local_set = set(int(i) for i in top_tf_local_idx)
    for pl in partners_by_tf.values():
        panel_local_set.update(int(p) for p in pl)
    panel_local = np.array(sorted(panel_local_set), dtype=int)
    if panel_local.size < 4:
        return pd.DataFrame()
    panel_global = hv_pool[panel_local]
    panel_names = [pool_gene_names[i] for i in panel_local]
    n_panel = len(panel_names)
    panel_is_tf = np.array([g in tf_set for g in panel_names], dtype=bool)

    # ---- Step 5: control + per-pert correlation matrices on the panel. ----
    ctrl_panel_expr = _to_dense(adata.X[ctrl_idx_full, :][:, panel_global])\
        .astype(np.float32, copy=False)
    ctrl_corr = _pearson_corrmat(ctrl_panel_expr)
    ctrl_means = ctrl_panel_expr.mean(axis=0)

    # Gene module assignment (shared across all panels of a perturbation).
    n_clust = max(2, min(n_gene_modules, n_panel - 1))
    cluster_labels = _cluster_genes(ctrl_corr, n_clust)
    n_modules_real = int(cluster_labels.max())
    palette = sns.color_palette("Set2", min(n_modules_real, 8))
    name_to_color = {
        g: palette[(cluster_labels[i] - 1) % len(palette)]
        for i, g in enumerate(panel_names)
    }
    # Legend handles for the gene-module palette — passed to the drawing
    # helpers so every TF panel shows what each colour means.
    module_handles = [
        mpatches.Patch(color=palette[(c - 1) % len(palette)],
                       label=f"Module {c}")
        for c in range(1, n_modules_real + 1)
    ]

    # Build the control graph (sample-size aware threshold).
    n_ctrl_used = int(ctrl_idx_full.size)
    G_ctrl = _build_network(ctrl_corr, panel_names, corr_threshold,
                            n_cells=n_ctrl_used)
    eff_ctrl = G_ctrl.graph.get("effective_threshold", corr_threshold)

    # Pre-compute every per-pert graph and its means, so we can spring-layout
    # on the union of all edges (control + every perturbation).
    all_perts = [
        p for p in adata.obs["perturbation"].astype(str).unique()
        if str(p).lower() != "control"
    ]
    if perturbations is None:
        perts_to_run = all_perts
    else:
        perts_to_run = [p for p in perturbations if p in all_perts]

    pert_data = {}            # pert -> dict(G_pert, means, expr, n_cells)
    for pert in perts_to_run:
        pmask = adata.obs["perturbation"].values == pert
        n_pert_cells = int(pmask.sum())
        if n_pert_cells < 5:
            continue
        pexpr = _to_dense(adata.X[pmask, :][:, panel_global])\
            .astype(np.float32, copy=False)
        pcorr = _pearson_corrmat(pexpr)
        pert_data[pert] = {
            "G": _build_network(pcorr, panel_names, corr_threshold,
                                n_cells=n_pert_cells),
            "means": pexpr.mean(axis=0),
            "corr": pcorr,
            "n_cells": n_pert_cells,
        }

    import networkx as nx
    G_union = nx.Graph()
    G_union.add_nodes_from(panel_names)
    G_union.add_edges_from(G_ctrl.edges())
    for d in pert_data.values():
        G_union.add_edges_from(d["G"].edges())
    pos = _spring_layout(G_union)

    # ---- Step 6: write the seed-TF and partner CSVs. ----
    seed_rows = []
    for i_loc in top_tf_local_idx:
        gname = pool_gene_names[i_loc]
        v = float(ctrl_pool_x[:, i_loc].var())
        m = float(ctrl_pool_x[:, i_loc].mean())
        seed_rows.append({"tf": gname, "var_ctrl": v, "mean_ctrl": m})
    seed_df = pd.DataFrame(seed_rows)
    seed_df.to_csv(tables / "tfnet_seed_tfs.csv", index=False)

    partner_rows = []
    z_full = pool_z_ctrl
    for tf_loc, partner_locs in partners_by_tf.items():
        tf_name = pool_gene_names[tf_loc]
        tf_z_ctrl = z_full[:, tf_loc]
        for p_loc in partner_locs:
            partner_z_ctrl = z_full[:, p_loc]
            r_ctrl = float((partner_z_ctrl * tf_z_ctrl).mean())
            # Also report max |r| across analysed conditions (the criterion
            # actually used to pick the partner).
            r_max = abs(r_ctrl)
            for z_other in z_matrices[1:]:
                r_p = float((z_other[:, p_loc] * z_other[:, tf_loc]).mean())
                if abs(r_p) > r_max:
                    r_max = abs(r_p)
            partner_rows.append({
                "tf": tf_name,
                "partner": pool_gene_names[p_loc],
                "pearson_r_ctrl": r_ctrl,
                "max_abs_r_any_condition": r_max,
            })
    pd.DataFrame(partner_rows).to_csv(tables / "tfnet_partners.csv", index=False)

    # ---- Step 7: render control + per-pert + per-pert-diff PNGs. ----
    # Same encoding scheme as run_gene_network:
    #   node size = degree (number of significant connections),
    #   line width ∝ |r|,
    #   node colour = gene module (shared across all 3 panels of a pert),
    #   TFs are drawn as squares (handled inside _draw_tf_network).
    _STATIC_NODE_BASE = 200
    _STATIC_NODE_SCALE = 80
    _STATIC_NODE_MAX = 1800
    _DIFF_NODE_BASE = 160
    _DIFF_NODE_SCALE = 65
    _DIFF_NODE_MAX = 1500
    _PANEL_SIZE = 12.0
    _LBL_FS = 22       # static panels (TFs get +2 internally)
    _LBL_FS_DIFF = 20  # diff panel
    _TITLE_FS_TF = 24
    _LEGEND_FS_TF = 18

    # Control plot (one shared file).
    ctrl_deg = dict(G_ctrl.degree())
    ctrl_size = {
        g: min(_STATIC_NODE_MAX,
               _STATIC_NODE_BASE + ctrl_deg.get(g, 0) * _STATIC_NODE_SCALE)
        for g in panel_names
    }
    fig_c, ax_c = plt.subplots(figsize=(_PANEL_SIZE, _PANEL_SIZE))
    _draw_tf_network(
        G_ctrl, pos, ax_c,
        tf_set=set(top_tf_names),
        name_to_color=name_to_color,
        node_sizes=ctrl_size,
        title=f"Control — TF-anchored network "
              f"(top {len(top_tf_names)} variable TFs + {n_panel - len(top_tf_names)} partners, "
              f"|r| ≥ {eff_ctrl:.2f}, n={n_ctrl_used} cells)",
        show_labels=True,
        label_font_size=_LBL_FS,
        default_node_size=_STATIC_NODE_BASE,
        title_font_size=_TITLE_FS_TF,
        legend_font_size=_LEGEND_FS_TF,
        module_handles=module_handles,
    )
    fig_c.tight_layout()
    fig_c.savefig(plots / "tfnet_control_network.png", dpi=120, bbox_inches="tight")
    plt.close(fig_c)

    # Per-pert plots + JSON sidecar.
    for pert, d in pert_data.items():
        safe = pert.replace("/", "_").replace(" ", "_").replace("+", "_")
        pmeans = d["means"]
        eff_pert = d["G"].graph.get("effective_threshold", corr_threshold)

        # Differential graph (built first so its degree drives diff node sizes).
        ctrl_e = set(tuple(sorted((u, v))) for u, v in G_ctrl.edges())
        pert_e = set(tuple(sorted((u, v))) for u, v in d["G"].edges())
        G_diff = nx.Graph()
        G_diff.add_nodes_from(panel_names)
        for u, v in ctrl_e | pert_e:
            rc = abs(G_ctrl[u][v]["weight"]) if (u, v) in ctrl_e else 0.0
            rp = abs(d["G"][u][v]["weight"]) if (u, v) in pert_e else 0.0
            if (u, v) in ctrl_e and (u, v) in pert_e:
                status = "shared"; w = max(rc, rp)
            elif (u, v) in ctrl_e:
                status = "lost";   w = rc
            else:
                status = "gained"; w = rp
            G_diff.add_edge(u, v, status=status, weight=float(w))

        pert_deg = dict(d["G"].degree())
        diff_deg = dict(G_diff.degree())
        pert_size = {
            g: min(_STATIC_NODE_MAX,
                   _STATIC_NODE_BASE + pert_deg.get(g, 0) * _STATIC_NODE_SCALE)
            for g in panel_names
        }
        change_size = {
            g: min(_DIFF_NODE_MAX,
                   _DIFF_NODE_BASE + diff_deg.get(g, 0) * _DIFF_NODE_SCALE)
            for g in panel_names
        }

        fig_p, ax_p = plt.subplots(figsize=(_PANEL_SIZE, _PANEL_SIZE))
        _draw_tf_network(
            d["G"], pos, ax_p,
            tf_set=set(top_tf_names),
            name_to_color=name_to_color,
            node_sizes=pert_size,
            title=f"{pert} — TF-anchored network "
                  f"(|r| ≥ {eff_pert:.2f}, n={d['n_cells']} cells)",
            show_labels=True,
            label_font_size=_LBL_FS,
            default_node_size=_STATIC_NODE_BASE,
            title_font_size=_TITLE_FS_TF,
            legend_font_size=_LEGEND_FS_TF,
            module_handles=module_handles,
        )
        fig_p.tight_layout()
        fig_p.savefig(plots / f"tfnet_{safe}_pert_network.png",
                      dpi=120, bbox_inches="tight")
        plt.close(fig_p)

        fig_d, ax_d = plt.subplots(figsize=(_PANEL_SIZE, _PANEL_SIZE))
        _draw_tf_diff_network(
            G_diff, pos, ax_d,
            tf_set=set(top_tf_names),
            name_to_color=name_to_color,
            node_sizes=change_size,
            title=f"Differential TF-anchored network: {pert} vs control",
            show_labels=True,
            label_font_size=_LBL_FS_DIFF,
            default_node_size=_DIFF_NODE_BASE,
            title_font_size=_TITLE_FS_TF,
            legend_font_size=_LEGEND_FS_TF,
            module_handles=module_handles,
        )
        fig_d.tight_layout()
        fig_d.savefig(plots / f"tfnet_{safe}_diff_network.png",
                      dpi=120, bbox_inches="tight")
        plt.close(fig_d)

        # JSON sidecar — same schema as genenet_*_network.json plus an
        # `is_tf` flag per node, so the interactive report can style TFs.
        log2fc = ((pmeans - ctrl_means) / np.log(2)).tolist()
        ctrl_edges = [
            {"s": u, "t": v,
             "r": round(float(G_ctrl[u][v]["weight"]), 4),
             "pos": bool(G_ctrl[u][v]["positive"])}
            for u, v in G_ctrl.edges()
        ]
        pert_edges_list = [
            {"s": u, "t": v,
             "r": round(float(d["G"][u][v]["weight"]), 4),
             "pos": bool(d["G"][u][v]["positive"])}
            for u, v in d["G"].edges()
        ]
        payload = {
            "perturbation": pert,
            "genes": panel_names,
            "is_tf": [bool(b) for b in panel_is_tf],
            "cluster_ids": cluster_labels.tolist(),
            "ctrl_means": [round(float(m), 4) for m in ctrl_means],
            "pert_means": [round(float(m), 4) for m in pmeans],
            "log2fc": [round(float(f), 4) for f in log2fc],
            "positions": {
                g: [round(float(xy[0]), 4), round(float(xy[1]), 4)]
                for g, xy in pos.items()
            },
            "ctrl_edges": ctrl_edges,
            "pert_edges": pert_edges_list,
            "threshold": corr_threshold,
        }
        (tables / f"tfnet_{safe}_network.json").write_text(
            json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
        )

    return seed_df
