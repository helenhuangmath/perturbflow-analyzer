# =============================================================================
# perturbflow/analyzer/interactive.py  (v2)
#
# Self-contained interactive HTML report for Perturb-seq results.
#
# Design:
#   7 tabs: Home | QC | UMAP Explorer | Heatmaps | Perturbation |
#           Gene Expression | Parameters
#   * Scales to thousands of perturbations: top-N dropdown + "Show all".
#   * Home: dataset summary, pipeline steps, top-perturbation bar.
#   * QC: violins for expression metrics + cells-per-perturbation bar.
#   * UMAP Explorer: color by perturbation / cell_state / score / gene.
#   * Heatmaps: base64-embedded PNG images from EDA step.
#   * Perturbation: volcano, UMAP highlight, DEG bar, paginated table,
#     pathway enrichment bar.
#   * Gene Expression: box plots per perturbation + log2FC bar.
#   * Parameters: config table with descriptions.
#   * All charts Plotly.js, bundled next to the HTML for offline cluster use.
# =============================================================================

from __future__ import annotations

import base64
import json
import shutil
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from .utils import ensure_dir


def _gsea_running_curve(ranked_genes: list[str], scores: np.ndarray, gene_set: set[str]) -> tuple[float, list[int], list[float], list[int]]:
    hits = np.array([g in gene_set for g in ranked_genes], dtype=bool)
    n_hit = int(hits.sum())
    n_miss = len(hits) - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0, [], [], []
    weights = np.abs(scores)
    hit_weights = weights * hits
    hit_denom = float(hit_weights.sum()) or 1.0
    running = np.cumsum(np.where(hits, hit_weights / hit_denom, -1.0 / n_miss))
    max_es = float(running.max())
    min_es = float(running.min())
    es = max_es if abs(max_es) >= abs(min_es) else min_es
    hit_idx = np.where(hits)[0].astype(int).tolist()
    # Downsample the running curve for compact HTML while preserving shape.
    if len(running) > 250:
        keep = np.unique(np.linspace(0, len(running) - 1, 250).astype(int))
        curve = [round(float(running[i]), 4) for i in keep]
        curve_x = keep.astype(int).tolist()
    else:
        curve = [round(float(v), 4) for v in running]
        curve_x = list(range(len(running)))
    return es, curve_x, curve, hit_idx[:150]


def _compute_gsea_rows(deg_df: pd.DataFrame, top_n: int = 15) -> list[dict]:
    """Small report-side preranked GSEA summary using embedded gene sets.

    This is intentionally lightweight: it reports a weighted running-sum ES and
    leading-edge genes for quick interpretation in the browser. The heavier
    permutation/FDR workflow can be added later as a pipeline step.
    """
    try:
        from .pathways import GENE_SETS
    except Exception:
        return []
    if deg_df.empty or "gene" not in deg_df.columns or "log2fc" not in deg_df.columns:
        return []
    df = deg_df[["gene", "log2fc"]].dropna().copy()
    if df.empty:
        return []
    df["gene"] = df["gene"].astype(str)
    df["score"] = pd.to_numeric(df["log2fc"], errors="coerce").fillna(0.0)
    df = df.sort_values("score", ascending=False)
    ranked = df["gene"].tolist()
    scores = df["score"].to_numpy(dtype=float)
    universe = set(ranked)
    rows = []
    for term, genes in GENE_SETS.items():
        members = set(genes) & universe
        if len(members) < 5:
            continue
        es, curve_x, curve, hit_idx = _gsea_running_curve(ranked, scores, members)
        if es == 0:
            continue
        if es > 0:
            ordered = [g for g in ranked if g in members][:10]
        else:
            ordered = [g for g in reversed(ranked) if g in members][:10]
        rows.append({
            "term": term,
            "es": round(float(es), 4),
            "direction": "up" if es > 0 else "down",
            "n_genes": int(len(members)),
            "leading_edge": "|".join(ordered),
            "curve_x": curve_x,
            "curve": curve,
            "hit_idx": hit_idx,
        })
    rows.sort(key=lambda r: abs(r["es"]), reverse=True)
    return rows[:top_n]


def _to_dense(x):
    return x.A if hasattr(x, "A") else np.asarray(x)


def _safe_float(v):
    if v is None:
        return None
    try:
        f = float(v)
        return None if (np.isnan(f) or np.isinf(f)) else round(f, 6)
    except Exception:
        return None


def _clean_record(row: dict) -> dict:
    out = {}
    for k, v in row.items():
        if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
            out[k] = None
        elif isinstance(v, (np.integer,)):
            out[k] = int(v)
        elif isinstance(v, (np.floating,)):
            out[k] = _safe_float(v)
        elif isinstance(v, bool):
            out[k] = bool(v)
        else:
            out[k] = v
    return out


# ---------------------------------------------------------------------------
# Data extraction
# ---------------------------------------------------------------------------

def _extract_data(adata, csv_dir: Path, n_top_genes: int, max_cells_per_group: int) -> dict:
    data: dict = {}

    n_pert = int(adata.obs["perturbation"].nunique()) if "perturbation" in adata.obs.columns else 0
    data["summary"] = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_perturbations": n_pert,
    }

    if "perturbation" in adata.obs.columns:
        cpp_raw = adata.obs["perturbation"].astype(str).value_counts()
        all_perts = cpp_raw.index.tolist()
        if "control" in all_perts:
            all_perts = ["control"] + [p for p in all_perts if p != "control"]
        data["perturbations"] = all_perts
    else:
        data["perturbations"] = []

    cpp_path = csv_dir / "eda_cells_per_perturbation.csv"
    if cpp_path.exists():
        cpp_df = pd.read_csv(cpp_path)
    else:
        cpp_df = (
            adata.obs["perturbation"].astype(str)
            .value_counts().rename_axis("perturbation").reset_index(name="n_cells")
        )
    cpp_df["is_control"] = cpp_df["perturbation"].astype(str) == "control"
    data["cells_per_pert"] = [_clean_record(r) for r in cpp_df.to_dict(orient="records")]

    cell_path = csv_dir / "cell_level_summary.csv"
    qc_cols = ["perturbation_score", "guide_confidence_score",
               "target_expr_reduction", "perturbation_burden",
               "pseudotime", "escape_probability"]
    if cell_path.exists():
        cell_df = pd.read_csv(cell_path)
        for col in qc_cols:
            if col not in cell_df.columns:
                cell_df[col] = np.nan
    else:
        cell_df = adata.obs.reset_index(drop=True)
        for col in qc_cols:
            if col not in cell_df.columns:
                cell_df[col] = np.nan
    for qc_raw in ["n_genes_by_counts", "total_counts", "pct_counts_mt", "n_genes",
                    "percent_mito", "UMI_count", "core_adjusted_UMI_count",
                    "z_gemgroup_UMI", "core_scale_factor"]:
        if qc_raw not in cell_df.columns and qc_raw in adata.obs.columns:
            cell_df[qc_raw] = adata.obs[qc_raw].values
    if len(cell_df) > 3000:
        rng = np.random.default_rng(42)
        cell_df = cell_df.iloc[rng.choice(len(cell_df), 3000, replace=False)]
    qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "n_genes",
                  "percent_mito", "UMI_count", "core_adjusted_UMI_count",
                  "z_gemgroup_UMI", "core_scale_factor",
                  "perturbation_score", "guide_confidence_score", "target_expr_reduction",
                  "perturbation_burden", "pseudotime", "escape_probability"]
    keep_cols = [c for c in qc_metrics + ["perturbation", "cell_state"] if c in cell_df.columns]
    data["qc_cells"] = [_clean_record(r) for r in cell_df[keep_cols].to_dict(orient="records")]

    eff_path = csv_dir / "effect_decomposition.csv"
    data["effect_df"] = [_clean_record(r) for r in pd.read_csv(eff_path).to_dict(orient="records")] \
        if eff_path.exists() else []

    degs_path = csv_dir / "deg_summary.csv"
    data["deg_summary"] = [_clean_record(r) for r in pd.read_csv(degs_path).to_dict(orient="records")] \
        if degs_path.exists() else []

    data["deg"] = {}
    data["gsea"] = {}
    top_gene_union: List[str] = []
    top_gene_set: set = set()
    for f in sorted(csv_dir.glob("deg_*.csv")):
        if "deg_summary" in f.stem or "deg_enrichment" in f.stem:
            continue
        pert_name = f.stem[4:]
        df_full = pd.read_csv(f)
        # Per-pert DEG table is bundled into the HTML; trim to top 100 rows
        # to keep interactive_report.html small.
        df = df_full.head(100)
        data["deg"][pert_name] = [_clean_record(r) for r in df.to_dict(orient="records")]
        data["gsea"][pert_name] = _compute_gsea_rows(df_full)
        for g in df_full.head(n_top_genes)["gene"].tolist():
            if g not in top_gene_set:
                top_gene_set.add(g)
                top_gene_union.append(g)

    data["enrichment"] = {}
    for f in sorted(csv_dir.glob("deg_enrichment_*.csv")):
        pert_name = f.stem[len("deg_enrichment_"):]
        data["enrichment"][pert_name] = [_clean_record(r)
                                         for r in pd.read_csv(f).to_dict(orient="records")]

    data["umap"] = []
    if "X_umap" in adata.obsm:
        um = adata.obsm["X_umap"]
        n = min(5000, adata.n_obs)
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(adata.n_obs, size=n, replace=False))
        obs = adata.obs.iloc[idx]
        for i, cell_i in enumerate(idx):
            row_d: dict = {
                "x": round(float(um[cell_i, 0]), 4),
                "y": round(float(um[cell_i, 1]), 4),
                "pert": str(obs["perturbation"].iloc[i]) if "perturbation" in obs.columns else "",
            }
            if "cell_state" in obs.columns:
                row_d["state"] = str(obs["cell_state"].iloc[i])
            if "perturbation_score" in obs.columns:
                sc = _safe_float(obs["perturbation_score"].iloc[i])
                if sc is not None:
                    row_d["score"] = sc
            data["umap"].append(row_d)

    gene_idx_map = {g: int(i) for i, g in enumerate(adata.var_names)}
    valid_genes = [g for g in top_gene_union if g in gene_idx_map][: min(n_top_genes * 2, 300)]
    data["genes"] = valid_genes
    data["gene_expr"] = {}
    if valid_genes and "perturbation" in adata.obs.columns:
        perts_arr = adata.obs["perturbation"].astype(str).values
        unique_perts = sorted(set(perts_arr))
        for gene in valid_genes:
            gi = gene_idx_map[gene]
            entry: dict = {}
            for p in unique_perts:
                cidx = np.where(perts_arr == p)[0]
                if len(cidx) == 0:
                    continue
                if len(cidx) > max_cells_per_group:
                    rng2 = np.random.default_rng(abs(hash(gene + p)) % (2**31))
                    cidx = rng2.choice(cidx, max_cells_per_group, replace=False)
                vals = _to_dense(adata.X[cidx, gi]).flatten()
                entry[p] = [round(float(v), 4) for v in vals]
            data["gene_expr"][gene] = entry

    data["umap_genes"] = []
    data["umap_gene_expr"] = {}
    if data["umap"] and valid_genes:
        # Browser-side UMAP coloring needs per-cell values, not per-perturbation
        # medians. Keep this bounded so the HTML stays portable.
        umap_genes = valid_genes[:60]
        if "X_umap" in adata.obsm:
            n = min(5000, adata.n_obs)
            rng = np.random.default_rng(0)
            idx = np.sort(rng.choice(adata.n_obs, size=n, replace=False))
            for gene in umap_genes:
                gi = gene_idx_map.get(gene)
                if gi is None:
                    continue
                vals = _to_dense(adata.X[idx, gi]).flatten()
                data["umap_gene_expr"][gene] = [round(float(v), 4) for v in vals]
            data["umap_genes"] = [g for g in umap_genes if g in data["umap_gene_expr"]]

    hm_names = {
        "gene_by_cell": "eda_gene_by_cell_heatmap.png",
        "gene_by_pert": "eda_clustered_gene_by_pert_heatmap.png",
        "gene_corr":    "eda_gene_correlation_heatmap.png",
        "pert_sim":     "eda_perturbation_similarity.png",
        "umap_pert":    "umap_perturbation.png",
        "umap_state":   "umap_cell_state.png",
        "cluster_summary":        "eda_gene_pert_cluster_summary_heatmap.png",
    }
    # Static PNGs — store relative paths so the browser loads them lazily
    # from the sibling plots/ folder instead of base64-embedding ~30 large
    # images into the HTML. Field names are kept (heatmap_b64 / genenet_b64
    # / eda_corr_pert) for backward compat with the JS bindings.
    data["heatmap_b64"] = {}
    plots_dir = csv_dir.parent / "plots"
    for key, fname in hm_names.items():
        p = plots_dir / fname
        if p.exists():
            data["heatmap_b64"][key] = f"plots/{fname}"

    data["genenet_b64"] = {"control_heatmap": None, "control_network": None, "perturbations": {}}
    ctrl_hm = plots_dir / "genenet_control_heatmap.png"
    if ctrl_hm.exists():
        data["genenet_b64"]["control_heatmap"] = f"plots/{ctrl_hm.name}"
    ctrl_net = plots_dir / "genenet_control_network.png"
    if ctrl_net.exists():
        data["genenet_b64"]["control_network"] = f"plots/{ctrl_net.name}"
    # Per-pert single-axis network plots (preferred, square layout).
    for net_path in sorted(plots_dir.glob("genenet_*_pert_network.png")):
        pert = net_path.stem.replace("genenet_", "").replace("_pert_network", "")
        data["genenet_b64"]["perturbations"].setdefault(pert, {})
        data["genenet_b64"]["perturbations"][pert]["pert_network"] = f"plots/{net_path.name}"
    for net_path in sorted(plots_dir.glob("genenet_*_ctrl_network.png")):
        pert = net_path.stem.replace("genenet_", "").replace("_ctrl_network", "")
        data["genenet_b64"]["perturbations"].setdefault(pert, {})
        data["genenet_b64"]["perturbations"][pert]["ctrl_network"] = f"plots/{net_path.name}"
    # Legacy combined ctrl+pert network (kept as fallback for old runs).
    for net_path in sorted(plots_dir.glob("genenet_*_network.png")):
        if net_path.name.endswith("_pert_network.png") or net_path.name.endswith("_ctrl_network.png") \
                or net_path.name.endswith("_diff_network.png") or net_path.name == "genenet_control_network.png":
            continue
        pert = net_path.stem.replace("genenet_", "").replace("_network", "")
        data["genenet_b64"]["perturbations"].setdefault(pert, {})
        data["genenet_b64"]["perturbations"][pert]["network"] = f"plots/{net_path.name}"
    for diff_path in sorted(plots_dir.glob("genenet_*_diff_network.png")):
        pert = diff_path.stem.replace("genenet_", "").replace("_diff_network", "")
        data["genenet_b64"]["perturbations"].setdefault(pert, {})
        data["genenet_b64"]["perturbations"][pert]["diff"] = f"plots/{diff_path.name}"
    for hm_path in sorted(plots_dir.glob("genenet_*_heatmap_comparison.png")):
        pert = hm_path.stem.replace("genenet_", "").replace("_heatmap_comparison", "")
        data["genenet_b64"]["perturbations"].setdefault(pert, {})
        data["genenet_b64"]["perturbations"][pert]["heatmap"] = f"plots/{hm_path.name}"

    # Per-perturbation gene–gene correlation heatmaps from eda step
    data["eda_corr_pert"] = {}
    for cp in sorted(plots_dir.glob("eda_gene_corr_vs_*.png")):
        pert = cp.stem.replace("eda_gene_corr_vs_", "")
        data["eda_corr_pert"][pert] = f"plots/{cp.name}"

    # Interactive perturbation similarity matrix. Prefer the bundled cosine
    # similarity when available; otherwise the static PNG remains the fallback.
    data["pert_sim"] = {}
    sim_path = csv_dir.parent / "bundle" / "perturbation_similarity.parquet"
    if False and sim_path.exists():
        try:
            sim_long = pd.read_parquet(sim_path)
            sim_df = sim_long.pivot(
                index="perturbation_a",
                columns="perturbation_b",
                values="cosine",
            )
            labels = sim_df.index.astype(str).tolist()
            data["pert_sim"] = {
                "perts": labels,
                "matrix": [[_safe_float(v) for v in row] for row in sim_df.values],
                "label": "Cosine similarity",
            }
        except Exception:
            pass
    if False and not data["pert_sim"] and "perturbation" in adata.obs.columns:
        try:
            from .eda import _eda_gene_indices, _pseudobulk_mean, _zscore_columns

            gene_idx = _eda_gene_indices(adata, n_top_genes)
            pb = _pseudobulk_mean(adata, gene_idx)
            if not pb.empty and pb.shape[0] >= 3:
                pb_z = _zscore_columns(pb)
                mat = pb_z.values
                norms = np.linalg.norm(mat, axis=1, keepdims=True) + 1e-12
                sim = (mat / norms) @ (mat / norms).T
                np.fill_diagonal(sim, 1.0)
                labels = pb.index.astype(str).tolist()
                data["pert_sim"] = {
                    "perts": labels,
                    "matrix": [[_safe_float(v) for v in row] for row in sim],
                    "label": "Cosine similarity",
                }
        except Exception:
            pass

    # Interactive network JSON (nodes, edges, spring-layout positions).
    data["genenet_json"] = {}
    for json_path in sorted(csv_dir.glob("genenet_*_network.json")):
        safe_key = json_path.stem[len("genenet_"):-len("_network")]
        try:
            data["genenet_json"][safe_key] = json.loads(json_path.read_text())
        except Exception:
            pass

    # ── Cell state enrichment (state_enrich step) ─────────────────────────────
    data["state_enrich"] = {}
    se_path = csv_dir / "state_enrich_matrix.csv"
    if se_path.exists():
        try:
            se_df = pd.read_csv(se_path, index_col=0)
            data["state_enrich"] = {
                "perturbations": se_df.index.tolist(),
                "cell_states": se_df.columns.tolist(),
                "matrix": [[round(float(v), 3) for v in row] for row in se_df.values],
                "_is_fallback": False,
            }
        except Exception:
            pass
    # Fallback: cluster proportions log2FC vs control
    if not data["state_enrich"]:
        cp_path = csv_dir / "eda_cluster_proportions.csv"
        if cp_path.exists():
            try:
                cp_df = pd.read_csv(cp_path, index_col=0)
                if "control" in cp_df.index:
                    ctrl_row = cp_df.loc["control"]
                    pert_df = cp_df.drop("control")
                    ratio_df = np.log2((pert_df + 0.01) / (ctrl_row.values + 0.01))
                    data["state_enrich"] = {
                        "perturbations": ratio_df.index.tolist(),
                        "cell_states": ["Cluster " + str(c) for c in ratio_df.columns.tolist()],
                        "matrix": [[round(float(v), 3) for v in row] for row in ratio_df.values],
                        "_is_fallback": True,
                    }
            except Exception:
                pass

    # ── Perturbation effect correlation (regulatory step) ─────────────────────
    data["pert_corr_before"] = {}
    pcb_path = csv_dir / "pert_effect_corr_before.csv"
    if pcb_path.exists():
        try:
            pcb = pd.read_csv(pcb_path, index_col=0)
            data["pert_corr_before"] = {
                "perts": pcb.index.tolist(),
                "matrix": [[round(float(v), 4) for v in row] for row in pcb.values],
            }
        except Exception:
            pass

    data["pert_corr_after"] = {}
    pca_path = csv_dir / "pert_effect_corr_after.csv"
    if pca_path.exists():
        try:
            pca = pd.read_csv(pca_path, index_col=0)
            data["pert_corr_after"] = {
                "perts": pca.index.tolist(),
                "matrix": [[round(float(v), 4) for v in row] for row in pca.values],
            }
        except Exception:
            pass

    # ── TF regulatory heatmap + network (regulatory step) ─────────────────────
    data["tf_regulatory"] = {}
    tfh_path = csv_dir / "tf_regulatory_matrix.csv"
    if tfh_path.exists():
        try:
            tfh = pd.read_csv(tfh_path, index_col=0)
            data["tf_regulatory"] = {
                "guides": tfh.index.tolist(),
                "genes": tfh.columns.tolist(),
                "matrix": [[round(float(v), 3) for v in row] for row in tfh.values],
            }
        except Exception:
            pass

    data["tf_network"] = {}
    tfn_path = csv_dir / "tf_regulatory_network.json"
    if tfn_path.exists():
        try:
            data["tf_network"] = json.loads(tfn_path.read_text())
        except Exception:
            pass

    data["config"] = {}
    summary_path = csv_dir.parent / "summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text())
            data["config"] = summary.get("config", {})
        except Exception:
            pass
    checkpoint_path = csv_dir.parent / "checkpoint.json"
    if checkpoint_path.exists():
        try:
            checkpoint = json.loads(checkpoint_path.read_text())
            data["config"].setdefault("_completed_steps", checkpoint.get("completed_steps", []))
        except Exception:
            pass

    # ── C-score ───────────────────────────────────────────────────────────────
    data["cscore"] = {"summary": [], "edges": {}, "b64": {}}
    cscore_sum_path = csv_dir / "cscore_summary.csv"
    if cscore_sum_path.exists():
        try:
            cs_df = pd.read_csv(cscore_sum_path)
            data["cscore"]["summary"] = [_clean_record(r) for r in cs_df.to_dict(orient="records")]
        except Exception:
            pass
    for edge_path in sorted(csv_dir.glob("cscore_edges_*.csv")):
        safe_key = edge_path.stem[len("cscore_edges_"):]
        try:
            edf = pd.read_csv(edge_path)
            data["cscore"]["edges"][safe_key] = [_clean_record(r) for r in edf.to_dict(orient="records")]
        except Exception:
            pass
    cscore_png_map = {
        "ranked_bar":    "cscore_ranked_bar.png",
        "decomposition": "cscore_decomposition.png",
        "vs_deg":        "cscore_vs_deg.png",
        "hub_heatmap":   "cscore_gene_hub_heatmap.png",
        "module_rewiring": "cscore_module_rewiring.png",
    }
    for key, fname in cscore_png_map.items():
        p = plots_dir / fname
        if p.exists():
            data["cscore"]["b64"][key] = f"plots/{fname}"

    return data


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

def _html_template() -> str:
    return (
        '<!DOCTYPE html>\n'
        '<html lang="en">\n'
        '<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<title>PerturbFlow-Analyzer \u2014 Interactive Report</title>\n'
        '<style>\n'
        ':root{'
        # GEPIA-inspired dark-navy palette: deep navy primary, lighter blue
        # secondary, mostly-white panels, subtly bluish backgrounds.
        '--bg:#f3f6f9;--panel:#ffffff;--border:#c0d1de;'
        '--accent:#1f4e79;--accent2:#2a7fa6;--text:#1a2733;--muted:#5a6b7d;'
        '--danger:#b04a5a;--warn:#b58b1d;--card-bg:#ffffff;--tab-active:#e1ecf5;'
        '--nav:#1f4e79;--nav-fg:#ffffff;--soft:#ecf2f7;--font:Arial,Helvetica,sans-serif;}\n'
        '*,*::before,*::after{box-sizing:border-box;margin:0;padding:0;}\n'
        'body{background:var(--bg);color:var(--text);font-family:var(--font);'
        'font-size:15px;display:flex;flex-direction:column;min-height:100vh;}\n'
        'header{background:linear-gradient(#0c2742,#163d61);border-bottom:1px solid #08203a;'
        'padding:14px 28px;display:flex;align-items:center;gap:18px;}\n'
        'header h1{font-size:24px;font-weight:700;color:var(--nav-fg);}\n'
        'header span{color:rgba(255,255,255,0.82);font-size:14px;}\n'
        '.layout{display:flex;flex:1;overflow:hidden;height:calc(100vh - 57px);}\n'
        '.sidebar{width:240px;min-width:200px;background:var(--panel);'
        'border-right:1px solid var(--border);display:flex;flex-direction:column;'
        'padding:16px 12px;overflow-y:auto;flex-shrink:0;}\n'
        '.sidebar h3{font-size:12px;text-transform:uppercase;letter-spacing:.08em;'
        'color:var(--muted);margin-bottom:8px;margin-top:14px;}\n'
        '.sidebar h3:first-child{margin-top:0;}\n'
        '.sidebar input,.sidebar select,input,select{background:#fff;'
        'border:1px solid var(--border);border-radius:6px;color:var(--text);'
        'padding:10px 12px;font-size:15px;outline:none;margin-bottom:6px;min-height:40px;}\n'
        '.sidebar input,.sidebar select{width:100%;}\n'
        '.sidebar input:focus,.sidebar select:focus,input:focus,select:focus{border-color:var(--accent);box-shadow:0 0 0 3px rgba(31,78,121,.14);}\n'
        '.sidebar button{width:100%;background:var(--tab-active);border:1px solid var(--border);'
        'border-radius:6px;color:var(--text);padding:9px 10px;font-size:14px;'
        'cursor:pointer;margin-bottom:4px;text-align:left;}\n'
        '.sidebar button:hover{border-color:var(--accent);color:var(--accent);}\n'
        '.sidebar button.active{background:var(--accent);color:#fff;font-weight:600;}\n'
        '.content{flex:1;overflow-y:auto;padding:24px 28px;}\n'
        '.tabs{display:flex;gap:4px;border-bottom:1px solid var(--border);'
        'margin-bottom:22px;flex-wrap:wrap;}\n'
        '.tab-btn{background:none;border:none;border-bottom:3px solid transparent;'
        'color:var(--muted);padding:10px 22px;cursor:pointer;font-size:18px;font-weight:600;'
        'transition:color .15s,border-color .15s;}\n'
        '.tab-btn:hover{color:var(--text);}\n'
        '.tab-btn.active{color:var(--accent);border-bottom-color:var(--accent);}\n'
        '.tab-pane{display:none;}.tab-pane.active{display:block;}\n'
        '.cards{display:flex;gap:16px;flex-wrap:wrap;margin-bottom:22px;}\n'
        '.card{background:var(--card-bg);border:1px solid var(--border);border-radius:6px;'
        'padding:18px 22px;flex:1 1 160px;min-width:140px;}\n'
        '.card .label{font-size:12px;color:var(--muted);text-transform:uppercase;'
        'letter-spacing:.06em;margin-bottom:6px;}\n'
        '.card .value{font-size:28px;font-weight:700;color:var(--accent);}\n'
        '.card .value.green{color:var(--accent2);}\n'
        'h2{font-size:20px;font-weight:600;margin-bottom:16px;margin-top:4px;color:var(--text);}\n'
        'h3{font-size:16px;font-weight:600;margin-bottom:10px;margin-top:18px;color:var(--text);}\n'
        'p{color:var(--muted);font-size:14px;line-height:1.6;margin-bottom:10px;}\n'
        '.chart{background:#fff;border:1px solid var(--border);border-radius:6px;overflow:hidden;margin-bottom:18px;}\n'
        '.chart-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(440px,1fr));gap:18px;}\n'
        '.chart-full{width:100%;}\n'
        '.tbl-wrap{overflow-x:auto;border-radius:8px;border:1px solid var(--border);margin-bottom:14px;}\n'
        'table{width:100%;border-collapse:collapse;font-size:13px;}\n'
        'thead th{background:var(--tab-active);padding:9px 12px;text-align:left;'
        'font-size:13px;font-weight:600;color:var(--muted);cursor:pointer;white-space:nowrap;'
        'border-bottom:1px solid var(--border);user-select:none;}\n'
        'thead th:hover{color:var(--accent);}\n'
        'tbody tr:nth-child(even){background:var(--card-bg);}\n'
        'tbody tr:hover{background:var(--tab-active);}\n'
        'tbody td{padding:7px 12px;border-bottom:1px solid var(--border);}\n'
        'tbody td.up{color:#b04a5a;}tbody td.dn{color:#1f4e79;}tbody td.sig{color:var(--accent2);}\n'
        '.pagination{display:flex;gap:8px;align-items:center;margin-top:8px;}\n'
        '.pagination button{background:var(--tab-active);border:1px solid var(--border);'
        'border-radius:5px;color:var(--text);padding:4px 12px;cursor:pointer;font-size:13px;}\n'
        '.pagination button:hover{border-color:var(--accent);}\n'
        '.pagination span{font-size:13px;color:var(--muted);}\n'
        '.hm-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(520px,1fr));gap:22px;}\n'
        '.hm-card{background:var(--card-bg);border:1px solid var(--border);border-radius:6px;padding:16px;}\n'
        '.hm-card h3{margin-top:0;}\n'
        '.hm-card img{width:100%;border-radius:4px;margin-top:8px;}\n'
        '.hm-card.hm-wide img{height:560px;object-fit:contain;background:#fff;}\n'
        '.metric-help{background:#fff;border:1px solid var(--border);border-radius:6px;'
        'padding:14px 16px;margin:0 0 18px;}\n'
        '.metric-help dl{display:grid;grid-template-columns:minmax(170px,240px) 1fr;gap:8px 14px;}\n'
        '.metric-help dt{font-weight:700;color:var(--accent);font-size:13px;}\n'
        '.metric-help dd{color:var(--muted);font-size:13px;line-height:1.45;}\n'
        '.hm-viewer{background:#fff;border:1px solid var(--border);border-radius:6px;padding:16px;margin-bottom:22px;}\n'
        '.hm-canvas{height:650px;overflow:auto;border:1px solid var(--border);background:#fff;text-align:center;}\n'
        '.hm-canvas img{transform-origin:top left;max-width:none;image-rendering:auto;}\n'
        '.network-img{width:100%;border:1px solid var(--border);border-radius:4px;background:#fff;}\n'
        '.control-panel{display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;margin-bottom:14px;}\n'
        '.control-panel label{font-size:13px;color:var(--muted);display:block;margin-bottom:4px;}\n'
        '.badge{display:inline-block;padding:2px 8px;border-radius:12px;font-size:11px;'
        'font-weight:600;margin-left:6px;background:var(--tab-active);color:var(--accent);}\n'
        '.steps-row{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:18px;}\n'
        '.step-chip{background:var(--tab-active);border:1px solid var(--border);'
        'border-radius:20px;padding:4px 14px;font-size:13px;color:var(--muted);}\n'
        '.step-chip.done{border-color:var(--accent2);color:var(--accent2);}\n'
        '.param-table{width:100%;border-collapse:collapse;font-size:14px;margin-bottom:18px;}\n'
        '.param-table th{text-align:left;padding:8px 12px;background:var(--tab-active);'
        'color:var(--muted);font-size:13px;font-weight:600;}\n'
        '.param-table td{padding:7px 12px;border-bottom:1px solid var(--border);}\n'
        '.param-table td.pn{color:var(--accent);font-family:monospace;font-size:13px;}\n'
        '.param-table td.pv{color:var(--accent2);font-family:monospace;font-size:13px;}\n'
        '::-webkit-scrollbar{width:6px;height:6px;}\n'
        '::-webkit-scrollbar-track{background:transparent;}\n'
        '::-webkit-scrollbar-thumb{background:var(--border);border-radius:3px;}\n'
        '::-webkit-scrollbar-thumb:hover{background:var(--muted);}\n'
        '</style>\n'
        '</head>\n'
        '<body>\n'
        '<header>\n'
        '  <h1>PerturbFlow-Analyzer</h1>\n'
        '  <span id="hdr-summary">Loading\u2026</span>\n'
        '</header>\n'
        '<div class="layout">\n'
        '<aside class="sidebar">\n'
        '  <h3>Navigate</h3>\n'
        '  <button class="active" onclick="showTab(\'home\',this)">Summary</button>\n'
        '  <button onclick="showTab(\'qc\',this)">QC</button>\n'
        '  <button onclick="showTab(\'heatmaps\',this)">Heatmaps</button>\n'
        '  <button onclick="showTab(\'pert\',this)">Perturbation</button>\n'
        '  <button onclick="showTab(\'network\',this)">Gene Network</button>\n'
        '  <button onclick="showTab(\'gene\',this)">Gene Expression</button>\n'
        '  <button onclick="showTab(\'states\',this)">Cell States</button>\n'
        '  <button onclick="showTab(\'cscore\',this)">C-Score</button>\n'
        '  <button onclick="showTab(\'regulatory\',this)">Regulatory</button>\n'
        '  <button onclick="showTab(\'params\',this)">Parameters</button>\n'
        '  <h3>Perturbation</h3>\n'
        '  <input id="sb-pert-search" type="text" placeholder="Search perturbation\u2026"\n'
        '         oninput="filterSidebarPerts(this.value)">\n'
        '  <div id="sb-pert-list" style="max-height:220px;overflow-y:auto;"></div>\n'
        '  <button id="sb-pert-show-all" style="display:none;font-size:12px;color:var(--accent);"\n'
        '          onclick="showAllPerts()"></button>\n'
        '  <h3>Gene</h3>\n'
        '  <input id="sb-gene-search" type="text" placeholder="Search gene\u2026"\n'
        '         oninput="filterSidebarGenes(this.value)">\n'
        '  <div id="sb-gene-list" style="max-height:160px;overflow-y:auto;"></div>\n'
        '</aside>\n'
        '<main class="content">\n'
        '  <div class="tabs">\n'
        '    <button class="tab-btn active" onclick="showTab(\'home\',null,\'tab\')">Summary</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'qc\',null,\'tab\')">QC</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'heatmaps\',null,\'tab\')">Heatmaps</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'pert\',null,\'tab\')">Perturbation</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'network\',null,\'tab\')">Gene Network</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'gene\',null,\'tab\')">Gene Expression</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'states\',null,\'tab\')">Cell States</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'cscore\',null,\'tab\')">C-Score</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'regulatory\',null,\'tab\')">Regulatory</button>\n'
        '    <button class="tab-btn" onclick="showTab(\'params\',null,\'tab\')">Parameters</button>\n'
        '  </div>\n'
        + _TAB_HOME
        + _TAB_QC
        + _TAB_HEATMAPS
        + _TAB_PERT
        + _TAB_NETWORK
        + _TAB_GENE
        + _TAB_STATES
        + _TAB_CSCORE
        + _TAB_REGULATORY
        + _TAB_PARAMS
        + '</main>\n'
        '</div>\n'
        '<script>\n'
        'let D = {DATA_STUB};\n'
        + _JS_BODY
        + '</script>\n'
        '</body>\n'
        '</html>\n'
    )


def _copy_plotly_bundle(out: Path) -> None:
    """Place Plotly next to the report so it opens without internet access."""
    candidates: list[Path] = []
    try:
        import plotly  # type: ignore

        candidates.append(Path(plotly.__file__).resolve().parent / "package_data" / "plotly.min.js")
    except Exception:
        pass
    candidates.extend([
        Path("/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env/lib/python3.10/site-packages/plotly/package_data/plotly.min.js"),
        Path("/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env/lib/python3.12/site-packages/plotly/package_data/plotly.min.js"),
    ])
    for src in candidates:
        if src.exists():
            dst = out / "plotly.min.js"
            if not dst.exists() or dst.stat().st_size != src.stat().st_size:
                shutil.copy2(src, dst)
            return


_TAB_HOME = (
    '<div class="tab-pane active" id="tab-home">\n'
    '  <h2>Dataset Summary</h2>\n'
    '  <div class="cards" id="home-cards"></div>\n'
    '  <h3>Analysis Summary</h3>\n'
    '  <p id="home-analysis-text"></p>\n'
    '  <div class="steps-row" id="home-steps"></div>\n'
    '  <h3>Top Perturbations by Effect</h3>\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr 480px;gap:18px;align-items:start;">\n'
    '    <div class="chart" id="home-top-perts-chart" style="height:480px;"></div>\n'
    '    <div class="chart" id="home-deg-bar" style="height:480px;"></div>\n'
    '    <div style="width:480px;flex-shrink:0;">\n'
    '      <div id="home-effect-scatter" style="height:480px;"></div>\n'
    '      <div style="text-align:right;margin-top:6px;">\n'
    '        <button onclick="downloadEffectPng()"\n'
    '                style="background:var(--tab-active);border:1px solid var(--border);border-radius:6px;\n'
    '                       color:var(--accent);padding:6px 14px;cursor:pointer;font-size:12px;">\n'
    '          ⬇ Save Effect Decomposition as PNG\n'
    '        </button>\n'
    '      </div>\n'
    '    </div>\n'
    '  </div>\n'
    '  <hr style="border:none;border-top:1px solid var(--border);margin:32px 0 22px;">\n'
    '  <h3>Cell Embedding (UMAP)</h3>\n'
    '  <p style="color:var(--muted);font-size:14px;line-height:1.6;margin-bottom:16px;">\n'
    '    UMAP dimensionality reduction of all cells. Use the controls below to choose which view to display.\n'
    '  </p>\n'
    '  <div class="control-panel" style="margin-bottom:12px;">\n'
    '    <div>\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Color by</label>\n'
    '      <select id="umap-color-by" style="min-width:200px;" onchange="onUmapColorChange()">\n'
    '        <option value="perturbation">Perturbation</option>\n'
    '        <option value="state">Cell State</option>\n'
    '        <option value="score">Perturbation Score</option>\n'
    '        <option value="gene">Gene expression</option>\n'
    '      </select>\n'
    '    </div>\n'
    '    <div id="umap-gene-wrap" style="display:none;">\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Gene</label>\n'
    '      <input id="umap-gene-input" type="text" placeholder="Type gene name\u2026"\n'
    '             style="min-width:200px;" list="umap-gene-list" oninput="debounceRenderUMAP()">\n'
    '      <datalist id="umap-gene-list"></datalist>\n'
    '    </div>\n'
    '  </div>\n'
    '  <div class="chart chart-full" id="umap-chart" style="height:560px;aspect-ratio:1.15/1;max-width:780px;margin:0 auto;"></div>\n'
  '  <div style="text-align:right;max-width:780px;margin:8px auto 0;">\n'
  '    <button onclick="downloadUmapPng()"\n'
  '            style="background:var(--tab-active);border:1px solid var(--border);border-radius:6px;\n'
  '                   color:var(--accent);padding:7px 16px;cursor:pointer;font-size:13px;">\n'
  '      ⬇ Save UMAP as PNG\n'
  '    </button>\n'
  '  </div>\n'
    '</div>\n'
)

_TAB_QC = (
    '<div class="tab-pane" id="tab-qc">\n'
    '  <h2>Quality Control</h2>\n'
    '  <h3>QC Summary</h3>\n'
    '  <div class="cards" id="qc-summary-cards" style="margin-bottom:18px;"></div>\n'
    '  <div class="cards" id="qc-metric-cards" style="margin-bottom:18px;"></div>\n'
    '  <div class="metric-help">\n'
    '    <dl>\n'
    '      <dt>Guide confidence</dt><dd>Per-cell confidence that the assigned guide/perturbation is reliable, when provided by upstream guide calling or Mixscape-style scoring.</dd>\n'
    '      <dt>Target expr reduction</dt><dd>Estimated decrease of the target gene expression relative to control cells. Higher values indicate stronger knockdown for CRISPRi-style perturbations.</dd>\n'
    '      <dt>Perturbation burden</dt><dd>Aggregate strength of perturbation-associated transcriptional change for a cell; useful for spotting overly strong, stressed, or multi-hit cells.</dd>\n'
    '      <dt>Escape probability</dt><dd>Estimated probability that a cell escaped the intended perturbation response, usually because target knockdown or downstream signature is weak.</dd>\n'
    '    </dl>\n'
    '  </div>\n'
    '  <div class="chart" id="qc-cells-bar" style="height:280px;"></div>\n'
    '  <div class="chart-grid" id="qc-violin-grid"></div>\n'
    '</div>\n'
)

_TAB_UMAP = (
    '<div class="tab-pane" id="tab-umap">\n'
    '  <h2>UMAP Explorer</h2>\n'
    '  <div class="control-panel">\n'
    '    <div>\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Color by</label>\n'
    '      <select id="umap-color-by" style="min-width:200px;" onchange="onUmapColorChange()">\n'
    '        <option value="perturbation">Perturbation</option>\n'
    '        <option value="state">Cell State</option>\n'
    '        <option value="score">Perturbation Score</option>\n'
    '        <option value="gene">Gene expression</option>\n'
    '      </select>\n'
    '    </div>\n'
    '    <div id="umap-gene-wrap" style="display:none;">\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Gene</label>\n'
    '      <input id="umap-gene-input" type="text" placeholder="Type gene name\u2026"\n'
    '             style="min-width:200px;" list="umap-gene-list" oninput="renderUMAP()">\n'
    '      <datalist id="umap-gene-list"></datalist>\n'
    '    </div>\n'
    '  </div>\n'
    '  <div class="chart chart-full" id="umap-chart" style="height:760px;aspect-ratio:1.2/1;max-width:920px;margin:0 auto;"></div>\n'
  '  <div style="text-align:right;max-width:920px;margin:8px auto 0;">\n'
  '    <button onclick="downloadUmapPng()"\n'
  '            style="background:var(--tab-active);border:1px solid var(--border);border-radius:6px;\n'
  '                   color:var(--accent);padding:7px 16px;cursor:pointer;font-size:13px;">\n'
  '      ⬇ Save UMAP as PNG\n'
  '    </button>\n'
  '  </div>\n'
    '</div>\n'
)

_TAB_HEATMAPS = (
    '<div class="tab-pane" id="tab-heatmaps">\n'
    '  <h2>Heatmaps</h2>\n'
    '  <p>Pre-computed heatmaps from the EDA pipeline step. Heatmaps use the configured top HVGs where available; correlation matrices are Pearson; similarity is cosine.</p>\n'
    '  <div class="hm-grid" id="heatmaps-grid"></div>\n'
    '  <hr style="border:none;border-top:1px solid var(--border);margin:28px 0 20px;">\n'
    '  <h3 style="font-size:16px;margin-bottom:12px;">Gene\u2013Gene Correlation: Control vs Perturbation</h3>\n'
    '  <div class="control-panel" style="margin-bottom:14px;">\n'
    '    <div style="flex:1;min-width:260px;"><label>Perturbation</label>'
    '<select id="corr-pert-selector" onchange="renderCorrPert(this.value)"></select></div>\n'
    '  </div>\n'
    '  <div id="corr-pert-img-wrap" style="text-align:center;">\n'
    '    <img id="corr-pert-img" style="max-width:100%;border:1px solid var(--border);border-radius:4px;" alt="Gene correlation heatmap">\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_PERT = (
    '<div class="tab-pane" id="tab-pert">\n'
    '  <h2>Perturbation Analysis</h2>\n'
    '  <div class="control-panel">\n'
    '    <div style="flex:1;min-width:220px;">\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">\n'
    '        Select perturbation <span id="pert-count-badge" class="badge"></span>\n'
    '      </label>\n'
    '      <select id="pert-selector" style="width:100%;" onchange="renderPert(this.value)"></select>\n'
    '    </div>\n'
    '    <button id="pert-show-all-btn"\n'
    '            style="background:var(--tab-active);border:1px solid var(--border);border-radius:6px;\n'
    '                   color:var(--accent);padding:7px 14px;cursor:pointer;font-size:13px;"\n'
    '            onclick="expandPertDropdown()"></button>\n'
    '  </div>\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr;gap:18px;align-items:stretch;margin-bottom:14px;">\n'
    '    <div class="chart" id="pert-volcano" style="height:500px;min-width:0;"></div>\n'
    '    <div class="chart" id="pert-umap-hl" style="height:500px;min-width:0;"></div>\n'
    '  </div>\n'
    '  <div class="control-panel" style="margin-top:4px;">\n'
    '    <div><label>DEG bar filter</label><select id="deg-bar-filter" onchange="reRenderDegBar()" style="min-width:220px;">\n'
    '      <option value="sig">Significant only (top 20)</option>\n'
    '      <option value="top20lfc">Top 20 by |log\u2082FC|</option>\n'
    '      <option value="all">All DEGs (top 40)</option>\n'
    '    </select></div>\n'
    '  </div>\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr;gap:22px;align-items:start;margin-top:14px;">\n'
    '    <div>\n'
    '      <h3 style="margin-top:0;">Top DEGs</h3>\n'
    '      <div class="chart" id="pert-top-deg-bar" style="height:440px;"></div>\n'
    '    </div>\n'
    '    <div>\n'
    '      <h3 style="margin-top:0;">Pathway Enrichment</h3>\n'
    '      <div class="chart" id="pert-enrichment-bar" style="height:440px;"></div>\n'
    '    </div>\n'
    '  </div>\n'
    '  <h3>GSEA (preranked curated gene sets)</h3>\n'
    '  <p style="font-size:13px;color:var(--muted);margin-top:-4px;">Exploratory preranked GSEA using the embedded curated gene-set catalogue in <code>perturbflow/analyzer/pathways.py</code> (Hallmark-like, KEGG/Reactome-style, and TF-target sets). This is not a full GO database unless GO terms are explicitly added to that catalogue. Positive ES means enrichment among up-regulated genes; negative ES means enrichment among down-regulated genes. Click a term bar to show the running enrichment plot.</p>\n'
    '  <div class="chart" id="pert-gsea-bar" style="height:420px;"></div>\n'
    '  <div class="chart" id="pert-gsea-waterfall" style="height:360px;display:none;"></div>\n'
    '  <h3>DEG Table <span id="deg-table-info" class="badge"></span></h3>\n'
    '  <div class="control-panel">\n'
    '    <div style="flex:1;min-width:260px;"><label>Search genes and pathways</label><input id="deg-filter" type="search" placeholder="Filter DEG results..." oninput="filterDegTable(this.value)"></div>\n'
    '    <div><label>Rows per page</label><select id="deg-page-size" onchange="setDegPageSize(this.value)"><option value="10" selected>10</option><option value="25">25</option><option value="100">100</option><option value="200">200</option></select></div>\n'
    '    <button onclick="downloadDegTable()" style="background:var(--tab-active);border:1px solid var(--border);border-radius:6px;color:var(--accent);padding:9px 18px;cursor:pointer;font-size:14px;white-space:nowrap;">\u2b07 Download CSV</button>\n'
    '  </div>\n'
    '  <div class="tbl-wrap">\n'
    '    <table id="deg-table">\n'
    '      <thead><tr>\n'
    '        <th onclick="sortDegTable(0)">Gene</th>\n'
    '        <th onclick="sortDegTable(1)">log\u2082FC</th>\n'
    '        <th onclick="sortDegTable(2)">p-adj</th>\n'
    '        <th onclick="sortDegTable(3)">Sig</th>\n'
    '        <th onclick="sortDegTable(4)">Pathway</th>\n'
    '      </tr></thead>\n'
    '      <tbody id="deg-tbody"></tbody>\n'
    '    </table>\n'
    '  </div>\n'
    '  <div class="pagination">\n'
    '    <button onclick="degPage(-1)">\u2190 Prev</button>\n'
    '    <span id="deg-page-info">Page 1</span>\n'
    '    <button onclick="degPage(1)">Next \u2192</button>\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_NETWORK = (
    '<div class="tab-pane" id="tab-network">\n'
    '  <h2>Gene Regulatory Network</h2>\n'
    '  <p>Co-expression networks and differential gene networks from Pearson correlations.</p>\n'
    '  <div class="control-panel" style="margin-bottom:18px;">\n'
    '    <div style="flex:1;min-width:260px;"><label>Perturbation</label>'
    '<select id="network-selector" onchange="renderNetwork(this.value)"></select></div>\n'
    '  </div>\n'
    '  <!-- Row 1: three square networks side by side -->\n'
    '  <div style="display:grid;grid-template-columns:repeat(3,minmax(260px,1fr));gap:20px;margin-bottom:28px;">\n'
    '    <div style="aspect-ratio:1/1;display:flex;flex-direction:column;">\n'
    '      <h3 id="network-ctrl-label" style="font-size:15px;margin:0 0 8px;text-align:center;">Control Network</h3>\n'
    '      <img id="network-ctrl-net-img" style="width:100%;height:100%;object-fit:contain;border:1px solid var(--border);border-radius:4px;background:#fff;" alt="Control co-expression network">\n'
    '    </div>\n'
    '    <div style="aspect-ratio:1/1;display:flex;flex-direction:column;">\n'
    '      <h3 id="network-pert-label" style="font-size:15px;margin:0 0 8px;text-align:center;">Perturbation Network</h3>\n'
    '      <img id="network-pert-img" style="width:100%;height:100%;object-fit:contain;border:1px solid var(--border);border-radius:4px;background:#fff;" alt="Perturbation co-expression network">\n'
    '    </div>\n'
    '    <div style="aspect-ratio:1/1;display:flex;flex-direction:column;">\n'
    '      <h3 id="network-diff-label" style="font-size:15px;margin:0 0 8px;text-align:center;">Differential Network</h3>\n'
    '      <img id="network-diff-img" style="width:100%;height:100%;object-fit:contain;border:1px solid var(--border);border-radius:4px;background:#fff;" alt="Differential network">\n'
    '    </div>\n'
    '  </div>\n'
    '  <!-- Row 2: bottom heatmaps for the genes used in the networks -->\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:24px;">\n'
    '    <div>\n'
    '      <h3 style="font-size:15px;margin-bottom:8px;">Control Gene Co-expression Heatmap</h3>\n'
    '      <div class="chart" id="network-ctrl-heatmap-chart" style="height:620px;"></div>\n'
    '    </div>\n'
    '    <div>\n'
    '      <h3 id="network-heatmap-label" style="font-size:15px;margin-bottom:8px;">Gene Expression Heatmap</h3>\n'
    '      <img id="network-heatmap-img" style="width:100%;max-height:620px;object-fit:contain;border:1px solid var(--border);border-radius:4px;" alt="Gene expression heatmap">\n'
    '      <div style="display:flex;justify-content:flex-end;margin-top:8px;">\n'
    '        <button id="network-dl-btn" style="background:var(--tab-active);border:1px solid var(--border);border-radius:6px;color:var(--accent);padding:6px 14px;cursor:pointer;font-size:13px;">\u2b07 Download</button>\n'
    '      </div>\n'
    '    </div>\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_GENE = (
    '<div class="tab-pane" id="tab-gene">\n'
    '  <h2>Gene Expression</h2>\n'
    '  <p>For the selected gene, the left plot shows sampled single-cell expression grouped by perturbation. The right plot summarizes each perturbation as median log\u2082 fold-change relative to control, so this page is perturbation-based.</p>\n'
    '  <div style="display:flex;gap:12px;margin-bottom:18px;flex-wrap:wrap;align-items:flex-end;">\n'
    '    <div style="flex:1;min-width:220px;">\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Select gene</label>\n'
    '      <select id="gene-selector" style="width:100%;" onchange="renderGeneExpr(this.value)"></select>\n'
    '    </div>\n'
    '    <div>\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Custom gene</label>\n'
    '      <input id="gene-custom-input" type="text" placeholder="Type gene name\u2026"\n'
    '             style="min-width:180px;"\n'
    '             onkeydown="if(event.key===\'Enter\')renderGeneExpr(this.value.trim())">\n'
    '    </div>\n'
    '  </div>\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr;gap:18px;align-items:stretch;">\n'
    '    <div class="chart" id="gene-box" style="height:480px;min-width:0;"></div>\n'
    '    <div class="chart" id="gene-fc-bar" style="height:480px;min-width:0;"></div>\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_STATES = (
    '<div class="tab-pane" id="tab-states">\n'
    '  <h2>Cell State Enrichment / Depletion</h2>\n'
    '  <div style="background:#f0f7ff;border-left:4px solid var(--accent);padding:14px 18px;'
    'border-radius:4px;margin-bottom:18px;">\n'
    '    <p style="margin:0 0 8px;font-size:14px;font-weight:600;">How to read this heatmap</p>\n'
    '    <ul style="margin:0;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '      <li><span style="color:#B2182B;font-weight:700;">\u25a0 Red (positive value)</span>'
    ' \u2014 the perturbation <b>enriches</b> cells in this state: a larger fraction of perturbed cells'
    ' occupy this state compared to unperturbed controls.</li>\n'
    '      <li><span style="color:#2166AC;font-weight:700;">\u25a0 Blue (negative value)</span>'
    ' \u2014 the perturbation <b>depletes</b> cells from this state: fewer perturbed cells'
    ' reside in this state relative to controls.</li>\n'
    '      <li><b>Value</b> = signed \u2212log\u2081\u2080(q-value). |value|\u202f>\u202f1.3 corresponds to'
    ' q\u202f<\u202f0.05; |value|\u202f>\u202f2 corresponds to q\u202f<\u202f0.01.'
    ' Larger magnitude = stronger and more significant shift.</li>\n'
    '      <li>Each <b>row</b> is one perturbation; each <b>column</b> is one cell state (Leiden cluster).</li>\n'
    '      <li>Rows near zero indicate perturbations with no significant cell-state redistribution.</li>\n'
    '    </ul>\n'
    '  </div>\n'
    '  <div class="control-panel">\n'
    '    <div><label>Min abs. signal</label>'
    '<input id="states-min-q" type="number" value="0" min="0" step="0.5"'
    ' style="width:80px;" oninput="renderStateEnrich()"></div>\n'
    '  </div>\n'
    '  <div class="chart chart-full" id="states-heatmap-chart" style="height:640px;"></div>\n'
    '  <div style="background:#f8f9fb;border:1px solid var(--border);border-radius:6px;'
    'padding:14px 18px;margin-top:18px;">\n'
    '    <p style="margin:0;font-size:13px;color:var(--muted);line-height:1.7;">\n'
    '      <b>Statistical method:</b> Chi-square test of independence comparing the cell-state'
    ' composition of each perturbation group vs control cells.'
    ' The signed \u2212log\u2081\u2080(q) is positive when the perturbation increases the fraction of'
    ' cells in a given state (enrichment) and negative when it decreases it (depletion).'
    ' BH\u2013FDR correction is applied across all perturbation\u202f\u00d7\u202fcell-state pairs.'
    ' When the formal test cannot be run, log\u2082\u202fFC of cluster proportions vs control is'
    ' displayed instead.\n'
    '    </p>\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_REGULATORY = (
    '<div class="tab-pane" id="tab-regulatory">\n'
    '  <h2>Regulatory Analysis</h2>\n'
    '  <div style="display:flex;flex-direction:column;gap:10px;margin-bottom:18px;">\n'
    '    <div style="font-size:18px;font-weight:700;color:var(--accent);">Effect Correlation</div>\n'
    '    <div style="font-size:18px;font-weight:700;color:var(--accent);">TF Regulatory Heatmap</div>\n'
    '    <div style="font-size:18px;font-weight:700;color:var(--accent);">TF Regulatory Network</div>\n'
    '  </div>\n'

    '  <div id="reg-corr-pane">\n'
    '    <h3>Effect Correlation</h3>\n'
    '    <div style="background:#f0f7ff;border-left:4px solid var(--accent);padding:14px 18px;'
    'border-radius:4px;margin-bottom:16px;">\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">What this heatmap shows</p>\n'
    '      <p style="margin:0 0 10px;font-size:13px;color:#3a4a58;line-height:1.7;">'
    'Each cell is the <b>Pearson correlation coefficient</b> between the log2 fold-change'
    ' vectors of two perturbations, computed across all genes. In other words, it asks:'
    ' if I summarise each perturbation as its own gene-wide expression-change profile,'
    ' how similar are any two profiles to each other?'
    ' Both rows and columns are perturbations, so the matrix is symmetric and the diagonal'
    ' is always one (a perturbation is perfectly correlated with itself).</p>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">Colour legend</p>\n'
    '      <ul style="margin:0 0 10px;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li><span style="color:#B2182B;font-weight:700;">Red cells (positive correlation)</span>'
    ' indicate that the two perturbations push the transcriptome in the <b>same direction</b>.'
    ' Genes that go up under one perturbation also tend to go up under the other, and likewise'
    ' for genes that go down. Strong red blocks along the diagonal point to functionally'
    ' redundant transcription factors, members of the same complex, or perturbations that'
    ' converge on a shared pathway.</li>\n'
    '        <li><span style="color:#2166AC;font-weight:700;">Blue cells (negative correlation)</span>'
    ' indicate that the two perturbations push the transcriptome in <b>opposite directions</b>.'
    ' This often reflects antagonistic regulators, where one factor activates and the other'
    ' represses the same downstream programme.</li>\n'
    '        <li><span style="color:#888;font-weight:700;">Near-zero (white) cells</span>'
    ' mean the two perturbations affect largely independent gene sets and have no measurable'
    ' overlap in their transcriptional consequences.</li>\n'
    '      </ul>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">The two panels</p>\n'
    '      <ul style="margin:0 0 10px;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li><b>Left panel (Without modeling cell states):</b> raw correlation computed'
    ' on log2 fold-change vectors derived from all cells of each perturbation pooled together.'
    ' This view is the simplest summary but can be confounded if one perturbation pushes cells'
    ' into a different cell state, because the apparent gene-expression change may reflect a'
    ' shift in cell-type composition rather than a direct regulatory effect.</li>\n'
    '        <li><b>Right panel (After modeling cell states):</b> the log2 fold-change is'
    ' first recomputed within each cell state separately, then averaged using the control'
    ' population\u2019s state composition as fixed weights. This removes the composition'
    ' confound, so the right panel reflects <b>direct transcriptional effects</b> rather than'
    ' indirect cell-state redistribution. Differences between the two panels point to'
    ' perturbations whose apparent similarity is driven by composition rather than wiring.</li>\n'
    '      </ul>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">How to interpret it</p>\n'
    '      <ul style="margin:0;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li>Look for <b>red blocks along the diagonal</b> after hierarchical clustering:'
    ' these are groups of perturbations with shared downstream programmes and are good candidates'
    ' for being members of the same regulatory module or protein complex.</li>\n'
    '        <li>Look for <b>blue off-diagonal cells</b> linking two distant clusters: these'
    ' are candidate antagonistic pairs (one activator, one repressor) acting on the same set'
    ' of target genes.</li>\n'
    '        <li>If a perturbation has a strong red block in the left panel but not in the'
    ' right panel, its similarity to its neighbours is mostly explained by cell-state shifts'
    ' rather than direct regulation.</li>\n'
    '      </ul>\n'
    '    </div>\n'
    '    <div class="chart-grid">\n'
    '      <div class="chart" id="reg-corr-before-chart" style="height:500px;"></div>\n'
    '      <div class="chart" id="reg-corr-after-chart" style="height:500px;"></div>\n'
    '    </div>\n'
    '  </div>\n'

    '  <div id="reg-tf-pane">\n'
    '    <h3>TF Regulatory Heatmap</h3>\n'
    '    <div style="background:#f0f7ff;border-left:4px solid var(--accent);padding:14px 18px;'
    'border-radius:4px;margin-bottom:16px;">\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">What this heatmap shows</p>\n'
    '      <p style="margin:0 0 10px;font-size:13px;color:#3a4a58;line-height:1.7;">'
    'This panel is <b>not a correlation matrix</b>. Each cell summarises a single'
    ' differential-expression test: it asks how strongly, and in which direction, knocking'
    ' out the transcription factor in that row changes the expression of the target gene'
    ' in that column. Rows are perturbed transcription-factor guides. Columns are restricted'
    ' to the genes that were also perturbed in the screen, so the matrix captures'
    ' <b>regulation among the perturbed transcription factors themselves</b>,'
    ' that is, a transcription-factor-versus-transcription-factor wiring map.</p>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">How each cell is computed</p>\n'
    '      <p style="margin:0 0 10px;font-size:13px;color:#3a4a58;line-height:1.7;">'
    'For each (row, column) pair the pipeline reads the per-perturbation DEG table and'
    ' takes two numbers: the BH-adjusted p-value (the q-value) and the log2 fold-change.'
    ' The cell value is the <b>signed minus-log10 q-value</b>: the magnitude is the'
    ' negative of the base-10 logarithm of the q-value (so larger numbers mean stronger'
    ' statistical evidence), and the sign is taken from the direction of the log2 fold-change'
    ' (positive when the knockout raises the target, negative when it lowers it).</p>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">Colour legend</p>\n'
    '      <ul style="margin:0 0 10px;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li><span style="color:#B2182B;font-weight:700;">Red cells (positive value)</span>'
    ' mean the knockout <b>raises</b> the target gene\u2019s expression. Because losing the'
    ' factor leads to more of the target, the factor is interpreted as a <b>repressor</b>'
    ' of that gene under normal conditions.</li>\n'
    '        <li><span style="color:#2166AC;font-weight:700;">Blue cells (negative value)</span>'
    ' mean the knockout <b>lowers</b> the target gene\u2019s expression. Losing the factor'
    ' leads to less of the target, so the factor is interpreted as an <b>activator</b>'
    ' of that gene under normal conditions.</li>\n'
    '        <li><span style="color:#888;font-weight:700;">Near-white cells</span> indicate'
    ' that the knockout produced no statistically significant change in that target,'
    ' or that the change was small in magnitude and direction was uncertain.</li>\n'
    '      </ul>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">How to interpret it</p>\n'
    '      <ul style="margin:0;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li>Each <b>row</b> is the regulatory footprint of one transcription factor:'
    ' which other factors it normally represses (red entries) and which it normally activates'
    ' (blue entries) within the panel.</li>\n'
    '        <li>Each <b>column</b> shows which transcription factors converge on the same'
    ' target gene, so it gives a quick view of the upstream regulators of that target.</li>\n'
    '        <li>The <b>diagonal</b> serves as a sanity check: knocking out a transcription'
    ' factor should normally reduce its own transcript, producing a strong blue diagonal cell.</li>\n'
    '        <li>Pairs of cells that mirror each other across the diagonal (factor A on factor B,'
    ' and factor B on factor A) reveal feedback or feed-forward relationships: matching colours'
    ' suggest mutual reinforcement, opposite colours suggest a feedback loop.</li>\n'
    '        <li>Hierarchical clustering of rows and columns groups factors with similar'
    ' regulatory footprints, often recovering known co-regulators or members of the same'
    ' complex.</li>\n'
    '      </ul>\n'
    '    </div>\n'
    '    <div class="chart" id="reg-tf-heatmap-chart" style="height:820px;max-width:760px;margin:0 auto;"></div>\n'
    '  </div>\n'

    '  <div id="reg-net-pane">\n'
    '    <h3>TF Regulatory Network</h3>\n'
    '    <div style="background:#f0f7ff;border-left:4px solid var(--accent);padding:14px 18px;'
    'border-radius:4px;margin-bottom:16px;">\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">What this network shows</p>\n'
    '      <p style="margin:0 0 10px;font-size:13px;color:#3a4a58;line-height:1.7;">'
    'A directed graph where every node is a perturbed transcription factor and every edge'
    ' represents a statistically significant regulatory relationship inferred from the same'
    ' DEG tables that feed the heatmap above. An edge is drawn from factor A to factor B'
    ' when knocking out factor A produces a significant change in the expression of factor B'
    ' (BH-adjusted q-value below the FDR threshold and absolute log2 fold-change above the'
    ' magnitude threshold). The graph is the thresholded, directed version of the heatmap.</p>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">Edge legend</p>\n'
    '      <ul style="margin:0 0 10px;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li><span style="color:#1a9850;font-weight:700;">Green arrows (activating edges)</span>'
    ' indicate that knocking out the source factor <b>lowers</b> the target. The source factor'
    ' is therefore inferred to be a <b>positive regulator</b> of the target under normal'
    ' conditions.</li>\n'
    '        <li><span style="color:#B2182B;font-weight:700;">Red arrows (inhibiting edges)</span>'
    ' indicate that knocking out the source factor <b>raises</b> the target. The source factor'
    ' is therefore inferred to be a <b>negative regulator</b> (repressor) of the target under'
    ' normal conditions.</li>\n'
    '        <li>The <b>thickness of an arrow</b> is proportional to the absolute log2 fold-change'
    ' of the target gene when the source factor is knocked out, so thicker edges mean larger'
    ' regulatory effects.</li>\n'
    '      </ul>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">Node legend</p>\n'
    '      <ul style="margin:0 0 10px;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li>The <b>size of each node</b> reflects its total regulatory degree, that is,'
    ' the number of incoming edges plus the number of outgoing edges. Large nodes are highly'
    ' connected and are good candidate <b>master regulators</b> within the panel.</li>\n'
    '        <li>Nodes with many <b>outgoing</b> edges are upstream factors that influence many'
    ' others; nodes with many <b>incoming</b> edges are downstream targets converged on by'
    ' multiple regulators.</li>\n'
    '      </ul>\n'
    '      <p style="margin:0 0 8px;font-size:14px;font-weight:600;">How to interpret it</p>\n'
    '      <ul style="margin:0;padding-left:20px;font-size:13px;color:#3a4a58;line-height:1.9;">\n'
    '        <li>Look for <b>hub nodes</b> with many outgoing arrows: these are likely the'
    ' top of the regulatory hierarchy in the screened panel.</li>\n'
    '        <li><b>Reciprocal edges</b> between two nodes (an arrow in each direction) reveal'
    ' feedback loops. If both arrows are the same colour, the two factors mutually reinforce'
    ' each other; if the colours are opposite, they form a negative-feedback loop.</li>\n'
    '        <li>Use the threshold control below to make the network sparser (raise the cut-off)'
    ' to see only the strongest edges, or denser (lower the cut-off) to see weaker but'
    ' potentially biologically meaningful relationships.</li>\n'
    '      </ul>\n'
    '    </div>\n'
    '    <div class="control-panel">\n'
    '      <div><label>Minimum absolute log2 fold-change for edges</label>'
    '<input id="tfnet-lfc-thresh" type="number" value="0.3" min="0" step="0.1"'
    ' style="width:80px;" oninput="renderTFNetwork()"></div>\n'
    '    </div>\n'
    '    <div class="chart chart-full" id="reg-tfnet-chart" style="height:640px;"></div>\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_CSCORE = (
    '<div class="tab-pane" id="tab-cscore">\n'
    '  <h2>Connectivity Score (C-score)</h2>\n'
    '  <p>The Connectivity Score quantifies how much each perturbation rewires the gene co-expression network'
    ' relative to control. It decomposes into three orthogonal components: C_gain'
    ' (new edges gained), C_loss (edges lost), and C_shift'
    ' (weight change on preserved edges). C_total = C_gain + C_loss + C_shift.</p>\n'
    '  <div class="cards" id="cscore-cards" style="margin-bottom:22px;"></div>\n'
    '  <div style="display:grid;grid-template-columns:1fr;gap:22px;margin-bottom:22px;">\n'
    '    <div class="hm-card">\n'
    '      <h3 style="font-size:17px;margin-top:0;">Perturbations Ranked by C_total</h3>\n'
    '      <p style="font-size:13px;color:var(--muted);margin-bottom:8px;">'
    'Stacked bar: green\u202f=\u202fC_gain, red\u202f=\u202fC_loss, orange\u202f=\u202fC_shift.</p>\n'
    '      <img id="cscore-ranked-bar-img" style="width:100%;min-height:540px;object-fit:contain;border-radius:4px;" alt="C-score ranked bar">\n'
    '    </div>\n'
    '    <div class="hm-card">\n'
    '      <h3 style="font-size:17px;margin-top:0;">Edge Gain vs Loss</h3>\n'
    '      <p style="font-size:13px;color:var(--muted);margin-bottom:8px;">'
    'Grouped bars showing C_gain and C_loss per perturbation.</p>\n'
    '      <img id="cscore-module-img" style="width:100%;min-height:420px;object-fit:contain;border-radius:4px;" alt="Gain vs loss">\n'
    '    </div>\n'
    '  </div>\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr;gap:22px;margin-bottom:22px;">\n'
    '    <div class="hm-card">\n'
    '      <h3 style="font-size:17px;margin-top:0;">C-score Decomposition</h3>\n'
    '      <p style="font-size:13px;color:var(--muted);margin-bottom:8px;">'
    'Scatter: x\u202f=\u202fC_gain, y\u202f=\u202fC_loss. Dot size\u202f=\u202fC_shift. Color\u202f=\u202fC_total.</p>\n'
    '      <img id="cscore-decomp-img" style="width:100%;border-radius:4px;" alt="C-score decomposition">\n'
    '    </div>\n'
    '    <div class="hm-card">\n'
    '      <h3 style="font-size:17px;margin-top:0;">Network Rewiring vs DE Genes</h3>\n'
    '      <p style="font-size:13px;color:var(--muted);margin-bottom:8px;">'
    'C_total vs number of differentially expressed genes (Pearson r annotated).</p>\n'
    '      <img id="cscore-vs-deg-img" style="width:100%;border-radius:4px;" alt="C-score vs DEG">\n'
    '    </div>\n'
    '  </div>\n'
    '  <div style="display:grid;grid-template-columns:1fr 1fr;gap:22px;margin-bottom:22px;">\n'
    '    <div class="hm-card">\n'
    '      <h3 style="font-size:17px;margin-top:0;">Gene Hub Rewiring</h3>\n'
    '      <p style="font-size:13px;color:var(--muted);margin-bottom:8px;">'
    'Heatmap of \u0394degree (pert \u2212 ctrl) for top 25 hub genes. Red\u202f=\u202fgained connections, blue\u202f=\u202flost.</p>\n'
    '      <img id="cscore-hub-img" style="width:100%;min-height:620px;object-fit:contain;border-radius:4px;" alt="Hub heatmap">\n'
    '    </div>\n'
    '  </div>\n'
    '  <h3>Per-Perturbation Edge Explorer</h3>\n'
    '  <div class="control-panel" style="margin-bottom:10px;">\n'
    '    <div>\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Perturbation</label>\n'
    '      <select id="cscore-pert-selector" style="min-width:220px;"'
    ' onchange="renderCscoreEdges(this.value)"></select>\n'
    '    </div>\n'
    '    <div>\n'
    '      <label style="font-size:13px;color:var(--muted);display:block;margin-bottom:4px;">Gene</label>\n'
    '      <select id="cscore-gene-selector" style="min-width:220px;"'
    ' onchange="renderCscoreEdges(document.getElementById(\'cscore-pert-selector\').value,this.value)"></select>\n'
    '    </div>\n'
    '  </div>\n'
    '  <div class="chart" id="cscore-edge-scatter" style="height:420px;margin-bottom:16px;"></div>\n'
    '  <div class="tbl-wrap">\n'
    '    <table>\n'
    '      <thead><tr><th>Gene A</th><th>Gene B</th><th>Status</th>'
    '<th>ctrl r</th><th>pert r</th><th>\u0394r</th></tr></thead>\n'
    '      <tbody id="cscore-edge-tbody"></tbody>\n'
    '    </table>\n'
    '  </div>\n'
    '</div>\n'
)

_TAB_PARAMS = (
    '<div class="tab-pane" id="tab-params">\n'
    '  <h2>Analysis Parameters</h2>\n'
    '  <p>Pipeline settings are grouped by analysis step so the report reflects the multi-step PerturbFlow-Analyzer workflow.</p>\n'
    '  <h3>Pipeline Steps</h3>\n'
    '  <div class="tbl-wrap" style="margin-top:12px;">\n'
    '    <table class="param-table">\n'
    '      <thead><tr><th>Step</th><th>Status</th><th>Purpose</th><th>Key parameters</th></tr></thead>\n'
    '      <tbody id="step-param-tbody"></tbody>\n'
    '    </table>\n'
    '  </div>\n'
    '  <h3>All Run Parameters</h3>\n'
    '  <div class="tbl-wrap" style="margin-top:18px;">\n'
    '    <table class="param-table">\n'
    '      <thead><tr><th>Parameter</th><th>Value</th><th>Description</th></tr></thead>\n'
    '      <tbody id="param-tbody"></tbody>\n'
    '    </table>\n'
    '  </div>\n'
    '</div>\n'
)

_JS_BODY = r"""
const PARAM_DESC = {
  min_genes:"Minimum genes per cell for QC pass",
  max_pct_mt:"Maximum % mitochondrial reads for QC pass",
  min_cells_per_perturbation:"Minimum cells to retain a perturbation",
  random_state:"Random seed for UMAP/Leiden reproducibility",
  n_neighbors:"k-NN graph neighbors (UMAP/scoring)",
  leiden_resolution:"Leiden clustering resolution",
  eda_n_top_genes:"Top genes in EDA heatmaps",
  eda_max_cells_heatmap:"Max cells for gene\u00d7cell heatmap",
  deg_n_top_perturbations:"Perturbations analyzed in DEG step",
  deg_logfc_threshold:"|log\u2082FC| cutoff for significance",
  deg_pval_threshold:"Adjusted p-value cutoff",
  deg_n_top_heatmap:"Top DEGs per perturbation in heatmap",
  genenet_n_top_genes:"HVGs for co-expression network",
  genenet_n_gene_clusters:"Gene module clusters",
  genenet_corr_threshold:"Pearson correlation threshold for edges",
  state_enrich_min_cells:"Minimum cells in perturbation group to test",
  state_enrich_fdr_threshold:"BH-FDR threshold for state enrichment significance",
  regulatory_fdr_threshold:"FDR cutoff for TF regulatory edges",
  regulatory_lfc_threshold:"|log\u2082FC| cutoff for TF regulatory edges",
};
const STEP_PARAMS = [
  {step:'qc', purpose:'Filter low-quality cells and sparse perturbation groups.', params:['min_genes','max_pct_mt','min_cells_per_perturbation']},
  {step:'preprocess', purpose:'Normalize counts, select features, embed cells, and assign cell-state clusters.', params:['random_state','n_neighbors','leiden_resolution']},
  {step:'eda', purpose:'Create exploratory cell-count, gene-expression, correlation, and similarity heatmaps.', params:['eda_n_top_genes','eda_max_cells_heatmap']},
  {step:'score', purpose:'Estimate perturbation burden, guide confidence, target reduction, and perturbation score.', params:['n_neighbors','n_top_genes_signature']},
  {step:'effects', purpose:'Separate transcriptional response from cell-state shift effects.', params:['n_neighbors']},
  {step:'trajectory', purpose:'Quantify pseudotime and commitment-shift behavior after perturbation.', params:['random_state']},
  {step:'programs', purpose:'Infer perturbation-associated gene programs.', params:['n_top_genes_signature']},
  {step:'interaction', purpose:'Summarize interaction scores for combinatorial perturbation labels when present.', params:[]},
  {step:'deg', purpose:'Identify differential expression against control and produce volcano plots and DEG tables.', params:['deg_n_top_perturbations','deg_logfc_threshold','deg_pval_threshold','deg_n_top_heatmap']},
  {step:'genenet', purpose:'Compare control and perturbation gene co-expression network structure.', params:['genenet_n_top_genes','genenet_n_gene_clusters','genenet_corr_threshold']},
  {step:'state_enrich', purpose:'Chi-square test of cell state enrichment/depletion for each perturbation vs control.', params:['state_enrich_min_cells','state_enrich_fdr_threshold']},
  {step:'regulatory', purpose:'Perturbation effect correlation and TF-TF regulatory heatmap and network.', params:['regulatory_fdr_threshold','regulatory_lfc_threshold']},
  {step:'report', purpose:'Write static and interactive web reports.', params:['eda_n_top_genes']},
  {step:'bundle', purpose:'Emit viewer-ready parquet and JSON artifacts.', params:['bundle_top_de_per_pert','bundle_min_de_hits_per_gene','bundle_schema_version']},
];

const _tabInit = {};
let _dataReady=null;
function emptyDataStub(){
  return {
    summary:{n_cells:0,n_genes:0,n_perturbations:0},
    perturbations:[],cells_per_pert:[],qc_cells:[],effect_df:[],
    deg_summary:[],deg:{},gsea:{},enrichment:{},umap:[],genes:[],
    gene_expr:{},umap_genes:[],umap_gene_expr:{},heatmap_b64:{},
    genenet_b64:{},eda_corr_pert:{},pert_sim:{},genenet_json:{},
    state_enrich:{},pert_corr_before:{},pert_corr_after:{},
    tf_regulatory:{},tf_network:{},config:{_completed_steps:[]},
    cscore:{summary:[],edges:{},b64:{}}
  };
}
function loadReportData(){
  if(_dataReady) return _dataReady;
  const hdr=document.getElementById('hdr-summary');
  if(hdr) hdr.textContent='Loading results data...';
  _dataReady=fetch('interactive_data.json', {cache:'no-store'})
    .then(r=>{
      if(!r.ok) throw new Error('Could not load interactive_data.json: HTTP '+r.status);
      return r.json();
    })
    .then(data=>{
      D=Object.assign(emptyDataStub(), data||{});
      if(hdr) hdr.textContent='Results data loaded';
      return D;
    })
    .catch(err=>{
      showStartupError(new Error(String(err && err.message || err) + '. If you opened the HTML directly with file://, serve the results folder with a small web server, for example: python -m http.server 8899'));
      throw err;
    });
  return _dataReady;
}
let _plotlyLoading=false;
let _plotlyLoaded=false;
let _plotlyCallbacks=[];
function ensurePlotly(cb){
  if(typeof Plotly!=='undefined'){_plotlyLoaded=true;cb();return;}
  _plotlyCallbacks.push(cb);
  if(_plotlyLoading) return;
  _plotlyLoading=true;
  const s=document.createElement('script');
  s.src='plotly.min.js';
  s.async=true;
  s.onload=()=>{
    _plotlyLoaded=true;
    const q=_plotlyCallbacks.splice(0);
    q.forEach(fn=>{try{fn();}catch(e){showStartupError(e);}});
  };
  s.onerror=()=>showStartupError(new Error('Could not load plotly.min.js next to interactive_report.html'));
  document.head.appendChild(s);
}
function showStartupError(err){
  const msg=(err&&err.stack)||String(err||'Unknown error');
  const hdr=document.getElementById('hdr-summary');
  if(hdr) hdr.textContent='Report error';
  const target=document.querySelector('.content')||document.body;
  const old=document.getElementById('startup-error');
  if(old) old.remove();
  target.insertAdjacentHTML('afterbegin',
    '<div id="startup-error" class="metric-help" style="border-color:#b04a5a;background:#fff6f6;">'
    +'<b>Report JavaScript error</b><pre style="white-space:pre-wrap;margin-top:8px;font-size:12px;color:#7a1f2c;">'
    +msg.replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]))+'</pre></div>');
}
function showTab(name, btnEl, src) {
  document.querySelectorAll('.tab-pane').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.sidebar button').forEach(el => el.classList.remove('active'));
  const pane = document.getElementById('tab-' + name);
  if (pane) pane.classList.add('active');
  document.querySelectorAll('.tab-btn').forEach(b => {
    if ((b.dataset && b.dataset.tab === name) || (b.getAttribute('onclick')||'').includes("'"+name+"'")) b.classList.add('active');
  });
  document.querySelectorAll('.sidebar button').forEach(b => {
    if ((b.dataset && b.dataset.tab === name) || (b.getAttribute('onclick')||'').includes("'"+name+"'")) b.classList.add('active');
  });
  if (!_tabInit[name]) { _tabInit[name] = true; _initTab(name); }
}
function bindTabButtons(){
  document.querySelectorAll('button[onclick^="showTab("]').forEach(btn=>{
    const raw = btn.getAttribute('onclick') || '';
    const m = raw.match(/showTab\(['"]([^'"]+)['"]/);
    if(!m) return;
    btn.dataset.tab = m[1];
    if(btn.dataset.boundTab === '1') return;
    btn.dataset.boundTab = '1';
    btn.addEventListener('click', ev=>{
      ev.preventDefault();
      showTab(m[1], btn, btn.classList.contains('tab-btn') ? 'tab' : 'sidebar');
    });
  });
}
function _initTab(n) {
  if (n==='home')       renderHome();
  if (n==='qc')         ensurePlotly(renderQC);
  if (n==='heatmaps')   renderHeatmaps();
  if (n==='pert')       ensurePlotly(initPertTab);
  if (n==='network')    ensurePlotly(initNetworkTab);
  if (n==='gene')       ensurePlotly(initGeneTab);
  if (n==='states')     ensurePlotly(renderStateEnrich);
  if (n==='regulatory') ensurePlotly(initRegulatoryTab);
  if (n==='cscore')     ensurePlotly(initCscoreTab);
  if (n==='params')     renderParams();
}

const BL = {
  paper_bgcolor:'#ffffff', plot_bgcolor:'#ffffff',
  font:{color:'#203040',size:15,family:'Arial,Helvetica,sans-serif'},
  margin:{t:50,r:20,b:60,l:60},
  colorway:['#1f4e79','#169b8f','#8c6bb1','#c08a2b','#5f8f3f','#b04a5a','#4c78a8','#72b7b2','#9d755d'],
  hoverlabel:{font:{size:13}},
  xaxis:{zeroline:false,linecolor:'#4c5f70',mirror:false,ticks:'outside',showgrid:false},
  yaxis:{zeroline:false,linecolor:'#4c5f70',mirror:false,ticks:'outside',showgrid:false},
};
function mkL(o){return Object.assign({},BL,o);}
const PC = {responsive:true,displayModeBar:true,
  // Built-in modebar camera-icon download → 4× DPI for print-quality PNGs.
  toImageButtonOptions:{format:'png',scale:4,height:1100,width:1400}};
// Publication-quality diverging colorscale (ColorBrewer RdBu):
// 0 = blue (depleted / negative) → 0.5 = white (neutral) → 1 = red (enriched / positive)
const DIVG = [[0,'#2166AC'],[0.2,'#4393C3'],[0.35,'#92C5DE'],[0.45,'#D1E5F0'],[0.5,'#F7F7F7'],[0.55,'#FDDBC7'],[0.65,'#F4A582'],[0.8,'#D6604D'],[1,'#B2182B']];
const QC_METRICS = [
  {key:'n_genes_by_counts',label:'nFeature (genes/cell)'},
  {key:'n_genes',label:'n_genes'},
  {key:'total_counts',label:'nCount (total UMI)'},
  {key:'UMI_count',label:'UMI count (raw)'},
  {key:'core_adjusted_UMI_count',label:'Adjusted UMI'},
  {key:'z_gemgroup_UMI',label:'z-score UMI'},
  {key:'core_scale_factor',label:'Core scale factor'},
  {key:'pct_counts_mt',label:'% MT reads'},
  {key:'percent_mito',label:'% Mito'},
  {key:'perturbation_score',label:'Perturbation score'},
  {key:'guide_confidence_score',label:'Guide confidence'},
  {key:'target_expr_reduction',label:'Target expr reduction'},
  {key:'perturbation_burden',label:'Perturbation burden'},
  {key:'escape_probability',label:'Escape probability'}
];
function fmtNum(v){
  if(v==null || !isFinite(v)) return 'NA';
  const a=Math.abs(v);
  if(a>=1000) return Math.round(v).toLocaleString();
  if(a>=10) return v.toFixed(1);
  if(a>=1) return v.toFixed(2);
  return v.toFixed(3);
}
function median(vals){
  if(!vals.length) return null;
  const s=[...vals].sort((a,b)=>a-b);
  const m=Math.floor(s.length/2);
  return s.length%2?s[m]:(s[m-1]+s[m])/2;
}
function getQCSummaryCards(){
  const nPerts=D.perturbations?D.perturbations.length:D.summary.n_perturbations;
  const nCtrl=D.cells_per_pert.filter(r=>r.is_control).reduce((a,r)=>a+r.n_cells,0);
  const nPertCells=D.cells_per_pert.filter(r=>!r.is_control).reduce((a,r)=>a+r.n_cells,0);
  return [
    {l:'Total Cells',v:D.summary.n_cells.toLocaleString()},
    {l:'Genes Measured',v:D.summary.n_genes.toLocaleString()},
    {l:'Perturbations (n)',v:nPerts.toLocaleString()},
    {l:'Control Cells',v:nCtrl.toLocaleString()},
    {l:'Perturbed Cells',v:nPertCells.toLocaleString()}
  ];
}
function getQCMetricCards(){
  return QC_METRICS.map(({key,label})=>{
    const vals=D.qc_cells.map(r=>r[key]).filter(v=>v!=null&&isFinite(v));
    if(!vals.length) return null;
    return {l:label,v:fmtNum(median(vals)),c:'',sub:'median'};
  }).filter(Boolean);
}
function renderCardSet(id,cards){
  const el=document.getElementById(id);
  if(!el) return;
  el.innerHTML = cards.map(c=>
    '<div class="card"><div class="label">'+c.l+'</div><div class="value '+(c.c||'')+'">'+c.v+'</div>'
    +(c.sub?'<div class="label" style="margin-top:4px;">'+c.sub+'</div>':'')+'</div>').join('');
}

/* HOME */
function renderHome(){
  const s = D.summary;
  document.getElementById('hdr-summary').textContent =
    s.n_cells.toLocaleString()+' cells \u00b7 '+s.n_genes.toLocaleString()+' genes \u00b7 '+s.n_perturbations.toLocaleString()+' perturbations';
  const cards = [{l:'Cells',v:s.n_cells.toLocaleString(),c:''},{l:'Genes',v:s.n_genes.toLocaleString(),c:''},
    {l:'Perturbations',v:s.n_perturbations.toLocaleString(),c:''}];
  document.getElementById('home-cards').innerHTML = cards.map(c=>
    '<div class="card"><div class="label">'+c.l+'</div><div class="value '+c.c+'">'+c.v+'</div></div>').join('');
  document.getElementById('home-analysis-text').textContent =
    'This dataset contains '+s.n_cells.toLocaleString()+' single cells with '+s.n_perturbations.toLocaleString()+
    ' perturbations and '+s.n_genes.toLocaleString()+' measured genes. '+
    (D.deg_summary.length?'DEG analysis was performed on '+D.deg_summary.length+' perturbations. ':'')+
    'Use the tabs above to explore QC, UMAP, heatmaps, per-perturbation DEGs, and gene expression.';
  const ALL_STEPS=['qc','preprocess','eda','score','effects','trajectory','programs','interaction','state_enrich','deg','genenet','regulatory','report','bundle'];
  const done = new Set(D.config._completed_steps||[]);
  document.getElementById('home-steps').innerHTML = ALL_STEPS.map(st=>
    '<span class="step-chip '+(done.has(st)?'done':'')+'">'+( done.has(st)?'\u2713 ':'')+st+'</span>').join('');
  ['home-top-perts-chart','home-deg-bar','home-effect-scatter'].forEach(id=>{
    const el=document.getElementById(id);
    if(el) el.innerHTML='<p style="padding:18px;color:var(--muted)">Loading chart…</p>';
  });
  const gl=document.getElementById('umap-gene-list');
  if(gl) gl.innerHTML=(D.umap_genes||D.genes||[]).map(g=>'<option value="'+g+'"></option>').join('');
  renderStaticUmap();
  setTimeout(renderHomeTopPertChart, 50);
  setTimeout(renderHomeDegChart, 120);
  setTimeout(renderHomeEffectScatter, 190);
}

function renderHomeTopPertChart(){
  if(typeof Plotly==='undefined'){ensurePlotly(renderHomeTopPertChart);return;}
  const sd = D.deg_summary.length
    ? [...D.deg_summary].sort((a,b)=>(b.n_de_total||0)-(a.n_de_total||0)).slice(0,20)
    : [...D.cells_per_pert].filter(r=>!r.is_control).sort((a,b)=>b.n_cells-a.n_cells).slice(0,20);
  Plotly.newPlot('home-top-perts-chart',[{type:'bar',orientation:'h',width:0.38,
    x:sd.map(r=>r.n_de_total??r.n_cells),y:sd.map(r=>r.perturbation),
    marker:{color:'#1f4e79'},hovertemplate:'%{y}: %{x}<extra></extra>'}],
      mkL({title:{text:'Top Perturbations (Total)',font:{size:17}},
         xaxis:{title:{text:D.deg_summary.length?'DEG count':'Cell count',font:{size:14}}},
         yaxis:{automargin:true,tickfont:{size:13},autorange:'reversed'},
         margin:{t:50,r:20,b:60,l:200}}),PC);
}

function renderHomeDegChart(){
  if(typeof Plotly==='undefined'){ensurePlotly(renderHomeDegChart);return;}
  if(D.deg_summary.length){
    const td=[...D.deg_summary].sort((a,b)=>(b.n_de_total||0)-(a.n_de_total||0)).slice(0,15);
    Plotly.newPlot('home-deg-bar',[
      {name:'Up',type:'bar',orientation:'h',width:0.38,x:td.map(r=>r.n_de_up||0),y:td.map(r=>r.perturbation),marker:{color:'#b04a5a'}},
      {name:'Down',type:'bar',orientation:'h',width:0.38,x:td.map(r=>-(r.n_de_down||0)),y:td.map(r=>r.perturbation),marker:{color:'#1f4e79'}}],
      mkL({title:{text:'DEG Directional',font:{size:17}},barmode:'relative',
           xaxis:{title:{text:'Gene count',font:{size:14}}},yaxis:{automargin:true,tickfont:{size:13},autorange:'reversed'},
           margin:{t:50,r:20,b:60,l:200},showlegend:true}),PC);
  } else {
    const el=document.getElementById('home-deg-bar');
    if(el) el.innerHTML='<p style="padding:18px;color:var(--muted)">DEG summary was not found.</p>';
  }
}

function renderHomeEffectScatter(){
  if(typeof Plotly==='undefined'){ensurePlotly(renderHomeEffectScatter);return;}
  if(D.effect_df.length){
    Plotly.newPlot('home-effect-scatter',[{type:'scatter',mode:'markers+text',
      x:D.effect_df.map(r=>r.transcriptional_score),y:D.effect_df.map(r=>r.state_shift_score),
      text:D.effect_df.map(r=>r.perturbation),textposition:'top center',textfont:{size:12},
      hovertemplate:'<b>%{text}</b><br>Trans: %{x:.3f}<br>State: %{y:.3f}<extra></extra>',
      marker:{color:'#8c6bb1',size:9,opacity:0.85}}],
      mkL({title:{text:'Effect Decomposition',font:{size:17}},
           xaxis:{title:{text:'Transcriptional score',font:{size:14}}},
           yaxis:{title:{text:'State shift score',font:{size:14}}},
           width:460,height:460,autosize:false}),PC);
  }
}

/* QC */
function renderQC(){
  renderCardSet('qc-summary-cards', getQCSummaryCards());
  renderCardSet('qc-metric-cards', getQCMetricCards());
  const cpp=[...D.cells_per_pert].sort((a,b)=>b.n_cells-a.n_cells).slice(0,40);
  Plotly.newPlot('qc-cells-bar',[{type:'bar',width:0.4,x:cpp.map(r=>r.perturbation),y:cpp.map(r=>r.n_cells),
    marker:{color:cpp.map(r=>r.is_control?'#169b8f':'#1f4e79')},
    hovertemplate:'%{x}: %{y} cells<extra></extra>'}],
    mkL({title:{text:'Cells per Perturbation (top 40)',font:{size:17}},
         xaxis:{tickangle:-35,tickfont:{size:11}},yaxis:{title:{text:'Cell count',font:{size:14}}}}),PC);
  const grid=document.getElementById('qc-violin-grid');
  grid.innerHTML='';
  QC_METRICS.forEach(({key,label})=>{
    const vals=D.qc_cells.map(r=>r[key]).filter(v=>v!=null&&isFinite(v));
    if(!vals.length) return;
    const id='qc-v-'+key;
    const div=document.createElement('div');
    div.className='chart';div.id=id;div.style.height='280px';
    grid.appendChild(div);
    Plotly.newPlot(id,[{type:'violin',y:vals,box:{visible:true},meanline:{visible:true},
      fillcolor:'rgba(31,78,121,.22)',line:{color:'#1f4e79'},name:label,points:false}],
      mkL({title:{text:label,font:{size:16}},yaxis:{title:{text:label,font:{size:13}}},
           showlegend:false,margin:{t:45,r:15,b:40,l:65}}),PC);
  });
}

/* UMAP */
function onUmapColorChange(){
  const v=document.getElementById('umap-color-by').value;
  document.getElementById('umap-gene-wrap').style.display=v==='gene'?'block':'none';
  const gi=document.getElementById('umap-gene-input');
  if(v==='gene' && gi && !gi.value && (D.umap_genes||[]).length) gi.value=D.umap_genes[0];
  renderUMAP();
}
let _umapDebounce=null;
function debounceRenderUMAP(){
  clearTimeout(_umapDebounce);
  _umapDebounce=setTimeout(renderUMAP, 250);
}
function renderStaticUmap(){
  const el=document.getElementById('umap-chart');
  if(!el) return;
  const img=D.heatmap_b64.umap_pert||D.heatmap_b64.umap_state;
  const btn='<div style="text-align:center;margin-top:10px;">'
    +'<button onclick="renderUMAP()" style="background:var(--accent);border:1px solid var(--accent);border-radius:6px;color:#fff;padding:8px 18px;cursor:pointer;font-size:13px;">Load interactive UMAP</button>'
    +'</div>';
  if(img){
    el.innerHTML='<div style="height:100%;display:flex;flex-direction:column;align-items:center;justify-content:center;padding:12px;">'
      +'<img src="'+img+'" alt="UMAP" style="max-width:100%;max-height:490px;object-fit:contain;">'
      +btn+'</div>';
  } else {
    el.innerHTML='<p style="padding:20px;color:var(--muted)">Interactive UMAP is ready but not loaded yet.</p>'+btn;
  }
}
function renderUMAP(){
  if(typeof Plotly==='undefined'){ensurePlotly(renderUMAP);return;}
  const pts=D.umap;
  if(!pts.length){
    const el=document.getElementById('umap-chart');
    const colorBy=document.getElementById('umap-color-by').value;
    const b64=colorBy==='state'?(D.heatmap_b64.umap_state||D.heatmap_b64.umap_pert):D.heatmap_b64.umap_pert;
    if(b64){el.innerHTML='<img src="'+b64+'" style="max-width:100%;max-height:550px;display:block;margin:auto;" alt="UMAP">';}
    else{el.innerHTML='<p style="padding:20px;color:var(--muted)">No UMAP data (X_umap not found in h5ad).</p>';}
    return;
  }
  const colorBy=document.getElementById('umap-color-by').value;
  const hlPert=(document.getElementById('umap-pert-highlight')?.value||'').trim().toLowerCase();
  let traces=[];
  if(colorBy==='score'){
    traces=[{type:'scattergl',mode:'markers',x:pts.map(p=>p.x),y:pts.map(p=>p.y),
      text:pts.map(p=>p.pert),hovertemplate:'%{text}<br>(%{x:.2f},%{y:.2f})<extra></extra>',
      name:'Perturbation score',showlegend:true,
      marker:{color:pts.map(p=>p.score??0),colorscale:'Viridis',showscale:true,size:8,opacity:0.85,
              colorbar:{title:{text:'Perturbation score',font:{size:14}},tickfont:{size:12},x:1.02}}}}];
  } else if(colorBy==='gene'){
    let gName=document.getElementById('umap-gene-input').value.trim();
    if(!gName && (D.umap_genes||[]).length){
      gName=D.umap_genes[0];
      document.getElementById('umap-gene-input').value=gName;
    }
    const gVals=(D.umap_gene_expr&&D.umap_gene_expr[gName])||pts.map(()=>0);
    traces=[{type:'scattergl',mode:'markers',x:pts.map(p=>p.x),y:pts.map(p=>p.y),
      text:pts.map(p=>p.pert),hovertemplate:gName+': %{marker.color:.3f}<br>%{text}<extra></extra>',
      name:gName||'Gene',showlegend:false,
      marker:{color:gVals,colorscale:'Viridis',showscale:true,size:7,opacity:0.85,
              colorbar:{title:{text:gName||'Expression',font:{size:14}},tickfont:{size:12},x:1.02}}}];
  } else {
    const colorKey=colorBy==='state'?'state':'pert';
    const groups={};
    pts.forEach(p=>{const g=p[colorKey]??'unknown';if(!groups[g])groups[g]={x:[],y:[],t:[]};
      groups[g].x.push(p.x);groups[g].y.push(p.y);groups[g].t.push(p.pert);});
    // Plot order: control first so non-control sits on top.
    const ctrlKeys={control:1,ctrl:1,nontargeting:1,'non-targeting':1,nt:1,scramble:1,'safe-targeting':1,'safe_targeting':1};
    const grpEntries=Object.entries(groups).sort((a,b)=>{
      const aCtrl=ctrlKeys[String(a[0]).toLowerCase()]?1:0;
      const bCtrl=ctrlKeys[String(b[0]).toLowerCase()]?1:0;
      return bCtrl-aCtrl;
    });
    const C=['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'];
    grpEntries.forEach(([grp,g],i)=>{
      const isHL=hlPert&&grp.toLowerCase().includes(hlPert);
      const notHL=hlPert&&!isHL;
      const isCtrl=!!ctrlKeys[String(grp).toLowerCase()];
      const baseColor=isCtrl?'#9aa0a6':C[i%C.length];
      traces.push({type:'scattergl',mode:'markers',name:grp,x:g.x,y:g.y,text:g.t,
        hovertemplate:'%{text}<br>(%{x:.2f},%{y:.2f})<extra></extra>',
        marker:{color:notHL?'rgba(160,170,180,.45)':baseColor,
                size:notHL?6:10, opacity:notHL?0.22:0.85,
                line:{width:0.4,color:'#1d3557'}},
        showlegend:true});
    });
  }
  const plotTitle = '';
  Plotly.react('umap-chart',traces,mkL({
    title:{text:plotTitle,font:{size:24}},
    xaxis:{title:{text:'UMAP 1',font:{size:18}},zeroline:false,showgrid:false,linecolor:'#4c5f70',ticks:'outside',tickfont:{size:14},
           scaleanchor:'y',scaleratio:1},
    yaxis:{title:{text:'UMAP 2',font:{size:18}},zeroline:false,showgrid:false,linecolor:'#4c5f70',ticks:'outside',tickfont:{size:14}},
    showlegend:true,
    legend:{font:{size:12},tracegroupgap:4,
            x:1.02,y:1.0,xanchor:'left',yanchor:'top',
            bgcolor:'rgba(255,255,255,0.85)',
            bordercolor:'#888',borderwidth:1},
    margin:{t:45,r:colorBy==='gene'?90:170,b:60,l:70}}),PC);
}

/* HEATMAPS */
function renderHeatmaps(){
  const titles={gene_by_cell:'Gene \u00d7 Cell Expression',gene_by_pert:'All Genes \u00d7 Perturbations Clusters',
    gene_corr:'Gene\u2013Gene Correlation (all cells)',pert_sim:'Perturbation\u2013Perturbation Similarity',
    cluster_summary:'Gene Expression Cluster \u00d7 Perturbation Cluster Summary'};
  const grid=document.getElementById('heatmaps-grid');
  grid.innerHTML='';
  if(!D.heatmap_b64||!Object.keys(D.heatmap_b64).length){
    grid.innerHTML='<p style="color:var(--muted);padding:16px;">Heatmap images not found. Run the eda step.</p>';return;
  }
  Object.entries(titles).forEach(([key,title])=>{
    const b64=D.heatmap_b64[key];if(!b64)return;
    const d=document.createElement('div');d.className=(key==='gene_by_cell'||key==='gene_by_pert')?'hm-card hm-wide':'hm-card';
    if(key==='pert_sim') d.style.cssText='max-width:640px;justify-self:start;';
    d.innerHTML='<h3 style="font-size:17px;">'+title+'</h3>'
      +'<button onclick="dlHeatmapImg(\''+key+'\',\''+title+'\')" '
      +'style="font-size:13px;background:var(--tab-active);border:1px solid var(--border);'
      +'border-radius:5px;color:var(--accent);padding:5px 14px;cursor:pointer;margin-bottom:10px;">'
      +'\u2b07 Download PNG</button>'
      +'<img src="'+b64+'" alt="'+title+'" style="width:100%;border-radius:4px;">';
    grid.appendChild(d);
  });
  // Per-perturbation correlation selector
  const corrSel=document.getElementById('corr-pert-selector');
  const corrPerts=Object.keys(D.eda_corr_pert||{}).sort();
  corrSel.innerHTML='';
  if(corrPerts.length){
    corrPerts.forEach(p=>{const o=document.createElement('option');o.value=p;o.textContent=p;corrSel.appendChild(o);});
    renderCorrPert(corrPerts[0]);
  } else {
    document.getElementById('corr-pert-img-wrap').innerHTML='<p style="color:var(--muted);padding:12px;">Per-perturbation correlation plots not found. Run the full pipeline.</p>';
  }
}
function getPertSimData(){
  if(D.pert_sim&&D.pert_sim.perts&&D.pert_sim.matrix) return D.pert_sim;
  if(D.pert_corr_after&&D.pert_corr_after.perts&&D.pert_corr_after.matrix){
    return {perts:D.pert_corr_after.perts,matrix:D.pert_corr_after.matrix,label:'Effect correlation'};
  }
  return null;
}
function initPertSimSelector(){
  const sim=getPertSimData(),sel=document.getElementById('pert-sim-selector');
  if(!sim||!sel)return;
  sel.innerHTML='';
  sim.perts.forEach(p=>{const o=document.createElement('option');o.value=p;o.textContent=p;sel.appendChild(o);});
  renderPertSim(sim.perts[0]);
}
function renderPertSim(selected){
  const sim=getPertSimData();if(!sim)return;
  const labels=sim.perts,idx=labels.indexOf(selected);
  const shapes=idx>=0?[
    {type:'rect',xref:'x',yref:'y',x0:idx-0.5,x1:idx+0.5,y0:-0.5,y1:labels.length-0.5,line:{color:'#b04a5a',width:3},fillcolor:'rgba(0,0,0,0)'},
    {type:'rect',xref:'x',yref:'y',x0:-0.5,x1:labels.length-0.5,y0:idx-0.5,y1:idx+0.5,line:{color:'#b04a5a',width:3},fillcolor:'rgba(0,0,0,0)'}
  ]:[];
  Plotly.react('pert-sim-plot',[{type:'heatmap',z:sim.matrix,x:labels,y:labels,zmin:-1,zmax:1,
    colorscale:'RdBu',reversescale:false,colorbar:{title:{text:sim.label||'Similarity'}},
    hovertemplate:'%{y} vs %{x}<br>'+(sim.label||'Similarity')+': %{z:.3f}<extra></extra>'}],
    mkL({title:{text:'Perturbation Similarity'+(selected?' — '+selected:'')},
         xaxis:{tickangle:-45,automargin:true,scaleanchor:'y',scaleratio:1},
         yaxis:{automargin:true,autorange:'reversed'},
         shapes:shapes,width:720,height:720,autosize:false,margin:{t:60,r:90,b:130,l:130}}),PC);
}
function renderCorrPert(pert){
  const b64=(D.eda_corr_pert||{})[pert];
  const img=document.getElementById('corr-pert-img');
  if(img)img.src=b64||'';
}
function dlHeatmapImg(key,title){
  const a=document.createElement('a');a.href=D.heatmap_b64[key];
  a.download=title.replace(/[\s\u00d7\u2013]/g,'_')+'.png';document.body.appendChild(a);a.click();document.body.removeChild(a);
}

/* PERTURBATION TAB */
const TOP_N_PERT=10;
let _pertExp=false,_degData=[],_degFilter='',_degSortCol=1,_degSortAsc=false,_degPage=0,_degPageSize=10,_curPert='';

function initPertTab(){
  const allPerts=Object.keys(D.deg);
  const sm={};D.deg_summary.forEach(r=>{sm[r.perturbation]=r.n_de_total||0;});
  const sp=allPerts.slice().sort((a,b)=>(sm[b]||0)-(sm[a]||0));
  buildPertDropdown(sp,_pertExp);
  document.getElementById('pert-count-badge').textContent=allPerts.length+' total';
  if(sp.length)renderPert(sp[0]);
}

function buildPertDropdown(sp,showAll){
  const sel=document.getElementById('pert-selector');
  const btn=document.getElementById('pert-show-all-btn');
  const toShow=showAll?sp:sp.slice(0,TOP_N_PERT);
  const prev=sel.value;
  sel.innerHTML='';
  const grp=document.createElement('optgroup');
  grp.label=showAll?'All perturbations':'Top '+Math.min(TOP_N_PERT,sp.length)+' by DEG count';
  toShow.forEach(p=>{const o=document.createElement('option');o.value=p;o.textContent=p;if(p===prev)o.selected=true;grp.appendChild(o);});
  sel.appendChild(grp);
  if(!showAll&&sp.length>TOP_N_PERT){btn.style.display='';btn.textContent='Show all ('+sp.length+') \u25bc';}
  else btn.style.display='none';
  if(!sel.value&&toShow.length)sel.value=toShow[0];
}

function expandPertDropdown(){
  _pertExp=true;
  const sp=Object.keys(D.deg).sort((a,b)=>{const sm={};D.deg_summary.forEach(r=>{sm[r.perturbation]=r.n_de_total||0;});return(sm[b]||0)-(sm[a]||0);});
  buildPertDropdown(sp,true);
}

function renderPert(pert){
  if(!pert)return;
  _curPert=pert;
  const rows=D.deg[pert]||[];
  // volcano
  const vx=rows.map(r=>r.log2fc??0),vy=rows.map(r=>r.neg_log10_padj??0);
  const vc=rows.map(r=>r.significant?(r.log2fc??0)>0?'#b04a5a':'#1f4e79':'#9aa7b4');
  const fmean=vx.reduce((a,b)=>a+b,0)/(vx.length||1);
  const vs=[...vx].sort((a,b)=>a-b);
  const fmed=vs[Math.floor(vs.length/2)]??0;
  const ymax=Math.max(...vy,1)*1.05;
  Plotly.react('pert-volcano',[
    {type:'scatter',mode:'markers',name:'Genes',x:vx,y:vy,text:rows.map(r=>r.gene),
     hovertemplate:'<b>%{text}</b><br>log\u2082FC: %{x:.3f}<br>-log\u2081\u2080(padj): %{y:.3f}<extra></extra>',
     marker:{color:vc,size:6,opacity:0.75}},
    {type:'scatter',mode:'lines',name:'Mean FC ('+fmean.toFixed(2)+')',x:[fmean,fmean],y:[0,ymax],
     line:{color:'#c08a2b',dash:'dash',width:1.5},hoverinfo:'none'},
    {type:'scatter',mode:'lines',name:'Median FC ('+fmed.toFixed(2)+')',x:[fmed,fmed],y:[0,ymax],
     line:{color:'#169b8f',dash:'dot',width:1.5},hoverinfo:'none'}],
    mkL({title:{text:'Volcano \u2014 '+pert,font:{size:18}},
         xaxis:{title:{text:'log\u2082 Fold Change',font:{size:15}}},
         yaxis:{title:{text:'-log\u2081\u2080(adj. p-value)',font:{size:15}}},
         height:500,autosize:true,
         showlegend:true,legend:{font:{size:13}},margin:{t:55,r:25,b:65,l:80}}),PC);
  // UMAP highlight
  if(D.umap.length){
    const ga=D.umap.filter(p=>p.pert===pert),gb=D.umap.filter(p=>p.pert!==pert);
    Plotly.react('pert-umap-hl',[
      {type:'scattergl',mode:'markers',name:'Other',x:gb.map(p=>p.x),y:gb.map(p=>p.y),
       marker:{color:'#b8c3cf',size:3,opacity:0.35},hoverinfo:'skip'},
      {type:'scattergl',mode:'markers',name:pert,x:ga.map(p=>p.x),y:ga.map(p=>p.y),
       text:ga.map(p=>p.pert),hovertemplate:'%{text}<extra></extra>',
       marker:{color:'#b04a5a',size:6,opacity:0.9}}],
      mkL({title:{text:'UMAP \u2014 '+pert+' highlighted',font:{size:26}},
           xaxis:{title:{text:'UMAP 1',font:{size:22}},zeroline:false,showgrid:false,linecolor:'#4c5f70',ticks:'outside',tickfont:{size:18}},
           yaxis:{title:{text:'UMAP 2',font:{size:22}},zeroline:false,showgrid:false,linecolor:'#4c5f70',ticks:'outside',tickfont:{size:18}},
           height:500,autosize:true,
           showlegend:true,legend:{font:{size:14},bgcolor:'rgba(255,255,255,0.85)',bordercolor:'#888',borderwidth:1}}),PC);
  }
  // top DEG bar
  reRenderDegBar();
  // pathway enrichment
  const allEnrich=(D.enrichment[pert]||[]).filter(r=>(r.n_overlap||0)>2).slice();
  const enrichSig=allEnrich.filter(r=>r.significant);
  const enrichTop=allEnrich.filter(r=>!r.significant).sort((a,b)=>(a.fdr??1)-(b.fdr??1));
  const enrich=[...enrichSig,...enrichTop].slice(0,50);
  const enrichEl=document.getElementById('pert-enrichment-bar');
  // Always purge any existing Plotly chart before toggling content,
  // otherwise Plotly.react fails after innerHTML was set to a message.
  try{Plotly.purge('pert-enrichment-bar');}catch(e){}
  if(enrich.length){
    const eht=Math.max(300,enrich.length*26+100);
    enrichEl.style.height=eht+'px';
    enrichEl.innerHTML=''; // clear any previous message
    const sigBadge=enrichSig.length?' <span style="color:#169b8f;font-size:13px;">('+enrichSig.length+' significant)</span>':'';
    Plotly.newPlot('pert-enrichment-bar',[{type:'bar',orientation:'h',width:0.7,
      x:enrich.map(r=>r.fdr>0?-Math.log10(r.fdr):0),
      y:enrich.map(r=>r.pathway),
      text:enrich.map(r=>'n='+r.n_overlap+'/'+r.n_pathway+', OR='+(r.odds_ratio!=null?Number(r.odds_ratio).toFixed(2):'')),
      hovertemplate:'<b>%{y}</b><br>-log\u2081\u2080(FDR): %{x:.2f}<br>%{text}<extra></extra>',
      marker:{color:enrich.map(r=>r.significant?'#169b8f':'#c08a2b')},
      customdata:enrich.map(r=>r.query_type||'sig_deg')}],
      mkL({title:{text:'Pathway Enrichment \u2014 '+pert,font:{size:16}},
           xaxis:{title:{text:'-log\u2081\u2080(FDR)',font:{size:14}}},
           yaxis:{automargin:true,tickfont:{size:12}},
           shapes:[{type:'line',
             x0:-Math.log10(0.05),x1:-Math.log10(0.05),
             y0:-0.5,y1:enrich.length-0.5,
             line:{color:'#b04a5a',dash:'dash',width:1.5}}],
           margin:{t:50,r:20,b:60,l:300}}),PC);
  } else {
    enrichEl.style.height='80px';
    enrichEl.textContent='No pathway enrichment terms with more than 2 overlapping genes for '+pert+'. When fewer than 5 significant DEGs are available, the pipeline falls back to top-200 genes by adjusted p-value.';
    enrichEl.style.padding='16px';
    enrichEl.style.color='var(--muted)';
  }
  renderGsea(pert);
  // DEG table
  _degData=rows;_degPage=0;_degSortCol=1;_degSortAsc=false;_degFilter='';
  const filt=document.getElementById('deg-filter'); if(filt)filt.value='';
  renderDegTable();
}

function renderGsea(pert){
  const el=document.getElementById('pert-gsea-bar');
  if(!el)return;
  const wf=document.getElementById('pert-gsea-waterfall');
  const rows=(D.gsea&&D.gsea[pert])?D.gsea[pert].slice(0,20):[];
  try{Plotly.purge('pert-gsea-bar');}catch(e){}
  if(wf){try{Plotly.purge('pert-gsea-waterfall');}catch(e){} wf.style.display='none';}
  if(!rows.length){
    el.style.height='80px';
    el.style.padding='16px';
    el.style.color='var(--muted)';
    el.textContent='No preranked GSEA terms available for '+pert+'.';
    return;
  }
  el.innerHTML='';
  const ht=Math.max(360,rows.length*24+110);
  el.style.height=ht+'px';
  Plotly.react('pert-gsea-bar',[{type:'bar',orientation:'h',width:0.7,
    x:rows.map(r=>r.es),y:rows.map(r=>r.term),
    text:rows.map(r=>'n='+r.n_genes+'; leading edge: '+(r.leading_edge||'').split('|').slice(0,5).join(', ')),
    customdata:rows.map((r,i)=>i),
    hovertemplate:'<b>%{y}</b><br>ES: %{x:.3f}<br>%{text}<extra></extra>',
    marker:{color:rows.map(r=>r.es>=0?'#B2182B':'#2166AC')}}],
    mkL({title:{text:'Preranked GSEA \u2014 '+pert,font:{size:16}},
         xaxis:{title:{text:'Enrichment score (ES)',font:{size:14}},zeroline:true,zerolinecolor:'#555'},
         yaxis:{automargin:true,tickfont:{size:12}},
         margin:{t:50,r:30,b:60,l:260}}),PC).then(gd=>{
    gd.on('plotly_click',ev=>{
      const i=ev.points&&ev.points[0]?ev.points[0].customdata:0;
      renderGseaWaterfall(rows[i]);
    });
    renderGseaWaterfall(rows[0]);
  });
}

function renderGseaWaterfall(row){
  const el=document.getElementById('pert-gsea-waterfall');
  if(!el||!row||!row.curve||!row.curve.length)return;
  el.style.display='';
  const hitX=(row.hit_idx||[]), hitY=hitX.map(()=>0);
  Plotly.react('pert-gsea-waterfall',[
    {type:'scatter',mode:'lines',name:'Running ES',x:row.curve_x,y:row.curve,
     line:{color:row.es>=0?'#B2182B':'#2166AC',width:2.5},
     hovertemplate:'Rank: %{x}<br>Running ES: %{y:.3f}<extra></extra>'},
    {type:'scatter',mode:'markers',name:'Gene-set hits',x:hitX,y:hitY,
     marker:{symbol:'line-ns-open',size:12,color:'#222',line:{width:1.5}},
     hovertemplate:'Hit rank: %{x}<extra></extra>'}
  ],mkL({title:{text:'GSEA running enrichment plot \u2014 '+row.term,font:{size:15}},
         xaxis:{title:{text:'Rank in DEG list',font:{size:13}}},
         yaxis:{title:{text:'Running enrichment score',font:{size:13}},zeroline:true,zerolinecolor:'#777'},
         margin:{t:55,r:25,b:55,l:75}}),PC);
}

function reRenderDegBar(){
  const rows=D.deg[_curPert]||[];
  const barFilt=document.getElementById('deg-bar-filter')?.value||'sig';
  let shown=[];
  if(barFilt==='sig') shown=rows.filter(r=>r.significant).slice(0,20);
  else if(barFilt==='top20lfc') shown=[...rows].sort((a,b)=>Math.abs(b.log2fc??0)-Math.abs(a.log2fc??0)).slice(0,20);
  else shown=rows.slice(0,40);
  const degBarHt=Math.max(340,shown.length*26+80);
  document.getElementById('pert-top-deg-bar').style.height=degBarHt+'px';
  Plotly.react('pert-top-deg-bar',[{type:'bar',orientation:'h',width:0.7,
    x:shown.map(r=>r.log2fc),y:shown.map(r=>r.gene),
    marker:{color:shown.map(r=>(r.log2fc??0)>0?'#b04a5a':'#1f4e79')},
    hovertemplate:'%{y}: %{x:.3f}<extra></extra>'}],
    mkL({title:{text:'Top Significant DEGs \u2014 '+_curPert,font:{size:17}},
         xaxis:{title:{text:'log\u2082FC',font:{size:14}}},
         yaxis:{automargin:true,tickfont:{size:14}},
         margin:{t:50,r:20,b:60,l:220}}),PC);
}

function renderDegTable(){
  const cols=['gene','log2fc','padj','significant','pathway_membership'];
  const q=(_degFilter||'').toLowerCase();
  const filtered=q?_degData.filter(r=>String(r.gene||'').toLowerCase().includes(q)||
    String(r.pathway_membership||r.top_pathway||'').toLowerCase().includes(q)):_degData;
  const sorted=[...filtered].sort((a,b)=>{
    const c=cols[_degSortCol];let av=a[c],bv=b[c];
    if(c==='log2fc'){av=Math.abs(av??0);bv=Math.abs(bv??0);}
    if(av==null)av=_degSortAsc?Infinity:-Infinity;
    if(bv==null)bv=_degSortAsc?Infinity:-Infinity;
    if(typeof av==='string')av=av.toLowerCase();
    if(typeof bv==='string')bv=bv.toLowerCase();
    if(av<bv)return _degSortAsc?-1:1;if(av>bv)return _degSortAsc?1:-1;return 0;
  });
  const total=sorted.length,start=_degPage*_degPageSize,page=sorted.slice(start,start+_degPageSize);
  document.getElementById('deg-table-info').textContent=total+' genes';
  document.getElementById('deg-page-info').textContent='Page '+(_degPage+1)+' / '+(Math.ceil(total/_degPageSize)||1);
  document.getElementById('deg-tbody').innerHTML=page.map(r=>{
    const fc=r.log2fc!=null?r.log2fc.toFixed(3):'\u2014';
    const pj=r.padj!=null?r.padj.toExponential(2):'\u2014';
    const fcCls=(r.log2fc??0)>0?'up':'dn';
    const sg=r.significant?'<span style="color:var(--accent2)">\u25cf</span>':'';
    const pw=(r.pathway_membership||r.top_pathway||'').split('|').slice(0,2).join(', ');
    return '<tr><td><b>'+(r.gene||'')+'</b></td><td class="'+fcCls+'">'+fc+'</td><td>'+pj+'</td><td>'+sg+'</td>'
      +'<td title="'+(r.pathway_membership||'')+'" style="font-size:12px;color:var(--muted);">'+pw+'</td></tr>';
  }).join('');
}

function sortDegTable(col){
  if(_degSortCol===col)_degSortAsc=!_degSortAsc;else{_degSortCol=col;_degSortAsc=col===0;}
  _degPage=0;renderDegTable();
}
function degPage(d){
  const q=(_degFilter||'').toLowerCase();
  const n=q?_degData.filter(r=>String(r.gene||'').toLowerCase().includes(q)||
    String(r.pathway_membership||r.top_pathway||'').toLowerCase().includes(q)).length:_degData.length;
  _degPage=Math.max(0,Math.min(Math.ceil(n/_degPageSize)-1,_degPage+d));
  renderDegTable();
}
function filterDegTable(q){_degFilter=q||'';_degPage=0;renderDegTable();}
function setDegPageSize(v){_degPageSize=parseInt(v,10)||10;_degPage=0;renderDegTable();}
function downloadDegTable(){
  if(!_degData.length)return;
  const cols=['gene','log2fc','padj','significant','pathway_membership'];
  const header=cols.join(',');
  const rows=_degData.map(r=>cols.map(c=>{
    const v=r[c]??'';return (typeof v==='string'&&v.includes(','))?'"'+v+'"':v;
  }).join(','));
  const csv=[header,...rows].join('\n');
  const a=document.createElement('a');a.href='data:text/csv;charset=utf-8,'+encodeURIComponent(csv);
  a.download=(_curPert||'degs')+'_DEG_table.csv';document.body.appendChild(a);a.click();document.body.removeChild(a);
}

/* GENE NETWORK */
function initNetworkTab(){
  const sel=document.getElementById('network-selector');
  const b64=D.genenet_b64||{};
  // Filter out spurious diff keys
  const pertKeys=Object.keys(b64.perturbations||{}).filter(k=>!k.endsWith('_diff'));
  sel.innerHTML='';
  pertKeys.sort().forEach(p=>{const o=document.createElement('option');o.value=p;o.textContent=p;sel.appendChild(o);});
  if(sel.options.length) renderNetwork(sel.value);
  else document.getElementById('tab-network').querySelector('p').textContent='No gene-network plots found. Run the genenet step.';
}
function renderNetwork(pert){
  const b64=D.genenet_b64||{};
  const entry=(b64.perturbations||{})[pert]||{};
  // Row 1: three square networks (control, perturbation, differential).
  const ctrlNetImg=document.getElementById('network-ctrl-net-img');
  // Prefer the per-pert control network (same gene set as the pert plot);
  // fall back to the global control network if not available.
  const ctrlNetSrc=entry.ctrl_network||b64.control_network;
  if(ctrlNetImg)ctrlNetImg.src=ctrlNetSrc||'';
  const cnLbl=document.getElementById('network-ctrl-label');
  if(cnLbl)cnLbl.textContent='Control Network';
  const pertImg=document.getElementById('network-pert-img');
  if(pertImg)pertImg.src=entry.pert_network||entry.network||'';
  const lbl=document.getElementById('network-pert-label');
  if(lbl)lbl.textContent=pert+' Network';
  const diffImg=document.getElementById('network-diff-img');
  if(diffImg)diffImg.src=entry.diff||'';
  const dlbl=document.getElementById('network-diff-label');
  if(dlbl)dlbl.textContent=pert+' Differential Network';
  // Row 2: bottom heatmaps for the genes shown in the networks.
  const hmImg=document.getElementById('network-heatmap-img');
  if(hmImg)hmImg.src=entry.heatmap||'';
  const hmLbl=document.getElementById('network-heatmap-label');
  if(hmLbl)hmLbl.textContent=pert+' Gene Expression Heatmap';
  renderNetworkCtrlHeatmap(pert);
  const dlBtn=document.getElementById('network-dl-btn');
  if(dlBtn)dlBtn.onclick=function(){
    if(!entry.heatmap)return;
    const a=document.createElement('a');a.href=entry.heatmap;
    a.download=pert+'_heatmap_comparison.png';document.body.appendChild(a);a.click();document.body.removeChild(a);
  };
}

function renderNetworkCtrlHeatmap(pert){
  const el=document.getElementById('network-ctrl-heatmap-chart');
  if(!el)return;
  const net=(D.genenet_json||{})[pert];
  if(!net||!net.genes||!net.genes.length){
    el.innerHTML='<p style="padding:16px;color:var(--muted);">No interactive control correlation matrix for '+pert+'.</p>';
    return;
  }
  const genes=net.genes;
  const idx={}; genes.forEach((g,i)=>{idx[g]=i;});
  const mat=genes.map((_,i)=>genes.map((__,j)=>i===j?1:0));
  (net.ctrl_edges||[]).forEach(e=>{
    const i=idx[e.s],j=idx[e.t];
    if(i==null||j==null)return;
    mat[i][j]=e.r;mat[j][i]=e.r;
  });
  const side=Math.max(500,Math.min(720,genes.length*12+160));
  el.style.height=side+'px';
  Plotly.react('network-ctrl-heatmap-chart',[{type:'heatmap',z:mat,x:genes,y:genes,
    zmin:-1,zmax:1,zmid:0,colorscale:'RdBu',reversescale:false,
    colorbar:{title:{text:'Pearson r'}},
    hovertemplate:'%{y} vs %{x}<br>r: %{z:.3f}<extra></extra>'}],
    mkL({title:{text:'Control co-expression \u2014 '+pert,font:{size:16}},
         xaxis:{tickangle:-60,automargin:true,scaleanchor:'y',scaleratio:1},
         yaxis:{automargin:true,autorange:'reversed'},
         width:side,height:side,autosize:false,margin:{t:55,r:90,b:150,l:150}}),PC);
}

/* GENE EXPRESSION */
function initGeneTab(){
  const sel=document.getElementById('gene-selector');
  sel.innerHTML='';
  D.genes.forEach(g=>{const o=document.createElement('option');o.value=g;o.textContent=g;sel.appendChild(o);});
  if(D.genes.length)renderGeneExpr(D.genes[0]);
}

function renderGeneExpr(gene){
  if(!gene)return;
  const expr=D.gene_expr[gene];
  if(!expr){document.getElementById('gene-box').innerHTML='<p style="padding:16px;color:var(--muted);">No data for: '+gene+'</p>';return;}
  const sel=document.getElementById('gene-selector');
  if(D.genes.includes(gene)&&sel.value!==gene)sel.value=gene;
  // Put control first
  let perts=Object.keys(expr);
  if(perts.includes('control'))perts=['control',...perts.filter(p=>p!=='control')];
  // theme_classic: white background, no grid, axis lines only
  const classicLayout={
    title:{text:'Expression \u2014 '+gene,font:{size:17,color:'#203040',family:'Arial,Helvetica,sans-serif'}},
    paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',
    yaxis:{title:{text:'log-norm expression',font:{size:14}},zeroline:false,showgrid:false,
           linecolor:'#203040',linewidth:1.5,ticks:'outside',mirror:false},
    xaxis:{tickangle:-35,tickfont:{size:12},zeroline:false,showgrid:false,
           linecolor:'#203040',linewidth:1.5,ticks:'outside',mirror:false},
    showlegend:false,margin:{t:55,r:20,b:100,l:75},
    hoverlabel:{font:{size:13}},
  };
  Plotly.react('gene-box',perts.map(p=>({type:'box',name:p,y:expr[p],
    boxmean:false,boxpoints:'outliers',whiskerwidth:0.5,
    marker:{color:p==='control'?'#169b8f':'#4c78a8',size:4,opacity:0.7,
            outliercolor:p==='control'?'#169b8f':'#4c78a8'},
    line:{color:p==='control'?'#169b8f':'#4c78a8',width:1.5},
    fillcolor:p==='control'?'rgba(22,155,143,.05)':'rgba(76,120,168,.05)',
    hovertemplate:p+'<br>%{y:.3f}<extra></extra>'})),
    classicLayout,PC);
  const ctrl=expr['control']||[];
  const cmed=ctrl.length?[...ctrl].sort((a,b)=>a-b)[Math.floor(ctrl.length/2)]:0;
  const fcData=perts.filter(p=>p!=='control').map(p=>{
    const vals=expr[p];const med=vals.length?[...vals].sort((a,b)=>a-b)[Math.floor(vals.length/2)]:0;
    return {pert:p,fc:(med-cmed)/Math.LN2*Math.LOG2E};
  }).sort((a,b)=>b.fc-a.fc);
  Plotly.react('gene-fc-bar',[{type:'bar',orientation:'h',width:0.7,
    x:fcData.map(r=>r.fc),y:fcData.map(r=>r.pert),
    marker:{color:fcData.map(r=>r.fc>0?'#b04a5a':'#1f4e79')},
    hovertemplate:'%{y}: %{x:.3f}<extra></extra>'}],
    mkL({title:{text:'log\u2082FC vs Control \u2014 '+gene,font:{size:17}},
         xaxis:{title:{text:'Median log\u2082FC',font:{size:14}}},
         yaxis:{automargin:true,tickfont:{size:12}}}),PC);
}

/* PARAMS */
function renderParams(){
  const cfg=D.config||{};
  const done = new Set(cfg._completed_steps||[]);
  const stepBody=document.getElementById('step-param-tbody');
  stepBody.innerHTML=STEP_PARAMS.map(s=>{
    const pv=s.params.length?s.params.map(k=>'<span class="pn">'+k+'</span>: <span class="pv">'+(cfg[k]!==undefined?JSON.stringify(cfg[k]):'default')+'</span>').join('<br>'):'No tunable parameters';
    return '<tr><td class="pn">'+s.step+'</td><td>'+(done.has(s.step)?'Completed':'Not recorded')+'</td><td style="color:var(--muted);font-size:13px;">'+s.purpose+'</td><td>'+pv+'</td></tr>';
  }).join('');
  const keys=Object.keys(cfg).filter(k=>!k.startsWith('_'));
  const tbody=document.getElementById('param-tbody');
  if(!keys.length){tbody.innerHTML='<tr><td colspan="3" style="color:var(--muted);padding:14px;">No config data available.</td></tr>';return;}
  tbody.innerHTML=keys.map(k=>'<tr><td class="pn">'+k+'</td><td class="pv">'+JSON.stringify(cfg[k])+'</td><td style="color:var(--muted);font-size:13px;">'+(PARAM_DESC[k]||'')+'</td></tr>').join('');
}

/* Save the live Plotly UMAP as a high-resolution PNG. Filename keys off the
   "Color by" select so e.g. score / gene / state runs land on disk under
   distinguishable names. */
function downloadUmapPng() {
  const el = document.getElementById('umap-chart');
  if (!el) return;
  const colorBy = (document.getElementById('umap-color-by')||{}).value || 'umap';
  const gene = (document.getElementById('umap-gene-input')||{}).value || '';
  const stem = colorBy === 'gene' && gene ? 'umap_gene_' + gene : 'umap_' + colorBy;
  Plotly.downloadImage(el, {
    format: 'png', filename: stem,
    height: 1600, width: 1920, scale: 4,
  });
}

/* Save the home-tab "Effect Decomposition" Plotly scatter as a PNG.
   Compact square geometry (600×600 layout) but rendered at scale 6 →
   3600×3600 px output, so the file is small in figure-units yet sharp. */
function downloadEffectPng() {
  const el = document.getElementById('home-effect-scatter');
  if (!el) return;
  Plotly.downloadImage(el, {
    format: 'png', filename: 'effect_decomposition',
    height: 600, width: 600, scale: 6,
  });
}

/* CELL STATES tab */
function renderStateEnrich() {
  const se = D.state_enrich;
  const el = document.getElementById('states-heatmap-chart');
  if (!se || !se.perturbations || !se.perturbations.length) {
    el.innerHTML = '<p style="padding:20px;color:var(--muted);">No cell state enrichment data. Run the state_enrich step.</p>';
    return;
  }
  const isFallback = se._is_fallback;
  const minQ = parseFloat(document.getElementById('states-min-q')?.value || 0);
  const idxKeep = se.perturbations.map((_,i) =>
    Math.max(...se.matrix[i].map(Math.abs)) >= minQ ? i : -1
  ).filter(i => i >= 0);
  const perts = idxKeep.map(i => se.perturbations[i]);
  const zRows = idxKeep.map(i => se.matrix[i]);
  const vmax = Math.min(Math.max(...zRows.flat().map(Math.abs), 1), 6);
  const ht = Math.max(380, perts.length * 18 + 80);
  el.style.height = ht + 'px';
  Plotly.react('states-heatmap-chart', [{
    type: 'heatmap',
    x: se.cell_states, y: perts, z: zRows,
    colorscale: DIVG, reversescale: false,
    zmid: 0, zmin: -vmax, zmax: vmax,
    colorbar: {
      title: {text: isFallback ? 'log\u2082 ratio vs ctrl' : 'signed \u2212log\u2081\u2080(q)', font: {size: 13}},
      tickfont: {size: 11}, thickness: 14,
    },
    hovertemplate: 'Perturbation: <b>%{y}</b><br>Cell state: %{x}<br>' + (isFallback ? 'log\u2082 ratio' : 'signed \u2212log\u2081\u2080(q)') + ': %{z:.2f}<extra></extra>',
  }], mkL({
    title: {text: isFallback ? 'Cluster proportion (log\u2082FC vs control)' : 'Cell state enrichment or depletion', font: {size: 17}},
    xaxis: {title: {text: 'Cell state', font: {size: 14}}, tickfont: {size: 12}},
    yaxis: {
      title: {text: 'Perturbed genes', font: {size: 14}},
      tickfont: {size: 11}, automargin: true,
    },
    margin: {t: 60, r: 140, b: 60, l: 160},
  }), PC);
}

/* REGULATORY tab */
function initRegulatoryTab() {
  renderPertCorr(D.pert_corr_before, 'reg-corr-before-chart', 'Without modeling cell states');
  renderPertCorr(D.pert_corr_after,  'reg-corr-after-chart',  'After modeling cell states');
  renderTFHeatmap();
  renderTFNetwork();
}

function renderPertCorr(corr_data, elem_id, title) {
  const el = document.getElementById(elem_id);
  if (!corr_data || !corr_data.perts || corr_data.perts.length < 2) {
    el.innerHTML = '<p style="padding:16px;color:var(--muted);">No correlation data available.</p>';
    return;
  }
  const n = corr_data.perts.length;
  const ht = Math.max(360, n * 22 + 100);
  el.style.height = ht + 'px';
  Plotly.react(elem_id, [{
    type: 'heatmap',
    x: corr_data.perts, y: corr_data.perts, z: corr_data.matrix,
    colorscale: DIVG, reversescale: false,
    zmid: 0, zmin: -1, zmax: 1,
    colorbar: {title: {text: 'Pearson r', font: {size: 15}}, tickfont: {size: 13}, thickness: 14},
    hovertemplate: '%{y} \u00d7 %{x}: <b>%{z:.3f}</b><extra></extra>',
  }], mkL({
    title: {text: title, font: {size: 22}},
    xaxis: {tickfont: {size: 16}, tickangle: -45},
    yaxis: {
      title: {text: 'Perturbation', font: {size: 17}},
      tickfont: {size: 16}, automargin: true,
    },
    margin: {t: 55, r: 30, b: 130, l: 160},
  }), PC);
}

function renderTFHeatmap() {
  const el = document.getElementById('reg-tf-heatmap-chart');
  const d = D.tf_regulatory;
  if (!d || !d.guides || d.guides.length < 2) {
    el.innerHTML = '<p style="padding:16px;color:var(--muted);">No TF regulatory data. The perturbed gene names must overlap with gene names in the dataset for this analysis.</p>';
    return;
  }
  const vmax = Math.min(Math.max(...d.matrix.flat().map(Math.abs), 1), 8);
  const ht = Math.max(620, d.guides.length * 34 + 180);
  el.style.height = ht + 'px';
  el.style.maxWidth = '760px';
  el.style.margin = '0 auto';
  Plotly.react('reg-tf-heatmap-chart', [{
    type: 'heatmap',
    x: d.genes, y: d.guides, z: d.matrix,
    colorscale: DIVG, reversescale: false,
    zmid: 0, zmin: -vmax, zmax: vmax,
    colorbar: {
      title: {text: 'signed \u2212log\u2081\u2080(q-value)', font: {size: 16}},
      tickfont: {size: 13}, thickness: 14,
      tickvals: [-vmax, -vmax/2, 0, vmax/2, vmax],
      ticktext: ['Activating<br>(KO\u2192\u2193)', '', '0', '', 'Inhibiting<br>(KO\u2192\u2191)'],
    },
    hovertemplate: 'Guide: <b>%{y}</b><br>Gene: <b>%{x}</b><br>signed \u2212log\u2081\u2080(q): %{z:.2f}<extra></extra>',
  }], mkL({
    title: {text: 'TF regulatory relationships', font: {size: 22}},
    xaxis: {title: {text: 'Genes', font: {size: 18}}, tickfont: {size: 16}, tickangle: -60},
    yaxis: {title: {text: 'Guides', font: {size: 18}}, tickfont: {size: 16}, automargin: true},
    margin: {t: 70, r: 130, b: 180, l: 150},
  }), PC);
}

function renderTFNetwork() {
  const el = document.getElementById('reg-tfnet-chart');
  const net = D.tf_network;
  if (!net || !net.nodes || !net.nodes.length) {
    el.innerHTML = '<p style="padding:16px;color:var(--muted);">No TF network data. The perturbed gene names must overlap with gene names in the dataset.</p>';
    return;
  }
  const lfc_thresh = parseFloat(document.getElementById('tfnet-lfc-thresh')?.value || 0.3);
  const pos = net.positions || {};
  const traces = [];

  const _addEdgeTrace = (edges, color, name) => {
    const ex = [], ey = [];
    const anns = [];
    for (const e of edges) {
      if (e.weight < lfc_thresh) continue;
      const s = pos[e.source], t = pos[e.target];
      if (!s || !t) continue;
      ex.push(s[0], t[0], null);
      ey.push(s[1], t[1], null);
      anns.push({
        ax: s[0], ay: s[1], axref: 'x', ayref: 'y',
        x: t[0], y: t[1], xref: 'x', yref: 'y',
        showarrow: true, arrowhead: 2, arrowsize: 1.2,
        arrowwidth: Math.max(1, e.weight * 2),
        arrowcolor: color, opacity: 0.7,
      });
    }
    if (ex.length) {
      traces.push({
        type: 'scatter', mode: 'lines', name,
        x: ex, y: ey, line: {color, width: 1.2}, opacity: 0.6, hoverinfo: 'skip',
      });
    }
    return anns;
  };

  const anns_act = _addEdgeTrace(net.edges_activating || [], '#2166AC', 'Activating');
  const anns_inh = _addEdgeTrace(net.edges_inhibiting || [], '#B2182B', 'Inhibiting');

  // Degree for node sizing
  const deg = {};
  for (const e of [...(net.edges_activating||[]), ...(net.edges_inhibiting||[])]) {
    deg[e.source] = (deg[e.source]||0) + 1;
    deg[e.target] = (deg[e.target]||0) + 1;
  }
  const nx_arr = net.nodes.map(n => pos[n]?.[0] ?? 0);
  const ny_arr = net.nodes.map(n => pos[n]?.[1] ?? 0);
  traces.push({
    type: 'scatter', mode: 'markers+text', name: 'TF',
    x: nx_arr, y: ny_arr, text: net.nodes,
    textposition: 'top center', textfont: {size: 16, color: '#203040'},
    hovertemplate: '<b>%{text}</b><br>Degree: %{marker.size}<extra></extra>',
    marker: {
      size: net.nodes.map(n => Math.max(10, (deg[n]||0) * 4 + 10)),
      color: '#2a9d8f', line: {color: '#fff', width: 1.5}, opacity: 0.9,
    },
  });

  Plotly.react('reg-tfnet-chart', traces, mkL({
    title: {text: 'TF Regulatory Network', font: {size: 22}},
    xaxis: {
      title: {text: 'Layout X', font: {size: 16}},
      visible: true, zeroline: true, zerolinecolor: '#b8c2cc',
      showgrid: true, gridcolor: '#e7edf2', tickfont: {size: 13},
    },
    yaxis: {
      title: {text: 'Layout Y', font: {size: 16}},
      visible: true, zeroline: true, zerolinecolor: '#b8c2cc',
      showgrid: true, gridcolor: '#e7edf2', tickfont: {size: 13},
      scaleanchor: 'x', scaleratio: 1,
    },
    showlegend: true,
    legend: {font: {size: 17}},
    annotations: [...anns_act, ...anns_inh],
    margin: {t: 70, r: 35, b: 70, l: 75},
  }), PC);
}

/* C-SCORE tab */
function initCscoreTab() {
  const cs = D.cscore || {};
  const summary = cs.summary || [];
  const b64 = cs.b64 || {};
  const noData = !summary.length;

  // Summary cards: max C_total, mean C_gain, mean C_loss
  const cards = document.getElementById('cscore-cards');
  if (cards) {
    if (noData) {
      cards.innerHTML = '<p style="color:var(--muted);padding:12px;">No C-score data found. Run the cscore step.</p>';
    } else {
      const maxT = Math.max(...summary.map(r => r.c_total || 0)).toFixed(3);
      const meanG = (summary.reduce((a, r) => a + (r.c_gain || 0), 0) / summary.length).toFixed(3);
      const meanL = (summary.reduce((a, r) => a + (r.c_loss || 0), 0) / summary.length).toFixed(3);
      const meanS = (summary.reduce((a, r) => a + (r.c_shift || 0), 0) / summary.length).toFixed(3);
      cards.innerHTML = [
        {l: 'Perturbations scored', v: summary.length, c: ''},
        {l: 'Max C_total', v: maxT, c: ''},
        {l: 'Mean C_gain', v: meanG, c: 'green'},
        {l: 'Mean C_loss', v: meanL, c: ''},
        {l: 'Mean C_shift', v: meanS, c: ''},
      ].map(c => '<div class="card"><div class="label">' + c.l + '</div><div class="value ' + c.c + '">' + c.v + '</div></div>').join('');
    }
  }

  // Static PNG images
  const imgs = {
    'cscore-ranked-bar-img': b64.ranked_bar,
    'cscore-decomp-img': b64.decomposition,
    'cscore-vs-deg-img': b64.vs_deg,
    'cscore-hub-img': b64.hub_heatmap,
    'cscore-module-img': b64.module_rewiring,
  };
  Object.entries(imgs).forEach(([id, src]) => {
    const el = document.getElementById(id);
    if (el) el.src = src || '';
  });

  // Edge explorer selector
  const sel = document.getElementById('cscore-pert-selector');
  if (sel) {
    sel.innerHTML = '';
    const perts = Object.keys(cs.edges || {}).sort();
    perts.forEach(p => { const o = document.createElement('option'); o.value = p; o.textContent = p; sel.appendChild(o); });
    if (perts.length) renderCscoreEdges(perts[0]);
  }
}

function renderCscoreEdges(pert, geneFilter) {
  const cs = D.cscore || {};
  const edges = (cs.edges || {})[pert] || [];
  const summary = (cs.summary || []).find(r => r.perturbation === pert) || {};
  const geneSel = document.getElementById('cscore-gene-selector');
  const genes = Array.from(new Set(edges.flatMap(e => [e.gene_A, e.gene_B]).filter(Boolean))).sort();
  if (geneSel) {
    const prior = geneFilter || geneSel.value || '__all__';
    geneSel.innerHTML = '<option value="__all__">All genes</option>' +
      genes.map(g => '<option value="' + g.replace(/"/g, '&quot;') + '">' + g + '</option>').join('');
    geneSel.value = genes.includes(prior) ? prior : '__all__';
    geneFilter = geneSel.value;
  }
  const shownEdges = geneFilter && geneFilter !== '__all__'
    ? edges.filter(e => e.gene_A === geneFilter || e.gene_B === geneFilter)
    : edges;

  // Scatter: ctrl_r vs pert_r, coloured by status
  const colMap = {gained: '#2dc653', lost: '#e63946', shared: '#9aa7b4'};
  const groups = {};
  shownEdges.forEach(e => {
    const g = e.status || 'shared';
    if (!groups[g]) groups[g] = {x: [], y: [], t: []};
    groups[g].x.push(e.ctrl_r);
    groups[g].y.push(e.pert_r);
    groups[g].t.push(e.gene_A + ' \u2013 ' + e.gene_B);
  });
  const traces = Object.entries(groups).map(([status, g]) => ({
    type: 'scatter', mode: 'markers', name: status,
    x: g.x, y: g.y, text: g.t,
    hovertemplate: '<b>%{text}</b><br>ctrl r: %{x:.3f}<br>pert r: %{y:.3f}<extra></extra>',
    marker: {color: colMap[status] || '#888', size: 7, opacity: 0.75},
  }));
  // Diagonal guide
  const allR = shownEdges.flatMap(e => [e.ctrl_r, e.pert_r]).filter(v => v != null);
  const rmin = Math.min(...allR, -0.2), rmax = Math.max(...allR, 0.2);
  traces.push({
    type: 'scatter', mode: 'lines', name: 'No change',
    x: [rmin, rmax], y: [rmin, rmax],
    line: {color: '#aaa', dash: 'dash', width: 1.5}, hoverinfo: 'skip',
  });
  const cStr = summary.c_total != null ? (' \u2014 C_total\u202f=\u202f' + (summary.c_total).toFixed(3)) : '';
  const gStr = geneFilter && geneFilter !== '__all__' ? (' \u2014 ' + geneFilter) : '';
  Plotly.react('cscore-edge-scatter', traces, mkL({
    title: {text: 'Edge weights: ctrl vs pert \u2014 ' + pert + gStr + cStr, font: {size: 17}},
    xaxis: {title: {text: 'Control Pearson r', font: {size: 14}}, zeroline: true, zerolinecolor: '#ccc'},
    yaxis: {title: {text: 'Perturbation Pearson r', font: {size: 14}}, zeroline: true, zerolinecolor: '#ccc'},
    showlegend: true, legend: {font: {size: 13}},
    margin: {t: 55, r: 20, b: 65, l: 80},
  }), PC);

  // Edge table (top 200 by |delta_r|)
  const sorted = [...shownEdges].sort((a, b) => Math.abs(b.delta_r || 0) - Math.abs(a.delta_r || 0)).slice(0, 200);
  const colClass = {gained: 'sig', lost: 'up', shared: ''};
  document.getElementById('cscore-edge-tbody').innerHTML = sorted.map(e => {
    const dr = e.delta_r != null ? e.delta_r.toFixed(3) : '\u2014';
    const cr = e.ctrl_r != null ? e.ctrl_r.toFixed(3) : '\u2014';
    const pr = e.pert_r != null ? e.pert_r.toFixed(3) : '\u2014';
    const cls = colClass[e.status] || '';
    return '<tr><td>' + (e.gene_A || '') + '</td><td>' + (e.gene_B || '') + '</td>'
      + '<td class="' + cls + '">' + (e.status || '') + '</td>'
      + '<td>' + cr + '</td><td>' + pr + '</td><td class="' + cls + '">' + dr + '</td></tr>';
  }).join('');
}

/* SIDEBAR */
const SB_TOP=10;let _sbAll=false;
function buildSidebarPerts(filter){
  const all=D.perturbations||[];
  const filt=filter?all.filter(p=>p.toLowerCase().includes(filter.toLowerCase())):all;
  const show=(_sbAll||filter)?filt:filt.slice(0,SB_TOP);
  document.getElementById('sb-pert-list').innerHTML=show.map(p=>
    '<button style="margin-bottom:2px;font-size:12px;padding:4px 8px;" onclick="selectSidebarPert(\''+p.replace(/'/g,"\\'")+'\')">'
    +p+'</button>').join('');
  const btn=document.getElementById('sb-pert-show-all');
  if(!filter&&!_sbAll&&filt.length>SB_TOP){btn.style.display='';btn.textContent='Show all '+filt.length+' \u25bc';}
  else btn.style.display='none';
}
function filterSidebarPerts(q){buildSidebarPerts(q);}
function showAllPerts(){_sbAll=true;buildSidebarPerts('');}
function selectSidebarPert(pert){
  showTab('pert',null,'tab');
  const sel=document.getElementById('pert-selector');
  if(!Array.from(sel.options).some(o=>o.value===pert))expandPertDropdown();
  sel.value=pert;renderPert(pert);
}
function buildSidebarGenes(filter){
  const genes=D.genes||[];
  const filt=filter?genes.filter(g=>g.toLowerCase().includes(filter.toLowerCase())):genes;
  document.getElementById('sb-gene-list').innerHTML=filt.slice(0,30).map(g=>
    '<button style="margin-bottom:2px;font-size:12px;padding:4px 8px;" onclick="selectSidebarGene(\''+g.replace(/'/g,"\\'")+'\')">'
    +g+'</button>').join('');
}
function filterSidebarGenes(q){buildSidebarGenes(q);}
function selectSidebarGene(gene){showTab('gene',null,'tab');renderGeneExpr(gene);}

document.addEventListener('DOMContentLoaded',()=>{
  try{
    bindTabButtons();
    document.getElementById('home-cards').innerHTML='<div class="card"><div class="label">Status</div><div class="value">Loading</div></div>';
    document.getElementById('home-analysis-text').textContent='Loading interactive report data...';
    loadReportData().then(()=>{
      buildSidebarPerts('');buildSidebarGenes('');
      _tabInit['home']=true;
      requestAnimationFrame(()=>{try{renderHome();}catch(e){showStartupError(e);}});
    }).catch(()=>{});
  }catch(e){showStartupError(e);}
});
"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def build_interactive_report(
    adata,
    output_dir: str,
    n_top_genes: int = 50,
    max_cells_per_group: int = 30,
) -> Path:
    """Build a self-contained interactive HTML report (v2).

    Args:
        adata               -- AnnData after QC + normalisation.
        output_dir          -- Root output directory (same as pipeline output).
        n_top_genes         -- Top DEG genes to include in expression view.
        max_cells_per_group -- Cells sampled per perturbation for box plots.

    Returns:
        Path to the generated HTML file.
    """
    out = ensure_dir(output_dir)
    csv_dir = out / "csv"

    data = _extract_data(adata, csv_dir, n_top_genes, max_cells_per_group)
    template = _html_template()
    json_str = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    data_path = out / "interactive_data.json"
    data_path.write_text(json_str, encoding="utf-8")
    stub = {
        "summary": data.get("summary", {}),
        "config": data.get("config", {}),
    }
    stub_str = json.dumps(stub, ensure_ascii=False, separators=(",", ":"))
    html = template.replace("{DATA_STUB}", stub_str)

    out_path = out / "interactive_report.html"
    out_path.write_text(html, encoding="utf-8")
    _copy_plotly_bundle(out)
    return out_path
