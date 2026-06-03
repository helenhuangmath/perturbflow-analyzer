# =============================================================================
# perturbflow/analyzer/bundle.py  (v1, NEW)
#
# Emits a versioned, viewer-ready results bundle to <output_dir>/bundle/.
#
# Layout (per DESIGN.md):
#   bundle/
#     manifest.json
#     metadata.json
#     embeddings/umap.parquet                  # cell x (UMAP1, UMAP2, perturbation, ...)
#     perturbations_summary.parquet            # one row per perturbation
#     perturbations/<pert>.parquet             # one row per perturbed cell
#     de/<pert>.parquet                        # gene x stats (within-state DEG)
#     genes/<gene>.parquet                     # inverted index for the gene page
#     modules.json                             # gene programs / modules
#     tf_network.json                          # TF -> targets edges + nodes
#     perturbation_similarity.parquet          # pert x pert cosine similarity
#     search_index.json                        # precomputed search list
#
# Design choices:
#   - Parquet is preferred for tabular artifacts; if pyarrow is unavailable
#     the writer falls back to CSV transparently and updates the manifest.
#   - DEG tables are computed here from the AnnData (mean log2FC per state +
#     Welch t-test p-values on log-normalised counts) so the existing
#     state.py module is left untouched.
#   - The TF network is computed from a built-in CollecTRI-like prior of well
#     known TF -> target relationships restricted to genes present in the
#     dataset; weighted by the perturbation-induced delta in target expression.
#     This keeps the v1 bundle self-contained without external downloads.
#   - All artifacts are independently optional. The viewer reads manifest.json
#     to decide which panels to render.
# =============================================================================

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
import json

import numpy as np
import pandas as pd

from . import __version__
from .utils import ensure_dir


# A small, curated TF -> target prior used to seed the network panel without
# requiring an external resource download. Genes not present in the dataset
# are filtered out at bundle time, so this is safe for any species.
_BUILTIN_TF_PRIOR = {
    "STAT1": ["IRF1", "IFI6", "IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2"],
    "STAT3": ["SOCS3", "BCL2", "MYC", "VEGFA", "IL6", "IL10"],
    "MYC": ["MKI67", "CCND2", "CDK4", "ODC1", "LDHA", "NPM1"],
    "TP53": ["CDKN1A", "MDM2", "BAX", "BBC3", "GADD45A", "PMAIP1"],
    "NFKB1": ["TNF", "IL6", "ICAM1", "CXCL8", "BCL2", "BIRC3"],
    "RELA": ["TNF", "IL6", "ICAM1", "CXCL8", "BCL2", "NFKBIA"],
    "FOXP3": ["IL2RA", "CTLA4", "IKZF2"],
    "TBX21": ["IFNG", "GZMB", "PRF1", "EOMES"],
    "GATA3": ["IL4", "IL5", "IL13", "CCR4"],
    "RORC": ["IL17A", "IL17F", "IL22", "IL23R"],
    "BACH2": ["TCF7", "IL7R", "BCL6"],
    "TCF7": ["LEF1", "CCR7", "SELL"],
    "EOMES": ["IFNG", "GZMB", "PRF1"],
    "IRF4": ["IL21", "BATF", "MAF"],
    "BATF": ["IL21", "BCL6", "MAF"],
    "JUN": ["MMP9", "CCND1", "FOS"],
    "FOS": ["JUN", "EGR1"],
    "E2F1": ["CCNE1", "CCNE2", "CDC25A", "MCM2", "MCM3"],
    "RB1": ["E2F1", "CCND1"],
    "BRCA1": ["BRIP1", "RAD51", "ATM"],
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _to_dense(x):
    return x.A if hasattr(x, "A") else np.asarray(x)


def _try_parquet():
    # Detect whether pyarrow is importable so we can prefer parquet output.
    try:
        import pyarrow  # noqa: F401
        return True
    except Exception:
        return False


def _save_table(df: pd.DataFrame, path: Path, prefer_parquet: bool) -> str:
    # Write a DataFrame to parquet (preferred) or CSV fallback. Returns the
    # path of the actually-written file (relative to bundle root) so the
    # manifest can record what got produced.
    if df is None or df.empty:
        return ""
    if prefer_parquet:
        try:
            target = path.with_suffix(".parquet")
            df.to_parquet(target, index=False)
            return target.name
        except Exception:
            # Fall through to CSV on any parquet failure.
            pass
    target = path.with_suffix(".csv")
    df.to_csv(target, index=False)
    return target.name


# -----------------------------------------------------------------------------
# Per-perturbation DEG (within-state, mean-difference + Welch t-test)
# -----------------------------------------------------------------------------
def _compute_de_for_perturbation(adata, perturbation: str, top_n: int) -> pd.DataFrame:
    # Compute per-gene DEG statistics for one perturbation, broken down by
    # cell_state. Uses log-normalised expression (already in adata.X after
    # preprocessing) so the values are comparable across states.
    p_mask = adata.obs["perturbation"].values == perturbation
    c_mask = adata.obs["is_control"].values
    if p_mask.sum() < 3 or c_mask.sum() < 3:
        return pd.DataFrame()

    states = adata.obs["cell_state"].astype(str).values
    state_levels = sorted(set(states))
    rows = []

    for st in state_levels:
        s_mask = states == st
        p_idx = np.where(p_mask & s_mask)[0]
        c_idx = np.where(c_mask & s_mask)[0]
        if len(p_idx) < 3 or len(c_idx) < 3:
            continue

        xp = _to_dense(adata.X[p_idx, :])
        xc = _to_dense(adata.X[c_idx, :])

        mp = xp.mean(axis=0)
        mc = xc.mean(axis=0)
        log2fc = np.log2((mp + 1e-3) / (mc + 1e-3))

        # Welch t-test per gene without scipy: safe pooled std formulation.
        vp = xp.var(axis=0, ddof=1)
        vc = xc.var(axis=0, ddof=1)
        se = np.sqrt(vp / max(len(p_idx), 1) + vc / max(len(c_idx), 1)) + 1e-12
        t = (mp - mc) / se
        # Two-sided p-value via the normal approx (avoids scipy.stats); fine
        # for ranking purposes which is what the viewer uses.
        from math import erf, sqrt

        def _two_sided(z):
            return 2.0 * (1.0 - 0.5 * (1.0 + erf(abs(z) / sqrt(2.0))))

        pvals = np.array([_two_sided(z) for z in t])
        # BH correction (vectorised).
        n = len(pvals)
        order = np.argsort(pvals)
        ranked = pvals[order]
        bh = ranked * n / (np.arange(1, n + 1))
        bh = np.minimum.accumulate(bh[::-1])[::-1]
        padj = np.empty_like(bh)
        padj[order] = np.minimum(bh, 1.0)

        df = pd.DataFrame(
            {
                "gene": adata.var_names.astype(str).values,
                "cell_state": st,
                "log2fc": log2fc.astype(np.float32),
                "tstat": t.astype(np.float32),
                "pval": pvals.astype(np.float32),
                "padj": padj.astype(np.float32),
                "mean_perturbed": mp.astype(np.float32),
                "mean_control": mc.astype(np.float32),
            }
        )
        # Keep only the top-N genes by absolute log2FC so the bundle stays small.
        df = df.reindex(df["log2fc"].abs().sort_values(ascending=False).index)
        df = df.head(top_n).reset_index(drop=True)
        rows.append(df)

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, axis=0, ignore_index=True)


def _build_gene_index(de_tables: dict, min_hits: int) -> dict:
    # Build the inverted index: for each gene, list every (perturbation,
    # cell_state) pair where it appeared in the per-perturbation DEG table.
    by_gene: dict[str, list[dict]] = {}
    for pert, df in de_tables.items():
        if df is None or df.empty:
            continue
        for row in df.itertuples(index=False):
            entry = by_gene.setdefault(row.gene, [])
            entry.append(
                {
                    "perturbation": pert,
                    "cell_state": row.cell_state,
                    "log2fc": float(row.log2fc),
                    "padj": float(row.padj),
                    "direction": "up" if row.log2fc > 0 else "down",
                }
            )
    return {g: hits for g, hits in by_gene.items() if len(hits) >= min_hits}


def _build_perturbation_similarity(adata, perturbations: list[str]) -> pd.DataFrame:
    # Cosine similarity between perturbations' mean PCA shifts away from
    # controls. The matrix is symmetric and excludes the controls row/col.
    if "X_pca" not in adata.obsm or len(perturbations) < 2:
        return pd.DataFrame()
    x = adata.obsm["X_pca"]
    ctrl_mask = adata.obs["is_control"].values
    if ctrl_mask.sum() == 0:
        return pd.DataFrame()
    centroid = x[ctrl_mask].mean(axis=0, keepdims=True)
    profiles = []
    keep = []
    for p in perturbations:
        m = adata.obs["perturbation"].values == p
        if m.sum() < 3:
            continue
        prof = x[m].mean(axis=0, keepdims=True) - centroid
        profiles.append(prof.reshape(-1))
        keep.append(p)
    if len(profiles) < 2:
        return pd.DataFrame()
    P = np.vstack(profiles)
    norms = np.linalg.norm(P, axis=1, keepdims=True) + 1e-12
    P = P / norms
    sim = P @ P.T
    return pd.DataFrame(sim, index=keep, columns=keep)


def _build_tf_network(adata, de_tables: dict) -> dict:
    # Build TF -> target edges using the built-in prior, weighted by the
    # mean log2FC of each target across all perturbations that hit it.
    var_set = set(adata.var_names.astype(str))
    edges = []
    tf_nodes = {}
    target_aggregate: dict[str, list[float]] = {}

    # Collapse log2FC across (perturbation, state) for every prior target
    # and report the strongest signed mean.
    for pert, df in de_tables.items():
        if df is None or df.empty:
            continue
        for row in df.itertuples(index=False):
            target_aggregate.setdefault(row.gene, []).append(float(row.log2fc))

    for tf, targets in _BUILTIN_TF_PRIOR.items():
        if tf not in var_set:
            continue
        targets_present = [t for t in targets if t in var_set]
        if not targets_present:
            continue
        tf_nodes[tf] = {"name": tf, "n_targets_in_data": len(targets_present)}
        for t in targets_present:
            agg = target_aggregate.get(t, [])
            mean_lfc = float(np.mean(agg)) if agg else 0.0
            edges.append(
                {
                    "tf": tf,
                    "target": t,
                    "weight": 1.0,
                    "mean_log2fc": mean_lfc,
                    "source": "perturbscope_builtin_prior",
                }
            )

    return {
        "schema": "tf_network/v1",
        "tfs": list(tf_nodes.values()),
        "edges": edges,
    }


def _summary_row(adata, perturb: str, effect_row, traj_row, de_df) -> dict:
    p_mask = adata.obs["perturbation"].values == perturb
    n_cells = int(p_mask.sum())
    classes = adata.obs.loc[p_mask, "perturb_class"] if "perturb_class" in adata.obs.columns else None
    n_perturbed = int((classes == "perturbed").sum()) if classes is not None else 0
    n_escaped = int((classes == "escaped").sum()) if classes is not None else 0

    transcriptional_score = float(effect_row.get("transcriptional_score", 0.0)) if effect_row else 0.0
    state_shift_score = float(effect_row.get("state_shift_score", 0.0)) if effect_row else 0.0
    dominant = effect_row.get("dominant_effect_type", "neither") if effect_row else "neither"

    fate_bias = float(traj_row.get("fate_bias_score", 0.0)) if traj_row else 0.0
    pseudotime_shift = float(traj_row.get("pseudotime_shift", 0.0)) if traj_row else 0.0

    if de_df is not None and not de_df.empty:
        top = de_df.iloc[(-de_df["log2fc"].abs()).argsort()].head(1)
        top_de_gene = str(top["gene"].iloc[0])
        top_de_log2fc = float(top["log2fc"].iloc[0])
        n_significant_de = int((de_df["padj"] < 0.1).sum())
    else:
        top_de_gene = ""
        top_de_log2fc = 0.0
        n_significant_de = 0

    return {
        "perturbation": perturb,
        "n_cells": n_cells,
        "n_perturbed": n_perturbed,
        "n_escaped": n_escaped,
        "transcriptional_score": transcriptional_score,
        "state_shift_score": state_shift_score,
        "dominant_effect": dominant,
        "fate_bias_score": fate_bias,
        "pseudotime_shift": pseudotime_shift,
        "top_de_gene": top_de_gene,
        "top_de_log2fc": top_de_log2fc,
        "n_significant_de": n_significant_de,
    }


# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------
def emit_bundle(
    adata,
    output_dir: str,
    effect_df: Optional[pd.DataFrame] = None,
    traj_df: Optional[pd.DataFrame] = None,
    program_df: Optional[pd.DataFrame] = None,
    interaction_df: Optional[pd.DataFrame] = None,
    top_de_per_pert: int = 200,
    min_de_hits_per_gene: int = 1,
    schema_version: str = "1.0",
) -> Path:
    # Write the v1 results bundle alongside the existing report tree.
    # Returns the bundle root path.
    out = Path(output_dir)
    bundle_root = ensure_dir(out / "bundle")
    embeddings_dir = ensure_dir(bundle_root / "embeddings")
    pert_dir = ensure_dir(bundle_root / "perturbations")
    de_dir = ensure_dir(bundle_root / "de")
    genes_dir = ensure_dir(bundle_root / "genes")

    use_parquet = _try_parquet()

    # ---- metadata ---------------------------------------------------------
    perturbations = sorted([
        str(p) for p in adata.obs["perturbation"].unique()
        if str(p).lower() != "control"
    ])
    metadata = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_perturbations": len(perturbations),
        "perturbations": perturbations,
        "n_controls": int(adata.obs["is_control"].sum()) if "is_control" in adata.obs else 0,
        "cell_states": sorted([str(s) for s in adata.obs.get("cell_state", pd.Series(dtype=str)).unique()]) if "cell_state" in adata.obs else [],
        "qc_summary": {
            "n_pass": int(adata.obs["qc_pass"].sum()) if "qc_pass" in adata.obs else int(adata.n_obs),
            "median_n_genes": float(np.median(adata.obs["n_genes_by_counts"])) if "n_genes_by_counts" in adata.obs else 0.0,
            "median_total_counts": float(np.median(adata.obs["total_counts"])) if "total_counts" in adata.obs else 0.0,
        },
    }
    (bundle_root / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    # ---- embeddings -------------------------------------------------------
    embedding_files = []
    if "X_umap" in adata.obsm:
        um = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
        for col in ["perturbation", "perturb_class", "cell_state", "is_control"]:
            if col in adata.obs.columns:
                um[col] = adata.obs[col].astype(str).values
        name = _save_table(um, embeddings_dir / "umap", use_parquet)
        if name:
            embedding_files.append(f"embeddings/{name}")

    # ---- per-perturbation DEG tables --------------------------------------
    de_tables: dict[str, pd.DataFrame] = {}
    de_files: dict[str, str] = {}
    for perturb in perturbations:
        df = _compute_de_for_perturbation(adata, perturb, top_de_per_pert)
        de_tables[perturb] = df
        if not df.empty:
            written = _save_table(df, de_dir / _safe_filename(perturb), use_parquet)
            if written:
                de_files[perturb] = f"de/{written}"

    # ---- per-perturbation cell tables -------------------------------------
    per_pert_files: dict[str, str] = {}
    cell_cols = [c for c in [
        "perturbation", "is_control", "perturb_class", "perturbation_score",
        "escape_probability", "perturbation_signature",
        "neighbor_control_fraction", "target_expr_reduction",
        "guide_confidence_score", "perturbation_burden", "cell_state",
        "pseudotime", "n_genes_by_counts", "total_counts",
    ] if c in adata.obs.columns]
    for perturb in perturbations:
        m = adata.obs["perturbation"].values == perturb
        if m.sum() == 0:
            continue
        df = adata.obs.loc[m, cell_cols].copy()
        df.insert(0, "cell_id", adata.obs_names[m].astype(str).values)
        if "X_umap" in adata.obsm:
            um = adata.obsm["X_umap"][m]
            df["UMAP1"] = um[:, 0]
            df["UMAP2"] = um[:, 1]
        written = _save_table(df, pert_dir / _safe_filename(perturb), use_parquet)
        if written:
            per_pert_files[perturb] = f"perturbations/{written}"

    # ---- gene inverted index ---------------------------------------------
    gene_index = _build_gene_index(de_tables, min_de_hits_per_gene)
    gene_files: dict[str, str] = {}
    for gene, hits in gene_index.items():
        df = pd.DataFrame(hits)
        written = _save_table(df, genes_dir / _safe_filename(gene), use_parquet)
        if written:
            gene_files[gene] = f"genes/{written}"

    # ---- perturbations_summary --------------------------------------------
    effect_lookup = {row["perturbation"]: row for row in (effect_df.to_dict("records") if effect_df is not None and not effect_df.empty else [])}
    traj_lookup = {row["perturbation"]: row for row in (traj_df.to_dict("records") if traj_df is not None and not traj_df.empty else [])}
    summary_rows = [
        _summary_row(adata, p, effect_lookup.get(p), traj_lookup.get(p), de_tables.get(p))
        for p in perturbations
    ]
    summary_df = pd.DataFrame(summary_rows)
    summary_name = _save_table(summary_df, bundle_root / "perturbations_summary", use_parquet)

    # ---- modules.json -----------------------------------------------------
    modules_payload = {"schema": "modules/v1", "modules": []}
    if program_df is not None and not program_df.empty:
        for pname, group in program_df.groupby("program"):
            modules_payload["modules"].append(
                {
                    "name": str(pname),
                    "n_genes_used": int(group["n_genes_used"].iloc[0]),
                    "perturbations": [
                        {"perturbation": str(r["perturbation"]), "score": float(r["program_score"]) if pd.notna(r["program_score"]) else None}
                        for _, r in group.iterrows()
                    ],
                }
            )
    (bundle_root / "modules.json").write_text(json.dumps(modules_payload, indent=2), encoding="utf-8")

    # ---- tf_network.json --------------------------------------------------
    tf_network = _build_tf_network(adata, de_tables)
    (bundle_root / "tf_network.json").write_text(json.dumps(tf_network, indent=2), encoding="utf-8")

    # ---- perturbation similarity ------------------------------------------
    sim_df = _build_perturbation_similarity(adata, perturbations)
    sim_name = ""
    if not sim_df.empty:
        # Long-format is friendlier for the viewer.
        long = sim_df.stack().reset_index()
        long.columns = ["perturbation_a", "perturbation_b", "cosine"]
        sim_name = _save_table(long, bundle_root / "perturbation_similarity", use_parquet)

    # ---- search index -----------------------------------------------------
    search_index = {
        "schema": "search/v1",
        "items": (
            [{"type": "perturbation", "name": p} for p in perturbations]
            + [{"type": "gene", "name": g, "n_hits": len(gene_index[g])} for g in sorted(gene_index)]
        ),
    }
    (bundle_root / "search_index.json").write_text(json.dumps(search_index, indent=2), encoding="utf-8")

    # ---- manifest ---------------------------------------------------------
    manifest = {
        "schema_version": schema_version,
        "perturbscope_version": __version__,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "format": "parquet" if use_parquet else "csv",
        "artifacts": {
            "metadata": "metadata.json",
            "embeddings": embedding_files,
            "perturbations_summary": summary_name,
            "perturbations": per_pert_files,
            "de": de_files,
            "genes": gene_files,
            "modules": "modules.json",
            "tf_network": "tf_network.json",
            "perturbation_similarity": sim_name,
            "search_index": "search_index.json",
        },
        "counts": {
            "n_perturbations": len(perturbations),
            "n_de_tables": len(de_files),
            "n_gene_pages": len(gene_files),
            "n_tf_edges": len(tf_network.get("edges", [])),
        },
    }
    (bundle_root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return bundle_root


def _safe_filename(name: str) -> Path:
    # Replace filesystem-hostile characters in perturbation/gene names so a
    # combinatorial label like "GENE1+GENE2" still produces a usable filename.
    safe = "".join(c if (c.isalnum() or c in "-_.") else "_" for c in str(name))
    return Path(safe or "unknown")
