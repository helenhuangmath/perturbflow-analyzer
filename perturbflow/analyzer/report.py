# =============================================================================
# perturbflow/analyzer/report.py
#
# Output report generation module.
#
# Collects all analysis results and writes a structured output folder:
#
#   <output_dir>/
#     plots/
#       umap_perturbation.png        -- UMAP coloured by perturbation label
#       perturbation_score_bar.png   -- bar chart of mean perturbation scores
#       effect_decomposition.png     -- scatter: transcriptional vs state-shift
#       trajectory_shift.png         -- bar chart of commitment shift index
#       interaction_scores.png       -- combinatorial interaction bar chart
#     csv/
#       cell_level_summary.csv       -- per-cell QC and scoring metadata
#       effect_decomposition.csv     -- per-perturbation effect scores
#       trajectory_effects.csv       -- per-perturbation trajectory metrics
#       program_scores.csv           -- per-perturbation gene program scores
#       interaction_scores.csv       -- combinatorial interaction results
#     summary.json                   -- dataset-level metadata
#     report.html                    -- minimal HTML index linking to outputs
#
# Plots are only generated when the relevant analysis step has been run.
# =============================================================================

from __future__ import annotations

from pathlib import Path
import json
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .utils import ensure_dir


def _save_df(df: pd.DataFrame, path: Path):
    # Write a DataFrame to CSV, silently skipping if the DataFrame is None or empty.
    if df is None or df.empty:
        return
    df.to_csv(path, index=False)


def write_report(adata, output_dir: str, effect_df=None, traj_df=None, program_df=None, interaction_df=None):
    # Write all plots, CSV tables, summary JSON, and HTML index to output_dir.
    #
    # Args:
    #   adata          -- final annotated AnnData (after all pipeline steps)
    #   output_dir     -- root directory to write outputs into
    #   effect_df      -- DataFrame from compute_effect_decomposition()
    #   traj_df        -- DataFrame from compute_trajectory_effects()
    #   program_df     -- DataFrame from infer_programs()
    #   interaction_df -- DataFrame from analyze_interactions()
    out = ensure_dir(output_dir)
    plots = ensure_dir(out / "plots")
    tables = ensure_dir(out / "csv")

    # -- CSV: per-cell metadata table --
    # Collect whichever QC/scoring columns are present to avoid KeyErrors.
    cols = [
        c
        for c in [
            "perturbation",
            "is_control",
            "qc_pass",
            "perturbation_burden",
            "target_expr_reduction",
            "guide_confidence_score",
            "perturbation_score",
            "escape_probability",
            "perturb_class",
            "cell_state",
            "pseudotime",
        ]
        if c in adata.obs.columns
    ]
    adata.obs[cols].to_csv(tables / "cell_level_summary.csv")

    # -- CSV: analysis result tables --
    _save_df(effect_df, tables / "effect_decomposition.csv")
    _save_df(traj_df, tables / "trajectory_effects.csv")
    _save_df(program_df, tables / "program_scores.csv")
    _save_df(interaction_df, tables / "interaction_scores.csv")

    # -- Plot: UMAP coloured by perturbation label --
    if "X_umap" in adata.obsm:
        um = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
        um["perturbation"] = adata.obs["perturbation"].astype(str).values
        perts = sorted(um["perturbation"].unique())
        # Reuse the EDA UMAP palette (Spectral, control=grey) so the report
        # UMAP matches eda_umap_perturbation byte-for-byte in style.
        from .eda import _spectral_palette
        palette = _spectral_palette(perts)
        ctrl_keys = {"control", "ctrl", "nontargeting", "non-targeting", "nt",
                     "scramble", "safe-targeting", "safe_targeting"}
        perts_ordered = [p for p in perts if p.lower() in ctrl_keys] + \
                        [p for p in perts if p.lower() not in ctrl_keys]
        n_p = len(perts)
        # Fixed legend font so this UMAP matches eda_umap_perturbation +
        # umap_cell_state byte-for-byte in style.
        leg_fs = 24
        n_leg_cols = max(1, min(2, n_p // 25 + 1))

        fig, ax = plt.subplots(figsize=(11, 11))
        for p in perts_ordered:
            mask = um["perturbation"] == p
            ax.scatter(
                um.loc[mask, "UMAP1"], um.loc[mask, "UMAP2"],
                s=3, alpha=0.65, color=palette[p], label=p,
                edgecolors="none", rasterized=True,
            )
        leg = ax.legend(
            title="perturbation", title_fontsize=leg_fs + 1,
            loc="upper left", bbox_to_anchor=(1.02, 1.0),
            fontsize=leg_fs, markerscale=3, ncol=n_leg_cols,
            frameon=True, framealpha=0.9,
        )
        leg.get_frame().set_edgecolor("#888888")
        ax.set_xlabel("UMAP 1", fontsize=36)
        ax.set_ylabel("UMAP 2", fontsize=36)
        ax.tick_params(labelsize=30)
        ax.set_aspect("equal", adjustable="box")  # square plot box
        ax.set_title("UMAP colored by perturbation", fontsize=34)
        fig.tight_layout(rect=[0, 0, 0.76, 1])
        fig.savefig(plots / "umap_perturbation.png", dpi=180)
        plt.close(fig)

    # -- Plot: mean perturbation score per perturbation (top 25) --
    if "perturbation_score" in adata.obs.columns:
        base = pd.DataFrame({
            "perturbation": adata.obs["perturbation"].astype(str).values,
            "perturbation_score": adata.obs["perturbation_score"].astype(float).values,
        })
        agg = base.groupby("perturbation", as_index=False)["perturbation_score"].mean().sort_values("perturbation_score", ascending=False).head(25)
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.barplot(data=agg, x="perturbation", y="perturbation_score", ax=ax, color="#2a9d8f")
        # Anchor tick labels at their right edge so they line up under each
        # bar instead of drifting to the right of it.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")
        ax.set_xlabel("")
        ax.set_ylabel("Mean perturbation score")
        ax.set_title("Mean perturbation score")
        fig.tight_layout()
        fig.savefig(plots / "perturbation_score_bar.png", dpi=180)
        plt.close(fig)

    # -- Plot: effect decomposition scatter (transcriptional vs state-shift) --
    # ggplot2 theme_classic() equivalent: white background, no grid lines of
    # any kind, only bottom + left axis spines, ticks pointing outward, no
    # legend frame. Wrapped in a seaborn axes_style context so the surrounding
    # plotting code is unaffected.
    if effect_df is not None and not effect_df.empty:
        with sns.axes_style("ticks", {"axes.grid": False, "grid.color": "white",
                                        "axes.edgecolor": "black"}):
            fig, ax = plt.subplots(figsize=(7, 6))
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")
            sns.scatterplot(
                data=effect_df,
                x="transcriptional_score",
                y="state_shift_score",
                hue="dominant_effect_type",
                ax=ax,
            )
            ax.set_title("Effect decomposition")
            ax.grid(False)
            ax.set_axisbelow(False)
            sns.despine(ax=ax)  # remove top + right spines
            ax.spines["left"].set_color("black")
            ax.spines["bottom"].set_color("black")
            ax.tick_params(direction="out", colors="black")
            leg = ax.get_legend()
            if leg is not None:
                leg.set_frame_on(False)
            fig.tight_layout()
            fig.savefig(plots / "effect_decomposition.png", dpi=180,
                        facecolor="white")
            plt.close(fig)

    # -- Plot: commitment shift index bar chart (top 25 perturbations) --
    if traj_df is not None and not traj_df.empty:
        fig, ax = plt.subplots(figsize=(10, 5))
        top = traj_df.sort_values("commitment_shift_index", ascending=False).head(25)
        sns.barplot(data=top, x="perturbation", y="commitment_shift_index", ax=ax, color="#e76f51")
        ax.tick_params(axis="x", rotation=60)
        ax.set_title("Commitment shift index")
        fig.tight_layout()
        fig.savefig(plots / "trajectory_shift.png", dpi=180)
        plt.close(fig)

    # -- Plot: combinatorial interaction scores --
    if interaction_df is not None and not interaction_df.empty:
        fig, ax = plt.subplots(figsize=(9, 5))
        sns.barplot(data=interaction_df, x="perturbation", y="interaction_score", hue="interaction_class", ax=ax)
        ax.tick_params(axis="x", rotation=60)
        ax.set_title("Combinatorial interaction scores")
        fig.tight_layout()
        fig.savefig(plots / "interaction_scores.png", dpi=180)
        plt.close(fig)

    # -- summary.json: dataset-level counts and file manifest --
    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_perturbations": int(adata.obs["perturbation"].nunique()) if "perturbation" in adata.obs else 0,
        "available_tables": sorted([p.name for p in tables.glob("*.csv")]),
        "available_plots": sorted([p.name for p in plots.glob("*.png")]),
    }
    with open(out / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    # -- report.html: minimal landing page linking to all outputs --
    html = out / "report.html"

    # Build dynamic plot sections so the page adapts to which steps ran.
    def _img_tag(name: str, caption: str) -> str:
        p = plots / name
        if not p.exists():
            return ""
        return f'<figure><img src="plots/{name}" style="max-width:100%"><figcaption>{caption}</figcaption></figure>'

    eda_section = f"""
<h2>Exploratory Data Analysis</h2>
{_img_tag("eda_cells_per_perturbation.png", "Cells per perturbation group")}
{_img_tag("eda_cluster_proportions.png", "Cell-state cluster proportions per perturbation")}
{_img_tag("eda_gene_by_cell_heatmap.png", "Gene × Cell expression heatmap (top HVGs)")}
{_img_tag("eda_clustered_gene_by_pert_heatmap.png", "All Genes × Perturbations Clusters")}
{_img_tag("eda_gene_correlation_heatmap.png", "Gene–gene Pearson correlation heatmap")}
{_img_tag("eda_perturbation_similarity.png", "Perturbation–perturbation cosine similarity (functionally similar groups cluster together)")}
{_img_tag("eda_gene_pert_cluster_summary_heatmap.png", "Gene expression cluster × Perturbation cluster summary heatmap (mean z-score per cluster pair)")}
<p>Cluster assignments: <code>csv/eda_gene_clusters.csv</code> &nbsp;|&nbsp; <code>csv/eda_pert_clusters.csv</code></p>
{_img_tag("umap_cell_state.png", "UMAP colored by cell-state cluster")}
{_img_tag("eda_umap_perturbation.png", "UMAP colored by perturbation")}
"""
    # Gene network section (shown only when genenet step ran).
    genenet_ctrl = _img_tag("genenet_control_heatmap.png", "Gene co-expression clusters \u2014 control cells")
    genenet_nets = "".join(
        _img_tag(p.name, p.stem.replace("genenet_", "").replace("_", " "))
        for p in sorted(plots.glob("genenet_*.png"))
        if p.name != "genenet_control_heatmap.png"
    )
    genenet_section = f"""
<h2>Gene Co-expression Networks</h2>
{genenet_ctrl}
{genenet_nets}
<p>Gene cluster assignments: <code>csv/genenet_gene_clusters.csv</code></p>
""" if genenet_ctrl.strip() or genenet_nets.strip() else ""
    # Collect DEG volcano plots dynamically.
    volcano_imgs = "".join(
        _img_tag(p.name, f"Volcano: {p.stem.replace('deg_volcano_', '')}")
        for p in sorted(plots.glob("deg_volcano_*.png"))
    )
    deg_section = f"""
<h2>Differential Expression Analysis</h2>
{_img_tag("deg_top_perturbations_heatmap.png", "log₂FC heatmap: top perturbations vs control")}
{volcano_imgs}
<p>DEG tables: <code>csv/deg_*.csv</code> &nbsp;|&nbsp; Summary: <code>csv/deg_summary.csv</code></p>
"""

    effects_section = f"""
<h2>Effect Decomposition &amp; Downstream Analysis</h2>
{_img_tag("umap_perturbation.png", "UMAP colored by perturbation")}
{_img_tag("perturbation_score_bar.png", "Mean perturbation score (top 25)")}
{_img_tag("effect_decomposition.png", "Effect decomposition: transcriptional vs state-shift")}
{_img_tag("trajectory_shift.png", "Commitment shift index (top 25)")}
{_img_tag("interaction_scores.png", "Combinatorial interaction scores")}

<h2>TF Regulatory Analysis</h2>
{_img_tag("tf_regulatory_heatmap.png", "TF regulatory heatmap: signed −log₁₀(q) for each (perturbed TF, target gene). Red = KO ↑ target (TF normally inhibits); Green = KO ↓ target (TF normally activates).")}
{_img_tag("tf_regulatory_network.png", "TF regulatory network: directed edges from each perturbed TF to its strongest targets. Node size = total in+out degree.")}
"""

    html.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>PerturbFlow-Analyzer Report</title>
  <style>
    body {{ font-family: sans-serif; max-width: 1100px; margin: 2rem auto; color: #222; }}
    h1   {{ border-bottom: 2px solid #2a9d8f; padding-bottom: .4rem; }}
    h2   {{ color: #2a9d8f; margin-top: 2rem; }}
    figure {{ margin: 1rem 0; }}
    figcaption {{ font-size: .85rem; color: #555; margin-top: .3rem; }}
    code {{ background: #f4f4f4; padding: 2px 5px; border-radius: 3px; }}
    .interactive-link {{ display:inline-block; margin:.5rem 0; padding:.5rem 1rem;
      background:#2a9d8f; color:#fff; border-radius:6px; text-decoration:none;
      font-weight:600; font-size:.9rem; }}
    .interactive-link:hover {{ background:#1e7268; }}
  </style>
</head>
<body>
<h1>PerturbFlow-Analyzer Analysis Report</h1>
<p>
  <strong>Cells:</strong> {int(adata.n_obs):,} &nbsp;|&nbsp;
  <strong>Genes:</strong> {int(adata.n_vars):,} &nbsp;|&nbsp;
  <strong>Perturbations:</strong> {int(adata.obs["perturbation"].nunique()) if "perturbation" in adata.obs else 0}
</p>
<p><a class="interactive-link" href="interactive_report.html">Open Interactive Report</a></p>
<p>Full artifacts: <code>plots/</code>, <code>csv/</code>, <code>bundle/</code>, <code>final_adata.h5ad</code></p>

{eda_section}
{deg_section}
{genenet_section}
{effects_section}
</body>
</html>""",
        encoding="utf-8",
    )

    return out
