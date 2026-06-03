# =============================================================================
# perturbflow/analyzer/pipeline.py  (v1)
#
# Top-level orchestrator. v1 keeps the original step graph and appends a new
# terminal `bundle` step that emits a versioned, viewer-ready results bundle.
#
# Step dependency order:
#   ingest -> qc -> preprocess -> score -> effects
#                                       -> trajectory
#                                       -> programs
#                                       -> interaction -> report -> bundle
# =============================================================================

from __future__ import annotations

import json
import datetime
from pathlib import Path
from typing import Iterable, Set

from .config import PipelineConfig
from .io import read_data
from .qc import run_qc, plot_qc_summary
from .preprocessing import normalize_and_embed
from .perturbation import score_effective_perturbation
from .state import compute_effect_decomposition
from .trajectory import compute_trajectory_effects
from .programs import infer_programs
from .interaction import analyze_interactions
from .eda import run_eda
from .deg import run_deg_analysis, identify_top_perturbations
from .genenet import run_gene_network, run_tf_gene_network
from .cscore import run_cscore
from .state_enrich import run_state_enrichment
from .regulatory import run_regulatory_analysis
from .report import write_report
from .interactive import build_interactive_report
from .bundle import emit_bundle
from .utils import ensure_dir


def _prepare_for_h5ad_write(adata):
    # HDF5 cannot store python lists in .obs or arbitrary dicts/lists in .uns.
    # Coerce them to strings so the write succeeds without losing information.
    for col in adata.obs.columns:
        series = adata.obs[col]
        if series.dtype != object:
            continue
        sample = next((value for value in series if value is not None), None)
        if isinstance(sample, list):
            adata.obs[col] = series.apply(
                lambda value: "|".join(map(str, value)) if isinstance(value, list) else value
            )
    for key, value in list(adata.uns.items()):
        if isinstance(value, (list, dict)):
            adata.uns[key] = json.dumps(value, default=str)
    return adata


def _load_checkpoint(out: Path) -> Set[str]:
    """Return the set of step names that have already completed."""
    ckpt = out / "checkpoint.json"
    if ckpt.exists():
        try:
            return set(json.loads(ckpt.read_text()).get("completed_steps", []))
        except Exception:
            pass
    return set()


def _save_checkpoint(out: Path, completed: Set[str]) -> None:
    """Persist the current set of completed steps to checkpoint.json."""
    (out / "checkpoint.json").write_text(
        json.dumps(
            {
                "completed_steps": sorted(completed),
                "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
            },
            indent=2,
        )
    )


def run_pipeline(
    input_path: str,
    output_dir: str,
    config: PipelineConfig | None = None,
    steps: Iterable[str] | None = None,
    perturbation_col: str | None = None,
    resume: bool = False,
    force_steps: Iterable[str] | None = None,
    clear_from: str | None = None,
):
    cfg = config or PipelineConfig()
    active_steps = list(steps) if steps else cfg.default_steps
    forced: Set[str] = set(force_steps) if force_steps else set()

    out = ensure_dir(output_dir)
    ensure_dir(Path(out) / "plots")
    ensure_dir(Path(out) / "csv")

    # ── Checkpoint handling ────────────────────────────────────────────────
    # clear_from: invalidate checkpoint from that step onward.
    completed: Set[str] = _load_checkpoint(out) if resume else set()
    if clear_from and clear_from in active_steps:
        clear_idx = active_steps.index(clear_from)
        to_clear = set(active_steps[clear_idx:])
        completed -= to_clear
        _save_checkpoint(out, completed)

    def _should_run(step: str) -> bool:
        if step not in active_steps:
            return False
        if step in forced:
            return True
        if step in completed:
            print(f"[checkpoint] skipping '{step}' (already completed)")
            return False
        return True

    adata = read_data(input_path, perturbation_col=perturbation_col)

    effect_df = None
    traj_df = None
    program_df = None
    interaction_df = None
    deg_results = None
    top_perts = None

    if _should_run("qc"):
        adata = run_qc(
            adata,
            min_genes=cfg.min_genes,
            max_pct_mt=cfg.max_pct_mt,
            min_cells_per_perturbation=cfg.min_cells_per_perturbation,
        )
        plot_qc_summary(adata, str(Path(out) / "plots" / "qc_summary.png"))
        completed.add("qc"); _save_checkpoint(out, completed)

    if _should_run("preprocess"):
        adata = normalize_and_embed(
            adata,
            random_state=cfg.random_state,
            leiden_resolution=cfg.leiden_resolution,
        )
        completed.add("preprocess"); _save_checkpoint(out, completed)

    # EDA runs right after normalisation so analysts can inspect the data
    # (cell counts, heatmaps, correlation) before the heavier scoring steps.
    if _should_run("eda"):
        run_eda(
            adata,
            output_dir=str(out),
            n_top_genes=cfg.eda_n_top_genes,
            max_cells_heatmap=cfg.eda_max_cells_heatmap,
        )
        completed.add("eda"); _save_checkpoint(out, completed)

    if _should_run("score"):
        adata = score_effective_perturbation(adata, n_neighbors=cfg.n_neighbors)
        completed.add("score"); _save_checkpoint(out, completed)

    if _should_run("effects"):
        adata, effect_df = compute_effect_decomposition(adata)
        completed.add("effects"); _save_checkpoint(out, completed)

    if _should_run("trajectory"):
        adata, traj_df = compute_trajectory_effects(adata)
        completed.add("trajectory"); _save_checkpoint(out, completed)

    if _should_run("programs"):
        adata, program_df = infer_programs(adata)
        completed.add("programs"); _save_checkpoint(out, completed)

    if _should_run("interaction"):
        adata, interaction_df = analyze_interactions(adata)
        completed.add("interaction"); _save_checkpoint(out, completed)

    if _should_run("state_enrich"):
        run_state_enrichment(
            adata,
            output_dir=str(out),
            min_cells=cfg.state_enrich_min_cells,
            fdr_threshold=cfg.state_enrich_fdr_threshold,
        )
        completed.add("state_enrich"); _save_checkpoint(out, completed)

    # DEG analysis runs after effects so that identify_top_perturbations can
    # use effect_df for smarter ranking when available.
    if _should_run("deg"):
        deg_results = run_deg_analysis(
            adata,
            output_dir=str(out),
            effect_df=effect_df,
            n_top_perturbations=cfg.deg_n_top_perturbations,
            logfc_threshold=cfg.deg_logfc_threshold,
            pval_threshold=cfg.deg_pval_threshold,
            n_top_deg_heatmap=cfg.deg_n_top_heatmap,
        )
        completed.add("deg"); _save_checkpoint(out, completed)

    # Gene network runs after DEG so we can reuse the same top-perturbation
    # ranking.  Falls back gracefully when effect_df is not available.
    if _should_run("genenet"):
        top_perts = identify_top_perturbations(
            adata,
            effect_df=effect_df,
            n_top=cfg.deg_n_top_perturbations,
        )
        run_gene_network(
            adata,
            output_dir=str(out),
            perturbations=top_perts,
            n_top_genes=cfg.genenet_n_top_genes,
            n_gene_clusters=cfg.genenet_n_gene_clusters,
            corr_threshold=cfg.genenet_corr_threshold,
        )
        completed.add("genenet"); _save_checkpoint(out, completed)

    # TF-anchored gene network: top-N highly-variable TFs (HumanTFs) plus
    # their top correlated partner genes in control cells. Reuses the same
    # top-perturbation list so panels are comparable to genenet output.
    if _should_run("tf_genenet"):
        # Always recompute top_perts so this step also works when genenet
        # was already checkpointed and skipped.
        top_perts = identify_top_perturbations(
            adata,
            effect_df=effect_df,
            n_top=cfg.deg_n_top_perturbations,
        )
        run_tf_gene_network(
            adata,
            output_dir=str(out),
            perturbations=top_perts,
            n_top_tfs=cfg.tf_genenet_n_top_tfs,
            n_partners_per_tf=cfg.tf_genenet_n_partners_per_tf,
            corr_threshold=cfg.tf_genenet_corr_threshold,
            tf_list_path=cfg.tf_genenet_tf_list_path,
            candidate_gene_pool_size=cfg.tf_genenet_candidate_pool_size,
            n_gene_modules=cfg.tf_genenet_n_modules,
            max_ctrl_cells=cfg.tf_genenet_max_ctrl_cells,
        )
        completed.add("tf_genenet"); _save_checkpoint(out, completed)

    # C-score: connectivity rewiring metric. Uses its own (larger) HVG panel
    # so the network has enough edges for C_gain/C_loss/C_shift to be stable.
    if _should_run("cscore"):
        run_cscore(
            adata,
            output_dir=str(out),
            perturbations=top_perts,
            n_top_genes=cfg.cscore_n_top_genes,
            corr_threshold=cfg.cscore_corr_threshold,
        )
        completed.add("cscore"); _save_checkpoint(out, completed)

    if _should_run("regulatory"):
        run_regulatory_analysis(
            adata,
            output_dir=str(out),
            fdr_threshold=cfg.regulatory_fdr_threshold,
            lfc_threshold=cfg.regulatory_lfc_threshold,
        )
        completed.add("regulatory"); _save_checkpoint(out, completed)

    if _should_run("report"):
        write_report(
            adata,
            output_dir=str(out),
            effect_df=effect_df,
            traj_df=traj_df,
            program_df=program_df,
            interaction_df=interaction_df,
        )
        # Always build the interactive report alongside the static one.
        build_interactive_report(
            adata,
            output_dir=str(out),
            n_top_genes=cfg.eda_n_top_genes,
        )
        completed.add("report"); _save_checkpoint(out, completed)

    # v1: emit the viewer-ready bundle. Always last so it can include every
    # artifact produced by earlier steps; gracefully skips per-artifact when
    # an upstream step was disabled.
    if _should_run("bundle"):
        emit_bundle(
            adata,
            output_dir=str(out),
            effect_df=effect_df,
            traj_df=traj_df,
            program_df=program_df,
            interaction_df=interaction_df,
            top_de_per_pert=cfg.bundle_top_de_per_pert,
            min_de_hits_per_gene=cfg.bundle_min_de_hits_per_gene,
            schema_version=cfg.bundle_schema_version,
        )
        completed.add("bundle"); _save_checkpoint(out, completed)

    adata = _prepare_for_h5ad_write(adata)
    adata.write_h5ad(Path(out) / "final_adata.h5ad")
    return adata
