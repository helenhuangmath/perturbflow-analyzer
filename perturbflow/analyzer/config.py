# =============================================================================
# perturbflow/analyzer/config.py (v1)
#
# Single source of truth for tunable parameters. v1 adds the `bundle` step,
# which writes a versioned, viewer-ready results bundle (parquet + JSON) to
# the output directory in addition to the legacy plots/CSV/HTML.
# =============================================================================

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import List


@dataclass
class PipelineConfig:
    min_genes: int = 200
    max_pct_mt: float = 20.0
    min_cells_per_perturbation: int = 20
    n_neighbors: int = 20
    n_top_genes_signature: int = 50
    random_state: int = 0
    leiden_resolution: float = 0.2

    # Bundle-specific knobs (v1).
    # Maximum number of DEGs to keep per (perturbation, cell_state) table.
    bundle_top_de_per_pert: int = 200
    # Genes appearing in fewer than this many DEG tables are pruned from the
    # gene index (keeps the gene browse list focused on responsive genes).
    bundle_min_de_hits_per_gene: int = 1
    # Schema version stamped into manifest.json for forward compatibility.
    bundle_schema_version: str = "1.0"

    # EDA knobs — controls the exploratory-analysis step that runs after
    # normalisation and produces cell-count bar charts, three heatmaps,
    # and cluster proportion plots.
    eda_n_top_genes: int = 1000
    eda_max_cells_heatmap: int = 500

    # DEG knobs — controls the differential-expression step that runs after
    # the effects step and produces per-perturbation volcano plots, DEG tables,
    # and a top-perturbations summary heatmap.
    deg_n_top_perturbations: int = 5
    deg_logfc_threshold: float = 0.5
    deg_pval_threshold: float = 0.05
    deg_n_top_heatmap: int = 20

    # Gene network knobs — controls the genenet step that produces co-expression
    # network graphs (control vs perturbed) and gene-cluster heatmaps.
    genenet_n_top_genes: int = 50
    genenet_n_gene_clusters: int = 5
    genenet_corr_threshold: float = 0.2

    # TF-anchored gene network knobs — top-N highly-variable TFs (HumanTFs)
    # plus their top correlated partner genes in control cells. Same edge
    # rule (|r| ≥ tf_genenet_corr_threshold) as genenet.
    tf_genenet_n_top_tfs: int = 10
    tf_genenet_n_partners_per_tf: int = 8
    tf_genenet_corr_threshold: float = 0.5
    tf_genenet_candidate_pool_size: int = 2000
    tf_genenet_n_modules: int = 5
    tf_genenet_max_ctrl_cells: int = 20000
    tf_genenet_tf_list_path: str | None = None

    # C-score knobs — connectivity rewiring score derived from the gene
    # co-expression networks. Uses its own (larger) HVG panel rather than the
    # genenet one, because cscore is a metric that needs hundreds-to-thousands
    # of edges for stable C_gain / C_loss / C_shift estimates, whereas genenet
    # plots want a small readable network.
    cscore_n_top_genes: int = 1000
    cscore_corr_threshold: float = 0.4

    # Cell state enrichment knobs — chi-square test of perturbation × cell state
    # composition relative to control cells.
    state_enrich_min_cells: int = 5
    state_enrich_fdr_threshold: float = 0.05

    # Regulatory network knobs — TF-TF regulatory heatmap and directed network
    # derived from DEG log2FC / padj tables.
    regulatory_fdr_threshold: float = 0.05
    regulatory_lfc_threshold: float = 0.3

    default_steps: List[str] = field(
        default_factory=lambda: [
            "ingest",
            "qc",
            "preprocess",
            "eda",
            "score",
            "effects",
            "trajectory",
            "programs",
            "interaction",
            "state_enrich",
            "deg",
            "genenet",
            "tf_genenet",
            "cscore",
            "regulatory",
            "report",
            "bundle",
        ]
    )

    @classmethod
    def from_json(cls, path: str | Path) -> "PipelineConfig":
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        valid = {k: v for k, v in payload.items() if k in cls.__dataclass_fields__}
        return cls(**valid)

    def to_json(self, path: str | Path) -> None:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(self.__dict__, handle, indent=2)
