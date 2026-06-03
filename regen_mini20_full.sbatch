#!/bin/bash
# =============================================================================
# regen_mini20_full.sbatch
#
# Regenerate figures + reports for results_replogle_weissman_k562_essential_mini20
# after the round of plotting/threshold patches:
#   - cscore: 1000 HVGs (was 50), threshold 0.4 (unchanged), runs on ALL perts
#   - DEG  : runs on ALL non-control perts (was top-5)
#   - genenet: |r| >= 0.2 (was 0.5), node filter EITHER ctrl OR pert
#   - UMAP : square + Spectral palette
#   - similarity heatmap: bigger axis text
#   - effect_decomposition: theme_classic style
#   - cscore_vs_deg: square + no Pearson r + all 20 perts
#   - live Plotly UMAP: Spectral colors + scaleanchor + bigger dots
#
# Existing per-pert volcano/enrichment files for the original 5 perts will be
# overwritten; new files for the remaining 15 perts will be created.
#
# Usage:
#   cd /vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1
#   sbatch regen_mini20_full.sbatch
# =============================================================================
#SBATCH --job-name=regen_mini20
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/regen_mini20_%j.out
#SBATCH --error=logs/regen_mini20_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
RESULTS_DIR="${CODE_DIR}/results_replogle_weissman_k562_essential_mini20"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/xdg-cache-${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

PYTHON_BIN="${ENV_DIR}/bin/python"
cd "${CODE_DIR}"

# Re-install local package (idempotent) so the patched modules are picked up.
"${PYTHON_BIN}" -m pip install -q -e "${CODE_DIR}" --no-deps

echo "[$(date)] Regenerating mini20 figures + reports"
echo "[$(date)] results_dir = ${RESULTS_DIR}"

"${PYTHON_BIN}" - <<'PYEOF'
import sys
from pathlib import Path

sys.path.insert(0, "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1")

import anndata
import numpy as np

from perturbscope.eda import run_eda
from perturbscope.genenet import run_gene_network
from perturbscope.deg import run_deg_analysis
from perturbscope.cscore import run_cscore
from perturbscope.interactive import build_interactive_report
from perturbscope.report import write_report

results_dir = "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20"
h5ad_path   = results_dir + "/final_adata.h5ad"

print(f"[load] {h5ad_path}")
adata = anndata.read_h5ad(h5ad_path)
print(f"[load] n_obs={adata.n_obs}  n_vars={adata.n_vars}")

# All non-control perts -> DEG + cscore should run on ALL of them so the
# cscore_vs_deg plot has every perturbation.
_CTRL_KW = {"control", "ctrl", "nontargeting", "non-targeting", "nt",
            "scramble", "safe-targeting", "safe_targeting"}
all_perts = sorted({p for p in adata.obs["perturbation"].astype(str)
                    if p.lower() not in _CTRL_KW})
print(f"[plan] {len(all_perts)} non-control perturbations -> "
      f"DEG + cscore on all of them")

# ------------------------------------------------------------------
# 1) EDA  -> regenerates UMAPs (square + Spectral) + similarity heatmap
# ------------------------------------------------------------------
print("[step 1/5] run_eda")
run_eda(
    adata,
    output_dir=results_dir,
    n_top_genes=1000,
    max_cells_heatmap=500,
    n_gene_clusters=10,
    n_pert_clusters=10,
)

# ------------------------------------------------------------------
# 2) genenet  -> regenerates ctrl / pert / diff network panels with
#                |r| >= 0.2 and the EITHER-condition node filter.
# ------------------------------------------------------------------
print("[step 2/5] run_gene_network")
run_gene_network(
    adata,
    results_dir,
    n_top_genes=50,
    n_gene_clusters=5,
    corr_threshold=0.2,
)

# ------------------------------------------------------------------
# 3) DEG  -> per-pert volcanoes + per-pert enrichment + deg_summary.csv
#    on ALL non-control perts.
# ------------------------------------------------------------------
print("[step 3/5] run_deg_analysis on all perts")
run_deg_analysis(
    adata,
    output_dir=results_dir,
    perturbations=all_perts,           # explicit -> bypasses top-N selection
)

# ------------------------------------------------------------------
# 4) cscore  -> 1000 HVGs (the new default), all perts, regenerates
#    cscore_summary.csv + cscore_vs_deg.png + ranked_bar / decomposition.
# ------------------------------------------------------------------
print("[step 4/5] run_cscore on all perts (1000 HVGs)")
run_cscore(
    adata,
    output_dir=results_dir,
    perturbations=all_perts,
    n_top_genes=1000,
    corr_threshold=0.4,
)

# ------------------------------------------------------------------
# 5) Interactive + static reports
# ------------------------------------------------------------------
print("[step 5/5] interactive + static reports")
out = build_interactive_report(adata, results_dir,
                               n_top_genes=1000, max_cells_per_group=50)
print(f"  interactive: {out}")

# Pass effect_df so write_report regenerates effect_decomposition.png with
# the patched theme_classic styling (otherwise the plot is skipped).
import pandas as _pd
_eff_path = results_dir + "/csv/effect_decomposition.csv"
_effect_df = _pd.read_csv(_eff_path) if __import__("os").path.exists(_eff_path) else None
write_report(adata, results_dir, effect_df=_effect_df)
print("  static: report.html")

import os as _os
print(f"  interactive size: {_os.path.getsize(out)/1e6:.1f} MB")

print("[done]")
PYEOF

echo "[$(date)] Done."
