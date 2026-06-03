#!/bin/bash
# =============================================================================
# regen_report.sbatch
#
# Re-generate the interactive HTML report from existing results (no re-analysis).
# This only runs the "report" step using the already-computed final_adata.h5ad.
#
# Usage:
#   cd /vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1
#   sbatch regen_report.sbatch
# =============================================================================
#SBATCH --job-name=regen_report
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/regen_report_%j.out
#SBATCH --error=logs/regen_report_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
RESULTS_DIR="${CODE_DIR}/results_test"
INPUT_H5AD="${RESULTS_DIR}/final_adata.h5ad"
CONFIG_JSON="${CODE_DIR}/configs/test.json"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/xdg-cache-${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

PYTHON_BIN="${ENV_DIR}/bin/python"

cd "${CODE_DIR}"

# Re-install local package (idempotent)
"${PYTHON_BIN}" -m pip install -q -e "${CODE_DIR}" --no-deps

echo "[$(date)] Regenerating interactive report from: ${RESULTS_DIR}"
echo "[$(date)] Input h5ad: ${INPUT_H5AD}"

"${PYTHON_BIN}" - <<'PYEOF'
import sys, json
from pathlib import Path
sys.path.insert(0, "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1")

import anndata
from perturbscope.interactive import build_interactive_report

results_dir = "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_test"
h5ad_path   = results_dir + "/final_adata.h5ad"

print(f"Loading {h5ad_path} ...")
adata = anndata.read_h5ad(h5ad_path)
print(f"  {adata.n_obs} cells, {adata.n_vars} genes")
print(f"  obsm keys: {list(adata.obsm.keys())}")
print(f"  obs cols:  {list(adata.obs.columns)[:15]}")

# Re-run EDA to regenerate all plots including the two new clustered heatmaps:
#   eda_clustered_gene_by_pert_heatmap.png
#   eda_gene_pert_cluster_summary_heatmap.png
# and save cluster assignments to csv/eda_gene_clusters.csv + csv/eda_pert_clusters.csv
print("Running EDA step (regenerates all EDA plots + cluster CSVs) ...")
from perturbscope.eda import run_eda
run_eda(
    adata,
    output_dir=results_dir,
    n_top_genes=1000,
    max_cells_heatmap=500,
    n_gene_clusters=10,
    n_pert_clusters=10,
)
print("EDA step done.")

# Re-run genenet so the report uses the current network plot sizing/layout.
print("Running genenet step ...")
from perturbscope.genenet import run_gene_network
run_gene_network(
    adata,
    results_dir,
    n_top_genes=50,
    n_gene_clusters=5,
    corr_threshold=0.4,
)
print("genenet step done.")

# Run cscore step (generates PNGs + CSVs from adata in-memory; reads no prior step output)
print("Running cscore step ...")
from perturbscope.cscore import run_cscore
run_cscore(adata, results_dir, n_top_genes=50, corr_threshold=0.4)
print("cscore step done.")

out = build_interactive_report(adata, results_dir, n_top_genes=1000, max_cells_per_group=50)
print(f"Report written to: {out}")

# Regenerate the static HTML report as well
print("Regenerating static HTML report ...")
from perturbscope.report import write_report
write_report(adata, results_dir)
print("Static report regenerated.")
PYEOF

echo "[$(date)] Done."
