#!/bin/bash
# One-off: rebuild ONLY the interactive_report.html for mini20 with the
# slimmed-down bundler (relative-path PNGs + DEG head(100)). No re-analysis.
#SBATCH --job-name=regen_html
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=logs/regen_html_%j.out
#SBATCH --error=logs/regen_html_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
RESULTS_DIR="${CODE_DIR}/results_replogle_weissman_k562_essential_mini20"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/xdg-cache-${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

cd "${CODE_DIR}"
"${ENV_DIR}/bin/python" -m pip install -q -e "${CODE_DIR}" --no-deps

"${ENV_DIR}/bin/python" - <<'PYEOF'
import sys
sys.path.insert(0, "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1")

import anndata
from perturbscope.interactive import build_interactive_report

results_dir = "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20"
adata = anndata.read_h5ad(results_dir + "/final_adata.h5ad")
print(f"[load] n_obs={adata.n_obs}  n_vars={adata.n_vars}")

out = build_interactive_report(adata, results_dir,
                               n_top_genes=1000, max_cells_per_group=50)
print(f"[done] {out}")
import os
print(f"[size] {os.path.getsize(out)/1e6:.1f} MB")
PYEOF
