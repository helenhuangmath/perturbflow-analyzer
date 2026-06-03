#!/bin/bash
# One-off: re-render effect_decomposition.png with the patched theme_classic
# styling, using the existing csv/effect_decomposition.csv.
#SBATCH --job-name=regen_eff
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=logs/regen_eff_%j.out
#SBATCH --error=logs/regen_eff_%j.err

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
import pandas as pd

from perturbscope.report import write_report

results_dir = "/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20"
adata = anndata.read_h5ad(results_dir + "/final_adata.h5ad")
effect_df = pd.read_csv(results_dir + "/csv/effect_decomposition.csv")
print(f"[load] adata n_obs={adata.n_obs}; effect_df rows={len(effect_df)}")

# Single call: write_report re-emits effect_decomposition.png + report.html
write_report(adata, results_dir, effect_df=effect_df)
print("[done]")
PYEOF
