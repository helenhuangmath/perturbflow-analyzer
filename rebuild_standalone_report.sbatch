#!/bin/bash
# Rebuild interactive_report.html for mini20 as a fully standalone file.
# Uses JSON.parse() instead of a raw object literal — 5-15x faster to parse
# in browsers, fixes the "Loading…" hang when opening via file://.
#SBATCH --job-name=regen_standalone_html
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=logs/regen_standalone_html_%j.out
#SBATCH --error=logs/regen_standalone_html_%j.err

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

"${ENV_DIR}/bin/python" make_standalone_report.py "${RESULTS_DIR}"
