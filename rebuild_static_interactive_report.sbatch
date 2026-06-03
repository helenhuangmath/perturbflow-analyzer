#!/bin/bash
#SBATCH --job-name=regen_static_html
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=logs/regen_static_html_%j.out
#SBATCH --error=logs/regen_static_html_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
RESULTS_DIR="${CODE_DIR}/results_replogle_weissman_k562_essential_mini20"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

cd "${CODE_DIR}"
"${ENV_DIR}/bin/python" make_static_interactive_report.py "${RESULTS_DIR}"
