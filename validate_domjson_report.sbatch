#!/bin/bash
#SBATCH --job-name=validate_domjson_html
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=logs/validate_domjson_html_%j.out
#SBATCH --error=logs/validate_domjson_html_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
REPORT="${CODE_DIR}/results_replogle_weissman_k562_essential_mini20/interactive_report.html"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

cd "${CODE_DIR}"
"${ENV_DIR}/bin/python" validate_domjson_report.py "${REPORT}"
