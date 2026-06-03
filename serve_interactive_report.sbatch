#!/bin/bash
# Regenerate the small fetch-based interactive report and serve the results
# directory over HTTP from a compute node.
#SBATCH --job-name=serve_ps_report
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/serve_ps_report_%j.out
#SBATCH --error=logs/serve_ps_report_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
RESULTS_DIR="${CODE_DIR}/results_replogle_weissman_k562_essential_mini20"
PORT="${PORT:-8899}"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

cd "${CODE_DIR}"
"${ENV_DIR}/bin/python" regen_interactive_report.py "${RESULTS_DIR}"

cd "${RESULTS_DIR}"
HOST="$(hostname -f 2>/dev/null || hostname)"
echo "Serving ${RESULTS_DIR}"
echo "Node: ${HOST}"
echo "Port: ${PORT}"
echo "URL on the compute node: http://${HOST}:${PORT}/interactive_report.html"
echo
echo "Forward this port from your workstation, then open:"
echo "  http://localhost:${PORT}/interactive_report.html"
echo
"${ENV_DIR}/bin/python" -m http.server "${PORT}" --bind 0.0.0.0
