#!/bin/bash
# =============================================================================
# run_perturbflow_analyzer.sbatch
#
# SLURM batch script for PerturbFlow Analyzer (with viewer-ready bundle step).
#
# Usage:
#   sbatch run_perturbflow_analyzer.sbatch <input_h5ad> <output_dir> <config_json> [steps_csv] [perturbation_col]
#
# Reuses the existing conda env at /vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env
# and pip-installs pyarrow on first run (idempotent — fast no-op on subsequent runs).
# =============================================================================
#SBATCH --job-name=perturbflow_analyzer
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/perturbflow_analyzer_%j.out
#SBATCH --error=logs/perturbflow_analyzer_%j.err

set -euo pipefail

# Per-job cache directories so parallel jobs do not clash on shared NFS.
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/xdg-cache-${SLURM_JOB_ID:-local}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

if [[ $# -lt 3 ]]; then
  echo "Usage: sbatch run_perturbflow_analyzer.sbatch <input_h5ad> <output_dir> <config_json> [steps_csv] [perturbation_col] [--resume] [force_steps_csv]"
  exit 1
fi

INPUT_H5AD="$1"
OUTPUT_DIR="$2"
CONFIG_JSON="$3"
STEPS_CSV="${4:-qc,preprocess,eda,score,effects,trajectory,programs,interaction,state_enrich,deg,genenet,regulatory,report,bundle}"
PERT_COL="${5:-}"
RESUME_FLAG="${6:-}"        # pass --resume to skip already-completed steps
FORCE_STEPS="${7:-}"        # comma-separated steps to re-run even if checkpointed

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbFlow"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"

# Resolve relative paths against CODE_DIR for inputs / configs and against
# the submit directory for outputs.
resolve_path() {
  local raw="$1"; local base="$2"
  if [[ "$raw" = /* ]]; then printf '%s\n' "$raw"; else printf '%s\n' "$base/$raw"; fi
}
INPUT_H5AD="$(resolve_path "$INPUT_H5AD" "$CODE_DIR")"
CONFIG_JSON="$(resolve_path "$CONFIG_JSON" "$CODE_DIR")"
OUTPUT_DIR="$(resolve_path "$OUTPUT_DIR" "$CODE_DIR")"

mkdir -p "$OUTPUT_DIR" "$CODE_DIR/logs"
cd "$CODE_DIR"

# Activate conda env.
source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "$ENV_DIR"

PYTHON_BIN="${ENV_DIR}/bin/python"

# Ensure pyarrow is installed for parquet bundle output. Idempotent on retry.
if ! "$PYTHON_BIN" -c "import pyarrow" 2>/dev/null; then
  echo "[setup] installing pyarrow into env"
  "$PYTHON_BIN" -m pip install --quiet pyarrow
fi

export PYTHONPATH="$CODE_DIR${PYTHONPATH:+:$PYTHONPATH}"

CMD=("$PYTHON_BIN" -m perturbflow.analyzer.cli run --input "$INPUT_H5AD" --output "$OUTPUT_DIR" --config "$CONFIG_JSON" --steps "$STEPS_CSV")
if [[ -n "$PERT_COL" ]];    then CMD+=(--perturbation-col "$PERT_COL"); fi
if [[ -n "$RESUME_FLAG" ]]; then CMD+=($RESUME_FLAG); fi
if [[ -n "$FORCE_STEPS" ]]; then CMD+=(--force-steps "$FORCE_STEPS"); fi

echo "============================================================"
echo "PerturbFlow Analyzer run"
echo "Job ID    : ${SLURM_JOB_ID:-local}"
echo "Node      : $(hostname)"
echo "Started   : $(date)"
echo "Input     : $INPUT_H5AD"
echo "Output    : $OUTPUT_DIR"
echo "Config    : $CONFIG_JSON"
echo "Steps     : $STEPS_CSV"
echo "Resume    : ${RESUME_FLAG:-no}"
echo "ForceSteps: ${FORCE_STEPS:-none}"
echo "Python    : $PYTHON_BIN"
echo "Command   : ${CMD[*]}"
echo "============================================================"

"${CMD[@]}"

echo "============================================================"
echo "Finished  : $(date)"
echo "Bundle    : $OUTPUT_DIR/bundle"
echo "Manifest  : $OUTPUT_DIR/bundle/manifest.json"
echo "============================================================"
