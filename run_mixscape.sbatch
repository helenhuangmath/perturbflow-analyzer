#!/bin/bash
# =============================================================================
# run_mixscape.sbatch
#
# Seurat/Mixscape Perturb-seq analysis wrapper.
#
# Usage:
#   sbatch run_mixscape.sbatch <input> <output_dir> [config_json] [perturbation_col]
#
# Input may be a Seurat .rds, .h5Seurat, .h5ad file, or a 10x matrix directory.
# =============================================================================
#SBATCH --job-name=mixscape
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mixscape_%j.out
#SBATCH --error=logs/mixscape_%j.err

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch run_mixscape.sbatch <input> <output_dir> [config_json] [perturbation_col]"
  exit 1
fi

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
R_SCRIPT="${R_SCRIPT:-/vast/projects/wherry/foundation-models-immuno/hhua/tools/r_env/bin/Rscript}"

INPUT="$1"
OUTPUT_DIR="$2"
CONFIG_JSON="${3:-$CODE_DIR/configs/mixscape_default.json}"
PERT_COL="${4:-}"

resolve_path() {
  local raw="$1"; local base="$2"
  if [[ "$raw" = /* ]]; then printf '%s\n' "$raw"; else printf '%s\n' "$base/$raw"; fi
}

INPUT="$(resolve_path "$INPUT" "$CODE_DIR")"
CONFIG_JSON="$(resolve_path "$CONFIG_JSON" "$CODE_DIR")"
OUTPUT_DIR="$(resolve_path "$OUTPUT_DIR" "$CODE_DIR")"

mkdir -p "$OUTPUT_DIR" "$CODE_DIR/logs"
cd "$CODE_DIR"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export R_LIBS_USER="/vast/projects/wherry/foundation-models-immuno/hhua/tools/R/library${R_LIBS_USER:+:$R_LIBS_USER}"
export RETICULATE_PYTHON="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env/bin/python"

CMD=("$R_SCRIPT" "$CODE_DIR/scripts/run_mixscape.R" --input "$INPUT" --output "$OUTPUT_DIR" --config "$CONFIG_JSON")
if [[ -n "$PERT_COL" ]]; then CMD+=(--perturbation-col "$PERT_COL"); fi

echo "============================================================"
echo "Seurat Mixscape run"
echo "Job ID    : ${SLURM_JOB_ID:-local}"
echo "Node      : $(hostname)"
echo "Started   : $(date)"
echo "Input     : $INPUT"
echo "Output    : $OUTPUT_DIR"
echo "Config    : $CONFIG_JSON"
echo "Rscript   : $R_SCRIPT"
echo "Command   : ${CMD[*]}"
echo "============================================================"

"${CMD[@]}"

echo "============================================================"
echo "Finished  : $(date)"
echo "Report    : $OUTPUT_DIR/report.html"
echo "Manifest  : $OUTPUT_DIR/manifest.json"
echo "============================================================"
