#!/bin/bash
#SBATCH --job-name=regen_cluster_props
#SBATCH --partition=genoa-std-mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=logs/regen_cluster_props_%j.out
#SBATCH --error=logs/regen_cluster_props_%j.err

set -euo pipefail

CODE_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1"
ENV_DIR="/vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env"
RESULTS_DIR="${CODE_DIR}/results_replogle_weissman_k562_essential_mini20"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/xdg-cache-${SLURM_JOB_ID:-local}"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}" "${CODE_DIR}/logs"

source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate "${ENV_DIR}"

cd "${CODE_DIR}"

"${ENV_DIR}/bin/python" - <<'PYEOF'
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

results_dir = Path("/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20")
csv_path = results_dir / "csv" / "eda_cluster_proportions.csv"
plots_dir = results_dir / "plots"

ct = pd.read_csv(csv_path, index_col=0)
state_cols = sorted(
    [str(c) for c in ct.columns],
    key=lambda x: int(x) if x.isdigit() else x,
)
ct = ct.loc[:, state_cols]

states = [str(c) for c in ct.columns]
state_palette = dict(zip(states, sns.color_palette("Spectral", max(2, len(states)))))

n_perts = len(ct)
fig_width = max(10, min(0.45 * n_perts, 50))
fig, ax = plt.subplots(figsize=(fig_width, 5))
bottom = np.zeros(n_perts)

for col in ct.columns:
    state = str(col)
    ax.bar(
        range(n_perts),
        ct[col].values,
        bottom=bottom,
        label=state,
        color=state_palette[state],
        edgecolor="none",
    )
    bottom += ct[col].values

ax.set_xticks(range(n_perts))
ax.set_xticklabels(ct.index.tolist(), rotation=90, fontsize=9)
ax.tick_params(axis="y", labelsize=9)
ax.set_ylabel("Proportion", fontsize=10)
ax.set_title("Cell-state cluster proportions per perturbation", fontsize=12)
ax.legend(
    title="state",
    title_fontsize=9,
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    fontsize=9,
)
fig.tight_layout()

for name in ("eda_cluster_proportions.png", "eda_cluster_proportions_states_spectral.png"):
    out = plots_dir / name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"[write] {out}")

plt.close(fig)
PYEOF
