# PerturbFlow-Analyzer Quick Start

PerturbFlow-Analyzer takes a Perturb-seq `h5ad` file, standardizes perturbation labels,
runs QC, scoring, DEG, network, regulatory, C-score, and report steps, then
writes an interactive HTML report plus a viewer-ready results bundle.

## 1. Install

```bash
cd /vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbFlow
python -m pip install -e .
```

On the Wherry cluster, use the existing environment:

```bash
source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate /vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env
python -m pip install -e /vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbFlow --no-deps
```

## 2. Prepare Your Data

Your input must be an AnnData `.h5ad` file with cells in rows, genes in columns,
and an `obs` column containing perturbation labels.

```bash
perturbflow prepare \
  --input my_raw_data.h5ad \
  --output prepared/my_data.perturbflow.h5ad \
  --perturbation-col guide_gene \
  --control-labels control,non-targeting,NT \
  --cell-state-col leiden
```

This writes:

- `obs["perturbation"]`: standardized perturbation labels
- `obs["perturbation_original"]`: original labels
- `obs["cell_state"]`: optional copied cell-state labels

## 3. Run The Pipeline

```bash
perturbflow run \
  --input prepared/my_data.perturbflow.h5ad \
  --output results/my_run \
  --config configs/cluster_default.json \
  --resume
```

Useful step-only reruns:

```bash
perturbflow run --input prepared/my_data.perturbflow.h5ad --output results/my_run --steps deg,report,bundle
perturbflow run --input prepared/my_data.perturbflow.h5ad --output results/my_run --force-steps report --resume
perturbflow list-steps
```

## 4. Run On Slurm

Use the mini20 scripts as templates:

```bash
sbatch regen_mini20_full.sbatch
sbatch regen_mini20_html_only.sbatch
```

For a new dataset, copy one of these scripts and update:

- `CODE_DIR`
- `ENV_DIR`
- `RESULTS_DIR`
- input `.h5ad`
- requested resources and wall time

## 5. Main Outputs

After the run, open:

```text
results/my_run/interactive_report.html
```

Other useful outputs:

- `results/my_run/report.html`: static report
- `results/my_run/csv/`: DEG, enrichment, network, regulatory, C-score tables
- `results/my_run/plots/`: static PNG figures
- `results/my_run/bundle/`: compact viewer-ready data bundle
- `results/my_run/checkpoint.json`: completed steps for resume

## 6. Expected Input Columns

Minimum:

- perturbation label column, passed with `--perturbation-col`

Optional but useful:

- cell-state or cluster column, passed with `--cell-state-col`
- QC metrics such as mitochondrial fraction or total counts

PerturbFlow-Analyzer will compute core QC/preprocessing features when they are missing.
