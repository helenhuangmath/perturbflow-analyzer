# PerturbFlow Quick Start

PerturbFlow is open, AnnData-native infrastructure for perturbation biology.
The current `analyzer` subpackage standardizes perturbation labels, runs QC, scoring, DEG, network,
regulatory, C-score, and report steps, then writes an interactive HTML report
plus a viewer-ready results bundle.

## 1. Install

```bash
cd PerturbFlow
python -m pip install -e ".[bundle]"
```

On the Wherry cluster, use the existing environment:

```bash
source /vast/parcc/spack/sw/apps/linux-sapphirerapids/anaconda3-2023.09-0-ieilyrkph5mewqcum3ajc4odlt2vakri/etc/profile.d/conda.sh
conda activate /vast/projects/wherry/foundation-models-immuno/hhua/tools/perturbscope_env
python -m pip install -e /vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbFlow --no-deps
```

## 2. Run The Real Example Dataset

The repository includes a small real-data subset for testing:

```text
examples/data/adamson_2016_upr_360x1000.h5ad
```

It is a 360-cell by 1,000-gene subset of the scPerturb Adamson/Weissman 2016
pilot Perturb-seq dataset, balanced across one non-targeting/control-like guide
barcode and five UPR perturbations. The subset is intentionally small so the
quick workflow can run on a laptop or login node.

Prepare the example AnnData:

```bash
perturbflow prepare \
  --input examples/data/adamson_2016_upr_360x1000.h5ad \
  --output prepared/adamson_2016_upr.perturbflow.h5ad \
  --perturbation-col guide_gene \
  --control-labels non-targeting \
  --cell-state-col cell_state_hint
```

Run a fast analysis path:

```bash
perturbflow analyzer \
  --input prepared/adamson_2016_upr.perturbflow.h5ad \
  --output results/adamson_2016_upr_quickstart \
  --config configs/quickstart.json \
  --no-resume
```

This quick config runs:

```text
qc, preprocess, eda, score, effects, deg, report, bundle
```

Open the main result:

```text
results/adamson_2016_upr_quickstart/interactive_report.html
```

Other useful files:

- `results/adamson_2016_upr_quickstart/report.html`
- `results/adamson_2016_upr_quickstart/csv/deg_summary.csv`
- `results/adamson_2016_upr_quickstart/csv/effect_decomposition.csv`
- `results/adamson_2016_upr_quickstart/plots/qc_summary.png`
- `results/adamson_2016_upr_quickstart/bundle/manifest.json`

Rebuild the example subset from the public scPerturb source:

```bash
python scripts/make_example_subset.py
```

## 3. Prepare Your Own Data

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

## 4. Run The Full Pipeline

```bash
perturbflow analyzer \
  --input prepared/my_data.perturbflow.h5ad \
  --output results/my_run \
  --config configs/cluster_default.json \
  --resume
```

Useful step-only reruns:

```bash
perturbflow analyzer --input prepared/my_data.perturbflow.h5ad --output results/my_run --steps deg,report,bundle
perturbflow analyzer --input prepared/my_data.perturbflow.h5ad --output results/my_run --force-steps report --resume
perturbflow list-steps
```

## 5. Run On A Cluster

On a cluster, run the same `perturbflow analyzer` command inside your usual
job environment. For each dataset, set the input `.h5ad`, output directory,
config JSON, requested resources, and wall time according to your local
scheduler.

## 6. Use From Python

Use `PerturbFlowAPI` when another program needs to prepare data or launch an
analysis run directly:

```python
from perturbflow import PerturbFlowAPI

api = PerturbFlowAPI(config="configs/quickstart.json", perturbation_col="perturbation")

prepared = api.prepare(
    "raw/my_data.h5ad",
    "prepared/my_data.perturbflow.h5ad",
    control_labels="control,NT",
    cell_state_col="cell_type",
)

adata = api.analyze(
    prepared,
    "results/my_run",
    steps=["qc", "preprocess", "deg", "report", "bundle"],
)
```

## 7. Main Outputs

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

## 8. Expected Input Columns

Minimum:

- perturbation label column, passed with `--perturbation-col`

Optional but useful:

- cell-state or cluster column, passed with `--cell-state-col`
- QC metrics such as mitochondrial fraction or total counts

PerturbFlow will compute core QC/preprocessing features when they are missing.
