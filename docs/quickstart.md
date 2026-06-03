# Quick Start

## Install

```bash
python -m pip install -e ".[bundle]"
```

## Prepare Data

```bash
perturbflow prepare \
  --input my_raw_data.h5ad \
  --output prepared/my_data.perturbflow.h5ad \
  --perturbation-col guide_gene \
  --control-labels control,non-targeting,NT \
  --cell-state-col leiden
```

## Run Analysis

```bash
perturbflow run \
  --input prepared/my_data.perturbflow.h5ad \
  --output results/my_run \
  --config configs/cluster_default.json \
  --resume
```

## Review Results

Open:

```text
results/my_run/interactive_report.html
```

## Create Agent Handoff

```bash
perturbflow interpret \
  --results results/my_run \
  --project-name "My Perturb-seq screen"
```
