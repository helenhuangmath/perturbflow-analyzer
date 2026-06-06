# Quick Start

## Install

```bash
python -m pip install -e ".[bundle]"
```

## Run The Real Example Dataset

The repository includes a small real-data subset:

```text
examples/data/adamson_2016_upr_360x1000.h5ad
```

It contains 360 cells and 1,000 genes from the scPerturb Adamson/Weissman 2016
pilot Perturb-seq dataset.

```bash
perturbflow prepare \
  --input examples/data/adamson_2016_upr_360x1000.h5ad \
  --output prepared/adamson_2016_upr.perturbflow.h5ad \
  --perturbation-col guide_gene \
  --control-labels non-targeting \
  --cell-state-col cell_state_hint
```

```bash
perturbflow analyzer \
  --input prepared/adamson_2016_upr.perturbflow.h5ad \
  --output results/adamson_2016_upr_quickstart \
  --config configs/quickstart.json \
  --no-resume
```

Open:

```text
results/adamson_2016_upr_quickstart/interactive_report.html
```

Rebuild the subset from the public source:

```bash
python scripts/make_example_subset.py
```

## Prepare Your Own Data

```bash
perturbflow prepare \
  --input my_raw_data.h5ad \
  --output prepared/my_data.perturbflow.h5ad \
  --perturbation-col guide_gene \
  --control-labels control,non-targeting,NT \
  --cell-state-col leiden
```

## Run Full Analysis

```bash
perturbflow analyzer \
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
