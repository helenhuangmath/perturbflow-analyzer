# PerturbFlow

PerturbFlow is a modular Perturb-seq framework. It currently ships an analyzer workflow for AnnData input, interactive reports, and agent-ready interpretation files, with a predictor namespace reserved for future perturbation response prediction models.

## Core Workflow

```bash
perturbflow prepare --input raw.h5ad --output prepared/data.h5ad --perturbation-col guide_gene
perturbflow analyzer --input prepared/data.h5ad --output results/run1 --resume
perturbflow interpret --results results/run1 --project-name "My Perturb-seq screen"
```

## What You Get

- QC and preprocessing summaries.
- Differential expression and enrichment-style summaries.
- Perturbation effect decomposition.
- Cell-state and trajectory effects.
- Gene-network and TF-network rewiring.
- Static and interactive HTML reports.
- Agent handoff files for LLM-assisted interpretation.

PerturbFlow keeps LLM integration explicit: it prepares interpretation context, but it does not upload data automatically.
