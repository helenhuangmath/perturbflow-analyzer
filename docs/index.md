# PerturbFlow

**Open infrastructure for perturbation biology.**

PerturbFlow is an open-source, AnnData-native platform for perturbation
experiments — Perturb-seq, pooled CRISPR screens, and single-cell multi-omics.
It standardizes data representation, reproducible workflows, mechanistic
(rewiring-aware) interpretation, and structured outputs, and complements the
scverse ecosystem rather than replacing it.

The project is organized around three aims: community-standard infrastructure,
model-ready perturbation biology, and community benchmarks. See the
[Vision](vision.md) page for details.

The current release ships an analyzer workflow for AnnData input, interactive
reports, and structured interpretation files. The `predictor` and `benchmark`
namespaces are reserved for prediction and evaluation work.

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
- Structured handoff files for downstream interpretation.

PerturbFlow prepares interpretation context, but it does not upload data automatically.
