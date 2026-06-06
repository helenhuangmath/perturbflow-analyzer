# PerturbFlow

**Open infrastructure for perturbation biology.**

PerturbFlow is an open-source, AnnData-native platform for perturbation
experiments — Perturb-seq, pooled CRISPR screens, and single-cell multi-omics.
It standardizes data representation, reproducible workflows, mechanistic
(rewiring-aware) interpretation, and AI-ready outputs, and complements the
scverse ecosystem rather than replacing it.

The project is organized around three aims — community-standard infrastructure,
AI-ready perturbation biology, and community benchmarks. See the
[Vision and Roadmap](vision.md) for details.

The current release ships an analyzer workflow for AnnData input, interactive
reports, and agent-ready interpretation files. The `predictor` and `benchmark`
namespaces are reserved for the prediction and evaluation work on the roadmap.

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
