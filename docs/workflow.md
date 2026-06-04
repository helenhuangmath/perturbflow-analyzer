# Analyzer Workflow

PerturbFlow Analyzer runs a resumable pipeline with checkpointing.

```text
qc -> preprocess -> eda -> score -> effects -> trajectory -> programs
-> interaction -> state_enrich -> deg -> genenet -> tf_genenet
-> cscore -> regulatory -> report -> bundle
```

## Step Control

Run selected steps:

```bash
perturbflow run --input prepared.h5ad --output results/run1 --steps deg,report,bundle
```

Force report regeneration:

```bash
perturbflow run --input prepared.h5ad --output results/run1 --force-steps report --resume
```

Clear a checkpoint from a specific step:

```bash
perturbflow run --input prepared.h5ad --output results/run1 --clear-from deg --resume
```
