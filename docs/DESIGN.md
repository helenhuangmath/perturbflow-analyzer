# PerturbFlow Architecture

This document describes the software architecture behind PerturbFlow. It is
intended for developers and advanced users who want to understand how data moves
through the package, how outputs are organized, and where new functionality can
be added.

For step-level scientific details, see [`METHOD.md`](METHOD.md). For output
descriptions, see [`RESULT.md`](RESULT.md).

## Design Goals

PerturbFlow is designed around four goals:

- **AnnData-native input and output**: perturbation experiments are represented
  in `.h5ad` files with standardized `.obs`, `.var`, `.uns`, and `.obsm`
  fields.
- **Reproducible workflows**: analysis parameters are stored in JSON
  configuration files, and pipeline steps can be resumed or rerun selectively.
- **Structured results**: figures, tables, final AnnData, interactive reports,
  and viewer-ready bundles are written in predictable locations.
- **Extensible modes**: Analyzer, Predictor, Benchmarker, and Agent Interpretation share
  the same standardized data model.

## Package Layout

```text
perturbflow/
├── analyzer/    # QC, scoring, DEG, network, regulatory, C-score, reports
├── predictor/   # Prediction-facing namespace and model-ready interfaces
├── benchmark/   # Benchmarking-facing namespace and evaluation interfaces
├── data/        # Data preparation namespace
├── workflows/   # End-to-end workflow namespace
├── viz/         # Visualization/reporting namespace
├── api.py       # PerturbFlowAPI for programmatic access
├── ai.py        # Agent interpretation exporter
├── cli.py       # Main perturbflow command
└── workflow.py  # Python workflow helpers
```

The current release focuses on the Analyzer workflow. Predictor and Benchmarker
interfaces are present so downstream tools can build against stable package
namespaces as those modes expand.

## Processing Flow

```text
raw AnnData (.h5ad)
  |
  v
perturbflow prepare
  |
  v
standardized AnnData
  |
  v
perturbflow analyzer
  |
  v
tables + figures + final AnnData + reports + result bundle
  |
  v
perturbflow interpret
  |
  v
agent interpretation package
```

The preparation step standardizes perturbation labels, control labels, guide
fields, and optional cell-state annotations. The analyzer then runs a configured
sequence of modular steps.

## Analyzer Pipeline

The analyzer pipeline is implemented in
[`perturbflow/analyzer/pipeline.py`](https://github.com/helenhuangmath/PerturbFlow/blob/main/perturbflow/analyzer/pipeline.py).
Each major analysis step lives in its own module under `perturbflow/analyzer/`.

Default steps include:

```text
qc -> preprocess -> eda -> score -> effects -> trajectory -> programs
-> interaction -> state_enrich -> deg -> genenet -> tf_genenet
-> cscore -> regulatory -> report -> bundle
```

The pipeline supports:

- `--steps` for running a selected subset.
- `--force-steps` for regenerating selected outputs.
- `--resume` for continuing from completed checkpoints.
- JSON config files for reproducible parameter settings.

## Result Organization

A completed run writes a result directory with several output types:

```text
results/my_run/
├── csv/                    # Tables for DEG, C-score, networks, summaries
├── figures/                # Static plot files
├── bundle/                 # Viewer-ready JSON/CSV/Parquet-style artifacts
├── interactive_report.html # Browser-based interactive report
├── report.html             # Static summary report
├── final.h5ad              # Final AnnData output when enabled
├── summary.json            # Run-level summary
└── checkpoint.json         # Completed pipeline-step state
```

The `bundle/` directory is the main boundary between the analysis engine and
interactive viewers. Viewers should read bundle metadata first and hide panels
whose artifacts are not available.

## Bundle Schema

The bundle is designed to be stable enough for downstream web apps and external
programs. The manifest records the schema version, PerturbFlow version, creation
time, artifact paths, and available result types.

Example:

```json
{
  "schema_version": "1.0",
  "perturbflow_version": "0.1.0",
  "created_utc": "2026-04-26T00:00:00Z",
  "format": "parquet",
  "artifacts": {
    "metadata": "metadata.json",
    "embeddings": ["embeddings/umap.parquet"],
    "perturbations_summary": "perturbations_summary.parquet",
    "search_index": "search_index.json"
  }
}
```

## Modes

### Analyzer

Analyzer is the current core mode. It standardizes input data, runs the analysis
pipeline, and writes tables, figures, reports, final AnnData, and result
bundles.

### Predictor

Predictor is the prediction-facing mode. It is intended to consume standardized
PerturbFlow inputs and model-ready exports, then write prediction results in a
format that can be compared against Analyzer outputs.

### Benchmarker

Benchmarker is the evaluation-facing mode. It is intended to compare models or
analysis methods using reproducible metrics, including baseline-aware,
distance-aware, and rewiring-aware summaries.

### AI

Agent Interpretation is implemented through `perturbflow interpret` and
`perturbflow.ai.write_agent_handoff`. It exports compact Markdown and JSON files
that summarize a run without including raw count matrices.

## Extension Points

New analysis functionality should usually be added as a module under
`perturbflow/analyzer/` and wired into the pipeline as a named step. A good step
should:

- Read standardized fields from AnnData.
- Write explicit CSV/JSON/figure outputs.
- Store compact summary values in `.uns` when later steps need them.
- Be optional, resumable, and safe to rerun.
- Add bundle artifacts when the output should be used by the interactive report
  or external viewers.

New external integrations should prefer `PerturbFlowAPI` or the result bundle
instead of reading internal intermediate objects directly.

## Public Documentation Role

This file is not required for running PerturbFlow. It is a developer-facing
overview of the architecture. Most users should start with the README and Quick
Start, then use `METHOD.md` and `RESULT.md` when they need method or output
details.
